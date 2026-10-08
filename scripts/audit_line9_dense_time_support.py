"""Independent raw scalar-frequency DFT and inverse audit of time-support accounting."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    r = json.loads((a.public / 'analysis.json').read_text('utf-8'))
    assert r['numerical_sha256'] == sha(a.numerical)
    assert r['script_sha256'] == sha(Path(__file__).with_name('diagnose_line9_dense_time_support.py'))
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]; names = list(h['ids'].asstr()[:])
        stored = h['partition_response'][:]; full = h['full'][:]
    np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
    direct = np.empty_like(stored); errors = []
    hashes = {row['id']: row['native_sha256'] for row in r['native']}
    for j, sid in enumerate(names):
        path = a.source / (sid + '.h5'); assert sha(path) == hashes[sid]
        with h5py.File(path) as h:
            assert h['rxs/rx1/Ez'].dtype == h['srcs/src1/excitation/samples'].dtype == np.float64
            x = h['rxs/rx1/Ez'][:]; s = h['srcs/src1/excitation/samples'][:]
            tr = np.arange(len(x)) * h.attrs['dt'] + h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
            ts = np.arange(len(s)) * h.attrs['dt'] + h['srcs/src1/excitation'].attrs['TimeSampleOffset']
        for c, center in enumerate(r['split_centers_native_ns']):
            lower = center - r['transition_width_ns'] / 2
            upper = center + r['transition_width_ns'] / 2
            late = np.zeros(len(x)); late[tr * 1e9 >= upper] = 1
            crossing = (tr * 1e9 > lower) & (tr * 1e9 < upper)
            late[crossing] = np.sin(np.pi / 2 * (tr[crossing] * 1e9 - lower) / (upper - lower)) ** 2
            # Boundaries must leave strictly early/late samples exactly intact.
            assert np.all(late[tr * 1e9 <= lower] == 0) and np.all(late[tr * 1e9 >= upper] == 1)
            for first in range(0, 501, 11):
                ff = f[first:first + 11, None]
                receiver = np.exp(-2j * np.pi * ff * tr)
                source = np.exp(-2j * np.pi * ff * ts) @ s
                values = receiver @ np.column_stack([x * (1 - late), x * late]) / source[:, None] / .025
                direct[first:first + len(ff), j, c] = values
            for k in range(2):
                err = float(np.linalg.norm(direct[:, j, c, k] - stored[:, j, c, k]) / np.linalg.norm(stored[:, j, c, k]))
                assert err < 1e-8; errors.append(err)
    t = np.arange(4008) / (4008 * 300000.); lo, hi = r['basal_common_gate_ns']
    keep = (t * 1e9 >= lo) & (t * 1e9 <= hi); kernel = np.exp(2j * np.pi * t[keep, None] * f)
    checks = {}; worst = 0
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= np.mean(w); rows = {(q['id'], q['split_center_native_ns']): q for q in r['metrics'][window]}
        count = 0
        for j, sid in enumerate(names):
            q = kernel @ (full[:, j] * w) / 501
            for c, center in enumerate(r['split_centers_native_ns']):
                z = kernel @ (direct[:, j, c] * w[:, None]) / 501
                e, late = z.T; norm = np.linalg.norm
                values = dict(early_over_total_L2=float(norm(e) / norm(q)), late_over_total_L2=float(norm(late) / norm(q)),
                    early_late_normalized_real_inner_product=float(np.vdot(e, late).real / (norm(e) * norm(late))))
                for name, value in values.items():
                    err = abs(value - rows[sid, center][name]) / max(1, abs(value))
                    assert err < 1e-8; worst = max(worst, err); count += 1
                assert norm(e + late - q) / norm(q) < 1e-8
        checks[window] = count
    out = dict(status='PASS_INDEPENDENT_RAW_TIME_PARTITION_DFT_INVERSE', auditor_sha256=sha(__file__),
        analysis_sha256=sha(a.public / 'analysis.json'), numerical_sha256=sha(a.numerical),
        raw_partition_DFT_relative_L2_max=max(errors), metric_scaled_error_max=worst,
        scalar_checks=checks, strict_support_checks=True,
        limits='Verifies native-time linear accounting and exact frequency/time conventions, not path identity or a production suppression method.')
    a.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(out))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
