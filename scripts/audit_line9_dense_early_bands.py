"""Independently verify raw four-band transforms and all cutoff sensitivity metrics."""
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
    assert r['script_sha256'] == sha(Path(__file__).with_name('diagnose_line9_dense_early_bands.py'))
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]; ids = list(h['ids'].asstr()[:])
        cumulative = h['cumulative_early_response'][:]; full = h['full'][:]
    np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
    centers = r['cumulative_centers_native_ns']; raw = np.empty_like(cumulative); worst_dft = 0.
    hashes = {q['id']: q['native_sha256'] for q in r['native']}
    for j, sid in enumerate(ids):
        path = a.source / (sid + '.h5'); assert sha(path) == hashes[sid]
        with h5py.File(path) as h:
            x = h['rxs/rx1/Ez'][:]; s = h['srcs/src1/excitation/samples'][:]
            assert x.dtype == s.dtype == np.float64
            tr = np.arange(len(x)) * h.attrs['dt'] + h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
            ts = np.arange(len(s)) * h.attrs['dt'] + h['srcs/src1/excitation'].attrs['TimeSampleOffset']
        weights = []
        for center in centers:
            u = np.clip((tr * 1e9 - (center - 10)) / 20, 0, 1)
            weights.append(1 - np.sin(np.pi * u / 2) ** 2)
        weighted = x[:, None] * np.array(weights).T
        for first in range(0, 501, 13):
            ff = f[first:first + 13, None]
            raw[first:first + len(ff), j] = (np.exp(-2j * np.pi * ff * tr) @ weighted) / (np.exp(-2j * np.pi * ff * ts) @ s)[:, None] / .025
        for c in range(len(centers)):
            err = np.linalg.norm(raw[:, j, c] - cumulative[:, j, c]) / np.linalg.norm(cumulative[:, j, c])
            assert err < 1e-8; worst_dft = max(worst_dft, float(err))
    t = np.arange(4008) / (4008 * 300000.); lo, hi = r['basal_common_gate_ns']
    mask = (t * 1e9 >= lo) & (t * 1e9 <= hi); e = np.exp(2j * np.pi * t[mask, None] * f)
    worst_metric = 0.; checks = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean(); count = 0
        records = {(q['id'], tuple(q['native_cutoff_centers_ns'])): q for q in r['metrics'][window]}
        for j, sid in enumerate(ids):
            q = e @ (full[:, j] * w) / 501
            cum = e @ (raw[:, j] * w[:, None]) / 501
            for config in r['configurations']:
                c0, c1, c2 = [cum[:, centers.index(c)] for c in config]
                bands = [c0, c1 - c0, c2 - c1, q - c2]
                actual = [float(np.linalg.norm(b) / np.linalg.norm(q)) for b in bands]
                recorded = records[sid, tuple(config)]['band_L2_over_full']
                error = max(abs(np.array(actual) - recorded) / np.maximum(1, abs(np.array(recorded))))
                assert error < 1e-8; worst_metric = max(worst_metric, float(error)); count += 4
                assert np.linalg.norm(sum(bands) - q) / np.linalg.norm(q) < 1e-10
        assert len(records) == len(ids) * len(r['configurations']); checks[window] = count
    out = dict(status='PASS_RAW_DFT_DIRECT_INVERSE_FOUR_BANDS_AND_CUTOFF_SENSITIVITY',
        auditor_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        native_DFT_relative_L2_max=worst_dft, scaled_metric_error_max=worst_metric,
        scalar_checks=checks, limits='Arithmetic and native support only; no independent physical path separation or field calibration.')
    a.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8'); print(json.dumps(out))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
