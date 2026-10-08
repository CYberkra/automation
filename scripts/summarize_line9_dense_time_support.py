"""Report both basal and wider windows without equating L2 ratios to energy shares."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    r = json.loads((a.public / 'analysis.json').read_text('utf-8'))
    audit = json.loads((a.public / 'independent_audit.json').read_text('utf-8'))
    assert audit['analysis_sha256'] == sha(a.public / 'analysis.json') and audit['status'].startswith('PASS')
    assert r['numerical_sha256'] == sha(a.numerical)
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]; full = h['full'][:]
        cumulative = h['cumulative_early_response'][:]; ids = list(h['ids'].asstr()[:])
    t = np.arange(4008) / (4008 * 300000.); centers = r['cumulative_centers_native_ns']
    windows = dict(basal_common=r['basal_common_gate_ns'], wide=[300., 450.], late=[450., 1100.])
    values = {}; errors = []
    for name, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean(); values[name] = {}
        allz = np.concatenate([full, cumulative.reshape(501, -1)], axis=1)
        padded = np.zeros((4008, allz.shape[1]), complex); padded[:501] = allz * w[:, None]
        fft = np.fft.ifft(padded, axis=0) * 8 * np.exp(2j * np.pi * 20e6 * t[:, None])
        for tag, (lo, hi) in windows.items():
            k = (t * 1e9 >= lo) & (t * 1e9 <= hi)
            direct = np.exp(2j * np.pi * t[k, None] * f) @ (allz * w[:, None]) / 501
            err = float(np.linalg.norm(direct - fft[k]) / np.linalg.norm(direct)); assert err < 1e-9; errors.append(err)
            q = direct[:, :4]; cum = direct[:, 4:].reshape(k.sum(), 4, len(centers)); rows = []
            for j, sid in enumerate(ids):
                for config in r['configurations']:
                    c0, c1, c2 = [cum[:, j, centers.index(c)] for c in config]
                    bands = np.column_stack([c0, c1 - c0, c2 - c1, q[:, j] - c2])
                    norm = np.linalg.norm; denominator = norm(q[:, j])
                    rows.append(dict(id=sid, native_cutoff_centers_ns=config,
                        band_L2_over_full=[float(norm(bands[:, b]) / denominator) for b in range(4)],
                        early_combined_over_full_L2=float(norm(c2) / denominator),
                        full_L2=float(denominator), full_envelope_peak_ns=float(t[k][np.argmax(abs(q[:, j]))] * 1e9)))
            values[name][tag] = rows
    out = dict(status='AUDITED_SPECTRA_MULTIWINDOW_TEMPORAL_ACCOUNTING_NOT_PATH_CERTIFICATION',
        script_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        independent_audit_sha256=sha(a.public / 'independent_audit.json'), numerical_sha256=sha(a.numerical),
        windows_ns=windows, metrics=values, independent_direct_vs_padded_IFFT_relative_L2_max=max(errors),
        limits='Wide/late windows added as exploratory diagnostics; basal common window remains inherited predeclared gate. Time-band norm ratios are not energy or physical-path fractions. Strong late-native packets outside basal gate must not be described as entirely early sidelobes. No clean-label or processing change.')
    a.out.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=out['status'], inverse_error=max(errors),
        central={name: {tag: [q for q in rows if q['native_cutoff_centers_ns'] == [40., 130., 260.]] for tag, rows in bywindow.items()} for name, bywindow in values.items()})))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
