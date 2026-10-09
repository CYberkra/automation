"""Direct-sum recomputation of paired interference; imports no analysis helpers."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read = lambda p: json.loads(p.read_text('utf-8'))
    assert not a.out.exists()
    report = read(a.public / 'analysis.json')
    source = read(a.source_public / 'analysis.json')
    audit = read(a.source_public / 'independent_audit.json')
    m = read(a.package / 'manifest.json')
    assert source['numerical_sha256'] == report['numerical_sha256'] == sha(a.numerical)
    assert report['analysis_sha256'] == audit['analysis_sha256'] == sha(a.source_public / 'analysis.json')
    assert audit['status'].startswith('PASS') and report['manifest_sha256'] == sha(a.package / 'manifest.json')
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]
        np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
        ids = h['ids'].asstr()[:].tolist()
        z = h['response'][:]
    time = np.arange(4008) / (4008 * 300000.)
    tn = time * 1e9
    count = 0
    worst = 0.
    phase_comparisons = []
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean()
        assert len(report['metrics'][window]) == 8
        keys = set()
        for record in report['metrics'][window]:
            chain, variant = record['chainage_m'], record['variant']
            key = (chain, variant, record['window'])
            assert key not in keys
            keys.add(key)
            station = next(s for s in m['stations'] if s['chainage_m'] == chain)
            bounds = [300., 450.] if record['window'] == 'fixed_wide' else [min(station['templates'][v]['basal_gate_ns'][0] for v in ['high', 'low']), max(station['templates'][v]['basal_gate_ns'][1] for v in ['high', 'low'])]
            assert bounds == record['gate_ns']
            prefix = f'{variant}_x{round(chain * 100):05d}_'
            index = [ids.index(prefix + role) for role in ['H0', 'H1']]
            keep = (tn >= bounds[0]) & (tn <= bounds[1])
            x = np.exp(2j * np.pi * tn[keep, None] * 1e-9 * f) @ (z[:, index] * w[:, None]) / 501
            padded = np.zeros((4008, 2), complex)
            padded[:501] = z[:, index] * w[:, None]
            x_fft = (np.fft.ifft(padded, axis=0) * 8 * np.exp(2j * np.pi * 20e6 * time[:, None]))[keep]
            b, q = x[:, 0], x[:, 1]
            d = q - b
            # Explicit real/imag sums, rather than the generating vdot formula.
            nb2 = np.sum(b.real**2 + b.imag**2)
            nd2 = np.sum(d.real**2 + d.imag**2)
            nq2 = np.sum(q.real**2 + q.imag**2)
            real = np.sum(d.real * b.real + d.imag * b.imag)
            imag = np.sum(d.real * b.imag - d.imag * b.real)
            values = dict(H0_norm=np.sqrt(nb2), delta_norm=np.sqrt(nd2), H1_norm=np.sqrt(nq2),
                H0_over_delta=np.sqrt(nb2 / nd2), H0_delta_inner_product_real=real, H0_delta_inner_product_imag=imag,
                H0_vs_delta_inner_phase_deg=np.arctan2(imag, real) * 180 / np.pi,
                complex_coherence_magnitude=np.hypot(real, imag) / np.sqrt(nb2 * nd2),
                cross_term_over_incoherent_squared_norm=2 * real / (nb2 + nd2),
                squared_norm_identity_relative_error=abs(nq2 - nb2 - nd2 - 2 * real) / nq2,
                H1_peak_ns=tn[keep][np.argmax(abs(q))], delta_peak_ns=tn[keep][np.argmax(abs(d))])
            # Weak late arrays magnify FFT/direct-sum phase differences. Use the
            # observed inverse-array discrepancy and a length-dependent dot
            # rounding envelope,not an arbitrary relaxed phase tolerance.
            bf = x_fft[:, 0]
            df = x_fft[:, 1] - bf
            norm = np.linalg.norm
            delta_b, delta_d = norm(bf - b), norm(df - d)
            unit_roundoff = np.finfo(float).eps / 2
            operations = 8 * len(b) + 64
            gamma = operations * unit_roundoff / (1 - operations * unit_roundoff)
            inner_bound = norm(d) * delta_b + norm(b) * delta_d + delta_d * delta_b
            inner_bound += gamma * (norm(d) * norm(b) + norm(df) * norm(bf))
            assert inner_bound < np.hypot(real, imag), 'Phase not resolved by CPU-method agreement'
            phase_bound_deg = np.arcsin(inner_bound / np.hypot(real, imag)) * 180 / np.pi
            phase_error_deg = abs(np.angle(np.exp(1j * (values['H0_vs_delta_inner_phase_deg'] - record['H0_vs_delta_inner_phase_deg']) * np.pi / 180))) * 180 / np.pi
            assert phase_error_deg <= phase_bound_deg
            phase_comparisons.append(dict(window=window, chainage_m=chain, variant=variant,
                gate=record['window'], observed_phase_difference_deg=float(phase_error_deg),
                CPU_inverse_agreement_phase_envelope_deg=float(phase_bound_deg),
                inverse_H0_difference_norm=float(delta_b), inverse_delta_difference_norm=float(delta_d)))
            for name, value in values.items():
                if name == 'H0_vs_delta_inner_phase_deg':
                    count += 1
                    continue
                error = float(abs(value - record[name]) / max(1, abs(value)))
                assert error < 1e-8, (key, name, value, record[name])
                count += 1
                worst = max(worst, error)
        assert keys == {(chain, variant, gate) for chain in [190., 187.5] for variant in ['high', 'low'] for gate in ['basal_union', 'fixed_wide']}
    result = dict(status='PASS_DIRECT_INVERSE_AND_EXPLICIT_COMPLEX_INTERFERENCE_SUMS',
        auditor_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        source_native_DFT_audit_sha256=sha(a.source_public / 'independent_audit.json'),
        scalar_checks=count, max_scaled_error=worst,
        phase_method_agreement=phase_comparisons,
        limits='Reuses hash-bound independently raw-DFT-audited spectra; independently checks inverse/sums/windows. Phase envelope uses actual FFT/direct-sum array discrepancy plus dot-rounding envelope,not a FDTD error certificate or evaluation direction bound. No newnative3D/field or path isolation.')
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source-public', 'public', 'package', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
