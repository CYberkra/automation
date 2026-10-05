"""Independent saved-product acceptance: full tone sum, carrier and amplitude."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(results):
    summary_path = results/'summary.json'
    summary = json.loads(summary_path.read_text('utf-8'))
    if summary['status'] != 'PASS_SFCW_PROCESSING_NOT_FIELD_VALIDATION':
        raise ValueError('completed processing required')
    for name, identity in summary['products'].items():
        if sha(results/name) != identity:
            raise ValueError('saved product changed: '+name)
    generator = Path(__file__).with_name('reconstruct_slope_fine2m_sfcw_v0_1.py')
    if sha(generator) != summary['script_sha256']:
        raise ValueError('generator changed since reconstruction')
    checks = {}
    with np.load(results/'sfcw_profiles.npz') as z:
        f, t = z['frequency_Hz'], z['sfcw_time_ns']*1e-9
        np.testing.assert_array_equal(f, np.linspace(20e6, 170e6, 501))
        assert z['spectrum_ideal'].shape == (501, 13)
        # Sample the WHOLE periodic interval, including wrapped sidelobes.
        indices = np.unique(np.r_[np.arange(0, len(t), 43), np.arange(240)])
        tone_phase = np.exp(2j*np.pi*t[indices, None]*f)
        for window in ('rectangular', 'hann'):
            weights = np.ones(501) if window == 'rectangular' else np.hanning(501)
            weights /= weights.mean()
            np.testing.assert_allclose(z[window+'_weights'], weights, rtol=1e-14, atol=1e-14)
            for name, spectrum in [('slope_rough', z['spectrum_slope_rough']),
                                   ('slope_fullcover', z['spectrum_slope_fullcover']),
                                   ('ideal', z['spectrum_ideal'])]:
                # No FFT/CZT or official reconstructor: sum physical-frequency tones.
                expected_complex = tone_phase@(spectrum*weights[:, None])/501
                expected_real = 2*expected_complex.real
                saved = z[f'{window}_{name}_real'][indices]
                error = float(np.linalg.norm(saved-expected_real)/np.linalg.norm(expected_real))
                if error > 1e-9:
                    raise ValueError('physical carrier/amplitude reconstruction mismatch')
                saved_complex = z[f'{window}_{name}_complex_envelope'][indices]*np.exp(2j*np.pi*f[0]*t[indices, None])
                phase_error = float(np.linalg.norm(saved_complex-expected_complex)/np.linalg.norm(expected_complex))
                if phase_error > 1e-9:
                    raise ValueError('complex phase was not retained')
                checks[f'{window}_{name}'] = dict(real_tone_sum_relative_L2=error,
                                                  complex_tone_sum_relative_L2=phase_error)
        # A negative control must reject the historical 95 MHz re-modulation.
        expected = 2*np.real(tone_phase@(z['spectrum_ideal']*z['hann_weights'][:, None])/501)
        wrong = 2*np.real(z['hann_ideal_complex_envelope'][indices]*np.exp(2j*np.pi*95e6*t[indices, None]))
        wrong_error = float(np.linalg.norm(wrong-expected)/np.linalg.norm(expected))
        if wrong_error < .1:
            raise ValueError('tone-sum check cannot distinguish historical carrier bug')
        checks['negative_control_95MHz_carrier_relative_L2'] = wrong_error
    certificate = dict(status='PASS_SAVED_PRODUCT_INDEPENDENT_TONE_SUM',
        checks=checks, checker_sha256=sha(__file__), summary_sha256=sha(summary_path),
        all_13_stations_checked=True, tested_delay_samples_per_station=len(indices),
        source_code_sha256={name: sha(Path(__file__).with_name(name)) for name in
            ('hs4_slope_bscan_fine2m_v0_1.py', 'research_operator_contract.py')},
        scope='Processing arithmetic/phase/amplitude validation only, not FDTD convergence or field validation')
    path = results/'independent_verification.json'
    if path.exists():
        raise ValueError('fresh verification required')
    path.write_text(json.dumps(certificate, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(certificate, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    main(p.parse_args().results)
