"""Polarity/gain negative controls for phase-aware factor comparisons; no FDTD."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from analyze_hs4_v4_factor_controls import compare
from hs_capsule_identity import sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new array-check directory required')
    t = np.linspace(150, 240, 901)
    z = np.exp(-((t-190)/12)**2)*np.exp(1j*.55*t)
    reference = SimpleNamespace(complex_envelope=z, real_bandpass=2*z.real)
    inverse = SimpleNamespace(complex_envelope=-z, real_bandpass=-2*z.real)
    doubled = SimpleNamespace(complex_envelope=2*z, real_bandpass=4*z.real)
    reversed_metrics = compare(reference, inverse, t)['underground']
    gain_metrics = compare(reference, doubled, t)['underground']
    if not (abs(reversed_metrics['signed_waveform_relative_L2']-2) < 1e-12
            and abs(reversed_metrics['complex_envelope_relative_L2']-2) < 1e-12
            and reversed_metrics['envelope_magnitude_relative_L2'] == 0
            and abs(reversed_metrics['complex_coherence_real']+1) < 1e-12
            and abs(gain_metrics['amplitude_norm_change_dB']-20*np.log10(2)) < 1e-12):
        raise ValueError('phase or gain counterexample failed')
    result = {'status': 'PASS', 'checks': 2, 'code_sha256': sha256(__file__),
              'metric_code_sha256': sha256(Path(__file__).with_name('analyze_hs4_v4_factor_controls.py')),
              'polarity_reversal': reversed_metrics, 'doubled_amplitude': gain_metrics,
              'solver_called': False, 'scope': 'Array counterexamples only, no physical validation.'}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('PASS:2 phase/gain counterexamples')


if __name__ == '__main__':
    main()
