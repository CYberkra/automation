"""Compare archived native-clock spectra into a new directory; no solver."""
import argparse
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New directory; must not exist')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = ROOT/'artifacts/research_checks/2026-09-25_hsg_clock'
    rows = {}
    for r, path in [(1, root/'ratio1'), (3, ROOT/'artifacts/research_checks/2026-09-25_hsg_smoke_analysis_r2'), (5, root/'ratio5')]:
        with np.load(path/'arrays.npz') as a:
            f = a['frequency_Hz']
            phase = np.angle(a['hsg_inside']/a['reference_inside'], deg=True)
        metadata = json.loads((path/'results.json').read_text(encoding='utf-8'))
        dt = metadata['metadata']['reference_inside']['dt_s']
        lag = (r-1)/(2*r)*dt
        predicted = -360*f*lag
        rows[str(r)] = {'predicted_lag_ps': lag*1e12,
                       'measured_lag_at_170MHz_ps': float(-phase[-1]/(360*f[-1])*1e12),
                       'phase_at_170MHz_deg': float(phase[-1]),
                       'predicted_phase_at_170MHz_deg': float(predicted[-1]),
                       'max_abs_phase_prediction_residual_deg': float(max(abs(phase-predicted))),
                       'comparisons': metadata['comparisons']}
    result = {'ratios': rows, 'phase_correction_applied': False, 'solver_modified': False,
              'interpretation': 'Ratio dependence tests the lag hypothesis; not proof of a solver defect or validated timestamp correction. Ratio1 also disables filtering/PML/interpolation.'}
    with (args.output/'comparison.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
