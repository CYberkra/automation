"""CPU spatial-split/interval sensitivity from saved migration ridges."""
import argparse
import json
from pathlib import Path
import numpy as np
from hs4_analysis_metrics import matched_profile_metrics
from hs_capsule_identity import sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--arrays', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new output directory required')
    rows = []
    with np.load(args.arrays) as data:
        x, truth = data['metric_x_m'], data['truth_z_m']
        for name in data.files:
            if not name.endswith('_ridge_z_m'):
                continue
            for trim in (0, 10):
                idx = slice(trim, -trim if trim else None)
                xx, tt, rr = x[idx], truth[idx], data[name][idx]
                for scale in (1., 1.577855042105263, 2.):
                    rows.append({'case': name.removesuffix('_ridge_z_m'),
                                 'metric_x_range_m': [float(xx[0]), float(xx[-1])],
                                 'metrics': matched_profile_metrics(rr, tt, float(np.mean(np.diff(xx))), scale)})
    result = {'status': 'COMPLETED_METRIC_SENSITIVITY_DIAGNOSTIC',
              'input_arrays_sha256': sha256(args.arrays), 'code_sha256': sha256(__file__),
              'metrics_code_sha256': sha256(Path(__file__).with_name('hs4_analysis_metrics.py')),
              'rows': rows, 'physical_resolution_bound_certified': False}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for case in sorted({r['case'] for r in rows}):
        selected = [r['metrics']['highpass_matched'] for r in rows if r['case'] == case]
        print(case, 'highpass corr range', min(m['corr'] for m in selected), max(m['corr'] for m in selected),
              'relative error range', min(m['relative_L2_error'] for m in selected),
              max(m['relative_L2_error'] for m in selected))


if __name__ == '__main__':
    main()
