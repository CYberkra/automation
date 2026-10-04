"""Fixed-window 12m vs 16m domain contribution for HS4 3d-xwide-centre.

Recomputes, from the two archived result arrays, how much widening the X
domain 12m to 16m (interior geometry rigidly shifted +2m) changes the centre
station products in the fixed diagnostic windows (160-180, 180-220,
160-220 ns). Paired differences share the boundary geometry, so their change
isolates the residual boundary contribution; the raw response change includes
boundary-sensitive late arrivals and is reported separately.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CHECKS = ROOT / 'artifacts/research_checks'
BASE = CHECKS / 'target_hs4_3d-centre_results_r1'
XWIDE = CHECKS / 'target_hs4_3d-xwide-centre_results_r1'
OUT = XWIDE / 'window_boundary_contribution.json'
WINDOWS = {'early': (160.0, 180.0), 'later': (180.0, 220.0), 'full': (160.0, 220.0)}
PRODUCTS = ('centre_rough_minus_flat', 'centre_rough_minus_halfspace', 'centre_raw_rough')


def relative(a, b, mask):
    norm = np.linalg.norm(a[mask])
    return float(np.linalg.norm((b - a)[mask]) / norm) if norm > 0 else None


def main():
    if OUT.exists():
        raise ValueError('result file already exists; preserve history')
    a = np.load(BASE / 'arrays.npz')
    b = np.load(XWIDE / 'arrays.npz')
    if not np.array_equal(a['time_ns'], b['time_ns']):
        raise ValueError('time axes differ')
    t = a['time_ns']
    record = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
              'base_results': str(BASE.relative_to(ROOT)),
              'xwide_results': str(XWIDE.relative_to(ROOT)),
              'windows_ns': {k: list(v) for k, v in WINDOWS.items()},
              'products': {}, 'code_sha256': __import__('hashlib').sha256(
                  Path(__file__).read_bytes()).hexdigest()}
    for name in PRODUCTS:
        entry = {}
        for label, (lo, hi) in WINDOWS.items():
            mask = (t >= lo) & (t <= hi)
            entry[label] = {
                'signed_relative_L2': relative(a[name + '_signed'], b[name + '_signed'], mask),
                'complex_relative_L2': relative(a[name + '_complex'], b[name + '_complex'], mask)}
        record['products'][name] = entry
    OUT.write_text(json.dumps(record, ensure_ascii=False, indent=2, allow_nan=False) + '\n',
                   encoding='utf-8')
    print(json.dumps(record['products'], indent=2))


if __name__ == '__main__':
    main()
