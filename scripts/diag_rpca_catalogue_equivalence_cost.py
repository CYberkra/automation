"""Diagnostic (not archived): RPCA catalogue-path equivalence + cost timing.

1) Equivalence: apply_configuration(x, B7/B8/B9_G1_BG, catalogue_version='0.2')
   must reproduce study_t3_operator_effects_rpca.py bit-for-bit (same
   rpca_control call: lam = factor/sqrt(max(shape))).
2) Cost: wall time per family x lambda (median of 3), machine-dependent
   diagnostic for the addendum doc only - NOT a reward term, not archived.

Run: artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe
     scripts/diag_rpca_catalogue_equivalence_cost.py
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from research_operator_contract import apply_configuration, rpca_control
from study_t3_damage_ladder import MOTHERS, load_bscan

IDS = {0.5: 'B7_G1_BG', 1.0: 'B8_G1_BG', 2.0: 'B9_G1_BG'}


def main():
    for fam, mother, _geo in MOTHERS:
        s, _t = load_bscan(mother)
        x = np.ascontiguousarray(s.T)
        lam0 = 1.0 / (max(x.shape) ** 0.5)
        for f, cid in IDS.items():
            y_direct, _ = rpca_control(x, lam0 * f)
            r = apply_configuration(x, cid, catalogue_version='0.2')
            diff = float(np.max(np.abs(r['output'] - y_direct)))
            times = []
            for _ in range(3):
                t0 = time.perf_counter()
                rpca_control(x, lam0 * f)
                times.append(time.perf_counter() - t0)
            times.sort()
            print(f'{fam:6s} {cid} lam{f:g}x | max|catalogue-direct|={diff:.3e} '
                  f'| t_med={times[1]:.3f}s t_min={times[0]:.3f}s t_max={times[2]:.3f}s '
                  f'| shape={x.shape}')
            assert diff == 0.0, 'catalogue path diverges from direct rpca_control'
    print('equivalence: all 9 outputs bit-identical')


if __name__ == '__main__':
    main()
