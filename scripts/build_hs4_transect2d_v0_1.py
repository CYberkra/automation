"""HS4 transect 2D fast-study generator v0.1 (2026-10-02).

Kinematic companion to the 3D HS4 rough-interface standard piece.
- Invariant axis: y (single cell). Profile axis: x (maps to 3D line at bi=6, x=1.625 m).
- Interface: binned transect z0tab[6, :] from hs4_interface_binned_table.npz (0.25 m bins),
  PML ring (first/last 1 m = 4 bins) held constant at nearest interior value, mirroring
  the 3D edge-controlled continuation.
- Numerics inherited from HS4 header (dx 0.05, time_window 600e-9, impulse waveform,
  HORIPML, Debye cover dispersion); pml_cfs line dropped for 2D (defaults) -- documented
  deviation: kinematics unaffected, amplitudes not cross-comparable to 3D capsules.
"""
import argparse
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs/hs4_interface_binned_table.npz'

HEADER = (
    '#title: HS4T2D transect-2D fast study v0.1 (kinematic companion of HS4; invariant y; profile x = 3D bi=6 transect)\n'
    '#domain: 12 0.05 33\n'
    '#dx_dy_dz: 0.05 0.05 0.05\n'
    '#time_window: 600e-9\n'
    '#omp_threads: 8\n'
    '#pml_cells: 20 0 20 20 0 20\n'
    '#pml_formulation: HORIPML\n'
    '#waveform: impulse 1 1 impulse\n'
    '#hertzian_dipole: y {tx:.4g} 0.025 27 impulse\n'
    '#rx: {rx:.4g} 0.025 27 {tag} Ey\n'
    '#material: 9 0.001 1 0 rock\n'
    '#material: 18.017 0.003 1 0 cover\n'
    '#add_dispersion_debye: 1 7.878 6.4567e-09 cover\n'
)


def build_geometry():
    z0tab = np.load(TABLE)['z0_m']
    prof = z0tab[6, :].copy()           # 48 bins, 0.25 m, along y in 3D -> x here
    prof[:4] = prof[4]                  # PML ring continuation (1 m each side)
    prof[-4:] = prof[-5]
    lines = ['#box: 0 0 0 12 0.05 12 rock']
    for bj in range(48):
        x0, x1 = bj * 0.25, (bj + 1) * 0.25
        z0 = prof[bj]
        lines.append(f'#box: {x0:g} 0 {z0:g} {x1:g} 0.05 12 cover')
    lines.append('#geometry_view: 0 0 0 12 0.05 33 0.05 0.05 0.05 hs4t2d_geom n')
    return lines, prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--x0', type=float, default=2.60)
    ap.add_argument('--x1', type=float, default=8.60)
    ap.add_argument('--step', type=float, default=0.05)
    ap.add_argument('--offset', type=float, default=1.30)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    lines, prof = build_geometry()
    geom = '\n'.join(lines) + '\n'
    xs = np.round(np.arange(args.x0, args.x1 + 1e-9, args.step), 4)
    for k, x in enumerate(xs):
        tag = f't{k + 1:02d}'
        txt = HEADER.format(tx=x, rx=x + args.offset, tag=tag) + geom
        (out / f'hs4t2d_{tag}.in').write_text(txt, encoding='utf-8')
    np.savez(out / 'hs4t2d_transect_table.npz', prof_m=prof, bin_m=0.25,
             x0=args.x0, x1=args.x1, step=args.step, offset=args.offset,
             n_traces=len(xs))
    print(f'traces: {len(xs)}  x {args.x0}..{args.x1} step {args.step}  relief '
          f'{prof[10:35].max() - prof[10:35].min():.3f} m')


if __name__ == '__main__':
    main()
