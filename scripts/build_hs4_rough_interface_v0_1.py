"""Build HS4 rough-interface half-space standard piece (.in) v0.1.

HS4 = HS1 flat-interface anchor with the bedrock-cover interface replaced by
the frozen natural-morphology field interface_field.npz (sha256 eabcad0b16db97d1…,
RMS 0.16 m / CL 2.5 m / seed 20261002, zero-mean Gaussian-correlated, quantized
to the 5 cm grid, clipped to [-2.00, +1.95] m so cover stays in (0, 5 m) and
never touches the surface).

Differences vs the legacy benchmark3d_r2_co build (whose defects HS1 fixed):
  1. Materials extend THROUGH the PML to the domain boundary: rock box spans
     0..12 m in x and y (legacy: 1..11 m, leaving a 1 m air ring).
  2. Edge-controlled continuation: the interface field's PML ring (outer 20
     cells = 1 m each side) is overwritten by edge replication of the inner
     region (np.pad edge mode), so the interface meets the PML with zero
     cross-boundary slope — no random cliffs inside the absorbing layer.
  3. Rock spans z 0..12 (legacy started at z=1, leaving air at the bottom PML).

Numerics are verbatim identical to frozen t07 / HS1-HS3 (domain 12x12x33,
dx 0.05, 600 ns, PML 20 HORIPML, impulse, Tx(6,5.35,27)/Rx(6,6.65,27), same
pml_cfs, same materials rock(9,0.001) / cover(18.017,0.003)+Debye).

Binning: same 0.25 m (5-cell) column bins as the frozen benchmark, extended to
the full domain: 48x48 = 2304 cover boxes; per-bin z0 = 9 + mean(eta_bin)
rounded to the 5 cm grid. The binned table is saved alongside the .in so that
any future 2D pair MUST be extracted from this final binned table (never from
the raw field) — per the HS design draft.

Usage (repo root):
  python scripts/build_hs4_rough_interface_v0_1.py --out <dir>
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIELD = ROOT / 'configs/research/benchmark3d_r2_co/interface_field.npz'
FIELD_SHA256_16 = 'eabcad0b16db97d1'

DX = 0.05
Z_IF0 = 9.0
Z_TOP = 12.0
DOMAIN = (12.0, 12.0, 33.0)
PML_CELLS = 20
BIN_CELLS = 5           # 0.25 m bins, same as frozen benchmark
ETA_CLIP = (-2.00, 1.95)

HEADER = """#title: HS4 rough-interface half-space standard piece v0.1 (natural-morphology interface, edge-controlled continuation, materials through PML); frozen numerics identical to B3D5CM-t07 / HS1-HS3
#domain_mode: 3D
#domain: 12 12 33
#dx_dy_dz: 0.05 0.05 0.05
#time_window: 600e-9
#omp_threads: 8
#pml_cells: 20 20 20 20 20 20
#pml_formulation: HORIPML
#waveform: impulse 1 1 impulse
#hertzian_dipole: x 6 5.35 27 impulse
#rx: 6 6.65 27 hs4 Ex
#pml_cfs: constant forward 0 0 constant forward 1 1 quartic forward 0 0.21235349838321013
#material: 9 0.001 1 0 rock
#material: 18.017 0.003 1 0 cover
#add_dispersion_debye: 1 7.878 6.4567e-09 cover
"""


def load_field():
    raw = FIELD.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    assert sha.startswith(FIELD_SHA256_16), f'interface field sha mismatch: {sha[:16]}'
    eta = np.load(FIELD)['eta']
    assert eta.shape == (240, 240)
    return eta, sha


def edge_controlled(eta):
    """Overwrite the PML ring (outer 20 cells) by edge replication of the inner region."""
    p = PML_CELLS
    inner = eta[p:-p, p:-p]
    return np.pad(inner, p, mode='edge')


def build(eta_ext):
    """Return (in_lines, binned z0 table [48x48])."""
    n = eta_ext.shape[0]
    nb = n // BIN_CELLS          # 48
    z0tab = np.zeros((nb, nb))
    lines = ['#box: 0 0 0 12 12 12 rock']
    for bi in range(nb):
        for bj in range(nb):
            eta_bin = eta_ext[bi * BIN_CELLS:(bi + 1) * BIN_CELLS,
                              bj * BIN_CELLS:(bj + 1) * BIN_CELLS].mean()
            z0 = Z_IF0 + float(np.clip(np.round(eta_bin / DX) * DX, *ETA_CLIP))
            z0tab[bi, bj] = z0
            x0, x1 = bi * BIN_CELLS * DX, (bi + 1) * BIN_CELLS * DX
            y0, y1 = bj * BIN_CELLS * DX, (bj + 1) * BIN_CELLS * DX
            lines.append(f'#box: {x0:g} {y0:g} {z0:g} {x1:g} {y1:g} {Z_TOP:g} cover')
    return lines, z0tab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)

    eta, sha = load_field()
    eta_ext = edge_controlled(eta)
    # sanity: PML ring is constant along the outward direction
    assert np.all(eta_ext[:PML_CELLS, :] == eta_ext[PML_CELLS, :][None, :])
    assert np.all(eta_ext[:, :PML_CELLS] == eta_ext[:, PML_CELLS][:, None])

    lines, z0tab = build(eta_ext)
    in_path = out / 'hs4_rough_halfspace.in'
    in_path.write_text(
        HEADER
        + lines[0] + '\n'
        + '\n'.join(lines[1:]) + '\n'
        + '#geometry_view: 0 0 0 12 12 33 0.05 0.05 0.05 hs4_geom n\n',
        encoding='utf-8')
    np.savez(out / 'hs4_interface_binned_table.npz',
             z0_m=z0tab, bin_m=BIN_CELLS * DX, z_if0=Z_IF0,
             field_sha256=sha, eta_extended=eta_ext)
    stats = {
        'field_sha256': sha,
        'bin_m': BIN_CELLS * DX,
        'n_cover_boxes': len(lines) - 1,
        'z0_min_m': float(z0tab.min()), 'z0_max_m': float(z0tab.max()),
        'z0_mean_m': float(z0tab.mean()),
        'pml_ring_constant': True,
    }
    (out / 'hs4_build_stats.json').write_text(json.dumps(stats, indent=1), encoding='utf-8')
    print(json.dumps(stats, indent=1))
    print('in_lines:', 1 + 1 + len(lines) + 1)


if __name__ == '__main__':
    main()
