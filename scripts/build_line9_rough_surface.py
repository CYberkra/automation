"""Seeded fractal micro-relief on the frozen Line9 terrain voxels; never solve.

Deviation note (approved scope S5, 2026-10-11): official #fractal_box: +
#add_surface_roughness: only applies to axis-aligned flat box faces and cannot
follow the real undulating Line9 terrain. We therefore apply the same fractal
idea directly to the imported voxel columns (official #geometry_objects_read:
path stays unchanged), with fixed seed and bounded amplitude.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

from prepare_line9_agl5 import surface_y

DL = 0.025
AIR, COVER = 0, 1
MAX_SHIFT_CELLS = 8  # ±0.2 m
SEED = 20261011


def fractal_profile(n, hurst=0.5, seed=SEED):
    """1D fractional-Brownian height profile via spectral synthesis, zero mean."""
    rng = np.random.default_rng(seed)
    freqs = np.fft.rfftfreq(n)
    amp = np.zeros_like(freqs, dtype=complex)
    amp[1:] = freqs[1:] ** -(hurst + 0.5)
    amp *= np.exp(1j * rng.uniform(0, 2 * np.pi, len(freqs)))
    h = np.fft.irfft(amp, n)
    h -= h.mean()
    return h / h.std()


def build(src, dst, rms_m=0.1):
    with h5py.File(src) as h:
        data = h['data'][:]
        others = {key: h[key][:] for key in h.keys() if key != 'data'}
        keys = list(h.keys())
        attrs = dict(h.attrs)
    nx = data.shape[0]
    shifts = np.round(fractal_profile(nx) * rms_m / DL).astype(int)
    shifts = np.clip(shifts, -MAX_SHIFT_CELLS, MAX_SHIFT_CELLS)
    out = data.copy()
    changed = 0
    for i in range(nx):
        s = shifts[i]
        if s == 0:
            continue
        col = data[i, :, 0]
        surf_cells = np.flatnonzero(col != AIR)
        top = surf_cells[-1]  # highest non-air voxel
        if s > 0:  # raise surface: air -> cover
            out[i, top + 1: top + 1 + s, 0] = COVER
        else:  # lower surface: cover -> air
            if col[top + s + 1: top + 1].max() > COVER:
                raise ValueError(f'column {i}: lowering would cut below cover')
            out[i, top + s + 1: top + 1, 0] = AIR
        changed += abs(s)
    with h5py.File(dst, 'w') as h:
        for key in keys:
            if key == 'data':
                h.create_dataset('data', data=out, compression='gzip')
            else:
                h.create_dataset(key, data=others[key])
        for k, v in attrs.items():
            h.attrs[k] = v
    return dict(shifts_min_cells=int(shifts.min()), shifts_max_cells=int(shifts.max()),
                shifts_rms_m=float(np.abs(shifts).std() * DL), changed_voxels=int(changed),
                seed=SEED, rms_target_m=rms_m, max_shift_cells=MAX_SHIFT_CELLS)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometries', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    stats = {}
    for role in ['H0', 'H1']:
        stats[role] = build(a.geometries / f'{role}.h5', a.out / f'{role}_rough.h5')
        with h5py.File(a.out / f'{role}_rough.h5') as h:
            d = h['data'][:, :, 0]
            for x in [169.35, 170.65]:
                s = surface_y(d, x)
                assert s + 5 + 0.2 < 40.5, 'antenna must stay clear of top PML'
    (a.out / 'rough_stats.json').write_text(json.dumps(stats, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(stats))
