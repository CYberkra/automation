# -*- coding: utf-8 -*-
"""Generate the benchmark-v3 2D scan batch (BF series) gprMax .in files.

Official-chain constraint (gprMax.toolboxes.SFCW README + processing.py):
the direct method divides the receiver spectrum by the STORED source spectrum
(H(f) = Y(f)/X(f)). Therefore any solver excitation is deconvolved away and
the synthesized SFCW output is waveform-independent. All .in files here use
the documented impulse excitation; the measured antenna wavelet is applied
POST-chain (device-domain convolution, same stage as the v6r ringing),
recorded per run in manifest.json.

Phases:
  phase1: height {5,7,10,15} m, sigma=0, 3 traces (t05-t07), each processed
          with wavelet injection off/on -> 12 .in, 24 outputs
  phase2: h=7 m, wavelet on, clutter sigma {0,2,3}, 13 traces (t01-t13,
          5 cm step) -> 39 .in, 39 outputs

Geometry follows benchmark3d_r2_co 2D pair conventions (dx 0.025, 600 ns,
HORIPML, same pml_cfs, same rough-interface staircase y 1..11 m in 0.25 m
bins, rock eps_r=9, cover eps_r=18.017 + Debye), with:
  - flight height h (antenna z = 12 + h, air headroom +6 m)
  - Tx y=5.35 / Rx y=6.65 (1.3 m offset), trace k shifts both by (k-7)*0.05
  - per-bin cover er perturbation: er = 18.017 + sigma * g(y), g = smoothed
    unit Gaussian (correlation length ~0.5 m), seed frozen per config

Output: verify_runs/bf2d_scan_20261002/{phase1,phase2}/BF2D-*.in +
manifest.json. Importable: build_geometry(h, sigma, seed) for the plotter.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(r'E:\automation_djh')
OUT = ROOT / 'verify_runs' / 'bf2d_scan_20261002'

# ---- frozen r2 interface staircase (z of each 0.25 m y-bin, y=1..11) ----
ZI_R2 = [8.8, 8.8, 8.75, 8.75, 8.7, 8.7, 8.7, 8.7, 8.75, 8.75,
         8.75, 8.75, 8.8, 8.8, 8.85, 8.85, 8.9, 8.9, 8.95, 8.95,
         8.95, 8.95, 9.0, 9.0, 9.0, 9.0, 9.05, 9.05, 9.05, 9.05,
         9.0, 9.0, 9.0, 9.0, 8.95, 8.95, 8.95, 8.95, 8.9, 8.9]
Y0, YBIN = 1.0, 0.25
ER_COVER, SIGMA_COVER = 18.017, 0.003
DEBYE_DELTA, DEBYE_TAU = 7.878, 6.4567e-09
Z_SURF = 12.0
TX_Y0, RX_Y0, DY = 5.35, 6.65, 0.05
PML_CFS = ('#pml_cfs: constant forward 0 0 constant forward 1 1 '
           'quartic forward 0 0.21235349838321013')


def cfg_seed(s):
    return int(hashlib.sha256(s.encode()).hexdigest()[:8], 16)


def clutter_g(n, seed, corr_bins=2):
    """Smoothed unit-variance Gaussian along bins (~0.5 m correlation)."""
    rng = np.random.default_rng(seed)
    g = rng.standard_normal(n)
    k = np.ones(corr_bins) / corr_bins
    g = np.convolve(g, k, mode='same')
    g -= g.mean()
    g /= g.std() + 1e-30
    return g


def build_geometry(h, sigma, seed):
    """Return list of cover bins (y0, y1, z0, z_surf, er)."""
    g = clutter_g(len(ZI_R2), seed) if sigma > 0 else np.zeros(len(ZI_R2))
    bins = []
    for i, zi in enumerate(ZI_R2):
        er = ER_COVER + sigma * g[i]
        bins.append((Y0 + i * YBIN, Y0 + (i + 1) * YBIN, zi, Z_SURF, er))
    return bins


def material_block(bins):
    """Deduplicate materials by er; return (materials, bin->matname)."""
    mats, assign = [], []
    for (_, _, _, _, er) in bins:
        key = round(er, 3)
        name = next((m[0] for m in mats if abs(m[1] - key) < 1e-9), None)
        if name is None:
            name = f'cov{len(mats):02d}'
            mats.append((name, key))
        assign.append(name)
    return mats, assign


def render_in(h, sigma, seed, trace):
    """One .in file text (impulse excitation, official-chain compliant)."""
    ztop = Z_SURF + h + 6.0
    ztx = Z_SURF + h
    dy_off = (trace - 7) * DY
    tx_y, rx_y = TX_Y0 + dy_off, RX_Y0 + dy_off
    bins = build_geometry(h, sigma, seed)
    mats, assign = material_block(bins)

    lines = [
        f'#title: BF2D-h{h}-s{sigma:g}-t{trace:02d} benchmark v3 2D scan'
        f' (seed {seed}, field bf2d_scan_20261002)',
        '#domain_mode: TM',
        f'#domain: inf 12 {ztop:g}',
        '#dx_dy_dz: 0.025 0.025 0.025',
        '#time_window: 1200e-9',
        '#omp_threads: 8',
        '#pml_cells: 0 40 40 0 40 40',
        '#pml_formulation: HORIPML',
        PML_CFS,
        '#waveform: impulse 1 1 impulse',
        '#material: 9 0.001 1 0 rock',
    ]
    for name, er in mats:
        lines.append(f'#material: {er:.3f} {SIGMA_COVER} 1 0 {name}')
        delta = DEBYE_DELTA * er / ER_COVER
        lines.append(f'#add_dispersion_debye: 1 {delta:.4f} {DEBYE_TAU} {name}')
    lines.append('#box: inf 1 1 inf 11 12 rock')
    for (name, (y0, y1, z0, z1, _)) in zip(assign, bins):
        lines.append(f'#box: inf {y0:g} {z0:g} inf {y1:g} {z1:g} {name}')
    lines.append(f'#hertzian_dipole: x inf {tx_y:g} {ztx:g} impulse')
    lines.append(f'#rx: inf {rx_y:g} {ztx:g} t{trace:02d} Ex')
    return '\n'.join(lines) + '\n'


def main():
    manifest = {'field': 'bf2d_scan_20261002', 'wavelet':
                str(ROOT / 'artifacts_check' / 'ringing_field'
                    / 'line9_wavelet.npz'), 'runs': []}
    n_in = 0
    plans = [
        ('phase1', [(h, 0.0, (5, 6, 7)) for h in (5, 7, 10, 15)],
         [False, True]),                       # wavelet off/on
        ('phase2', [(7, s, tuple(range(1, 14))) for s in (0.0, 2.0, 3.0)],
         [True]),
    ]
    for phase, cfgs, wavelet_flags in plans:
        pdir = OUT / phase
        pdir.mkdir(parents=True, exist_ok=True)
        for h, sigma, traces in cfgs:
            seed = cfg_seed(f'bf2d-{phase}-h{h}-s{sigma}')
            for tr in traces:
                name = f'BF2D-h{h}-s{sigma:g}-t{tr:02d}'
                p = pdir / f'{name}.in'
                if not p.exists():
                    p.write_text(render_in(h, sigma, seed, tr),
                                 encoding='utf-8')
                    n_in += 1
                for w in wavelet_flags:
                    manifest['runs'].append(
                        {'run': f'{name}-w{int(w)}', 'in': str(p),
                         'h': h, 'sigma': sigma, 'trace': tr,
                         'wavelet': w, 'seed': seed})
    # drop stale excitation-based .in from the superseded draft
    for stale in OUT.rglob('*-meas-*.in'):
        stale.unlink()
        print('removed stale', stale.name)
    wav = OUT / 'bf2d_wavelet.dat'
    if wav.exists():
        wav.unlink()
        print('removed stale', wav.name)
    (OUT / 'manifest.json').write_text(
        json.dumps(manifest, indent=1, ensure_ascii=False), encoding='utf-8')
    print(f'{n_in} new .in files, {len(manifest["runs"])} runs '
          f'-> {OUT / "manifest.json"}')


if __name__ == '__main__':
    main()
