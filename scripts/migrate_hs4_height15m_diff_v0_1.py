"""Kirchhoff migration diagnostic on the verified 13-station 15 m HS4T2D data.

CPU-only; no solver. Measures ridge errors using the same spatial split on
the ridge and truth; the selected split scale is not a resolution bound.
Migrates the official-chain band-limited
traces (rough, halfspace, rough-halfspace diff) with a two-layer (air/cover)
Fermat travel time through the flat surface. A single constant time offset is
calibrated from the observed ground-bounce peak versus its predicted air-only
two-way time, then applied to all images. The diff-image ridge is compared
with the true relief (full / >1.6 m low-pass / <1.6 m high-pass components).
"""
import argparse
import json
from dataclasses import replace
from pathlib import Path
import sys

import h5py  # noqa: F401  (used by response chain)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hs4_analysis_metrics import matched_profile_metrics, split_profile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
PATCH_NPZ = ROOT/'artifacts/research_checks/2026-10-04_hs4_local_patch_results/patch_arrays.npz'
C_AIR = 0.299792458       # m/ns
V_COVER = C_AIR/np.sqrt(18.017)
SURFACE_Z = 12.0
ANTENNA_Z = 27.0
SRC_X0, RX_X0, STEP = 14.6, 15.9, 0.5
PROFILE_SPLIT_SCALE_M = C_AIR/(2*0.095)  # legacy diagnostic scale, not a measured bound
IMAGE_X = np.arange(13.5, 22.51, 0.05)
IMAGE_Z = np.arange(8.2, 12.01, 0.025)
METRIC_X = (14.75, 21.0)  # inside station midpoint coverage, away from taper edges


def interface_relief():
    xs, zs = [], []
    ref = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre/centre_rough/profile.in'
    for line in ref.read_text('utf-8').splitlines():
        if not line.startswith('#box:'):
            continue
        p = [float(v) for v in line.split()[1:7]]
        if p[5] == 12.0 and p[1] == 0.0 and p[2] > 5.0:
            xs.append(0.5*(p[0]+p[3])); zs.append(p[2])
    order = np.argsort(xs)
    return np.array(xs)[order], np.array(zs)[order]


def leg_times(px, pz, ant):
    """Minimum antenna->P time refracting through the flat surface (dense scan)."""
    sx = np.linspace(min(px.min(), ant[0])-1.5, max(px.max(), ant[0])+1.5, 600)
    d_air = np.hypot(ant[0]-sx, ant[1]-SURFACE_Z)
    d_cov = np.hypot(px.ravel()[None, :]-sx[:, None], pz.ravel()[None, :]-SURFACE_Z)
    t = d_air[:, None]/C_AIR + d_cov/V_COVER
    return t.min(axis=0).reshape(px.shape)


def migrate(traces, t_axis, px, pz):
    img = np.zeros(px.shape)
    for i in range(traces.shape[1]):
        tau = leg_times(px, pz, (SRC_X0+STEP*i, ANTENNA_Z)) + leg_times(px, pz, (RX_X0+STEP*i, ANTENNA_Z))
        img += np.interp(tau, t_axis, traces[:, i], left=0.0, right=0.0)
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('new result directory required')
    from analyze_hs4_interface_echo_budget import (shim_gprmax, high_path, manifest_sha,
                                                   sha256, HIGH_INDICES)
    shim_gprmax(Path(r'E:\gprMax-v.4.0.0\gprMax-v.4.0.0'))
    from sfcw_official_loader_v0_2 import verify_official_runtime
    from analyze_hs4_height_wavefield import response
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    verify_official_runtime()
    inputs = {}
    for role in ('rough', 'halfspace'):
        for idx in HIGH_INDICES:
            capsule, rel = high_path(role, idx)
            if sha256(capsule/rel) != manifest_sha(capsule, rel):
                raise ValueError(f'identity mismatch: {rel}')
            inputs[f'{role}_{idx}'] = capsule/rel
    spectra = {k: response(p) for k, p in inputs.items()}
    for k, r in spectra.items():
        if not r.source_valid.all():
            raise ValueError(f'invalid source bins: {k}')
    rough = np.stack([spectra[f'rough_{i}'].response for i in HIGH_INDICES], axis=1)
    halfspace = np.stack([spectra[f'halfspace_{i}'].response for i in HIGH_INDICES], axis=1)
    diff = rough - halfspace
    template = spectra['rough_61']

    def product(spec):
        return reconstruct_time_response(replace(template, response=spec),
                                         window='hann', zero_pad_factor=8)
    prod = {k: product(v) for k, v in [('rough', rough), ('halfspace', halfspace), ('diff', diff)]}
    with np.load(PATCH_NPZ) as old:
        anchor = float(np.linalg.norm(old['high_base_spectrum'] - diff)
                       / np.linalg.norm(old['high_base_spectrum']))
    if anchor > 1e-9:
        raise ValueError('chain does not reproduce archived baseline')
    t_ns = prod['diff'].time*1e9
    traces = {k: prod[k].real_bandpass for k in prod}

    # Time-zero calibration: observed ground-bounce envelope peak vs predicted
    # air-only two-way time to the flat surface point under the midpoint.
    dt_obs = []
    for i in range(13):
        tx, rx = SRC_X0+STEP*i, RX_X0+STEP*i
        xmid = 0.5*(tx+rx)
        tau_pred = (np.hypot(tx-xmid, ANTENNA_Z-SURFACE_Z)
                    + np.hypot(rx-xmid, ANTENNA_Z-SURFACE_Z))/C_AIR
        env = np.abs(prod['rough'].complex_envelope[:, i])
        win = (t_ns > 85) & (t_ns < 135)
        dt_obs.append(float(t_ns[win][np.argmax(env[win])] - tau_pred))
    time_offset_ns = float(np.median(dt_obs))
    t_cal = t_ns - time_offset_ns

    PX, PZ = np.meshgrid(IMAGE_X, IMAGE_Z)
    images = {k: migrate(traces[k], t_cal, PX, PZ) for k in traces}

    # Ridge of |diff image| per x column, restricted to below-surface depths.
    ix0, ix1 = np.searchsorted(IMAGE_X, METRIC_X[0]), np.searchsorted(IMAGE_X, METRIC_X[1])
    sub = IMAGE_Z < SURFACE_Z-0.05
    zx = IMAGE_Z[sub]
    diff_sub = np.abs(images['diff'][np.ix_(sub, np.arange(ix0, ix1+1))])
    # Robust ridge: per-column envelope centroid above half the column max.
    ridge = np.empty(diff_sub.shape[1])
    for j in range(diff_sub.shape[1]):
        col = diff_sub[:, j]
        m = col.max()
        w = np.where(col >= 0.5*m, col, 0.0)
        ridge[j] = float((zx*w).sum()/w.sum()) if w.sum() > 0 else float(zx[np.argmax(col)])
    rx_x = IMAGE_X[ix0:ix1+1]
    true = np.interp(rx_x, *interface_relief())
    spacing = float(np.mean(np.diff(rx_x)))
    low_centered, high = split_profile(true, spacing, PROFILE_SPLIT_SCALE_M)
    low = low_centered + true.mean()
    profile_metrics = matched_profile_metrics(ridge, true, spacing, PROFILE_SPLIT_SCALE_M)

    # Template-match score: 2D normalised correlation of the diff-image
    # envelope against a Gaussian corridor along a candidate relief.
    env = diff_sub
    sigma = 0.185  # lambda_cover/4 corridor width
    def template(relief_z):
        T = np.exp(-(zx[:, None]-relief_z[None, :])**2/(2*sigma**2))
        return float((env*T).sum()/np.sqrt((env**2).sum()*(T**2).sum()))
    templates = {
        'true': true,
        'true_lowpass_gt_resolution_bound': low,
        'flat_mean_depth': np.full_like(true, true.mean()),
        'true_shifted_plus1m': np.interp(rx_x-1.0, rx_x, true, left=true[0], right=true[-1]),
        'true_shifted_plus2m': np.interp(rx_x-2.0, rx_x, true, left=true[0], right=true[-1]),
        'flat_plus_highpass_only': true.mean()+high,
    }
    scores = {name: template(r) for name, r in templates.items()}

    # Control: halfspace image should have no interface ridge energy.
    band = (IMAGE_Z >= 8.4) & (IMAGE_Z <= 9.7)
    e = {k: float((images[k][np.ix_(band, np.arange(ix0, ix1+1))]**2).sum()) for k in images}
    summary = {
        'status': 'COMPLETED_MIGRATION_DIAGNOSTIC',
        'code_sha256': sha256(__file__),
        'metrics_code_sha256': sha256(ROOT/'scripts/hs4_analysis_metrics.py'),
        'time_offset_ns_median': time_offset_ns,
        'time_offset_per_station_ns': dt_obs,
        'chain_anchor_relative_L2': anchor,
        'diagnostic_split_scale_m': float(PROFILE_SPLIT_SCALE_M),
        'metric_x_range_m': list(METRIC_X),
        'ridge_vs_true_relief': profile_metrics,
        'template_scores': scores,
        'interface_band_image_energy': e,
        'diff_over_rough_band_energy': e['diff']/e['rough'],
        'halfspace_over_diff_band_energy': e['halfspace']/e['diff'],
        'notes': 'Kirchhoff stack, 13 common-offset traces, two-layer Fermat times; no amplitude weighting, no clean-claim.'}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n',
                                         encoding='utf-8')
    np.savez_compressed(args.out/'migration_arrays.npz', metric_x_m=rx_x,
                        truth_z_m=true, ridge_z_m=ridge, image_x_m=IMAGE_X,
                        image_z_m=IMAGE_Z, **images)
    # Figures.
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fx, fz = interface_relief()
    for ax, name in zip(axes, ('rough', 'halfspace', 'diff')):
        v = np.abs(images[name]).max() or 1.0
        ax.pcolormesh(IMAGE_X, IMAGE_Z, images[name], cmap='RdBu_r', vmin=-v, vmax=v, shading='auto')
        ax.plot(fx, fz, 'k-', lw=1)
        ax.axhline(SURFACE_Z, color='k', ls=':', lw=.6)
        ax.set_ylabel('z (m)'); ax.set_title(f'migrated {name}')
        ax.invert_yaxis()
    axes[-1].set_xlabel('x (m)')
    fig.suptitle(f'Kirchhoff migration, 13 stations, 15 m height; time offset {time_offset_ns:.2f} ns')
    fig.tight_layout()
    fig.savefig(args.out/'migrated_images.png', dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(rx_x, ridge, label='diff-image ridge')
    ax.plot(rx_x, true, 'k--', label='true relief')
    ax.plot(rx_x, low, 'g:', label=f'true, diagnostic lowpass > {PROFILE_SPLIT_SCALE_M:.2f} m')
    ax.invert_yaxis(); ax.grid(alpha=.3); ax.legend()
    ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
    ax.set_title(f'ridge vs relief: corr full {summary["ridge_vs_true_relief"]["corr_full"]:.3f}, '
                 f'matched low-pass {profile_metrics["lowpass_matched"]["corr"]:.3f}, '
                 f'matched high-pass {profile_metrics["highpass_matched"]["corr"]:.3f}')
    fig.tight_layout()
    fig.savefig(args.out/'ridge_vs_relief.png', dpi=110)
    plt.close(fig)
    print(json.dumps(summary, indent=1))


if __name__=='__main__':
    main()
