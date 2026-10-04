"""Kirchhoff migration comparison across cover permittivity (scan v0.1).

CPU-only; no solver. Repeats the 15 m migration diagnostic of
migrate_hs4_height15m_diff_v0_1.py for each cover permittivity group of the
frozen scan capsule (eps_r 6 / 9 / 18.017-non-dispersive) and for the archived
dispersive eps_r=18.017 baseline, using identical grids, metrics and the
official SFCW chain. Per-case cover velocity v=c/sqrt(eps_r); the lateral
profile split scale c/(2f) is held fixed as a diagnostic choice. This is not
an independently measured resolution bound or a test isolating critical angle.
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
SCAN = ROOT/'artifacts/research_checks/2026-10-04_hs4_permittivity_scan'
PATCH_NPZ = ROOT/'artifacts/research_checks/2026-10-04_hs4_local_patch_results/patch_arrays.npz'
C_AIR = 0.299792458       # m/ns
SURFACE_Z = 12.0
ANTENNA_Z = 27.0
SRC_X0, RX_X0, STEP = 14.6, 15.9, 0.5
F_GHZ = 0.095
PROFILE_SPLIT_SCALE_M = C_AIR/(2*F_GHZ)  # diagnostic scale, not a measured bound
IMAGE_X = np.arange(13.5, 22.51, 0.05)
IMAGE_Z = np.arange(8.2, 12.01, 0.025)
METRIC_X = (14.75, 21.0)
CASES = [('eps6', 6.0, 'scan'), ('eps9', 9.0, 'scan'),
         ('eps18nd', 18.017, 'scan'), ('eps18disp', 18.017, 'baseline')]


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


def leg_times(px, pz, ant, v_cover):
    sx = np.linspace(min(px.min(), ant[0])-1.5, max(px.max(), ant[0])+1.5, 600)
    d_air = np.hypot(ant[0]-sx, ant[1]-SURFACE_Z)
    d_cov = np.hypot(px.ravel()[None, :]-sx[:, None], pz.ravel()[None, :]-SURFACE_Z)
    t = d_air[:, None]/C_AIR + d_cov/v_cover
    return t.min(axis=0).reshape(px.shape)


def migrate(traces, t_axis, px, pz, v_cover):
    img = np.zeros(px.shape)
    for i in range(traces.shape[1]):
        tau = leg_times(px, pz, (SRC_X0+STEP*i, ANTENNA_Z), v_cover) \
            + leg_times(px, pz, (RX_X0+STEP*i, ANTENNA_Z), v_cover)
        img += np.interp(tau, t_axis, traces[:, i], left=0.0, right=0.0)
    return img


def scan_paths(eps_tag):
    paths = {'rough': [], 'halfspace': []}
    for segment in ('left', 'centre', 'right'):
        for role in ('rough', 'halfspace'):
            folder = SCAN/f'{eps_tag}_{segment}_{role}'
            audit = json.loads((folder/'audit.json').read_text('utf-8'))
            from analyze_hs4_interface_echo_budget import sha256
            for row in sorted(audit, key=lambda r: r['trace']):
                p = folder/row['file']
                if sha256(p) != row['sha256']:
                    raise ValueError(f'identity mismatch: {p}')
                paths[role].append(p)
    return paths


def case_metrics(tag, eps, inputs, response, reconstruct_time_response):
    """inputs: role -> list of 13 h5 paths. Returns (summary, images, ridge, t_cal)."""
    v_cover = C_AIR/np.sqrt(eps)
    spectra = {role: [response(p) for p in inputs[role]] for role in ('rough', 'halfspace')}
    for role, rs in spectra.items():
        for r in rs:
            if not r.source_valid.all():
                raise ValueError(f'invalid source bins: {tag}/{role}')
    rough = np.stack([r.response for r in spectra['rough']], axis=1)
    halfspace = np.stack([r.response for r in spectra['halfspace']], axis=1)
    diff = rough - halfspace
    template = spectra['rough'][0]

    def product(spec):
        return reconstruct_time_response(replace(template, response=spec), window='hann', zero_pad_factor=8)
    prod = {k: product(v) for k, v in [('rough', rough), ('halfspace', halfspace), ('diff', diff)]}
    t_ns = prod['diff'].time*1e9
    traces = {k: prod[k].real_bandpass for k in prod}

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
    images = {k: migrate(traces[k], t_cal, PX, PZ, v_cover) for k in traces}

    ix0, ix1 = np.searchsorted(IMAGE_X, METRIC_X[0]), np.searchsorted(IMAGE_X, METRIC_X[1])
    sub = IMAGE_Z < SURFACE_Z-0.05
    zx = IMAGE_Z[sub]
    diff_sub = np.abs(images['diff'][np.ix_(sub, np.arange(ix0, ix1+1))])
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

    env = diff_sub
    sigma = 0.185
    def template_score(relief_z):
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
    scores = {name: template_score(r) for name, r in templates.items()}

    band = (IMAGE_Z >= 8.4) & (IMAGE_Z <= 9.7)
    e = {k: float((images[k][np.ix_(band, np.arange(ix0, ix1+1))]**2).sum()) for k in images}

    # Echo observables from the calibrated time traces (rough role).
    env_rough = np.abs(prod['rough'].complex_envelope)
    env_diff = np.abs(prod['diff'].complex_envelope)
    direct = []
    echo_t, echo_a = [], []
    tau_echo = 2*(ANTENNA_Z-SURFACE_Z)/C_AIR + 2*(SURFACE_Z-9.05)/v_cover
    for i in range(13):
        wd = (t_cal > 0) & (t_cal < 20)
        direct.append(float(env_rough[wd, i].max()))
        we = (t_cal > tau_echo-30) & (t_cal < tau_echo+30)
        j = int(np.argmax(env_diff[we, i]))
        echo_t.append(float(t_cal[we][j]))
        echo_a.append(float(env_diff[we, i].max()))
    summary = {
        'cover_eps_r': eps, 'v_cover_m_per_ns': float(v_cover),
        'sin_critical': float(1/np.sqrt(eps)),
        'critical_angle_deg': float(np.degrees(np.arcsin(1/np.sqrt(eps)))),
        'diagnostic_split_scale_m': float(PROFILE_SPLIT_SCALE_M),
        'time_offset_ns_median': time_offset_ns,
        'ridge_vs_true_relief': profile_metrics,
        'template_scores': scores,
        'interface_band_image_energy': e,
        'halfspace_over_diff_band_energy': e['halfspace']/e['diff'],
        'echo_peak_time_ns_median': float(np.median(echo_t)),
        'echo_over_direct_dB_median': float(20*np.log10(np.median(echo_a)/np.median(direct)))}
    return summary, images, ridge, rx_x, true, low, t_cal


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
    c = json.loads((SCAN/'execution_contract.json').read_text('utf-8'))
    v = json.loads((SCAN/'completed_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(SCAN/'execution_contract.json'):
        raise ValueError('completed verified scan capsule required')
    baseline = {'rough': [], 'halfspace': []}
    for role in ('rough', 'halfspace'):
        for idx in HIGH_INDICES:
            capsule, rel = high_path(role, idx)
            if sha256(capsule/rel) != manifest_sha(capsule, rel):
                raise ValueError(f'identity mismatch: {rel}')
            baseline[role].append(capsule/rel)

    results, images, ridges = {}, {}, {}
    for tag, eps, kind in CASES:
        inputs = scan_paths(tag) if kind == 'scan' else baseline
        summary, img, ridge, rx_x, true, low, t_cal = case_metrics(
            tag, eps, inputs, response, reconstruct_time_response)
        if tag == 'eps18disp':
            with np.load(PATCH_NPZ) as old:
                spectra_r = np.stack([response(p).response for p in baseline['rough']], axis=1)
                spectra_h = np.stack([response(p).response for p in baseline['halfspace']], axis=1)
                anchor = float(np.linalg.norm(old['high_base_spectrum']-(spectra_r-spectra_h))
                               / np.linalg.norm(old['high_base_spectrum']))
            if anchor > 1e-9:
                raise ValueError('chain does not reproduce archived baseline')
            summary['chain_anchor_relative_L2'] = anchor
        results[tag], images[tag], ridges[tag] = summary, img, (ridge, rx_x, true, low)
        print(tag, 'corr_low', round(summary['ridge_vs_true_relief']['lowpass_matched']['corr'], 3),
              'tpl_true', round(summary['template_scores']['true'], 3),
              'echo dB', round(summary['echo_over_direct_dB_median'], 1), flush=True)

    out = {'status': 'COMPLETED_MIGRATION_DIAGNOSTIC', 'cases': results,
           'code_sha256': sha256(__file__),
           'metrics_code_sha256': sha256(ROOT/'scripts/hs4_analysis_metrics.py'),
           'scan_contract_sha256': sha256(SCAN/'execution_contract.json'),
           'diagnostic_split_scale_m': float(PROFILE_SPLIT_SCALE_M),
           'notes': 'Identical Kirchhoff stack per case; v=c/sqrt(eps_r); matched profile splits; finite periodic interval; no physical resolution or recoverability certification.'}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(out, indent=2, ensure_ascii=False, allow_nan=False)+'\n',
                                         encoding='utf-8')
    saved = {'metric_x_m': rx_x, 'truth_z_m': true, 'image_x_m': IMAGE_X, 'image_z_m': IMAGE_Z}
    for tag in results:
        saved[tag+'_ridge_z_m'] = ridges[tag][0]
        for role, img in images[tag].items():
            saved[tag+'_'+role+'_image'] = img
    np.savez_compressed(args.out/'migration_arrays.npz', **saved)

    fx, fz = interface_relief()
    fig, axes = plt.subplots(4, 1, figsize=(12, 13), sharex=True)
    for ax, (tag, eps, kind) in zip(axes, CASES):
        im = images[tag]['diff']
        vm = np.abs(im).max() or 1.0
        ax.pcolormesh(IMAGE_X, IMAGE_Z, im, cmap='RdBu_r', vmin=-vm, vmax=vm, shading='auto')
        ax.plot(fx, fz, 'k-', lw=1)
        ax.axhline(SURFACE_Z, color='k', ls=':', lw=.6)
        ax.set_ylabel('z (m)')
        ax.set_title(f'{tag}: migrated rough-halfspace diff (eps_r={eps}{"" if kind=="scan" else ", dispersive"})')
        ax.invert_yaxis()
    axes[-1].set_xlabel('x (m)')
    fig.tight_layout()
    fig.savefig(args.out/'migrated_diff_by_eps.png', dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    for tag, eps, kind in CASES:
        ridge, rx_x, true, low = ridges[tag]
        ax.plot(rx_x, ridge, label=tag, lw=1.2)
    ax.plot(rx_x, true, 'k--', lw=1.4, label='true relief')
    ax.plot(rx_x, low, 'k:', lw=1, label=f'true, diagnostic lowpass > {PROFILE_SPLIT_SCALE_M:.2f} m')
    ax.invert_yaxis(); ax.grid(alpha=.3); ax.legend(ncol=3, fontsize=8)
    ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
    ax.set_title('diff-image ridge vs true relief, by cover permittivity')
    fig.tight_layout()
    fig.savefig(args.out/'ridge_vs_relief_by_eps.png', dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    xs = [r['cover_eps_r'] for _, r in results.items()]
    labels = list(results)
    for ax, (title, getter) in zip(axes, [
            ('matched lowpass corr (diagnostic split)', lambda r: r['ridge_vs_true_relief']['lowpass_matched']['corr']),
            ('template score (true)', lambda r: r['template_scores']['true']),
            ('echo/direct (dB)', lambda r: r['echo_over_direct_dB_median'])]):
        ax.plot(xs, [getter(r) for r in results.values()], 'o-')
        for x, y, lb in zip(xs, [getter(r) for r in results.values()], labels):
            ax.annotate(lb, (x, y), textcoords='offset points', xytext=(4, 4), fontsize=7)
        ax.set_xscale('log'); ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs])
        ax.set_xlabel('cover eps_r'); ax.set_title(title); ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(args.out/'metrics_vs_eps.png', dpi=110)
    plt.close(fig)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
