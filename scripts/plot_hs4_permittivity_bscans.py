"""Per-case band-limited B-scans for the cover-permittivity scan (CPU-only).

Renders rough / halfspace / rough-halfspace B-scans for the four cases
(eps6, eps9, eps18nd from the scan capsule; eps18disp from the archived
dispersive 15 m baseline) through the identical official SFCW chain, with
per-case calibrated time axes and the true-relief two-way-time curve overlaid
(per-case cover velocity). Display norm is one scalar per row shared by all
four cases (no per-trace scaling); this is a view of verified evidence, not a
processing step.
"""
import argparse
import json
from dataclasses import replace
from pathlib import Path
import sys

import h5py  # noqa: F401
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
C_AIR = 0.299792458
SURFACE_Z = 12.0
ANTENNA_Z = 27.0
SRC_X0, RX_X0, STEP = 14.6, 15.9, 0.5
CASES = [('eps6', 6.0, 'scan'), ('eps9', 9.0, 'scan'),
         ('eps18nd', 18.017, 'scan'), ('eps18disp', 18.017, 'baseline')]
T_WIN = (80.0, 240.0)


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('new result directory required')
    from analyze_hs4_interface_echo_budget import shim_gprmax, high_path, manifest_sha, sha256, HIGH_INDICES
    shim_gprmax(Path(r'E:\gprMax-v.4.0.0\gprMax-v.4.0.0'))
    from sfcw_official_loader_v0_2 import verify_official_runtime
    from analyze_hs4_height_wavefield import response
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    from migrate_hs4_permittivity_scan_v0_1 import scan_paths, SCAN
    verify_official_runtime()
    v = json.loads((SCAN/'completed_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS' or not v['completed']:
        raise ValueError('completed verified scan capsule required')
    baseline = {'rough': [], 'halfspace': []}
    for role in ('rough', 'halfspace'):
        for idx in HIGH_INDICES:
            capsule, rel = high_path(role, idx)
            if sha256(capsule/rel) != manifest_sha(capsule, rel):
                raise ValueError(f'identity mismatch: {rel}')
            baseline[role].append(capsule/rel)
    fx, fz = interface_relief()

    panels, meta = {}, {}
    for tag, eps, kind in CASES:
        inputs = scan_paths(tag) if kind == 'scan' else baseline
        v_cover = C_AIR/np.sqrt(eps)
        spectra = {role: [response(p) for p in inputs[role]] for role in ('rough', 'halfspace')}
        rough = np.stack([r.response for r in spectra['rough']], axis=1)
        halfspace = np.stack([r.response for r in spectra['halfspace']], axis=1)
        template = spectra['rough'][0]

        def product(spec):
            return reconstruct_time_response(replace(template, response=spec),
                                             window='hann', zero_pad_factor=8)
        prod = {k: product(vv) for k, vv in [('rough', rough), ('halfspace', halfspace),
                                             ('diff', rough-halfspace)]}
        t_ns = prod['diff'].time*1e9
        dt_obs = []
        for i in range(13):
            tx, rx = SRC_X0+STEP*i, RX_X0+STEP*i
            xmid = 0.5*(tx+rx)
            tau_pred = (np.hypot(tx-xmid, ANTENNA_Z-SURFACE_Z)
                        + np.hypot(rx-xmid, ANTENNA_Z-SURFACE_Z))/C_AIR
            env = np.abs(prod['rough'].complex_envelope[:, i])
            win = (t_ns > 85) & (t_ns < 135)
            dt_obs.append(float(t_ns[win][np.argmax(env[win])] - tau_pred))
        t_cal = t_ns - float(np.median(dt_obs))
        panels[tag] = {k: prod[k].real_bandpass for k in prod}
        # True-relief two-way time per station midpoint.
        relief_t = []
        for i in range(13):
            tx, rx = SRC_X0+STEP*i, RX_X0+STEP*i
            xmid = 0.5*(tx+rx)
            z = float(np.interp(xmid, fx, fz))
            relief_t.append((np.hypot(tx-xmid, ANTENNA_Z-SURFACE_Z)
                             + np.hypot(rx-xmid, ANTENNA_Z-SURFACE_Z))/C_AIR
                            + 2*(SURFACE_Z-z)/v_cover)
        meta[tag] = {'cover_eps_r': eps, 'v_cover_m_per_ns': float(v_cover),
                     'time_offset_ns_median': float(np.median(dt_obs)),
                     'relief_two_way_ns': relief_t}
    norms = {row: float(max(np.abs(panels[tag][row]).max() for tag, _, _ in CASES))
             for row in ('rough', 'halfspace', 'diff')}
    args.out.mkdir(parents=True)
    fig, axes = plt.subplots(3, 4, figsize=(17, 9), sharex=True, sharey=True)
    xmids = [0.5*(SRC_X0+STEP*i+RX_X0+STEP*i) for i in range(13)]
    # Calibrated time axis per case (identical chain => identical grid).
    tcal = {}
    for tag, eps, kind in CASES:
        inputs = scan_paths(tag) if kind == 'scan' else baseline
        pr = reconstruct_time_response(response(inputs['rough'][0]), window='hann', zero_pad_factor=8)
        tcal[tag] = pr.time*1e9 - meta[tag]['time_offset_ns_median']
    for col, (tag, eps, kind) in enumerate(CASES):
        t = tcal[tag]
        w = (t >= T_WIN[0]) & (t <= T_WIN[1])
        for row, name in enumerate(('rough', 'halfspace', 'diff')):
            ax = axes[row, col]
            data = panels[tag][name][w, :]
            n = norms[name]
            ax.pcolormesh(xmids, t[w], data, cmap='RdBu_r', vmin=-n, vmax=n, shading='auto')
            ax.plot(xmids, meta[tag]['relief_two_way_ns'], 'k-', lw=1.0)
            ax.set_ylim(T_WIN[1], T_WIN[0])
            if row == 0:
                ax.set_title(f'{tag} (eps_r={eps}{"" if kind=="scan" else ", disp"})', fontsize=10)
            if col == 0:
                ax.set_ylabel(f'{name}\nt (ns)')
            if row == 2:
                ax.set_xlabel('midpoint x (m)')
    fig.suptitle('Band-limited B-scans (official SFCW chain), row-shared fixed norms: '
                 f'rough ±{norms["rough"]:.3g}, halfspace ±{norms["halfspace"]:.3g}, diff ±{norms["diff"]:.3g}; '
                 'black = true-relief two-way time')
    fig.tight_layout()
    fig.savefig(args.out/'bscans_by_eps.png', dpi=110)
    plt.close(fig)
    # Second view: per-panel self-normalized dB envelope (display only; each
    # panel annotated with its own peak so weak cases stay visible).
    fig, axes = plt.subplots(3, 4, figsize=(17, 9), sharex=True, sharey=True)
    for col, (tag, eps, kind) in enumerate(CASES):
        t = tcal[tag]
        w = (t >= T_WIN[0]) & (t <= T_WIN[1])
        for row, name in enumerate(('rough', 'halfspace', 'diff')):
            ax = axes[row, col]
            env = np.abs(np.abs(panels[tag][name][w, :]))
            peak = float(env.max()) or 1.0
            ax.pcolormesh(xmids, t[w], 20*np.log10(np.maximum(env, peak*1e-4)/peak),
                          cmap='viridis', vmin=-80, vmax=0, shading='auto')
            ax.plot(xmids, meta[tag]['relief_two_way_ns'], 'r-', lw=1.0)
            ax.set_ylim(T_WIN[1], T_WIN[0])
            if row == 0:
                ax.set_title(f'{tag} (eps_r={eps}{"" if kind=="scan" else ", disp"})', fontsize=10)
            if col == 0:
                ax.set_ylabel(f'{name}\nt (ns)')
            if row == 2:
                ax.set_xlabel('midpoint x (m)')
            ax.text(0.02, 0.03, f'peak {peak:.3g}', transform=ax.transAxes, fontsize=7,
                    color='w', va='bottom')
    fig.suptitle('Per-panel self-normalized envelope (dB re panel peak, display only; '
                 'peaks annotated); red = true-relief two-way time')
    fig.tight_layout()
    fig.savefig(args.out/'bscans_by_eps_selfnorm_dB.png', dpi=110)
    plt.close(fig)
    summary = {'status': 'PASS', 'cases': meta, 'row_norms': norms,
               'time_window_ns': list(T_WIN),
               'notes': 'Row-shared fixed display norms; per-case ground-bounce time calibration; relief overlay uses per-case v=c/sqrt(eps_r).'}
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n',
                                         encoding='utf-8')
    print(json.dumps({'row_norms': norms}, indent=1))


if __name__ == '__main__':
    main()
