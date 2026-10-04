"""Aperture-scan analysis: 30 m aperture B-scans and sub-aperture shape fidelity.

Reads the completed 2026-10-05_hs4_aperture_scan_r1 capsule (2 materials x 2 roles
x 4 segments = 244 traces, 108 m domain, 8 m altitude, stations src 38.6..68.6 m
step 0.5 m) and:
- assembles full 61-station total/contrast B-scans per material;
- per-station interface picks (contrast envelope peak, +/-25 ns gate on the true
  relief two-way time) and correlation/RMSE against the true relief;
- sub-aperture analysis: central 13 stations (6 m, the archived layout) vs central
  25 (12 m) vs all 61 (30 m) shape fidelity, plus flat-segment control (stations
  over the flat interface should show a flat pick).
Outputs summary.json + figures into a fresh result directory.
Diagnostic; not physical acceptance.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from hs_capsule_identity import sha256
from analyze_hs4_height_wavefield import response
from check_hs4_v4_factor_evidence import direct_response
from plot_hs4_permittivity_bscans import C_AIR, SURFACE_Z, interface_relief
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

CHECKS = ROOT / 'artifacts/research_checks'
CAPSULE = CHECKS / '2026-10-05_hs4_aperture_scan_r1'
TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
SRC_X0, STEP, RX_OFFSET, ANT_Z, SHIFT = 38.6, 0.5, 1.3, 20.0, 36.0
N_STATIONS = 61
RELIEF_X_108 = (12.0 + SHIFT, 24.0 + SHIFT)  # relief segment in 108 m coords


def relief_twoway_108(xmids, eps):
    fx, fz = interface_relief()          # 36 m coords, x in [12, 24]
    fx108 = fx + SHIFT
    v = C_AIR/np.sqrt(eps)
    z = np.interp(xmids, fx108, fz)      # clamps to flat endpoints outside relief
    return 2*np.hypot(0.65, ANT_Z-SURFACE_Z)/C_AIR + 2*(SURFACE_Z - z)/v


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    contract = json.loads((CAPSULE/'execution_contract.json').read_text('utf-8'))
    ver = json.loads((CAPSULE/'completed_verification.json').read_text('utf-8'))
    if ver['status'] != 'PASS' or ver['contract_sha256'] != sha256(CAPSULE/'execution_contract.json'):
        raise ValueError('completed verification PASS required')

    spectra, dft_err = {}, {}
    for g in contract['groups']:
        audit = json.loads((Path(g['input']).parent/'audit.json').read_text('utf-8'))
        for row in audit:
            raw = Path(g['input']).parent/row['file']
            if sha256(raw) != row['sha256']:
                raise ValueError('raw identity differs: '+raw.name)
            r = response(raw)
            spectra[g['material_tag'], g['role'], row['station_index']] = r
            exact = direct_response(raw, r.frequency[TAKE])
            err = float(np.linalg.norm(exact - r.response[TAKE])/np.linalg.norm(exact))
            if err > 1e-9:
                raise ValueError('independent DFT failed: '+raw.name)
            dft_err[f"{g['id']}/{row['file']}"] = err

    template = next(iter(spectra.values()))
    tags = sorted({g['material_tag'] for g in contract['groups']})
    eps_of = {g['material_tag']: g['cover_eps_r'] for g in contract['groups']}
    station_x = SRC_X0 + STEP*np.arange(N_STATIONS)
    xmids = station_x + RX_OFFSET/2

    products, per_material = {}, {}
    time = None
    for tag in tags:
        eps = eps_of[tag]
        relief = relief_twoway_108(xmids, eps)
        for si in range(1, N_STATIONS+1):
            rough = spectra[tag, 'rough', si]
            contrast_spec = rough.response - spectra[tag, 'halfspace', si].response
            products[tag, si, 'total'] = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            products[tag, si, 'contrast'] = reconstruct_time_response(
                replace(template, response=contrast_spec), window='hann', zero_pad_factor=8)
        if time is None:
            time = products[tag, 1, 'total'].time * 1e9
        picks, rows = [], []
        for i, si in enumerate(range(1, N_STATIONS+1)):
            total = products[tag, si, 'total']
            contrast = products[tag, si, 'contrast']
            mgate = np.abs(time - relief[i]) <= 25.0
            env_c = np.abs(contrast.complex_envelope)
            a_echo = float(np.max(env_c[mgate]))
            a_direct = float(np.max(np.abs(total.complex_envelope)))
            pick = float(time[mgate][np.argmax(env_c[mgate])])
            picks.append(pick)
            rows.append({'station_index': si, 'src_x_m': float(station_x[i]),
                         'echo_over_direct_dB': float(20*np.log10(a_echo/a_direct)),
                         'pick_ns': pick, 'true_relief_ns': float(relief[i])})
        picks = np.array(picks)
        in_relief = (station_x >= RELIEF_X_108[0]) & (station_x <= RELIEF_X_108[1])
        flat = ~in_relief
        def corr(mask):
            if mask.sum() < 3:
                return None
            return float(np.corrcoef(picks[mask], relief[mask])[0, 1])
        def rmse(mask):
            return float(np.sqrt(np.mean((picks[mask]-relief[mask])**2)))
        central13 = np.abs(station_x - 53.6) <= 3.0
        central25 = np.abs(station_x - 53.6) <= 6.0
        per_material[tag] = {
            'cover_eps_r': eps,
            'echo_over_direct_dB_median': float(np.median([r['echo_over_direct_dB'] for r in rows])),
            'pick_corr_full_30m': corr(np.ones(N_STATIONS, bool)),
            'pick_corr_central12m': corr(central25),
            'pick_corr_central6m': corr(central13),
            'pick_corr_relief_segment': corr(in_relief),
            'pick_rmse_ns_full': rmse(np.ones(N_STATIONS, bool)),
            'flat_segment_pick_span_ns': float(picks[flat].max()-picks[flat].min()),
            'per_station': rows}

    # ---- figure 1: full-aperture B-scans
    fig, axes = plt.subplots(2, len(tags), figsize=(7.5*len(tags), 8), sharey='row', layout='constrained')
    if len(tags) == 1:
        axes = axes[:, None]
    for col, tag in enumerate(tags):
        relief = relief_twoway_108(xmids, eps_of[tag])
        for row, kind in enumerate(('total', 'contrast')):
            ax = axes[row, col]
            mat = np.stack([np.abs(products[tag, si, kind].complex_envelope) for si in range(1, N_STATIONS+1)], axis=1)
            mdb = 20*np.log10(np.maximum(mat, mat.max()*1e-9)/mat.max())
            ax.imshow(mdb, aspect='auto', cmap='gray', vmin=-60, vmax=0,
                      extent=[station_x[0]-STEP/2, station_x[-1]+RX_OFFSET+STEP/2, time[-1], time[0]])
            ax.plot(xmids, relief, 'r-', lw=1.0)
            ax.axvline(RELIEF_X_108[0], color='b', ls=':', lw=0.7)
            ax.axvline(RELIEF_X_108[1], color='b', ls=':', lw=0.7)
            ax.set_ylim(260, 0)
            ax.set_title(f'{tag} {kind}', fontsize=9)
            ax.set_xlabel('src x, 108 m domain (m)')
            if col == 0:
                ax.set_ylabel('time (ns)')
    fig.suptitle('30 m aperture B-scans at 8 m (self-normalized envelope dB; red = true relief two-way; blue dotted = relief segment)')
    a.out.mkdir(parents=True)
    fig.savefig(a.out/'aperture_bscans.png', dpi=140)
    plt.close(fig)

    # ---- figure 2: picks vs true relief per material
    fig2, axes2 = plt.subplots(1, len(tags), figsize=(7.5*len(tags), 4.5), sharey=True, layout='constrained')
    if len(tags) == 1:
        axes2 = [axes2]
    for col, tag in enumerate(tags):
        ax = axes2[col]
        rows = per_material[tag]['per_station']
        p = np.array([r['pick_ns'] for r in rows]); t = np.array([r['true_relief_ns'] for r in rows])
        ax.plot(station_x, t, 'r-', lw=1.5, label='true relief two-way')
        ax.plot(station_x, p, 'k.', ms=4, label='contrast pick')
        ax.axvline(RELIEF_X_108[0], color='b', ls=':', lw=0.7); ax.axvline(RELIEF_X_108[1], color='b', ls=':', lw=0.7)
        ax.invert_yaxis(); ax.set_xlabel('src x (m)'); ax.set_title(tag, fontsize=9)
        ax.legend(fontsize=7); ax.grid(alpha=0.3)
        m = per_material[tag]
        ax.text(0.02, 0.02, f"corr full={m['pick_corr_full_30m']:.3f} 6m={m['pick_corr_central6m']:.3f}",
                transform=ax.transAxes, fontsize=8)
    axes2[0].set_ylabel('time (ns)')
    fig2.savefig(a.out/'aperture_picks_vs_relief.png', dpi=140)
    plt.close(fig2)

    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'code_sha256': sha256(__file__), 'capsule': CAPSULE.name,
               'n_stations': N_STATIONS, 'aperture_m': float(station_x[-1]-station_x[0]+RX_OFFSET),
               'independent_dft_max_rel_l2': max(dft_err.values()),
               'per_material': {t: {k: v for k, v in per_material[t].items() if k != 'per_station'} for t in tags},
               'limits': '2D line source; picks gated +/-25 ns on true relief; flat-segment pick span includes relief-edge transition stations; diagnostic, not physical acceptance.'}
    with open(a.out/'summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary['per_material'], ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
