"""Debye-parameter scan analysis: echo levels, in-band dispersion slope vs measured anchor.

Reads the completed 2026-10-05_hs4_debye_scan_r1 capsule (13 parameter sets x
2 roles x 13 stations, 8 m altitude) and measures per parameter set:
- echo/direct and echo/surface dB (gated contrast peak, 13 stations);
- analytic band attenuation/dispersion: eps'(f), alpha(f) at 20/95/170 MHz,
  two-way loss at model cover and extrapolated to real-depth candidates
  (DECLARED MODEL, band-point values);
- relief-pick correlation;
- comparison against the measured Line9 anchors.
Outputs summary.json + annotated figures into a fresh result directory.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

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
from analyze_hs4_material_scan import real_anchors, relief_twoway, ratio_db

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

CHECKS = ROOT / 'artifacts/research_checks'
CAPSULE = CHECKS / '2026-10-05_hs4_debye_scan_r1'
TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
STEP, RX_OFFSET, ANT_Z = 0.5, 1.3, 20.0
EPS0 = 8.8541878e-12
C_SI = 2.99792458e8
REAL_DEPTHS_M = [13.0, 18.5]


def band_props(eps_s, sig, de, tau, freqs):
    w = 2*np.pi*freqs
    if de:
        er1 = eps_s + de/(1+(w*tau)**2)
        er2 = de*(w*tau)/(1+(w*tau)**2) + sig/(w*EPS0)
    else:
        er1 = np.full_like(freqs, eps_s)
        er2 = np.full_like(freqs, sig/(w*EPS0))
    alpha = (w/C_SI)*er2/(2*np.sqrt(er1))
    return er1, alpha


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

    fx, fz = interface_relief()
    station_x = np.array([14.6 + STEP*k for k in range(13)])
    xmids = station_x + RX_OFFSET/2
    freqs = np.linspace(20e6, 170e6, 501)

    spectra, dft_err = {}, {}
    for g in contract['groups']:
        audit = json.loads((Path(g['input']).parent/'audit.json').read_text('utf-8'))
        for row in audit:
            raw = Path(g['input']).parent/row['file']
            if sha256(raw) != row['sha256']:
                raise ValueError('raw identity differs')
            r = response(raw)
            spectra[g['material_tag'], g['role'], row['station_index']] = r
            exact = direct_response(raw, r.frequency[TAKE])
            err = float(np.linalg.norm(exact - r.response[TAKE])/np.linalg.norm(exact))
            if err > 1e-9:
                raise ValueError('independent DFT failed')
            dft_err[f"{g['id']}/{row['file']}"] = err

    template = next(iter(spectra.values()))
    info = {g['material_tag']: g for g in contract['groups']}
    tags = list(dict.fromkeys(g['material_tag'] for g in contract['groups']))
    products = {}
    time = None
    for tag in tags:
        for si in range(1, 122, 10):
            rough = spectra[tag, 'rough', si]
            cs = replace(template, response=rough.response - spectra[tag, 'halfspace', si].response)
            products[tag, si, 'total'] = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            products[tag, si, 'contrast'] = reconstruct_time_response(cs, window='hann', zero_pad_factor=8)
        if time is None:
            time = products[tag, 1, 'total'].time * 1e9

    surface_t = 2*np.hypot(0.65, ANT_Z-SURFACE_Z)/C_AIR
    per = {}
    for tag in tags:
        g = info[tag]
        er1, alpha = band_props(g['cover_eps_r'], g['cover_sigma'], g['debye_dEps'], g['debye_tau_s'], freqs)
        eps95 = float(er1[250])
        relief = relief_twoway(ANT_Z, xmids, fx, fz, C_AIR/np.sqrt(eps95))
        rows = []
        for i, si in enumerate(range(1, 122, 10)):
            total = products[tag, si, 'total']
            contrast = products[tag, si, 'contrast']
            mgate = np.abs(time - relief[i]) <= 25.0
            msurf = np.abs(time - surface_t) <= 12.0
            env_c = np.abs(contrast.complex_envelope)
            a_echo = float(np.max(env_c[mgate]))
            a_direct = float(np.max(np.abs(total.complex_envelope)))
            a_surface = float(np.max(np.abs(total.complex_envelope[msurf])))
            rows.append({'station_index': si,
                         'echo_over_direct_dB': ratio_db(a_echo, a_direct),
                         'echo_over_surface_dB': ratio_db(a_echo, a_surface),
                         'pick_ns': float(time[mgate][np.argmax(env_c[mgate])])})
        picks = np.array([r['pick_ns'] for r in rows])
        med = float(np.median([r['echo_over_direct_dB'] for r in rows]))
        am = float(alpha.mean())
        model_cover = float(np.mean(SURFACE_Z - np.interp(xmids, fx, fz)))
        per[tag] = {
            'debye_dEps': g['debye_dEps'], 'debye_tau_ns': g['debye_tau_s']*1e9,
            'cover_sigma': g['cover_sigma'],
            'eps_r_at_MHz': {'20': float(er1[0]), '95': eps95, '170': float(er1[-1])},
            'alpha_Np_m_at_MHz': {'20': float(alpha[0]), '95': float(alpha[250]), '170': float(alpha[-1])},
            'two_way_dB_model_cover': 8.686*2*am*model_cover,
            'two_way_dB_at_real_depth': {f'{d:g}m': 8.686*2*am*d for d in REAL_DEPTHS_M},
            'echo_over_direct_dB_median': med,
            'echo_over_surface_dB_median': float(np.median([r['echo_over_surface_dB'] for r in rows])),
            'extrapolated_echo_over_direct_dB': {f'{d:g}m': med + 8.686*2*am*(model_cover-d) for d in REAL_DEPTHS_M},
            'relief_pick_corr': float(np.corrcoef(picks, relief)[0, 1])}

    real = real_anchors()
    real_band = real['interface_band_db_mean']

    # ---- annotated figure: echo/direct (model + extrapolated) vs real anchor ----
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), layout='constrained')
    labels = [t for t in tags]
    x = np.arange(len(labels))
    med3 = [per[t]['echo_over_direct_dB_median'] for t in labels]
    ex13 = [per[t]['extrapolated_echo_over_direct_dB']['13m'] for t in labels]
    ex18 = [per[t]['extrapolated_echo_over_direct_dB']['18.5m'] for t in labels]
    ax = axes[0]
    ax.plot(x, med3, 'o-', label='模型覆盖层 3 m（实测值）')
    ax.plot(x, ex13, 's--', label='外推到 13 m 覆盖层')
    ax.plot(x, ex18, '^--', label='外推到 18.5 m 覆盖层')
    ax.axhline(real_band, color='k', ls=':', lw=1.5, label=f'实测界面条带锚点 {real_band:.1f} dB')
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('界面回波/直达波 (dB)'); ax.legend(fontsize=9); ax.grid(alpha=0.3)
    ax.set_title('Debye 参数细调：回波/直达 vs 实测锚点\n（外推为声明模型：带内均值衰减，不含扩散）', fontsize=12)
    ax2 = axes[1]
    for t in labels:
        g = info[t]
        er1, alpha = band_props(g['cover_eps_r'], g['cover_sigma'], g['debye_dEps'], g['debye_tau_s'], freqs)
        ax2.plot(freqs/1e6, 8.686*alpha*2, lw=1.2, label=t)
    ax2.set_xlabel('频率 (MHz)'); ax2.set_ylabel('双程衰减 (dB/m)')
    ax2.set_title('各参数组合的带内衰减谱 α(f)\n（SFCW 逐频链路直接感受的色散形状）', fontsize=12)
    ax2.legend(fontsize=6, ncol=2); ax2.grid(alpha=0.3)
    a.out.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out/'debye_scan_vs_anchor.png', dpi=140)
    plt.close(fig)

    summary = {'status': 'COMPLETED_SITE_ADAPTATION_DIAGNOSTIC',
               'code_sha256': sha256(__file__), 'capsule': CAPSULE.name,
               'independent_dft_max_rel_l2': max(dft_err.values()),
               'real_anchors': {k: v for k, v in real.items() if k not in ('t_ns', 'mean_db')},
               'per_parameter_set': per,
               'limits': '2D line source; extrapolation uses band-mean analytic alpha (Debye + conductive), no spreading; real anchors from one CSV; site-adaptation diagnostic.'}
    with open(a.out/'summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps({t: {'e/d_3m': per[t]['echo_over_direct_dB_median'],
                          'extrap18.5': per[t]['extrapolated_echo_over_direct_dB']['18.5m'],
                          'tw18.5': per[t]['two_way_dB_at_real_depth']['18.5m']}
                      for t in tags}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
