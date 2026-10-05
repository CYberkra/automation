"""Debye scan B-scan report figure (user requirement: every solved batch ships annotated B-scans).

Reads the completed 2026-10-05_hs4_debye_scan_r1 capsule plus the accepted
analysis summary.json, and renders grey-scale envelope-dB B-scan panels for six
representative parameter sets (total / rough-minus-halfspace contrast rows),
with Chinese annotations, per-set echo/direct dB in titles, the true relief
two-way curve overlay, and the measured Line9 anchor in the suptitle.
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
from analyze_hs4_material_scan import relief_twoway
from analyze_hs4_debye_scan import CAPSULE, TAKE, STEP, RX_OFFSET, ANT_Z, band_props

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

CHECKS = ROOT / 'artifacts/research_checks'
ANALYSIS = CHECKS / '2026-10-05_hs4_debye_scan_analysis_r1'
SHOW = [
    ('de0.5_tau1', 'Δε=0.5, τ=1ns（弛豫峰159MHz）'),
    ('de0.5_tau6.4567', 'Δε=0.5, τ=6.46ns【与实测相容】'),
    ('de1_tau6.4567', 'Δε=1, τ=6.46ns'),
    ('de2_tau6.4567', 'Δε=2, τ=6.46ns（过度损耗）'),
    ('de4_tau6.4567', 'Δε=4, τ=6.46ns（过度损耗）'),
    ('de7.878_tau6.4567_s0.003', '归档配方 Δε=7.878（原失真来源）'),
]
ROW_LABEL = {'total': '原始总场', 'contrast': '配对差分（起伏−全覆盖层=纯界面响应）'}


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
    summary = json.loads((ANALYSIS/'summary.json').read_text('utf-8'))
    per = summary['per_parameter_set']
    anchor = summary['real_anchors']['interface_band_db_mean']

    keep = {t for t, _ in SHOW}
    fx, fz = interface_relief()
    station_x = np.array([14.6 + STEP*k for k in range(13)])
    xmids = station_x + RX_OFFSET/2
    freqs = np.linspace(20e6, 170e6, 501)

    spectra, dft_max = {}, 0.0
    for g in contract['groups']:
        if g['material_tag'] not in keep:
            continue
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
            dft_max = max(dft_max, err)

    template = next(iter(spectra.values()))
    products = {}
    time = None
    for tag, _ in SHOW:
        for si in range(1, 122, 10):
            rough = spectra[tag, 'rough', si]
            cs = replace(template, response=rough.response - spectra[tag, 'halfspace', si].response)
            products[tag, si, 'total'] = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            products[tag, si, 'contrast'] = reconstruct_time_response(cs, window='hann', zero_pad_factor=8)
        if time is None:
            time = products[tag, 1, 'total'].time * 1e9

    info = {g['material_tag']: g for g in contract['groups']}
    fig, axes = plt.subplots(2, len(SHOW), figsize=(3.0*len(SHOW), 8.0), sharey=True, layout='constrained')
    for col, (tag, label) in enumerate(SHOW):
        g = info[tag]
        er1, _ = band_props(g['cover_eps_r'], g['cover_sigma'], g['debye_dEps'], g['debye_tau_s'], freqs)
        relief = relief_twoway(ANT_Z, xmids, fx, fz, C_AIR/np.sqrt(float(er1[250])))
        med = per[tag]['echo_over_direct_dB_median']
        ext = per[tag]['extrapolated_echo_over_direct_dB']['18.5m']
        a_direct = max(float(np.max(np.abs(products[tag, si, 'total'].complex_envelope)))
                       for si in range(1, 122, 10))
        for row, kind in enumerate(('total', 'contrast')):
            ax = axes[row, col]
            mat = np.stack([np.abs(products[tag, si, kind].complex_envelope) for si in range(1, 122, 10)], axis=1)
            if kind == 'total':
                mdb = 20*np.log10(np.maximum(mat, mat.max()*1e-9)/mat.max())
                vmin, vmax = -60, 0
            else:
                # shared cross-column scale: contrast relative to the direct wave, so
                # high-loss recipes visually sink below the floor instead of being renormalized
                mdb = 20*np.log10(np.maximum(mat, a_direct*1e-12)/a_direct)
                vmin, vmax = -70, -30
            ax.imshow(mdb, aspect='auto', cmap='gray', vmin=vmin, vmax=vmax,
                      extent=[station_x[0]-STEP/2, station_x[-1]+RX_OFFSET+STEP/2, time[-1], time[0]])
            ax.plot(xmids, relief, 'r-', lw=0.9)
            ax.set_xlim(station_x[0]-STEP/2, station_x[-1]+RX_OFFSET+STEP/2)
            ax.set_ylim(260, 0)
            if row == 0:
                ax.set_title(f'{label}\n回波/直达@3m {med:.1f} dB，外推18.5m {ext:.1f} dB', fontsize=8)
            else:
                ax.set_xlabel('发射站位 x (m)', fontsize=8)
            if col == 0:
                ax.set_ylabel(f'{ROW_LABEL[kind]}\n时间 (ns)', fontsize=9)
            else:
                ax.text(0.02, 0.02, ROW_LABEL[kind], transform=ax.transAxes, fontsize=7,
                        color='yellow', va='bottom')
    fig.suptitle(
        'HS4 起伏界面 8 m 航高 Debye 细调批 B-scan（上排：原始总场灰度包络 dB 各自归一；下排：配对差分相对直达波 dB，'
        f'跨列同一标尺 −70~−30，红线=真实起伏双程时；实测 Line9 界面条带锚点 {anchor:.1f} dB）\n'
        'Δε=0.5/τ=6.46ns 列：界面条带沿红线明亮可见；Δε≥2 与归档配方列：条带沉到标尺下限附近，与实测无此现象矛盾',
        fontsize=11)
    a.out.mkdir(parents=True)
    fig.savefig(a.out/'debye_bscans.png', dpi=130)
    plt.close(fig)
    print(json.dumps({'figure': str(a.out/'debye_bscans.png'), 'independent_dft_max': dft_max},
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
