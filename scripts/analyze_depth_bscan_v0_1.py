"""Depth B-scan visualization: 13 stations x 4 cases, total field and isolated interface.

Row 1: total field, shared fixed colour scale (deep events invisible by dynamic range).
Row 2: isolated interface (flat minus shared fullcover, declared target/background pair),
shared scale across depths, with the 95 MHz ray-predicted constant arrival overlaid.
Row 3: same isolated data with per-panel self-normalized dB envelope (declared display
enhancement only, not an amplitude-fidelity product).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from analyze_uav_local_free_space import FREQ, transfer
from green_halfspace_hed_v0_1 import C_SI, epsr_cover
from hs_capsule_identity import sha256

CASES = ('fullcover', 'depth3', 'depth10', 'depth18_5')
DEPTHS = {'depth3': 3., 'depth10': 10., 'depth18_5': 18.5}
LABELS = {'fullcover': '全覆盖背景', 'depth3': '平界面3m', 'depth10': '平界面10m', 'depth18_5': '平界面18.5m'}
AIR_NS = 16./C_SI*1e9


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('fresh analysis output required')
    a.out.mkdir(parents=True)
    c = json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    audit = json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha256(a.study/'execution_contract.json'):
        raise ValueError('completed raw identity required')
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}
    stations = c['geometry']['stations_tx_x_m']
    mats, times = {}, None
    dft_max = 0.
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        value, source, receiver = transfer(raw, 1, 'Ey', g['spacing_m'])
        take = np.array([0, 83, 167, 250, 333, 417, 500])
        sy = source.samples.copy(); ry = receiver.samples.copy()
        n = round(20e-9/receiver.dt)
        ry[-n:] *= .5*(1+np.cos(np.linspace(0, np.pi, n)))
        sft = np.exp(-2j*np.pi*FREQ[take, None]*(source.time_offset+source.dt*np.arange(len(sy))))@sy
        rft = np.exp(-2j*np.pi*FREQ[take, None]*(receiver.time_offset+receiver.dt*np.arange(len(ry))))@ry
        independent = rft/sft/g['spacing_m']
        dft_max = max(dft_max, float(np.linalg.norm(value[take]-independent)/np.linalg.norm(independent)))
        t = receiver.dt*np.arange(len(receiver.samples))*1e9
        times = t
        mats.setdefault(g['case'], np.zeros((len(t), len(stations))))[:, g['station_index']] = receiver.samples
    if dft_max > 1e-9:
        raise ValueError('independent direct DFT mismatch')
    iso = {case: mats[case]-mats['fullcover'] for case in DEPTHS}
    n95 = float(np.sqrt(epsr_cover(np.array([95e6])).real[0]))
    rows = []
    for case, depth in DEPTHS.items():
        ratio = np.linalg.norm(iso[case], axis=0)/np.linalg.norm(mats['fullcover'], axis=0)
        pred = AIR_NS + 2*depth*n95/C_SI*1e9
        peak_t = times[np.argmax(np.abs(iso[case]), axis=0)]
        rows.append({'case': case, 'depth_m': depth,
                     'isolated_over_fullcover_per_station_median': float(np.median(ratio)),
                     'isolated_over_fullcover_per_station_minmax': [float(ratio.min()), float(ratio.max())],
                     'predicted_arrival_ns': float(pred),
                     'isolated_peak_ns_median': float(np.median(peak_t)),
                     'isolated_peak_ns_span': [float(peak_t.min()), float(peak_t.max())]})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(3, 4, figsize=(17, 11), layout='constrained')
    extent = [stations[0], stations[-1], times[-1], 0]
    v_total = abs(mats['fullcover']).max()*0.005
    v_iso = max(abs(iso[k]).max() for k in DEPTHS)
    for j, case in enumerate(CASES):
        ax = axes[0, j]
        im = ax.imshow(mats[case], aspect='auto', cmap='RdBu_r', vmin=-v_total, vmax=v_total, extent=extent)
        ax.set_title(LABELS[case]+'：总场', fontsize=10)
        fig.colorbar(im, ax=ax, shrink=.8, label='Ey(V/m)')
        ax.set_ylim(650, 0)
    for j, case in enumerate(DEPTHS):
        ax = axes[1, j]
        im = ax.imshow(iso[case], aspect='auto', cmap='RdBu_r', vmin=-v_iso, vmax=v_iso, extent=extent)
        pred = AIR_NS + 2*DEPTHS[case]*n95/C_SI*1e9
        ax.axhline(pred, color='k', ls='--', lw=1, label=f'95MHz射线预测 {pred:.0f}ns')
        ax.set_title(LABELS[case]+'：隔离界面（共享标尺）', fontsize=10)
        ax.legend(fontsize=7, loc='lower right')
        fig.colorbar(im, ax=ax, shrink=.8, label='Ey(V/m)')
        env = 20*np.log10(np.abs(iso[case])/np.abs(iso[case]).max()+1e-12)
        ax2 = axes[2, j]
        im2 = ax2.imshow(env, aspect='auto', cmap='gray', vmin=-40, vmax=0, extent=extent)
        ax2.axhline(pred, color='r', ls='--', lw=1)
        ax2.set_title(LABELS[case]+'：隔离界面（自归一dB，显示增强）', fontsize=10)
        fig.colorbar(im2, ax=ax2, shrink=.8, label='dB')
    axes[1, 3].axis('off')
    axes[2, 3].axis('off')
    for ax in axes[:, 0]:
        ax.set_ylabel('时间(ns)')
    for ax in axes[2, :3]:
        ax.set_xlabel('发射站位 x(m)，接收=站位+1.3m')
    fig.suptitle('深度 B-scan：13站位×4场景（总场共享0.5%标尺；隔离=平界面−全覆盖背景，声明目标对比；自归一dB仅为显示增强）')
    fig.savefig(a.out/'depth_bscan_panels.png', dpi=140)
    plt.close(fig)
    summary = {'status': 'PASS', 'contract_sha256': sha256(a.study/'execution_contract.json'),
               'analysis_code_sha256': sha256(__file__), 'independent_DFT_relative_L2_max': dft_max,
               'rows': rows,
               'isolation_definition': 'interface = flat-interface case minus shared full-cover background (declared target/background pair)',
               'display_note': 'row3 per-panel self-normalized dB envelope is declared display enhancement, not amplitude fidelity',
               'scope': '2D line-source B-scan mechanism screening; depths are mechanism points, not Line9 truth.'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'dft_max': dft_max, 'rows': rows}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
