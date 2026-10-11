"""Render Chinese-labeled ablation figures from the ROG-produced NPZ/JSON (local, no gprMax)."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

VARIANT_LABEL = {
    'base': '基线（现行材料）',
    'c1_sandstone_paper': 'C1 砂岩=论文常数 εr9/σ0.001',
    'c2_cover_paper_soil': 'C2 覆盖层=论文土 εr12/σ0.003',
    'c3_cover_er184': 'C3 覆盖层 εr=18.4（归因证伪）',
}
ORDER = ['base', 'c1_sandstone_paper', 'c2_cover_paper_soil', 'c3_cover_er184']


def load(npz_path):
    z = np.load(npz_path)
    out = {}
    for cid in z.files:
        t, x, env, resp = z[cid]
        out[cid] = dict(t=t, x=x, env=env, resp=resp)
    return out


def fig_traces(cases, role, out):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                             gridspec_kw=dict(height_ratios=[2, 1]))
    for v in ORDER:
        d = cases[f'abl_{v}_x190_{role}']
        axes[0].plot(d['t'] * 1e9, d['x'], lw=.6, label=VARIANT_LABEL[v])
        axes[1].plot(d['t'] * 1e9, d['env'], lw=.8, label=VARIANT_LABEL[v])
    axes[0].set_ylabel('Ez (V/m)')
    axes[1].set_ylabel('包络 (V/m)')
    axes[1].set_xlabel('时间 (ns)')
    axes[0].legend(fontsize=8, loc='upper right')
    axes[0].set_title(f'材料消融 {role}：时域迹线（x190 站位，AGL 5 m，官方 impulse）')
    axes[1].set_title('包络（Hilbert）')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.set_xlim(0, 1200)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_diff(cases, out):
    fig, ax = plt.subplots(figsize=(11, 4))
    for v in ORDER:
        d1 = cases[f'abl_{v}_x190_H1']
        d0 = cases[f'abl_{v}_x190_H0']
        ax.plot(d1['t'] * 1e9, d1['x'] - d0['x'], lw=.6, label=VARIANT_LABEL[v])
    ax.set_xlabel('时间 (ns)')
    ax.set_ylabel('ΔEz (V/m)')
    ax.set_title('材料消融：H1−H0 配对差分（差分内为滑坡体相关响应）')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    ax.set_xlim(0, 1200)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_c3_shift(cases, analysis, out):
    fig, ax = plt.subplots(figsize=(11, 4.5))
    b = cases['abl_base_x190_H1']
    c = cases['abl_c3_cover_er184_x190_H1']
    ax.plot(b['t'] * 1e9, b['env'], lw=.9, label='基线 H1 包络')
    ax.plot(c['t'] * 1e9, c['env'], lw=.9, label='C3（覆盖层 εr=18.4）H1 包络')
    test = analysis['c3_shift_test']
    ax.axvline(test['base_multiple_peak_ns'], color='C0', ls='--', lw=.8,
               label=f"基线多次波峰 {test['base_multiple_peak_ns']:.0f} ns")
    ax.axvline(test['c3_multiple_peak_ns'], color='C1', ls='--', lw=.8,
               label=f"C3 实测峰 {test['c3_multiple_peak_ns']:.0f} ns")
    ax.axvline(test['predicted_c3_ns'], color='C1', ls=':', lw=1.2,
               label=f"C3 预测峰 {test['predicted_c3_ns']:.0f} ns（按覆盖层多次波 √εr 缩放）")
    ax.set_xlabel('时间 (ns)')
    ax.set_ylabel('包络 (V/m)')
    ax.set_title('归因检验：若强带为覆盖层内多次返回，C3 峰应按 √(18.4/11)≈1.293 缩放地下路径')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    ax.set_xlim(120, 550)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_freq(cases, out):
    freq = 20 + np.arange(501) * .3
    fig, ax = plt.subplots(figsize=(11, 4))
    for v in ORDER:
        d = cases[f'abl_{v}_x190_H1']
        ax.plot(freq, 20 * np.log10(np.abs(d['resp']) / np.max(np.abs(d['resp']))), lw=.8,
                label=VARIANT_LABEL[v])
    ax.set_xlabel('频率 (MHz)')
    ax.set_ylabel('归一化 |H(f)| (dB)')
    ax.set_title('材料消融 H1：官方 direct SFCW 频响（20–170 MHz，501 音调）')
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    ax.set_ylim(-60, 0)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--npz', type=Path, required=True)
    p.add_argument('--analysis', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cases = load(a.npz)
    analysis = json.loads(a.analysis.read_text('utf-8'))
    fig_traces(cases, 'H1', a.out / 'fig1_H1_traces.png')
    fig_traces(cases, 'H0', a.out / 'fig2_H0_traces.png')
    fig_diff(cases, a.out / 'fig3_H1_minus_H0.png')
    fig_c3_shift(cases, analysis, a.out / 'fig4_c3_shift_test.png')
    fig_freq(cases, a.out / 'fig5_sfcw_response.png')
    print('figures written to', a.out)
