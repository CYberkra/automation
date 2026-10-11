"""Render Chinese-labeled rough-vs-smooth sensitivity figures (local, no gprMax)."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

SURF_LABEL = {'rough': '粗糙地表（fBm H=0.5，RMS 5.4cm）', 'smooth': '平滑地表（消融批基线）'}


def load(npz_path):
    z = np.load(npz_path)
    out = {}
    for key in z.files:
        cid, name = key.rsplit('__', 1)
        out.setdefault(cid, {})[name] = z[key]
    return out


def fig_traces(cases, role, out):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                             gridspec_kw=dict(height_ratios=[2, 1]))
    for surf in ['smooth', 'rough']:
        d = cases[f'{role}_{surf}']
        axes[0].plot(d['t'] * 1e9, d['x'], lw=.6, label=SURF_LABEL[surf])
        axes[1].plot(d['t'] * 1e9, d['env'], lw=.8, label=SURF_LABEL[surf])
    axes[0].set_ylabel('Ez (V/m)')
    axes[1].set_ylabel('包络 (V/m)')
    axes[1].set_xlabel('时间 (ns)')
    axes[0].legend(fontsize=8)
    axes[0].set_title(f'地表起伏敏感性 {role}：时域迹线（x190 站位，AGL 5 m，官方 impulse）')
    axes[1].set_title('包络（Hilbert）')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.set_xlim(0, 1200)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_diff(cases, out):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for surf in ['smooth', 'rough']:
        d1, d0 = cases[f'H1_{surf}'], cases[f'H0_{surf}']
        axes[0].plot(d1['t'] * 1e9, d1['x'] - d0['x'], lw=.6, label=SURF_LABEL[surf])
    axes[0].set_ylabel('ΔEz (V/m)')
    axes[0].set_title('H1−H0 配对差分（滑体相关响应）：粗糙 vs 平滑')
    axes[0].legend(fontsize=8)
    d = cases['H1_rough']
    ds = cases['H1_smooth']
    axes[1].plot(d['t'] * 1e9, d['x'] - ds['x'], lw=.6, color='C2')
    axes[1].set_ylabel('ΔEz (V/m)')
    axes[1].set_xlabel('时间 (ns)')
    axes[1].set_title('起伏引入的变化：H1 粗糙 − H1 平滑')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.set_xlim(0, 1200)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_freq(cases, out):
    freq = 20 + np.arange(501) * .3
    fig, ax = plt.subplots(figsize=(11, 4))
    for surf in ['smooth', 'rough']:
        d = cases[f'H1_{surf}']
        ax.plot(freq, 20 * np.log10(np.abs(d['response']) / np.max(np.abs(d['response']))), lw=.8,
                label=SURF_LABEL[surf])
    ax.set_xlabel('频率 (MHz)')
    ax.set_ylabel('归一化 |H(f)| (dB)')
    ax.set_title('地表起伏敏感性 H1：官方 direct SFCW 频响（20–170 MHz，501 音调）')
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
    fig_traces(cases, 'H1', a.out / 'fig1_rough_H1_traces.png')
    fig_traces(cases, 'H0', a.out / 'fig2_rough_H0_traces.png')
    fig_diff(cases, a.out / 'fig3_rough_diff.png')
    fig_freq(cases, a.out / 'fig4_rough_sfcw.png')
    print('figures written to', a.out)
