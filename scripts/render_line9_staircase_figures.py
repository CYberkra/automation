"""Render staircase ablation figures (Chinese labels) from retrieved NPZ + analysis JSONs."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

GRIDS = {'d25': '2.5 cm（生产网格）', 'd125': '1.25 cm', 'd625': '0.625 cm（代理光滑参考）'}
SURFACES = {'flat': '平地表（对照，零阶梯）', 'slope': '5° 斜坡（纯阶梯）', 'rough': 'fBm 微起伏（RMS≈5 cm）'}
COLORS = {'d25': 'tab:red', 'd125': 'tab:orange', 'd625': 'tab:blue'}


def main(a):
    z = np.load(a.npz)
    zb = np.load(a.out / 'staircase_bandlimited.npz')
    ana = json.loads((a.out / 'staircase_analysis.json').read_text('utf-8'))
    det = json.loads((a.out / 'staircase_freq_detail.json').read_text('utf-8'))
    freq = z['frequency_Hz'] / 1e6
    t_ns = zb['time_s'] * 1e9
    fig_dir = a.out / 'figures'
    fig_dir.mkdir(exist_ok=True)

    # fig1: per-surface response magnitude across grids
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    for ax, (s, title) in zip(axes, SURFACES.items()):
        for g in ['d625', 'd125', 'd25']:
            r = z[f'stair_{s}_{g}__response']
            ax.plot(freq, 20 * np.log10(np.abs(r) / np.abs(r).max()), lw=0.9,
                    color=COLORS[g], label=GRIDS[g])
        ax.set_title(f'{title}：SFCW 响应幅度（各自归一）')
        ax.set_xlabel('频率 (MHz)')
        ax.legend(fontsize=8)
    axes[0].set_ylabel('响应幅度 (dB)')
    fig.suptitle('阶梯消融批：三种地表 × 三种网格的官方 direct SFCW 复谱（20–170 MHz，501 音调）')
    fig.tight_layout()
    fig.savefig(fig_dir / 'fig1_sfcw_response_magnitude.png', dpi=160)
    plt.close(fig)

    # fig2: per-tone complex relative difference vs reference grid
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    for ax, (s, title) in zip(axes, SURFACES.items()):
        ref = z[f'stair_{s}_d625__response']
        for g in ['d25', 'd125']:
            r = z[f'stair_{s}_{g}__response']
            ax.plot(freq, 100 * np.abs(r - ref) / np.abs(ref), lw=0.9,
                    color=COLORS[g], label=f'{GRIDS[g]} 对参考')
        ax.set_title(f'{title}：逐音调复数相对差')
        ax.set_xlabel('频率 (MHz)')
        ax.set_yscale('log')
        ax.legend(fontsize=8)
    axes[0].set_ylabel('逐音调相对差 (%)')
    fig.suptitle('阶梯消融批：各网格对 0.625 cm 代理光滑参考的逐音调复数相对差（对数轴）')
    fig.tight_layout()
    fig.savefig(fig_dir / 'fig2_per_tone_relative_difference.png', dpi=160)
    plt.close(fig)

    # fig3: band-limited impulse-response envelopes on common time base
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    for ax, (s, title) in zip(axes, SURFACES.items()):
        for g in ['d625', 'd125', 'd25']:
            env = np.abs(zb[f'{s}_{g}__bl'])
            ax.plot(t_ns, 20 * np.log10(env / env.max()), lw=0.8, color=COLORS[g], label=GRIDS[g])
        ax.set_title(f'{title}：带限冲激响应包络（共同时基）')
        ax.set_xlabel('时间 (ns)')
        ax.set_xlim(0, 900)
        ax.set_ylim(-80, 0)
        ax.legend(fontsize=8)
    axes[0].set_ylabel('包络 (dB)')
    fig.suptitle('阶梯消融批：由 501 音调复谱反演的带限冲激响应（20–170 MHz，无插值、共同时基）')
    fig.tight_layout()
    fig.savefig(fig_dir / 'fig3_bandlimited_envelopes.png', dpi=160)
    plt.close(fig)

    # fig4: summary bars + S5 physical anchor
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    labels, vals, cols = [], [], []
    for s in ['flat', 'slope', 'rough']:
        for g in ['d25', 'd125']:
            labels.append(f'{SURFACES[s].split("（")[0]}\n{g.replace("d25", "2.5cm").replace("d125", "1.25cm")}')
            vals.append(100 * ana['surfaces'][s]['frequency'][g]['E_raw'])
            cols.append(COLORS[g])
    ax.bar(range(len(vals)), vals, color=cols)
    ax.axhline(42, color='k', ls='--', lw=1.2, label='S5 实测物理微起伏效应下限（42%，多次波带）')
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('带内复数 RMS 相对差 E (%)')
    ax.set_title('阶梯/离散化误差（本批，≤1.6%）对比真实物理起伏效应（S5，42–73%）：低约两个量级')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / 'fig4_error_vs_physical_anchor.png', dpi=160)
    plt.close(fig)
    print(json.dumps(dict(figures=[str(p) for p in sorted(fig_dir.glob('*.png'))])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--npz', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    main(a)
