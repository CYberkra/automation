"""Visualize one batch2d_v1_mt gather (TGT, family BG, paired difference)."""
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt

for _f in ('msyh.ttc', 'simhei.ttf'):
    try:
        font_manager.fontManager.addfont(str(Path(r'C:\Windows\Fonts') / _f))
        matplotlib.rcParams['font.family'] = font_manager.FontProperties(
            fname=str(Path(r'C:\Windows\Fonts') / _f)).get_name()
        break
    except Exception:
        continue
matplotlib.rcParams['axes.unicode_minus'] = False
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-MT33'
BG = 'B2D-C3m-BG-MT33'
DT_NS = 5.896635841874211e-11 * 1e9
OUT = ROOT / 'artifacts/simulations/2026-09-26_gather_preview_c3_a0.png'


def load(rid):
    p = ROOT / 'artifacts/research_checks' / f'2026-09-26_{rid}' / f'{rid}.h5'
    with h5py.File(p, 'r') as f:
        x = np.stack([f[f'rxs/rx{k}/Ex'][:] for k in range(1, 34)], axis=1)
    return x


def panel(ax, x, title, t_ns, vmax):
    # 逐道归一化显示（显示用），数据本身不做任何处理
    g = x / np.maximum(np.abs(x).max(axis=0, keepdims=True), 1e-30)
    ax.imshow(g, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax,
              extent=(0.5, 33.5, t_ns[-1], t_ns[0]), interpolation='nearest')
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('道号 mt01..mt33（y 12.65→20.65 m）')
    ax.set_ylabel('时间 (ns)')


def main():
    tgt = load(CASE)
    bg = load(BG)
    diff = tgt - bg
    t_ns = np.arange(tgt.shape[0]) * DT_NS
    vmax = 0.5
    fig, axes = plt.subplots(1, 3, figsize=(16, 7.5), sharey=True)
    panel(axes[0], tgt, f'目标道集（含空腔）\n{CASE}', t_ns, vmax)
    panel(axes[1], bg, '族 BG 道集（无空腔）\nB2D-C3m-BG-MT33', t_ns, vmax)
    # 差分面板：固定公共色标（显示真实相对幅度）
    gd = diff / np.abs(diff).max()
    im = axes[2].imshow(gd, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax,
                        extent=(0.5, 33.5, t_ns[-1], t_ns[0]), interpolation='nearest')
    axes[2].set_title('配对差分 TGT−BG（公共色标）\n空腔响应（约 320 ns 附近）', fontsize=11)
    axes[2].set_xlabel('道号 mt01..mt33（y 12.65→20.65 m）')
    axes[2].axhline(332, color='k', lw=0.8, ls='--', alpha=0.6)
    for ax in axes[1:2]:
        ax.set_ylabel('时间 (ns)')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label='归一化幅度')
    fig.suptitle('batch2d_v1_mt 变偏移距单发多收道集（CSG，非 B-scan）：C3 族锚点 A0，Tx y15.35，偏移距 2.70→1.30→5.30 m', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT, dpi=130)
    print(OUT)


if __name__ == '__main__':
    main()
