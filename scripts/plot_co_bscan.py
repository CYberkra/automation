"""B-scan style preview of the batch2d_v1_co common-offset diff (11 traces)."""
from pathlib import Path
import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np

for _f in ('msyh.ttc', 'simhei.ttf'):
    try:
        font_manager.fontManager.addfont(str(Path(r'C:\Windows\Fonts') / _f))
        matplotlib.rcParams['font.family'] = font_manager.FontProperties(
            fname=str(Path(r'C:\Windows\Fonts') / _f)).get_name()
        break
    except Exception:
        continue
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / 'artifacts/research_checks'
OUT = R / '2026-09-26_batch2d_v1_co_analysis/co_bscan_preview.png'
DT_NS = 5.896635841874211e-11 * 1e9
BG = 'B2D-C3m-BG-CO11-t{:02d}'
TG = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-CO11-t{:02d}'


def load(rid):
    with h5py.File(R / f'2026-09-26_{rid}' / f'{rid}.h5', 'r') as f:
        return f['rxs/rx1/Ex'][:]


def main():
    bg = np.stack([load(BG.format(k)) for k in range(1, 12)], axis=1)
    tg = np.stack([load(TG.format(k)) for k in range(1, 12)], axis=1)
    diff = tg - bg
    t_ns = np.arange(bg.shape[0]) * DT_NS
    vmax = 0.5
    fig, axes = plt.subplots(1, 3, figsize=(16, 7.5), sharey=True)
    for ax, x, title in (
            (axes[0], tg, '常偏移目标道集（含空腔 y14-18）\nB2D-C3m-D10m-W4m-T0.5m-E20-S0.02-CO11'),
            (axes[1], bg, '常偏移族 BG（无空腔）\nB2D-C3m-BG-CO11'),
            (axes[2], diff, '配对差分 TGT-BG（公共色标）\n空腔响应（约 320-336 ns）')):
        g = x / np.maximum(np.abs(x).max(axis=0, keepdims=True), 1e-30)
        im = ax.imshow(g, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax,
                       extent=(0.5, 11.5, t_ns[-1], t_ns[0]), interpolation='nearest')
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('道 t01..t11（Rx y 12.65-22.65 m，偏移距恒 1.30 m）')
    axes[0].set_ylabel('时间 (ns)')
    axes[0].set_ylim(600, 0)
    axes[2].axhline(320, color='k', lw=0.8, ls='--', alpha=0.6)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label='归一化幅度')
    fig.suptitle('batch2d_v1_co 常偏移 B-scan（11 道 x 1.0 m）：Tx=Rx-1.30 m，'
                 '镜像对称轴 y=16 m（模型整体对称 + 互易）', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT, dpi=130)
    print(OUT)


if __name__ == '__main__':
    main()
