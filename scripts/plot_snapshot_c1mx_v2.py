# -*- coding: utf-8 -*-
"""SNAP-C1mX-TGT 波场快照渲染 v2：空气/地下分区独立归一化。

物理事实：透入地下的场强约为空气侧的 0.1–0.4%（-50 dB 量级），
线性同尺度显示时地下区域看起来是空的。v2 对地表以下区域单独限幅，
并在每面板标注地下放大倍数。
"""
from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SNAP_DIR = Path(r"E:\automation_djh\automation_repo\artifacts\simulations\2026-10-01_SNAP-C1mX-TGT\SNAP-C1mX-TGT_snaps")
OUT = Path(r"E:\automation_djh\fig_snapshot_c1mx_wavefield_v2.png")
TIMES = [20, 50, 90, 130, 160, 190, 240, 320]  # ns
SURFACE_ROW = 600  # z=30 m / 0.05 m


def load_ex(ns):
    with h5py.File(SNAP_DIR / f"s{ns:03d}.vtkhdf", "r") as h:
        return h["VTKHDF"]["CellData"]["Ex"][:, :, 0]


def overlay(ax):
    ax.axhline(30, color="k", lw=1.2)
    ax.axhline(29, color="k", lw=0.8, ls="--", alpha=0.7)
    ax.add_patch(plt.Rectangle((14, 19.75), 4, 0.5, fill=False,
                               edgecolor="#d62728", lw=1.4))
    ax.plot(15.35, 45, marker="v", color="k", ms=6)
    ax.plot(16.65, 45, marker="^", color="k", ms=5, mfc="none")


def main():
    fig, axes = plt.subplots(2, 4, figsize=(19, 13))
    for ax, ns in zip(axes.flat, TIMES):
        ex = load_ex(ns)
        air, gnd = ex[SURFACE_ROW:, :], ex[:SURFACE_ROW, :]
        sa = float(np.percentile(np.abs(air), 99.5))
        sg = float(np.percentile(np.abs(gnd), 99.5))
        disp = np.empty_like(ex)
        disp[SURFACE_ROW:, :] = np.clip(air / (sa + 1e-30), -1, 1)
        if sg > 1e-12:
            disp[:SURFACE_ROW, :] = np.clip(gnd / sg, -1, 1)
        else:
            disp[:SURFACE_ROW, :] = 0.0
        im = ax.imshow(disp, cmap="RdBu_r", vmin=-1, vmax=1,
                       extent=[0, 32, 50, 0], aspect="auto",
                       interpolation="bilinear")
        overlay(ax)
        amp = sa / sg if sg > 1e-12 else None
        title = f"t = {ns} ns"
        if amp is not None:
            title += f"（地下 ×{amp:.0f}）"
        else:
            title += "（地下尚无透射）"
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("y（m）")
        ax.set_ylabel("z（m）")
        print(f"t={ns} ns  地下放大 {f'{amp:.0f}x' if amp is not None else '无透射'}")
    fig.suptitle("C1mX 含目标模型波场快照（Ex）：空气侧与地下分区独立限幅；红框=目标（4 m×0.5 m @ 深 10 m）",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT, dpi=150)
    print(f"saved: {OUT}")


if __name__ == "__main__":
    main()
