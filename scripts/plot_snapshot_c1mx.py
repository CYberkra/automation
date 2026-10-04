# -*- coding: utf-8 -*-
"""SNAP-C1mX-TGT 波场快照渲染（2026-10-01 诊断例）。

8 个时刻（20/50/90/130/160/190/240/320 ns）的 Ex 波场，
几何叠加：地表 z=30、覆土/基岩界面 z=29、目标框 y14-18 × z19.75-20.25、Tx/Rx。
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
OUT = Path(r"E:\automation_djh\fig_snapshot_c1mx_wavefield.png")
TIMES = [20, 50, 90, 130, 160, 190, 240, 320]  # ns


def load_ex(ns):
    with h5py.File(SNAP_DIR / f"s{ns:03d}.vtkhdf", "r") as h:
        ex = h["VTKHDF"]["CellData"]["Ex"][:, :, 0]  # (z, y)
        sp = h["VTKHDF"].attrs["Spacing"]
    return ex, sp


def overlay(ax):
    ax.axhline(30, color="k", lw=1.2)                       # 地表
    ax.axhline(29, color="k", lw=0.8, ls="--", alpha=0.7)   # 覆土/基岩界面
    # 目标框 y14-18, z19.75-20.25
    ax.add_patch(plt.Rectangle((14, 19.75), 4, 0.5, fill=False,
                               edgecolor="#d62728", lw=1.4))
    ax.plot(15.35, 45, marker="v", color="k", ms=6)         # Tx
    ax.plot(16.65, 45, marker="^", color="k", ms=5, mfc="none")  # Rx


def main():
    fig, axes = plt.subplots(2, 4, figsize=(19, 13))
    for ax, ns in zip(axes.flat, TIMES):
        ex, sp = load_ex(ns)
        v = float(np.percentile(np.abs(ex), 99.5))
        im = ax.imshow(ex, cmap="RdBu_r", vmin=-v, vmax=v,
                       extent=[0, 32, 50, 0], aspect="auto",
                       interpolation="bilinear")
        overlay(ax)
        ax.set_title(f"t = {ns} ns", fontsize=12)
        ax.set_xlabel("y（m）")
        ax.set_ylabel("z（m）")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Ex (V/m)")
        print(f"t={ns} ns  max|Ex|={np.abs(ex).max():.3e}  p99.5={v:.3e}")
    fig.suptitle("C1mX 含目标模型波场快照（Ex）：天线 15 m 架高，覆土 1 m εr=18，基岩 εr=9，目标 4 m×0.5 m @ 深 10 m（红框）",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT, dpi=150)
    print(f"saved: {OUT}")


if __name__ == "__main__":
    main()
