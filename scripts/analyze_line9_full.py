# -*- coding: utf-8 -*-
"""实测 Line9 全线 + 多线纹理结构分析（只读）。

回答三个问题：
  Q1 振铃/条带在全线上有多连贯？（邻道相关沿道号分布）
  Q2 深部(>150ns)的纹理到底是什么？（功率包络、是否有条带基底）
  Q3 直达波强度/振铃强度沿线的起伏有多大？

图：
  fig_line9_full_texture.png  上：全线 2378 道 B-scan(0-400ns)
                              中：三段 200 道局部放大
                              下：功率包络(dB) + 邻道相关系数
  fig_lines_compare.png       六线各取 200 道中段的纹理对照
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

DIRS = [Path(r"E:\automation_djh\real_data_yingshan\营山测线数据"),
        Path(r"E:\automation_djh\temporary product")]
N_SAMPLES, T_WINDOW_NS = 501, 700.0
T_NS = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
LINES = ["Line3", "Line6", "Line7", "Line9", "LineL1", "LineX1"]


def load(name):
    for d in DIRS:
        p = d / f"{name}origin(36).csv"
        if p.exists():
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                hdr = [f.readline() for _ in range(4)]
            n_traces = int(hdr[2].split("=")[1].strip().rstrip(","))
            raw = np.loadtxt(p, delimiter=",", skiprows=4, dtype=np.float64)
            return raw[:, 3].reshape(n_traces, N_SAMPLES)
    raise FileNotFoundError(name)


def clim(a, q=99.0):
    w = T_NS[: a.shape[1]] >= 30
    return float(np.percentile(np.abs(a[:, w]), q))


def main():
    L9 = load("Line9")
    n = L9.shape[0]
    crop = T_NS <= 400
    t = T_NS[crop]

    # ---- 统计量 ----
    early = T_NS <= 60
    i0 = np.argmax(np.abs(L9[:, early]), axis=1)
    peak = np.abs(L9)[np.arange(n), i0]
    # 邻道相关：振铃窗 30-100 ns
    w_ring = (T_NS >= 30) & (T_NS <= 100)
    seg = L9[:, w_ring]
    seg = seg - seg.mean(axis=1, keepdims=True)
    num = np.sum(seg[:-1] * seg[1:], axis=1)
    den = np.sqrt(np.sum(seg[:-1] ** 2, axis=1) * np.sum(seg[1:] ** 2, axis=1)) + 1e-30
    cc_ring = num / den
    # 深部窗 200-400 ns 邻道相关
    w_deep = (T_NS >= 200) & (T_NS <= 400)
    segd = L9[:, w_deep]
    segd = segd - segd.mean(axis=1, keepdims=True)
    numd = np.sum(segd[:-1] * segd[1:], axis=1)
    dend = np.sqrt(np.sum(segd[:-1] ** 2, axis=1) * np.sum(segd[1:] ** 2, axis=1)) + 1e-30
    cc_deep = numd / dend
    # 平均功率包络（dB，全体道）
    p = 10 * np.log10(np.mean(L9 ** 2, axis=0) + 1e-30)

    print(f"Line9 n={n}")
    print(f"直达峰值幅度: 中位 {np.median(peak):.4g}  p10-p90 {np.percentile(peak,10):.3g}-{np.percentile(peak,90):.3g}")
    print(f"振铃窗(30-100ns)邻道相关: 中位 {np.median(cc_ring):.3f}  p25-p75 {np.percentile(cc_ring,25):.2f}-{np.percentile(cc_ring,75):.2f}")
    print(f"深部窗(200-400ns)邻道相关: 中位 {np.median(cc_deep):.3f}")
    for tt in (50, 100, 150, 200, 300, 400, 600):
        k = int(np.searchsorted(T_NS, tt))
        print(f"  功率包络 {tt:>3d} ns: {p[k]:.1f} dB (相对)")

    # ---- 图 1 ----
    fig = plt.figure(figsize=(16, 13))
    ax0 = fig.add_subplot(3, 1, 1)
    v = clim(L9[:, crop])
    ax0.imshow(L9[:, crop].T, aspect="auto", cmap="gray", vmin=-v, vmax=v,
               extent=[1, n, t[-1], t[0]])
    ax0.set_title("Line9 全线 2378 道（0–400 ns，99 分位限幅）", fontsize=12, fontweight="bold")
    ax0.set_ylabel("时间 (ns)")

    starts = [100, 1100, 2100]
    for j, st in enumerate(starts):
        ax = fig.add_subplot(3, 3, 4 + j)
        w200 = L9[st:st + 200][:, crop]
        ax.imshow(w200.T, aspect="auto", cmap="gray", vmin=-v, vmax=v,
                  extent=[st + 1, st + 200, t[-1], t[0]])
        ax.set_title(f"道 {st+1}–{st+200}", fontsize=10)
        if j == 0:
            ax.set_ylabel("时间 (ns)")

    axp = fig.add_subplot(3, 2, 5)
    axp.plot(T_NS, p - p.max(), color="#1F6FEB", lw=1.2)
    axp.set_xlim(0, 700)
    axp.set_xlabel("时间 (ns)")
    axp.set_ylabel("平均功率 (dB, 归一)")
    axp.set_title("全线平均功率包络", fontsize=10)
    axp.grid(alpha=0.3)

    axc = fig.add_subplot(3, 2, 6)
    axc.plot(np.arange(2, n + 1), cc_ring, color="#C0392B", lw=0.5, label="振铃窗 30–100 ns")
    axc.plot(np.arange(2, n + 1), cc_deep, color="#1F6FEB", lw=0.5, label="深部窗 200–400 ns")
    axc.set_ylim(-0.2, 1.02)
    axc.set_xlabel("道号")
    axc.set_ylabel("邻道相关系数")
    axc.set_title("纹理的道间相干性", fontsize=10)
    axc.legend(fontsize=9)
    axc.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_line9_full_texture.png", dpi=150)
    plt.close(fig)

    # ---- 图 2：六线对照 ----
    fig, axes = plt.subplots(2, 3, figsize=(17, 8), sharey=True)
    for ax, name in zip(axes.flat, LINES):
        a = load(name)
        m = a.shape[0] // 2
        w200 = a[m - 100: m + 100][:, crop]
        ax.imshow(w200.T, aspect="auto", cmap="gray", vmin=-clim(w200), vmax=clim(w200),
                  extent=[m - 99, m + 100, t[-1], t[0]])
        ax.set_title(f"{name}（中段 200 道）", fontsize=10)
    axes[0, 0].set_ylabel("时间 (ns)")
    axes[1, 0].set_ylabel("时间 (ns)")
    fig.suptitle("六条测线中段纹理对照（各自 99 分位限幅）", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(r"E:\automation_djh\fig_lines_compare.png", dpi=150)
    plt.close(fig)
    print("figs saved")


if __name__ == "__main__":
    main()
