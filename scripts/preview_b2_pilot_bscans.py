# -*- coding: utf-8 -*-
"""B2 试点批（2026-10-01）B-scan 验收图（零求解器，直接读归档 h5 重建）。

图 1：三个 A 层开发族 BG | TGT | TGT−BG 差分（差分应只显目标响应）。
图 2：B 层 4 个探针母模型 BG（几何/阶梯 sanity 检查）。
风格沿用 preview_bscan_operators.py：符号振幅灰度、横道号纵时间（向下增大）、
剔除 t<30 ns 后 99 分位自身限幅。
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(r"E:\automation_djh\automation_repo")
sys.path.insert(0, str(REPO / "scripts"))

from run_reward_protocol_b2_pilot_v0_1 import load_bscan  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

T_CROP_NS = 400.0
DAM = "D10m-W4m-T0.5m-E20-S0.02"

# (族标签, BG 母模型, BG 日期, TGT 母模型, TGT 日期)
FAMS = [
    ("C1mX 平地·覆土1m", "B2D-C1mX-BG", "2026-10-01", f"B2D-C1mX-{DAM}", "2026-10-01"),
    ("C3mX 平地·覆土3m", "B2D-C3mX-BG", "2026-09-28", f"B2D-C3mX-{DAM}", "2026-10-01"),
    ("S2X 坡地·无过渡带", "B2D-C3mS2X-BG", "2026-09-28", f"B2D-C3mS2X-{DAM}", "2026-10-01"),
]

PROBES = [
    ("C1p5mS1X 覆土1.5m·坡1", "B2D-C1p5mS1X-BG"),
    ("C1p5mS3X 覆土1.5m·坡3", "B2D-C1p5mS3X-BG"),
    ("C2mS3X 覆土2m·坡3", "B2D-C2mS3X-BG"),
    ("C2mS3TZX 覆土2m·坡3·TZ", "B2D-C2mS3TZX-BG"),
]

OUT1 = Path(r"E:\automation_djh\fig_b2_pilot_families.png")
OUT2 = Path(r"E:\automation_djh\fig_b2_pilot_probes.png")


def prep(mother, date):
    sigs, t_ns = load_bscan(mother, date)
    crop = t_ns <= T_CROP_NS
    return sigs[:, crop], t_ns[crop]


def panel(ax, sigs, t_ns, cmap, vlim, title=""):
    extent = [1, sigs.shape[0], t_ns[-1], t_ns[0]]
    ax.imshow(sigs.T, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim,
              extent=extent)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("道号")


def clim(sigs, t_ns, q=99.0):
    w = t_ns >= 30.0
    return float(np.percentile(np.abs(sigs[:, w]), q))


def main():
    # ---------- 图 1：三族 BG | TGT | 差分 ----------
    fig, axes = plt.subplots(3, 3, figsize=(17, 12))
    for r, (label, bg, bgd, tgt, tgtd) in enumerate(FAMS):
        s_bg, t_ns = prep(bg, bgd)
        s_tg, _ = prep(tgt, tgtd)
        assert s_bg.shape == s_tg.shape
        diff = s_tg - s_bg
        panel(axes[r, 0], s_bg, t_ns, "gray", clim(s_bg, t_ns),
              "BG（背景）" if r == 0 else "")
        panel(axes[r, 1], s_tg, t_ns, "gray", clim(s_tg, t_ns),
              "TGT（含目标：10m 处 4m×0.5m 损伤体）" if r == 0 else "")
        v = float(np.percentile(np.abs(diff), 99.9))
        panel(axes[r, 2], diff, t_ns, "RdBu_r", v,
              "TGT − BG 差分（目标响应）" if r == 0 else "")
        axes[r, 0].set_ylabel(f"{label}\n时间 (ns)")
        dpower = float(np.mean(diff ** 2) / np.mean(s_bg ** 2))
        print(f"{label}: BG/TGT shape={s_bg.shape} 差分功率比={dpower:.3e} "
              f"差分峰值={np.abs(diff).max():.4e} BG峰值={np.abs(s_bg).max():.4e}")
    fig.suptitle("B2 试点批验收 · 三个开发族：BG | TGT | TGT−BG（2026-10-01 归档，符号振幅灰度）",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT1, dpi=150)
    print(f"saved: {OUT1}")

    # ---------- 图 2：B 层 4 探针 BG ----------
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.2))
    for c, (label, mother) in enumerate(PROBES):
        s, t_ns = prep(mother, "2026-10-01")
        panel(axes[c], s, t_ns, "gray", clim(s, t_ns), label)
        axes[c].set_ylabel("时间 (ns)" if c == 0 else "")
        print(f"{label}: shape={s.shape} 峰值={np.abs(s).max():.4e}")
    fig.suptitle("B2 试点批验收 · B 层探针母模型 BG（2026-10-01 归档）", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT2, dpi=150)
    print(f"saved: {OUT2}")


if __name__ == "__main__":
    main()
