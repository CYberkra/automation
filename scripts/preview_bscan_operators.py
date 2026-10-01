# -*- coding: utf-8 -*-
"""既有归档 BG B-scan + 协议赢家算子 预览图（零求解器、零新仿真）。

行: C3mX(平地)->B7_G1_BG / S2X(T2坡无TZ)->B9_G1_BG / S2TZX(T2坡有TZ)->B9_G1_BG
列: 原始BG(dB) | 赢家算子处理后(dB, 同一色标) | 被抑制背景(原始-处理) | 中间道A-scan对比
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(r"E:\automation_djh\automation_repo")
sys.path.insert(0, str(REPO / "scripts"))

from study_t3_damage_ladder import load_bscan  # noqa: E402
from research_operator_contract import apply_configuration  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROWS = [
    ("B2D-C3mX-BG", "B7_G1_BG", "C3mX 平地 → B7（RPCA λ=0.5）"),
    ("B2D-C3mS2X-BG", "B9_G1_BG", "S2X T2坡·无过渡带 → B9（RPCA λ=2.0）"),
    ("B2D-C3mS2TZX-BG", "B9_G1_BG", "S2TZX T2坡·有过渡带 → B9（RPCA λ=2.0）"),
]
DB_FLOOR = -60.0
T_CROP_NS = 400.0  # 能量集中在 0-200 ns，裁剪显示窗口
OUT = Path(r"E:\automation_djh\fig_bscan_preview.png")


def to_db(sig, ref):
    return 20.0 * np.log10(np.abs(sig) / ref + 1e-12)


def main():
    fig, axes = plt.subplots(3, 4, figsize=(22, 13))
    for r, (mother, cfg, label) in enumerate(ROWS):
        sigs, t_ns = load_bscan(mother)          # (33, n_samples)
        x = np.ascontiguousarray(sigs.T)         # -> [sample, trace]
        res = apply_configuration(x, cfg, catalogue_version="0.2")
        out = res["output"].T                    # back to (33, n_samples)

        crop = t_ns <= T_CROP_NS                 # 只显示直达波/浅层窗口
        sigs, out, t_ns = sigs[:, crop], out[:, crop], t_ns[crop]

        ref = float(np.abs(sigs).max())
        removed = sigs - out                     # 被抑制的背景成分
        vlim = float(np.percentile(np.abs(removed), 99.9))
        # 各面板按自身 99 分位限幅（剔除 t<30 ns 直达耦合尖峰，避免吃掉色标动态）
        # 地震变密度惯例：零=中灰、正=白、负=黑
        w = t_ns >= 30.0
        alim_raw = float(np.percentile(np.abs(sigs[:, w]), 99.0))
        alim_out = float(np.percentile(np.abs(out[:, w]), 99.0))

        # imshow 行=纵轴：转置为 (samples, traces)，纵轴时间向下增大，直达波呈水平带
        extent = [1, sigs.shape[0], t_ns[-1], t_ns[0]]
        ax = axes[r, 0]
        ax.imshow(sigs.T, aspect="auto", cmap="gray", vmin=-alim_raw,
                  vmax=alim_raw, extent=extent)
        ax.set_ylabel(f"{label}\n时间 (ns)")
        ax.set_title("① 原始 BG B-scan（符号振幅，自身限幅）" if r == 0 else "")

        ax = axes[r, 1]
        im = ax.imshow(out.T, aspect="auto", cmap="gray", vmin=-alim_out,
                       vmax=alim_out, extent=extent)
        ax.set_title("② 赢家算子处理后（符号振幅，自身限幅）" if r == 0 else "")
        if r == 0:
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="振幅")

        ax = axes[r, 2]
        im = ax.imshow(removed.T, aspect="auto", cmap="RdBu_r", vmin=-vlim,
                       vmax=vlim, extent=extent)
        ax.set_title("③ 被抑制的背景（原始 − 处理后）" if r == 0 else "")
        if r == 0:
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

        ax = axes[r, 3]
        mid = sigs.shape[0] // 2
        ax.plot(t_ns, sigs[mid] / ref, color="0.45", lw=0.9,
                label="原始（第17道）")
        ax.plot(t_ns, out[mid] / ref, color="#d62728", lw=0.9,
                label="处理后")
        ax.set_xlim(t_ns[0], t_ns[-1])
        ax.set_ylim(-0.35, 0.35)
        ax.grid(alpha=0.3)
        ax.set_title("④ 中间道 A-scan 对比" if r == 0 else "")
        if r == 0:
            ax.legend(loc="upper right", fontsize=9)
        ax.set_xlabel("时间 (ns)")

        for c in range(3):
            axes[r, c].set_xlabel("道号")

        # 控制台数值摘要
        nb = float(np.mean(sigs ** 2))
        ob = float(np.mean(out ** 2))
        rb = float(np.mean(removed ** 2))
        print(f"{mother} {cfg}: traces={sigs.shape[0]} samples={sigs.shape[1]} "
              f"t=[{t_ns[0]:.1f},{t_ns[-1]:.1f}]ns 功率比(处理/原始)={ob / nb:.6f} "
              f"被抑制能量占比={rb / nb:.6f}")

    fig.suptitle("既有归档数据预览：三个开发母模型 BG B-scan × 协议赢家算子（2026-09-28 CO33 道集）",
                 fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT, dpi=150)
    print(f"saved: {OUT}")


if __name__ == "__main__":
    main()
