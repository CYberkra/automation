# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.4 —— 振铃能量按实测口径匹配，直达波保真优先。

v0.3 问题（analyze_direct_fidelity.py 实测）：
  - 直达窗波形相关仅 0.58（卷积项扭曲直达波）
  - E_ring/E_direct = 0.282，而实测 Line9 仅 0.029（p25-75 0.021-0.038），
    振铃能量高出约 10 倍——α 幅度标定没考虑仿真直达脉冲比仪器脉冲更紧凑。

v0.4 改法：β 不再按 α 幅度标定，改为**按能量比直接匹配实测分布**：
  per-trace: beta_i = sqrt( r_i * E_d / E_c )
    E_d = 直达窗(t0-5~t0+25ns)能量, E_c = (x⊛g) 在振铃窗(t0+40~t0+100ns)能量
    r_i = 目标能量比，沿道缓变，中位 0.029（范围 0.02~0.04，来自实测 p25-75）
  卷积形式保留（相位锁定）；OP2 地表编辑沿用。

输出：fig_augment_v04_bscan.png（原始 | v0.3 | v0.4 | 实测 Line9，同口径逐道归一化）
      fig_augment_v04_waveform.png
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
from augment_sim2real_v0_3 import (  # noqa: E402
    F0, TAU, K_LO, K_HI, JITTER_NS, T_CROP_NS, N_SAMPLES, T_WINDOW_NS,
    smooth_series, op1_ringing_conv, op2_surface_edit, load_line9, panel, clim,
)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SEED = 20261001
R_LO, R_HI = 0.020, 0.040   # 目标 E_ring/E_direct（实测 Line9 p25-p75 覆盖）


def op1_ringing_matched(sigs, t_ns, rng):
    """卷积振铃，β 按振铃/直达能量比匹配实测分布。"""
    n, nt = sigs.shape
    dt = float(t_ns[1] - t_ns[0])
    tg = np.arange(0, 5 * TAU, dt)
    g = np.exp(-tg / TAU) * np.sin(2 * np.pi * F0 * tg)
    rs = smooth_series(rng, n, R_LO, R_HI, 5.0)
    out = sigs.copy()
    w_direct = (t_ns >= 0.0) & (t_ns <= 60.0)
    idx = np.where(w_direct)[0]
    for i in range(n):
        i0 = idx[int(np.argmax(np.abs(sigs[i][w_direct])))]
        t0 = t_ns[i0]
        c = np.convolve(sigs[i], g)[:nt] * dt
        wd = (t_ns >= t0 - 5) & (t_ns <= t0 + 25)
        wr = (t_ns >= t0 + 40) & (t_ns <= t0 + 100)

        # 精确匹配：F(beta) = E_ring(x+beta*c)/E_direct(x+beta*c) - r = 0，二分求解
        def ratio(beta):
            y = sigs[i] + beta * c
            e_r = float(np.sum(y[wr] ** 2))
            e_d = float(np.sum(y[wd] ** 2)) + 1e-30
            return e_r / e_d

        r_hi_asym = float(np.sum(c[wr] ** 2)) / (float(np.sum(c[wd] ** 2)) + 1e-30)
        target = min(rs[i], 0.9 * r_hi_asym)
        lo_b, hi_b = 0.0, 1.0
        while ratio(hi_b) < target and hi_b < 1e6:
            hi_b *= 2.0
        for _ in range(40):
            mid_b = 0.5 * (lo_b + hi_b)
            if ratio(mid_b) < target:
                lo_b = mid_b
            else:
                hi_b = mid_b
        out[i] = sigs[i] + 0.5 * (lo_b + hi_b) * c
    return out


def trnorm(a):
    return a / (np.abs(a).max(axis=1, keepdims=True) + 1e-30)


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0 = sigs[:, crop]
    t = t_ns[crop]

    rng = np.random.default_rng(SEED)
    s_v03 = op2_surface_edit(op1_ringing_conv(sigs, t_ns, rng), t_ns, rng)[:, crop]
    rng = np.random.default_rng(SEED)
    s_v04 = op2_surface_edit(op1_ringing_matched(sigs, t_ns, rng), t_ns, rng)[:, crop]

    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    rcrop = t_real <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real[rcrop]

    s0n, sv3n, sv4n, r33n = trnorm(s0), trnorm(s_v03), trnorm(s_v04), trnorm(r33)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0n, t, clim(s0n, t), "仿真原始")
    panel(axes[1], sv3n, t, clim(sv3n, t), "v0.3（振铃能量过高）")
    panel(axes[2], sv4n, t, clim(sv4n, t), "v0.4（能量比匹配实测）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.4 对比 · C1mX-BG · 四面板同口径（逐道归一化，种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v04_bscan.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(t, sv4n[16], color="#C0392B", lw=1.0, label="仿真 + 增广 v0.4（C1mX-BG 道17）")
    ax.plot(tr, r33n[16], color="#1F6FEB", lw=1.0, alpha=0.8, label="实测 Line9 中点道")
    ax.set_xlabel("时间 (ns)")
    ax.set_ylabel("归一化振幅")
    ax.set_title("增广 v0.4 仿真道 vs 实测道（各自归一化）")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v04_waveform.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
