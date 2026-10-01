# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.6（修正版）—— 直接匹配带通包络/峰值比 + 起始门 + 相位抖动。

v0.6 初版问题：δ 时移把振铃首瓣与直达峰同相叠加（峰值+38%、到时漂 5.8ns），
且"能量比"口径与"眼睛看到的包络/峰值比"之间始终隔着子波长度差异。

修正版思路 —— 直接匹配视觉量本身：
  目标量 q = 带通(90-125MHz)包络在 t0+45~70ns 窗的中位 / 直达峰值幅度
  - 实测 Line9 全道标定 q 分布（p25~p75 作为目标采样区间）；
  - 每道解 β 使仿真道的 q 命中目标（包络对 β 单调，网格+插值）；
  - 振铃 c 加起始门（t0+12 前为 0，t0+12~22 Hann 上升）：
    物理上振铃是直达脉冲通过后的天线拖尾，且从结构上保证直达波零失真；
  - δ ∈ ±3.5ns 相位抖动（沿道缓变 corr≈3）+ σ=0.16 幅度白抖动（串珠纹理）。

核验门：直达窗相关≥0.95；峰值变化<2%；到时漂移=0；
       包络剖面（t0 相对 5 窗）中位与实测比 0.8~1.25；diff_rms≈0.18。

输出：fig_augment_v06_bscan.png / fig_augment_v06_profile.png / augment_v06_metrics.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(r"E:\automation_djh\automation_repo")
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, r"E:\automation_djh")

from run_reward_protocol_b2_pilot_v0_1 import load_bscan  # noqa: E402
from augment_sim2real_v0_3 import (  # noqa: E402
    F0, TAU, T_CROP_NS, N_SAMPLES, T_WINDOW_NS,
    smooth_series, op2_surface_edit, load_line9, panel, clim,
)
from augment_sim2real_v0_4 import trnorm  # noqa: E402
from augment_sim2real_v0_5 import (  # noqa: E402
    t0_index, win_mask, direct_fidelity, env_profile, lat_env_diffrms,
    op1_ringing_v05,
)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SEED = 20261001
JITTER_STD = 0.16
DELTA_NS = 3.5
BP_LO, BP_HI = 0.090, 0.125  # cycles/ns
QW_LO, QW_HI = 45.0, 70.0    # q 匹配窗（相对 t0）
ONSET0, ONSET1 = 12.0, 22.0  # 振铃起始门（相对 t0）


def _sos(dt):
    from scipy.signal import butter
    return butter(4, [BP_LO, BP_HI], btype="band", fs=1.0 / dt, output="sos")


def q_of_trace(x, t_ns, sos):
    """带通包络 t0+45~70 窗中位 / 直达峰值。"""
    from scipy.signal import hilbert, sosfiltfilt
    i0 = t0_index(x, t_ns)
    t0 = t_ns[i0]
    peak = abs(x[i0]) + 1e-30
    env = np.abs(hilbert(sosfiltfilt(sos, x)))
    m = win_mask(t_ns, t0, QW_LO, QW_HI)
    return float(np.median(env[m]) / peak)


def calibrate_q(real, t_real):
    sos = _sos(float(t_real[1] - t_real[0]))
    qs = [q_of_trace(real[i], t_real, sos) for i in range(real.shape[0])]
    qs = np.array(qs)
    return {k: float(np.percentile(qs, q)) for k, q in
            [("p10", 10), ("p25", 25), ("p50", 50), ("p75", 75), ("p90", 90)]}


def onset_gate(t_ns, t0):
    g = np.ones_like(t_ns)
    trel = t_ns - t0
    g[trel < ONSET0] = 0.0
    ramp = (trel >= ONSET0) & (trel < ONSET1)
    u = (trel[ramp] - ONSET0) / (ONSET1 - ONSET0)
    g[ramp] = 0.5 - 0.5 * np.cos(np.pi * u)
    return g


def op1_ringing_v06(sigs, t_ns, rng, q_lo, q_hi):
    """振铃注入：起始门 + δ 相位抖动 + β 按 q 命中 + 幅度白抖动。"""
    from scipy.signal import hilbert, sosfiltfilt
    n, nt = sigs.shape
    dt = float(t_ns[1] - t_ns[0])
    tg = np.arange(0, 5 * TAU, dt)
    g = np.exp(-tg / TAU) * np.sin(2 * np.pi * F0 * tg)
    sos = _sos(dt)
    qs = smooth_series(rng, n, q_lo, q_hi, 5.0)
    deltas = smooth_series(rng, n, -DELTA_NS, DELTA_NS, 3.0)
    jitter = np.clip(1.0 + JITTER_STD * rng.standard_normal(n), 0.4, None)
    beta_grid = np.concatenate([[0.0], np.logspace(-3, np.log10(3.0), 96)])
    out = sigs.copy()
    n_clip = 0
    for i in range(n):
        i0 = t0_index(sigs[i], t_ns)
        t0 = t_ns[i0]
        peak = abs(sigs[i][i0]) + 1e-30
        c = np.convolve(sigs[i], g)[:nt] * dt
        c = np.interp(t_ns - deltas[i], t_ns, c, left=0.0, right=0.0)
        c = c * onset_gate(t_ns, t0)
        Ya = hilbert(sosfiltfilt(sos, sigs[i]))
        Ca = hilbert(sosfiltfilt(sos, c))
        m = win_mask(t_ns, t0, QW_LO, QW_HI)
        env_med = np.median(np.abs(Ya[m][None, :] + beta_grid[:, None] * Ca[m][None, :]), axis=1)
        q_curve = env_med / peak
        target = qs[i]
        hit = np.where(q_curve >= target)[0]
        if len(hit) == 0:
            beta = beta_grid[-1]
            n_clip += 1
        else:
            j = hit[0]
            if j == 0:
                beta = 0.0
            else:
                b0, b1 = beta_grid[j - 1], beta_grid[j]
                q0, q1 = q_curve[j - 1], q_curve[j]
                beta = b0 + (b1 - b0) * (target - q0) / max(q1 - q0, 1e-30)
        out[i] = sigs[i] + beta * jitter[i] * c
    return out, n_clip


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0, t = sigs[:, crop], t_ns[crop]

    real = load_line9()
    t_real_full = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    calib = calibrate_q(real, t_real_full)
    print("Line9 全道标定 q（带通包络/峰值）:", {k: round(v, 4) for k, v in calib.items()})

    rng = np.random.default_rng(SEED)
    s_edit = op2_surface_edit(sigs, t_ns, rng)
    s_v06, n_clip = op1_ringing_v06(s_edit, t_ns, rng, calib["p25"], calib["p75"])
    s_v06 = s_v06[:, crop]
    print(f"β 扫描触顶道数: {n_clip}/{sigs.shape[0]}")

    rng5 = np.random.default_rng(SEED)
    s_edit5 = op2_surface_edit(sigs, t_ns, rng5)
    s_v05, _ = op1_ringing_v05(s_edit5, t_ns, rng5, 0.0295, 0.0464)
    s_v05 = s_v05[:, crop]

    mid = real.shape[0] // 2
    rcrop = t_real_full <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real_full[rcrop]

    s0n, v5n, v6n, r33n = trnorm(s0), trnorm(s_v05), trnorm(s_v06), trnorm(r33)

    cor_med, cor_min, pk_max, sh_max = direct_fidelity(s0, s_v06, t)
    print(f"直达保真: corr_med={cor_med:.4f} corr_min={cor_min:.4f} "
          f"peak_chg_max={pk_max * 100:.2f}% shift_max={sh_max:.2f}ns")

    rel_wins = [(25, 45), (45, 70), (70, 100), (100, 140), (140, 200)]
    prof_v5 = env_profile(v5n, t, rel_wins)
    prof_v6 = env_profile(v6n, t, rel_wins)
    prof_l9 = env_profile(r33n, tr, rel_wins)
    drift = [f"{a / max(b, 1e-9):.2f}x" for a, b in zip(prof_v6, prof_l9)]
    print("包络剖面 v6/实测 比:", drift)

    dr_v6 = lat_env_diffrms(v6n, t)
    dr_l9 = lat_env_diffrms(r33n, tr)
    print(f"横向抖动 diff_rms: v0.6={dr_v6:.3f} 实测={dr_l9:.3f}")

    metrics = {
        "seed": SEED, "mother": "B2D-C1mX-BG", "n_traces": int(sigs.shape[0]),
        "method": "bandpass envelope/peak ratio q match + onset gate + phase jitter",
        "q_window_ns": [QW_LO, QW_HI], "onset_gate_ns": [ONSET0, ONSET1],
        "calib_line9_q": calib,
        "bandpass_mhz": [BP_LO * 1000, BP_HI * 1000],
        "jitter_std": JITTER_STD, "delta_ns": DELTA_NS, "beta_clip_count": n_clip,
        "direct_fidelity": {"corr_med": cor_med, "corr_min": cor_min,
                            "peak_chg_max": pk_max, "shift_max_ns": sh_max},
        "env_profile": {"wins": rel_wins, "v05": prof_v5.tolist(),
                        "v06": prof_v6.tolist(), "line9": prof_l9.tolist(),
                        "v06_over_line9": [float(a / max(b, 1e-9)) for a, b in zip(prof_v6, prof_l9)]},
        "lateral_diff_rms": {"v06": dr_v6, "line9": dr_l9},
    }
    with open(r"E:\automation_djh\augment_v06_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0n, t, clim(s0n, t), "仿真原始")
    panel(axes[1], v5n, t, clim(v5n, t), "v0.5")
    panel(axes[2], v6n, t, clim(v6n, t), "v0.6（包络匹配+相位抖动）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.6 对比 · C1mX-BG · 四面板同口径（逐道归一化，种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v06_bscan.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    xc = [f"t0+{lo}~{hi}" for lo, hi in rel_wins]
    xs = np.arange(len(xc))
    w = 0.26
    ax.bar(xs - w, prof_v5, width=w, label="v0.5", color="#9E9E9E")
    ax.bar(xs, prof_v6, width=w, label="v0.6", color="#C0392B")
    ax.bar(xs + w, prof_l9, width=w, label="实测 Line9", color="#1F6FEB")
    for x, v in zip(xs, prof_v6):
        ax.text(x, v * 1.05, f"{v:.3f}", ha="center", fontsize=8)
    for x, v in zip(xs + w, prof_l9):
        ax.text(x + 0.02, v * 1.05, f"{v:.3f}", ha="center", fontsize=8, color="#1F6FEB")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(xc)
    ax.set_ylabel("带通 90-125MHz 包络中位（归一化，log）")
    ax.set_title("振铃包络剖面：v0.6 vs v0.5 vs 实测（t0 相对窗）")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v06_profile.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
