# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.5 —— 振铃不被地表编辑截断 + 包络剖面匹配 + 横向抖动。

v0.4 残留差距（diag_v04_t0rel.py 实测，t0 相对窗、带通 90-125MHz 包络）：
  G1 振铃带在 t0+100~140ns 被 OP2 汉宁窗混合顺带压掉（v0.4 env 0.0086 vs 实测 0.036）
     —— 物理上振铃激励与地表反射同时进入天线，编辑应先于振铃注入；
  G2 振铃幅值整体约为实测的 2/3（能量比口径里实测直达窗含长拖尾，标定被稀释）；
  G3 振铃包络横向起伏 diff_rms=0.057，实测 0.179 —— 缺逐道白抖动分量。

v0.5 改法：
  1) 算子顺序：OP2 地表编辑 -> OP1 振铃注入（振铃由编辑后道激励，不再被截断）；
  2) 匹配窗重定义：E_d = t0-6~t0+12（收紧峰窗，剔除拖尾稀释），
     E_r = t0+40~t0+140（延展到实测振铃真实覆盖区），目标比用全 Line9 2378 道重新标定；
  3) β 求解改网格扫描取小根（对 ratio(β) 非单调稳健）；
  4) 逐道幅度白抖动 σ=0.12 叠加在原缓变序列上（目标包络 diff_rms≈0.18）。

核验门：直达窗波形相关≥0.95；峰值变化<2%；到时零漂移；
       振铃包络剖面（t0 相对 6 窗）中位与实测偏差≤±25%。

输出：fig_augment_v05_bscan.png（原始 | v0.4 | v0.5 | 实测 Line9，同口径）
      fig_augment_v05_profile.png（振铃包络剖面三线对比）
      augment_v05_metrics.json（标定与核验指标留档）
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
from augment_sim2real_v0_4 import op1_ringing_matched, trnorm  # noqa: E402
from augment_sim2real_v0_4 import SEED as SEED04  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SEED = 20261001
JITTER_STD = 0.12          # 逐道幅度白抖动（目标包络 diff_rms≈0.18）
WD_LO, WD_HI = -6.0, 12.0  # E_d 窗（相对 t0，ns）
WR_LO, WR_HI = 40.0, 140.0  # E_r 窗（相对 t0，ns）


def t0_index(x, t_ns):
    w = (t_ns >= 0.0) & (t_ns <= 60.0)
    idx = np.where(w)[0]
    return idx[int(np.argmax(np.abs(x[w])))]


def win_mask(t_ns, t0, lo, hi):
    return (t_ns >= t0 + lo) & (t_ns <= t0 + hi)


def ring_ratio(y, t_ns, t0):
    wd = win_mask(t_ns, t0, WD_LO, WD_HI)
    wr = win_mask(t_ns, t0, WR_LO, WR_HI)
    e_r = float(np.sum(y[wr] ** 2))
    e_d = float(np.sum(y[wd] ** 2)) + 1e-30
    return e_r / e_d


def calibrate_targets(real, t_real):
    """全 Line9 道集标定 E_r/E_d 分布（新窗口径）。"""
    rs = []
    for i in range(real.shape[0]):
        i0 = t0_index(real[i], t_real)
        rs.append(ring_ratio(real[i], t_real, t_real[i0]))
    rs = np.array(rs)
    return {k: float(np.percentile(rs, q)) for k, q in
            [("p10", 10), ("p25", 25), ("p50", 50), ("p75", 75), ("p90", 90)]}


def op1_ringing_v05(sigs, t_ns, rng, r_lo, r_hi):
    """振铃注入：β 网格扫描取小根，逐道白抖动。"""
    n, nt = sigs.shape
    dt = float(t_ns[1] - t_ns[0])
    tg = np.arange(0, 5 * TAU, dt)
    g = np.exp(-tg / TAU) * np.sin(2 * np.pi * F0 * tg)
    rs = smooth_series(rng, n, r_lo, r_hi, 5.0)
    jitter = np.clip(1.0 + JITTER_STD * rng.standard_normal(n), 0.4, None)
    out = sigs.copy()
    beta_grid = np.concatenate([[0.0], np.logspace(-3, np.log10(2.0), 80)])
    n_clip = 0
    for i in range(n):
        i0 = t0_index(sigs[i], t_ns)
        t0 = t_ns[i0]
        c = np.convolve(sigs[i], g)[:nt] * dt
        target = rs[i]
        prev_b, prev_r = 0.0, ring_ratio(sigs[i], t_ns, t0)
        beta = None
        for b in beta_grid[1:]:
            r = ring_ratio(sigs[i] + b * c, t_ns, t0)
            if r >= target:
                beta = prev_b + (b - prev_b) * (target - prev_r) / max(r - prev_r, 1e-30)
                break
            prev_b, prev_r = b, r
        if beta is None:  # 目标超峰值：取扫描末端并计数
            beta = beta_grid[-1]
            n_clip += 1
        out[i] = sigs[i] + beta * jitter[i] * c
    return out, n_clip


def direct_fidelity(s0, s1, t_ns):
    """直达窗波形相关 / 峰值变化 / 到时漂移。"""
    cors, pkchg, shifts = [], [], []
    for i in range(s0.shape[0]):
        i0 = t0_index(s0[i], t_ns)
        t0 = t_ns[i0]
        wd = win_mask(t_ns, t0, -5.0, 25.0)
        a, b = s0[i][wd], s1[i][wd]
        cors.append(float(np.corrcoef(a, b)[0, 1]))
        pkchg.append(abs(np.abs(b).max() / (np.abs(a).max() + 1e-30) - 1.0))
        i1 = t0_index(s1[i], t_ns)
        shifts.append(abs(t_ns[i1] - t0))
    return (float(np.median(cors)), float(np.min(cors)),
            float(np.max(pkchg)), float(np.max(shifts)))


def env_profile(x, t_ns, rel_wins):
    """带通 90-125MHz 包络在 t0 相对窗的中位（全道中位）。"""
    from scipy.signal import hilbert, butter, sosfiltfilt
    dt = float(t_ns[1] - t_ns[0])
    sos = butter(4, [0.090, 0.125], btype="band", fs=1.0 / dt, output="sos")
    rows = []
    for i in range(x.shape[0]):
        t0 = t_ns[t0_index(x[i], t_ns)]
        env = np.abs(hilbert(sosfiltfilt(sos, x[i])))
        rows.append([float(np.median(env[win_mask(t_ns, t0, lo, hi)]))
                     for lo, hi in rel_wins])
    return np.median(np.array(rows), axis=0)


def lat_env_diffrms(x, t_ns):
    from scipy.signal import hilbert
    envs = []
    for i in range(x.shape[0]):
        t0 = t_ns[t0_index(x[i], t_ns)]
        m = win_mask(t_ns, t0, 40.0, 100.0)
        envs.append(np.abs(hilbert(x[i][m])).mean())
    e = np.array(envs) / np.mean(envs)
    return float(np.sqrt(np.mean(np.diff(e) ** 2)))


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0, t = sigs[:, crop], t_ns[crop]

    real = load_line9()
    t_real_full = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    calib = calibrate_targets(real, t_real_full)
    print("Line9 全道标定（新窗 E_r/E_d）:", {k: round(v, 4) for k, v in calib.items()})

    # v0.5：先 OP2 编辑，再注入振铃
    rng = np.random.default_rng(SEED)
    s_edit = op2_surface_edit(sigs, t_ns, rng)
    s_v05, n_clip = op1_ringing_v05(s_edit, t_ns, rng, calib["p25"], calib["p75"])
    s_v05 = s_v05[:, crop]
    print(f"β 扫描触顶道数: {n_clip}/{sigs.shape[0]}")

    # v0.4 对照
    rng4 = np.random.default_rng(SEED04)
    s_v04 = op2_surface_edit(op1_ringing_matched(sigs, t_ns, rng4), t_ns, rng4)[:, crop]

    mid = real.shape[0] // 2
    rcrop = t_real_full <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real_full[rcrop]

    s0n, v4n, v5n, r33n = trnorm(s0), trnorm(s_v04), trnorm(s_v05), trnorm(r33)

    # ---- 核验 ----
    cor_med, cor_min, pk_max, sh_max = direct_fidelity(s0, s_v05, t)
    print(f"直达保真: corr_med={cor_med:.4f} corr_min={cor_min:.4f} "
          f"peak_chg_max={pk_max * 100:.2f}% shift_max={sh_max:.2f}ns")

    rel_wins = [(25, 45), (45, 70), (70, 100), (100, 140), (140, 200)]
    prof_v4 = env_profile(v4n, t, rel_wins)
    prof_v5 = env_profile(v5n, t, rel_wins)
    prof_l9 = env_profile(r33n, tr, rel_wins)
    drift = [f"{a / max(b, 1e-9):.2f}x" for a, b in zip(prof_v5, prof_l9)]
    print("包络剖面 v5/实测 比:", drift)

    dr_v5 = lat_env_diffrms(v5n, t)
    dr_l9 = lat_env_diffrms(r33n, tr)
    print(f"横向抖动 diff_rms: v0.5={dr_v5:.3f} 实测={dr_l9:.3f}")

    metrics = {
        "seed": SEED, "mother": "B2D-C1mX-BG", "n_traces": int(sigs.shape[0]),
        "calib_line9_new_windows": calib,
        "windows_ns": {"wd": [WD_LO, WD_HI], "wr": [WR_LO, WR_HI]},
        "jitter_std": JITTER_STD, "beta_clip_count": n_clip,
        "direct_fidelity": {"corr_med": cor_med, "corr_min": cor_min,
                            "peak_chg_max": pk_max, "shift_max_ns": sh_max},
        "env_profile": {"wins": rel_wins, "v04": prof_v4.tolist(),
                        "v05": prof_v5.tolist(), "line9": prof_l9.tolist(),
                        "v05_over_line9": [float(a / max(b, 1e-9)) for a, b in zip(prof_v5, prof_l9)]},
        "lateral_diff_rms": {"v05": dr_v5, "line9": dr_l9},
    }
    with open(r"E:\automation_djh\augment_v05_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    # ---- 图1：四面板 B-scan ----
    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0n, t, clim(s0n, t), "仿真原始")
    panel(axes[1], v4n, t, clim(v4n, t), "v0.4")
    panel(axes[2], v5n, t, clim(v5n, t), "v0.5（顺序+新窗+抖动）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.5 对比 · C1mX-BG · 四面板同口径（逐道归一化，种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v05_bscan.png", dpi=150)
    plt.close(fig)

    # ---- 图2：振铃包络剖面 ----
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    xc = [f"t0+{lo}~{hi}" for lo, hi in rel_wins]
    xs = np.arange(len(xc))
    w = 0.26
    ax.bar(xs - w, prof_v4, width=w, label="v0.4", color="#9E9E9E")
    ax.bar(xs, prof_v5, width=w, label="v0.5", color="#C0392B")
    ax.bar(xs + w, prof_l9, width=w, label="实测 Line9", color="#1F6FEB")
    for x, v in zip(xs, prof_v5):
        ax.text(x, v * 1.05, f"{v:.3f}", ha="center", fontsize=8)
    for x, v in zip(xs + w, prof_l9):
        ax.text(x + 0.02, v * 1.05, f"{v:.3f}", ha="center", fontsize=8, color="#1F6FEB")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(xc)
    ax.set_ylabel("带通 90-125MHz 包络中位（归一化，log）")
    ax.set_title("振铃包络剖面：v0.5 vs v0.4 vs 实测（t0 相对窗）")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v05_profile.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
