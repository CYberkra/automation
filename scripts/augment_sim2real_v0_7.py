# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.7 —— v0.6 算子链原样，输入迁移到官方 20 MHz 载波。

修复的问题（独立审查 2026-10-03，见
automation_repo/docs/research/2026-10-03_sfcw_chain_review_and_repairs.md）：
v0_3–v0_6 的仿真输入来自 run_reward_protocol_b2_pilot_v0_1.load_bscan，
该函数用 |env|*cos(2π·95MHz·t+∠env) 重建带符号道（legacy 载波 bug），
物理频带被移位到 95–245 MHz。后果：v0.6 的 90–125 MHz 标定带对 legacy
仿真道落在移位频带下沿、对实测数据（真 20–170 MHz）落在中带，q 标定
跨不匹配频带。

v0.7 改动（仅输入表示，算子逐行不动）：
  - 仿真道改经 sfcw_official_loader_v0_2.load_bscan_both（官方
    reconstruct_time_response：载波 frequencies[0]=20 MHz，
    real_bandpass=2·Re(包络·载波)），与实测 Line9 同处真 20–170 MHz 带；
  - 同一加载同时给出 legacy 表示，v0.6 路径原样重跑作为连续性对照；
  - 标定带 90–125 MHz、振铃 F0=107 MHz、其余参数与 v0.6 完全相同，
    种子 20261001 不变。

注意幅度口径：官方 real_bandpass 含因子 2，legacy 为 1×；本脚本全部比较
为逐道归一化/比值，因子可消，不做跨链绝对幅值比较。

核验门（与 v0.6 相同）：直达窗相关≥0.95；峰值变化<2%；到时漂移=0；
包络剖面（t0 相对 5 窗）中位与实测比 0.8~1.25；diff_rms≈0.18。

输出：fig_augment_v07_bscan.png / fig_augment_v07_profile.png / augment_v07_metrics.json
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

from sfcw_official_loader_v0_2 import load_bscan_both  # noqa: E402
from augment_sim2real_v0_3 import (  # noqa: E402
    T_CROP_NS, N_SAMPLES, T_WINDOW_NS,
    op2_surface_edit, load_line9, panel, clim,
)
from augment_sim2real_v0_4 import trnorm  # noqa: E402
from augment_sim2real_v0_5 import (  # noqa: E402
    direct_fidelity, env_profile, lat_env_diffrms,
)
from augment_sim2real_v0_6 import (  # noqa: E402
    SEED, QW_LO, QW_HI, ONSET0, ONSET1, JITTER_STD, DELTA_NS, BP_LO, BP_HI,
    calibrate_q, op1_ringing_v06,
)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main():
    t_ns, sigs_off, sigs_leg = load_bscan_both("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0, t = sigs_off[:, crop], t_ns[crop]

    real = load_line9()
    t_real_full = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    calib = calibrate_q(real, t_real_full)
    print("Line9 全道标定 q（带通包络/峰值）:", {k: round(v, 4) for k, v in calib.items()})

    # v0.7 主路径：官方 20 MHz 载波输入
    rng = np.random.default_rng(SEED)
    s_edit = op2_surface_edit(sigs_off, t_ns, rng)
    s_v07, n_clip = op1_ringing_v06(s_edit, t_ns, rng, calib["p25"], calib["p75"])
    s_v07 = s_v07[:, crop]
    print(f"β 扫描触顶道数: {n_clip}/{sigs_off.shape[0]}")

    # 连续性对照：v0.6 路径在 legacy 95 MHz 载波输入上原样重跑
    rng6 = np.random.default_rng(SEED)
    s_edit6 = op2_surface_edit(sigs_leg, t_ns, rng6)
    s_v06, n_clip6 = op1_ringing_v06(s_edit6, t_ns, rng6, calib["p25"], calib["p75"])
    s_v06 = s_v06[:, crop]
    print(f"[对照 v0.6 legacy] β 触顶道数: {n_clip6}/{sigs_leg.shape[0]}")

    mid = real.shape[0] // 2
    rcrop = t_real_full <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real_full[rcrop]

    s0n, v6n, v7n, r33n = trnorm(s0), trnorm(s_v06), trnorm(s_v07), trnorm(r33)

    cor_med, cor_min, pk_max, sh_max = direct_fidelity(s0, s_v07, t)
    print(f"直达保真: corr_med={cor_med:.4f} corr_min={cor_min:.4f} "
          f"peak_chg_max={pk_max * 100:.2f}% shift_max={sh_max:.2f}ns")

    rel_wins = [(25, 45), (45, 70), (70, 100), (100, 140), (140, 200)]
    prof_v6 = env_profile(v6n, t, rel_wins)
    prof_v7 = env_profile(v7n, t, rel_wins)
    prof_l9 = env_profile(r33n, tr, rel_wins)
    drift = [f"{a / max(b, 1e-9):.2f}x" for a, b in zip(prof_v7, prof_l9)]
    print("包络剖面 v7/实测 比:", drift)

    dr_v7 = lat_env_diffrms(v7n, t)
    dr_l9 = lat_env_diffrms(r33n, tr)
    print(f"横向抖动 diff_rms: v0.7={dr_v7:.3f} 实测={dr_l9:.3f}")

    metrics = {
        "seed": SEED, "mother": "B2D-C1mX-BG", "n_traces": int(sigs_off.shape[0]),
        "method": "v0.6 operator chain on OFFICIAL 20 MHz carrier input "
                  "(sfcw_official_loader_v0_2); legacy 95 MHz rerun as continuity panel",
        "carrier": {"v07": "official frequencies[0]=20 MHz, real_bandpass=2*Re",
                    "v06_rerun": "legacy 95 MHz (bug convention, comparison only)"},
        "q_window_ns": [QW_LO, QW_HI], "onset_gate_ns": [ONSET0, ONSET1],
        "calib_line9_q": calib,
        "bandpass_mhz": [BP_LO * 1000, BP_HI * 1000],
        "jitter_std": JITTER_STD, "delta_ns": DELTA_NS,
        "beta_clip_count": n_clip, "beta_clip_count_v06_rerun": n_clip6,
        "direct_fidelity": {"corr_med": cor_med, "corr_min": cor_min,
                            "peak_chg_max": pk_max, "shift_max_ns": sh_max},
        "env_profile": {"wins": rel_wins, "v06_rerun": prof_v6.tolist(),
                        "v07": prof_v7.tolist(), "line9": prof_l9.tolist(),
                        "v07_over_line9": [float(a / max(b, 1e-9)) for a, b in zip(prof_v7, prof_l9)]},
        "lateral_diff_rms": {"v07": dr_v7, "line9": dr_l9},
        "supersedes": "augment_v06_metrics.json 的仿真输入载波口径（legacy 95 MHz）；"
                      "v0.6 指标保留为历史记录，不与 v0.7 直接比绝对值",
    }
    with open(r"E:\automation_djh\augment_v07_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0n, t, clim(s0n, t), "仿真原始（官方载波）")
    panel(axes[1], v6n, t, clim(v6n, t), "v0.6 重跑（legacy 载波对照）")
    panel(axes[2], v7n, t, clim(v7n, t), "v0.7（官方载波输入）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.7 对比 · C1mX-BG · 四面板同口径（逐道归一化，种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v07_bscan.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    xc = [f"t0+{lo}~{hi}" for lo, hi in rel_wins]
    xs = np.arange(len(xc))
    w = 0.26
    ax.bar(xs - w, prof_v6, width=w, label="v0.6 重跑（legacy）", color="#9E9E9E")
    ax.bar(xs, prof_v7, width=w, label="v0.7（官方载波）", color="#C0392B")
    ax.bar(xs + w, prof_l9, width=w, label="实测 Line9", color="#1F6FEB")
    for x, v in zip(xs, prof_v7):
        ax.text(x, v * 1.05, f"{v:.3f}", ha="center", fontsize=8)
    for x, v in zip(xs + w, prof_l9):
        ax.text(x + 0.02, v * 1.05, f"{v:.3f}", ha="center", fontsize=8, color="#1F6FEB")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(xc)
    ax.set_ylabel("带通 90-125MHz 包络中位（归一化，log）")
    ax.set_title("振铃包络剖面：v0.7 vs v0.6(legacy) vs 实测（t0 相对窗）")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v07_profile.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
