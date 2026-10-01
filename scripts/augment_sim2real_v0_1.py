# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.1（P2-1 振铃 + P2-2 地表反射压削/抖动）与效果预览。

算子定义（全部后处理、零求解器）：
  OP1 振铃注入（P2-1，参数全实测锚定）：
      ring(t) = A * exp(-t'/73ns) * sin(2*pi*107MHz*t' + phi), t' = t - t_direct
      A = alpha * 直达波峰值, alpha ~ U(0.04, 0.19), phi ~ U(0, 2pi)
  OP2 地表反射编辑（P2-2）：
      在 70–180 ns 窗内找地表反射峰 t_s，取 [t_s-30, t_s+90] 段，
      幅度乘 k ~ U(0.2, 0.64)（0.64×0.14 ≈ 0.09 实测上限），
      段内循环时移 delta ~ U(-15, 15) ns（架高抖动等效）。

预览输出：
  fig_augment_bscan.png     C1mX-BG：原始 | +振铃 | +振铃&压削抖动 | 实测 Line9 对照
  fig_augment_waveform.png  仿真增广道 vs 实测道（各自归一化）
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

SEED = 20261001
F0 = 0.107          # cycles/ns
TAU = 73.0          # ns
ALPHA_LO, ALPHA_HI = 0.04, 0.19
K_LO, K_HI = 0.20, 0.64
JITTER_NS = 15.0
T_CROP_NS = 400.0

LINE9 = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
N_SAMPLES, T_WINDOW_NS = 501, 700.0


def op1_ringing(sigs, t_ns, rng):
    """逐道注入振铃，返回增广后信号。"""
    out = sigs.copy()
    w = (t_ns >= 0.0) & (t_ns <= 60.0)
    for i in range(sigs.shape[0]):
        i0 = np.where(w)[0][int(np.argmax(np.abs(sigs[i][w])))]
        t0 = t_ns[i0]
        alpha = rng.uniform(ALPHA_LO, ALPHA_HI)
        phi = rng.uniform(0, 2 * np.pi)
        amp = alpha * float(np.abs(sigs[i, i0]))
        tt = t_ns - t0
        out[i] += np.where(tt >= 0,
                           amp * np.exp(-tt / TAU) * np.sin(2 * np.pi * F0 * tt + phi),
                           0.0)
    return out


def op2_surface_edit(sigs, t_ns, rng):
    """逐道压削 + 抖动地表反射段。"""
    out = sigs.copy()
    dt = float(t_ns[1] - t_ns[0])
    w = (t_ns >= 70.0) & (t_ns <= 180.0)
    idx = np.where(w)[0]
    for i in range(sigs.shape[0]):
        is_ = idx[int(np.argmax(np.abs(sigs[i][w])))]
        t_s = t_ns[is_]
        lo = int(np.searchsorted(t_ns, t_s - 30.0))
        hi = int(np.searchsorted(t_ns, t_s + 90.0))
        hi = min(hi, sigs.shape[1])
        if hi - lo < 8:
            continue
        seg = out[i, lo:hi].copy()
        k = rng.uniform(K_LO, K_HI)
        shift = int(round(rng.uniform(-JITTER_NS, JITTER_NS) / dt))
        seg = np.roll(seg, shift)
        if shift > 0:
            seg[:shift] = 0.0
        elif shift < 0:
            seg[shift:] = 0.0
        out[i, lo:hi] = seg * k
    return out


def load_line9():
    with open(LINE9, "r", encoding="utf-8", errors="replace") as f:
        hdr = [f.readline() for _ in range(4)]
    n_traces = int(hdr[2].split("=")[1].strip().rstrip(","))
    raw = np.loadtxt(LINE9, delimiter=",", skiprows=4, dtype=np.float64)
    return raw[:, 3].reshape(n_traces, N_SAMPLES)


def panel(ax, sigs, t_ns, vlim, title=""):
    extent = [1, sigs.shape[0], t_ns[-1], t_ns[0]]
    ax.imshow(sigs.T, aspect="auto", cmap="gray", vmin=-vlim, vmax=vlim, extent=extent)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("道号")


def clim(sigs, t_ns, q=99.0):
    return float(np.percentile(np.abs(sigs[:, t_ns >= 30.0]), q))


def main():
    rng = np.random.default_rng(SEED)
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0 = sigs[:, crop]
    t = t_ns[crop]
    s1 = op1_ringing(s0, t, rng)            # 仅 P2-1
    s2 = op2_surface_edit(s1, t, rng)       # P2-1 + P2-2

    # 实测 Line9：取中间 33 道，逐道归一化后与仿真同档显示
    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    r33 = real[mid - 16: mid + 17]
    rcrop = t_real <= T_CROP_NS
    r33 = r33[:, rcrop]
    tr = t_real[rcrop]
    r33n = r33 / np.abs(r33).max(axis=1, keepdims=True)
    # 仿真道也逐道归一化到同口径
    s2n = s2 / np.abs(s2).max(axis=1, keepdims=True)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0, t, clim(s0, t), "仿真原始（无增广）")
    panel(axes[1], s1, t, clim(s1, t), "仿真 + P2-1 振铃")
    panel(axes[2], s2, t, clim(s2, t), "仿真 + P2-1 + P2-2 压削/抖动")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道，逐道归一化）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广效果预览 · C1mX-BG（种子 20261001）", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_bscan.png", dpi=150)
    plt.close(fig)

    # ---------- 波形对比：仿真增广道17 vs 实测 Line9 中点道 ----------
    a = s2n[16]
    b = r33n[16]
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(t, a, color="#C0392B", lw=1.0, label="仿真 + 增广（C1mX-BG 道17）")
    ax.plot(tr, b, color="#1F6FEB", lw=1.0, alpha=0.8, label="实测 Line9 中点道")
    ax.set_xlabel("时间 (ns)")
    ax.set_ylabel("归一化振幅")
    ax.set_title("增广后仿真道 vs 实测道（各自归一化）")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_waveform.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
