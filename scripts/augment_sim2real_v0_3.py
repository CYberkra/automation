# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.3 —— 振铃改为卷积注入（系统冲激响应口径）。

v0.2→v0.3 的依据（fig_line9_full_texture.png 全线分析）：
  1) 实测直达波本身即多条纹（0–30 ns 约 4–5 条），与 30–100 ns 振铃带是
     同一个系统冲激响应 —— 加性"峰后挂尾"造成质地断层，应改为卷积：
         y = x + beta_i * (x ⊛ g),  g(t) = exp(-t/73ns) * sin(2π*107MHz*t)
     beta_i 逐道标定使振铃窗内 |x⊛g| 有效幅度 = alpha_i * 直达峰值，
     alpha_i 沿用实测分布（中位 0.11）。
  2) 实测振铃窗邻道相关中位 0.81（深部窗也有 0.83）—— alpha_i 缓变
     （高斯平滑 σ=5 道，±25%），相位锁定随直达波（卷积天然保证），
     不再随机相位。
  3) OP2 地表反射编辑沿用 v0.2（Hann 混合 + 缓变 k/δ）。

输出：fig_augment_v03_bscan.png（原始 | v0.2 | v0.3 | 实测 Line9）
      fig_augment_v03_waveform.png
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
ALPHA_MED = 0.11    # 实测中位
ALPHA_SPREAD = 0.25 # 相对离散（±25%）
SMOOTH_SIGMA = 5.0  # 道
K_LO, K_HI = 0.20, 0.64
JITTER_NS = 15.0
T_CROP_NS = 400.0

LINE9 = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
N_SAMPLES, T_WINDOW_NS = 501, 700.0


def smooth_series(rng, n, lo, hi, sigma):
    """沿道缓变随机序列。修复：reflect 填充 + valid 卷积，保证输出长度恰为 n
    （旧版 same 模式在核宽 > n 时输出长度错误且引入直流偏置，导致 tanh 打满上界）。"""
    x = rng.standard_normal(n)
    r = int(4 * sigma)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    xs = np.convolve(np.pad(x, r, mode="reflect"), k, mode="valid")[:n]
    xs = (xs - xs.mean()) / (xs.std() + 1e-12)
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    return mid + half * np.tanh(0.9 * xs)


def op1_ringing_conv(sigs, t_ns, rng):
    """卷积振铃：y = x + beta*(x⊛g)，beta 逐道标定到 alpha_i 口径。"""
    n, nt = sigs.shape
    dt = float(t_ns[1] - t_ns[0])
    tg = np.arange(0, 5 * TAU, dt)
    g = np.exp(-tg / TAU) * np.sin(2 * np.pi * F0 * tg)
    alphas = ALPHA_MED * (1.0 + ALPHA_SPREAD *
                          (smooth_series(rng, n, -1.0, 1.0, SMOOTH_SIGMA)))
    out = sigs.copy()
    w_direct = (t_ns >= 0.0) & (t_ns <= 60.0)
    idx = np.where(w_direct)[0]
    for i in range(n):
        i0 = idx[int(np.argmax(np.abs(sigs[i][w_direct])))]
        peak = float(np.abs(sigs[i, i0]))
        t0 = t_ns[i0]
        c = np.convolve(sigs[i], g)[:nt] * dt
        # 振铃窗：峰后 40–120 ns 的有效幅度作为标定基准
        w = (t_ns >= t0 + 40) & (t_ns <= t0 + 120)
        ref = float(np.sqrt(np.mean(c[w] ** 2))) + 1e-30
        beta = alphas[i] * peak / ref
        out[i] = sigs[i] + beta * c
    return out


def op2_surface_edit(sigs, t_ns, rng):
    """沿用 v0.2：k、δ 沿道缓变 + Hann 窗平滑混合。"""
    n = sigs.shape[0]
    ks = smooth_series(rng, n, K_LO, K_HI, 3.0)
    dts = smooth_series(rng, n, -JITTER_NS, JITTER_NS, 3.0)
    out = sigs.copy()
    dt = float(t_ns[1] - t_ns[0])
    w = (t_ns >= 70.0) & (t_ns <= 180.0)
    idx = np.where(w)[0]
    for i in range(n):
        is_ = idx[int(np.argmax(np.abs(sigs[i][w])))]
        t_s = t_ns[is_]
        lo = int(np.searchsorted(t_ns, t_s - 30.0))
        hi = min(int(np.searchsorted(t_ns, t_s + 90.0)), sigs.shape[1])
        m = hi - lo
        if m < 8:
            continue
        seg = out[i, lo:hi].copy()
        shift = int(round(dts[i] / dt))
        seg = np.roll(seg, shift)
        if shift > 0:
            seg[:shift] = 0.0
        elif shift < 0:
            seg[shift:] = 0.0
        win = 0.5 * (1 - np.cos(2 * np.pi * np.arange(m) / (m - 1)))
        out[i, lo:hi] = out[i, lo:hi] * (1 - win) + seg * ks[i] * win
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
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0 = sigs[:, crop]
    t = t_ns[crop]

    # v0.2 对照（加性、缓变相位）
    rng = np.random.default_rng(SEED)
    from augment_sim2real_v0_2 import op1_ringing as op1_v02, op2_surface_edit as op2_v02
    s_v02 = op2_v02(op1_v02(s0, t, rng), t, rng)

    # v0.3：卷积振铃
    rng = np.random.default_rng(SEED)
    s1 = op1_ringing_conv(s0, t, rng)
    s2 = op2_surface_edit(s1, t, rng)

    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    rcrop = t_real <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real[rcrop]

    # 同口径显示：四个面板全部逐道归一化，限幅取归一化后 t>30ns 的 99 分位
    def trnorm(a):
        return a / (np.abs(a).max(axis=1, keepdims=True) + 1e-30)

    s0n, sv02n, s2n, r33n = trnorm(s0), trnorm(s_v02), trnorm(s2), trnorm(r33)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0n, t, clim(s0n, t), "仿真原始")
    panel(axes[1], sv02n, t, clim(sv02n, t), "v0.2（加性振铃）")
    panel(axes[2], s2n, t, clim(s2n, t), "v0.3（卷积振铃+相位锁定）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.3 对比 · C1mX-BG · 四面板同口径（逐道归一化，种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v03_bscan.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(t, s2n[16], color="#C0392B", lw=1.0, label="仿真 + 增广 v0.3（C1mX-BG 道17）")
    ax.plot(tr, r33n[16], color="#1F6FEB", lw=1.0, alpha=0.8, label="实测 Line9 中点道")
    ax.set_xlabel("时间 (ns)")
    ax.set_ylabel("归一化振幅")
    ax.set_title("增广 v0.3 仿真道 vs 实测道（各自归一化）")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v03_waveform.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
