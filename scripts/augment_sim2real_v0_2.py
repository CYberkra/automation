# -*- coding: utf-8 -*-
"""sim2real 增广算子 v0.2 —— 针对 v0.1 两处差距的修订：

  修订 1（道间相干性）：振铃是系统属性（实测道间相关 0.72–0.96），
      α、φ 以及 P2-2 的 k、δ 改为沿道缓变（高斯平滑随机游走，σ=3 道），
      消除 v0.1 振铃带过碎与深部椒盐纹理。
  修订 2（窗段平滑）：P2-2 编辑段用 Hann 窗与原信号混合，
      消除硬切造成的地表反射碎段。

算子参数不变（全部实测锚定）：
  OP1 ring(t) = A*exp(-t'/73ns)*sin(2π*107MHz*t' + φ), α~[0.04,0.19], t0=直达峰
  OP2 地表段 [t_s-30, t_s+90]，k~[0.2,0.64]，δ~[-15,15] ns

输出：
  fig_augment_v02_bscan.png     原始 | v0.1 | v0.2 | 实测 Line9
  fig_augment_v02_waveform.png  增广道 vs 实测道
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
SMOOTH_SIGMA_TRACES = 3.0
T_CROP_NS = 400.0

LINE9 = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
N_SAMPLES, T_WINDOW_NS = 501, 700.0


def smooth_series(rng, n, lo, hi, sigma=SMOOTH_SIGMA_TRACES):
    """沿道缓变的随机序列：白噪声 → 高斯平滑 → tanh 压缩到 [lo, hi]。"""
    x = rng.standard_normal(n)
    r = int(4 * sigma)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    xs = np.convolve(np.pad(x, r, mode="reflect"), k, mode="same")[r:-r] if False else \
        np.convolve(x, k, mode="same")
    xs = xs / (xs.std() + 1e-12)
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    return mid + half * np.tanh(0.9 * xs)


def op1_ringing(sigs, t_ns, rng):
    """v0.2：α、φ 沿道缓变。"""
    n = sigs.shape[0]
    alphas = smooth_series(rng, n, ALPHA_LO, ALPHA_HI)
    phis = smooth_series(rng, n, 0.0, 2 * np.pi)
    out = sigs.copy()
    w = (t_ns >= 0.0) & (t_ns <= 60.0)
    idx = np.where(w)[0]
    for i in range(n):
        i0 = idx[int(np.argmax(np.abs(sigs[i][w])))]
        tt = t_ns - t_ns[i0]
        amp = alphas[i] * float(np.abs(sigs[i, i0]))
        out[i] += np.where(tt >= 0,
                           amp * np.exp(-tt / TAU) * np.sin(2 * np.pi * F0 * tt + phis[i]),
                           0.0)
    return out


def op2_surface_edit(sigs, t_ns, rng):
    """v0.2：k、δ 沿道缓变 + Hann 窗平滑混合。"""
    n = sigs.shape[0]
    ks = smooth_series(rng, n, K_LO, K_HI)
    dts = smooth_series(rng, n, -JITTER_NS, JITTER_NS)
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
        win = 0.5 * (1 - np.cos(2 * np.pi * np.arange(m) / (m - 1)))  # Hann
        out[i, lo:hi] = out[i, lo:hi] * (1 - win) + seg * ks[i] * win
    return out


def op1_v01(sigs, t_ns, rng):
    """v0.1 独立随机版（仅用于对照列）。"""
    out = sigs.copy()
    w = (t_ns >= 0.0) & (t_ns <= 60.0)
    idx = np.where(w)[0]
    for i in range(sigs.shape[0]):
        i0 = idx[int(np.argmax(np.abs(sigs[i][w])))]
        alpha = rng.uniform(ALPHA_LO, ALPHA_HI)
        phi = rng.uniform(0, 2 * np.pi)
        amp = alpha * float(np.abs(sigs[i, i0]))
        tt = t_ns - t_ns[i0]
        out[i] += np.where(tt >= 0,
                           amp * np.exp(-tt / TAU) * np.sin(2 * np.pi * F0 * tt + phi), 0.0)
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

    rng = np.random.default_rng(SEED)
    s_v01 = op1_v01(s0, t, rng)             # v0.1 对照
    rng = np.random.default_rng(SEED)
    s1 = op1_ringing(s0, t, rng)            # v0.2
    s2 = op2_surface_edit(s1, t, rng)

    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    rcrop = t_real <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real[rcrop]
    r33n = r33 / np.abs(r33).max(axis=1, keepdims=True)
    s2n = s2 / np.abs(s2).max(axis=1, keepdims=True)

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    panel(axes[0], s0, t, clim(s0, t), "仿真原始")
    panel(axes[1], s_v01, t, clim(s_v01, t), "v0.1（逐道独立随机）")
    panel(axes[2], s2, t, clim(s2, t), "v0.2（道间缓变+窗平滑）")
    panel(axes[3], r33n, tr, clim(r33n, tr), "实测 Line9（33 道，逐道归一化）")
    axes[0].set_ylabel("时间 (ns)")
    fig.suptitle("sim2real 增广 v0.2 对比 · C1mX-BG（种子 20261001）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(r"E:\automation_djh\fig_augment_v02_bscan.png", dpi=150)
    plt.close(fig)

    a = s2n[16]
    b = r33n[16]
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(t, a, color="#C0392B", lw=1.0, label="仿真 + 增广 v0.2（C1mX-BG 道17）")
    ax.plot(tr, b, color="#1F6FEB", lw=1.0, alpha=0.8, label="实测 Line9 中点道")
    ax.set_xlabel("时间 (ns)")
    ax.set_ylabel("归一化振幅")
    ax.set_title("增广 v0.2 仿真道 vs 实测道（各自归一化）")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_augment_v02_waveform.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
