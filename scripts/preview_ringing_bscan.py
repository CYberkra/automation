# -*- coding: utf-8 -*-
"""振铃注入演示：在已有试点批 B-scan 上叠加实测振铃模板（107 MHz, tau=73 ns）。

振铃模型（标准加性天线混响模型）：
  ring(t) = A * exp(-(t-t0)/tau) * sin(2*pi*f0*(t-t0)),  t >= t0
  - t0：每道直达波峰值时刻（15~70 ns 窗口内 |sig| 最大处）
  - A = alpha * 直达波峰值幅度，alpha=0.25
输出：
  fig_ringing_bscan.png     三族 BG：原始 | 加振铃 | 振铃成分
  fig_ringing_waveform.png  C1mX 第 17 道时域波形对比
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

F0_HZ = 107e6       # 振铃中心频率（实测标定）
TAU_NS = 73.0       # 衰减时间常数（实测标定）
ALPHA = 0.11        # 振铃幅度 / 直达波峰值（实测标定中位，见 calibrate_ringing_alpha.py）
T_CROP_NS = 400.0

FAMS = [
    ("C1mX 平地·覆土1m", "B2D-C1mX-BG", "2026-10-01"),
    ("C3mX 平地·覆土3m", "B2D-C3mX-BG", "2026-09-28"),
    ("S2X 坡地·覆土3m", "B2D-C3mS2X-BG", "2026-09-28"),
]


def add_ringing(sigs, t_ns, alpha=ALPHA):
    """对每道注入振铃，返回 (加振铃信号, 振铃成分)。"""
    out = sigs.copy()
    ring_all = np.zeros_like(sigs)
    w = (t_ns >= 15.0) & (t_ns <= 70.0)
    for i in range(sigs.shape[0]):
        seg = np.abs(sigs[i][w])
        i0 = np.where(w)[0][int(np.argmax(seg))]
        t0 = t_ns[i0]
        amp = alpha * float(np.abs(sigs[i, i0]))
        tt = t_ns - t0
        ring = np.where(tt >= 0,
                        amp * np.exp(-tt / TAU_NS) * np.sin(2 * np.pi * (F0_HZ * 1e-9) * tt),
                        0.0)
        out[i] += ring
        ring_all[i] = ring
    return out, ring_all


def panel(ax, sigs, t_ns, cmap, vlim, title=""):
    extent = [1, sigs.shape[0], t_ns[-1], t_ns[0]]
    ax.imshow(sigs.T, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim, extent=extent)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("道号")


def clim(sigs, t_ns, q=99.0):
    return float(np.percentile(np.abs(sigs[:, t_ns >= 30.0]), q))


def main():
    # ---------- 图 1：三族 原始 | 加振铃 | 振铃成分 ----------
    fig, axes = plt.subplots(3, 3, figsize=(17, 12))
    for r, (label, mother, date) in enumerate(FAMS):
        sigs, t_ns = load_bscan(mother, date)
        crop = t_ns <= T_CROP_NS
        s = sigs[:, crop]
        t = t_ns[crop]
        s_ring, ring = add_ringing(s, t)
        panel(axes[r, 0], s, t, "gray", clim(s, t),
              "原始（无振铃）" if r == 0 else "")
        panel(axes[r, 1], s_ring, t, "gray", clim(s_ring, t),
              f"加振铃（f0=107 MHz, τ=73 ns, α={ALPHA}）" if r == 0 else "")
        v = float(np.percentile(np.abs(ring), 99.9))
        panel(axes[r, 2], ring, t, "RdBu_r", v,
              "振铃成分（注入部分）" if r == 0 else "")
        axes[r, 0].set_ylabel(f"{label}\n时间 (ns)")
    fig.suptitle("振铃注入演示 · B-scan（符号振幅灰度，t<30 ns 剔除后 99 分位限幅）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(r"E:\automation_djh\fig_ringing_bscan.png", dpi=150)
    plt.close(fig)

    # ---------- 图 2：C1mX 第 17 道波形对比 ----------
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= 300.0
    s = sigs[16, crop]
    t = t_ns[crop]
    s_ring, ring = add_ringing(s[None, :], t)
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(t, s, color="#9AA5B1", lw=1.0, label="原始")
    ax.plot(t, s_ring[0], color="#C0392B", lw=1.0, label=f"加振铃（α={ALPHA}）")
    ax.plot(t, ring[0], color="#1F6FEB", lw=1.0, label="振铃成分")
    ax.set_xlabel("时间 (ns)")
    ax.set_ylabel("振幅")
    ax.set_title("C1mX-BG 第 17 道 · 振铃注入前后波形（107 MHz / τ=73 ns 指数衰减正弦）")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_ringing_waveform.png", dpi=150)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
