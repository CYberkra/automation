# -*- coding: utf-8 -*-
"""9号测线 AGC 显示 + PPT 地面真值标注（黏土/基岩界面12-16m、21m异常、真高8m）。"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

CSV = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
OUT = Path(r"E:\automation_djh\fig_line9_annotated.png")
N_SAMPLES, N_TRACES, T_NS = 501, 2378, 700.0

raw = np.loadtxt(CSV, delimiter=",", skiprows=4, dtype=np.float64)
traces = raw[:, 3].reshape(N_TRACES, N_SAMPLES)
t_ns = np.linspace(0, T_NS, N_SAMPLES)
dist_m = np.arange(N_TRACES) * 0.09093

env = np.abs(hilbert(traces, axis=1))
k = 25
kernel = np.ones(k) / k
rms = np.sqrt(np.apply_along_axis(
    lambda r: np.convolve(r, kernel, mode="same"), 1, env ** 2))
agc = traces / (rms + 1e-6)
alim = float(np.percentile(np.abs(agc[:, t_ns >= 25.0]), 98.0))

# 速度标定：界面 14 m ↔ 实测强带 ~445 ns => v ≈ 0.063 m/ns (er≈22.7)
V_M_NS = 2 * 14.0 / 445.0
ER_EST = (0.299792458 / V_M_NS) ** 2
print(f"标定速度 v={V_M_NS:.4f} m/ns, 等效 er={ER_EST:.1f}")


def depth(t):  # 双程走时转深度
    return V_M_NS * t / 2.0


fig, ax = plt.subplots(figsize=(20, 9))
ax.imshow(agc.T, aspect="auto", cmap="gray", vmin=-alim, vmax=alim,
          extent=[0, dist_m[-1], t_ns[-1], t_ns[0]])

# 地面真值标注
t_ifz = 2 * np.array([12.0, 16.0]) / V_M_NS
ax.axhspan(t_ifz[0], t_ifz[1], color="#ff7f0e", alpha=0.15)
ax.axhline(t_ifz[0], color="#ff7f0e", ls="--", lw=1.2)
ax.axhline(t_ifz[1], color="#ff7f0e", ls="--", lw=1.2)
ax.text(dist_m[-1] * 0.62, t_ifz[0] - 18,
        "黏土/基岩界面 12–16 m（ZK08 钻孔证实 ≈14 m）", color="#ff7f0e",
        fontsize=13, fontweight="bold")

t_anom = 2 * 21.0 / V_M_NS
ax.axhline(t_anom, color="#d62728", ls="--", lw=1.2)
ax.text(dist_m[-1] * 0.62, t_anom + 12, "地下异常 ≈21 m（PPT 标注）",
        color="#d62728", fontsize=13, fontweight="bold")

t_surf = 2 * 8.0 / 0.299792458
ax.axhline(t_surf, color="#1f77b4", ls=":", lw=1.2)
ax.text(dist_m[-1] * 0.02, t_surf + 12, "地表反射（真高8 m → ≈53 ns）",
        color="#1f77b4", fontsize=12)

ax.set_xlabel("测线距离 (m)")
ax.set_ylabel("时间 (ns)")
ax.set_title(f"9号测线 AGC 显示 + 地面真值标注（速度标定 v≈{V_M_NS:.3f} m/ns，εr≈{ER_EST:.0f}；"
             f"四川营山大秧坪滑坡，无人机载 SFCW 雷达）")

# 右侧深度轴
ax2 = ax.secondary_yaxis("right", functions=(depth, lambda d: 2 * d / V_M_NS))
ax2.set_ylabel(f"估算深度 (m，按 v={V_M_NS:.3f} m/ns)")

fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"saved: {OUT}")
