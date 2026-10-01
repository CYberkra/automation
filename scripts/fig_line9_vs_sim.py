# -*- coding: utf-8 -*-
"""真实9号测线单道 vs 仿真C3mX单道：早时波形对比（解释一道 vs 两道直达波）。"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

REPO = Path(r"E:\automation_djh\automation_repo")
sys.path.insert(0, str(REPO / "scripts"))
from study_t3_damage_ladder import load_bscan  # noqa: E402

CSV = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
OUT = Path(r"E:\automation_djh\fig_line9_vs_sim_early.png")

# 真实数据：取测线中部一道
raw = np.loadtxt(CSV, delimiter=",", skiprows=4, dtype=np.float64)
traces = raw[:, 3].reshape(2378, 501)
t_real = np.linspace(0, 700.0, 501)
tr = traces[1200]
h_real = raw[1200 * 501, 4]
t_surf_real = 2 * h_real / 0.299792458

# 仿真：C3mX 第17道
sigs, t_sim = load_bscan("B2D-C3mX-BG")
sm = sigs[16]

fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

ax = axes[0]
ax.plot(t_real, tr / np.abs(tr).max(), lw=0.9, color="0.2")
env = np.abs(hilbert(tr))
ax.plot(t_real, env / env.max(), lw=1.0, color="#d62728", alpha=0.7,
        label="包络")
ax.axvline(t_surf_real, color="#1f77b4", ls="--", lw=1)
ax.text(t_surf_real + 3, 0.8, f"地表反射理论到时 ≈{t_surf_real:.0f} ns\n(本道真高 {h_real:.1f} m)",
        color="#1f77b4", fontsize=10)
ax.set_ylim(-1.05, 1.05)
ax.set_xlim(0, 250)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
ax.set_ylabel("归一化振幅")
ax.set_title(f"真实：9号测线第1201道 —— 直达波后拖着 ~150 ns 系统振铃，"
             f"地表反射被淹没在振铃里")

ax = axes[1]
ax.plot(t_sim, sm / np.abs(sm).max(), lw=0.9, color="0.2")
ax.axvline(4.5, color="#2ca02c", ls="--", lw=1)
ax.text(8, 0.8, "直达耦合 4.5 ns", color="#2ca02c", fontsize=10)
ax.axvline(100.5, color="#1f77b4", ls="--", lw=1)
ax.text(104, 0.6, "地表反射 100.5 ns（架高15 m恒定）", color="#1f77b4",
        fontsize=10)
ax.set_ylim(-1.05, 1.05)
ax.grid(alpha=0.3)
ax.set_ylabel("归一化振幅")
ax.set_xlabel("时间 (ns)")
ax.set_title("仿真：C3mX 第17道 —— SFCW+Hann 重建脉冲紧凑，"
             "直达波与地表反射干净分离成两道")

fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"saved: {OUT}")
print(f"真实道真高={h_real:.2f} m, 地表反射理论={t_surf_real:.1f} ns")
