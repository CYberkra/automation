# -*- coding: utf-8 -*-
"""9号测线真实数据解析与概览图（只读分析，不改原始文件）。"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

CSV = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
OUT = Path(r"E:\automation_djh\fig_line9_overview.png")

N_SAMPLES = 501
T_WINDOW_NS = 700.0
N_TRACES = 2378

print("读取 CSV ...")
raw = np.loadtxt(CSV, delimiter=",", skiprows=4, dtype=np.float64)
print(f"原始行数={raw.shape[0]}, 列数={raw.shape[1]}")
assert raw.shape[0] == N_TRACES * N_SAMPLES, raw.shape

lon = raw[:, 0].reshape(N_TRACES, N_SAMPLES)
lat = raw[:, 1].reshape(N_TRACES, N_SAMPLES)
elev = raw[:, 2].reshape(N_TRACES, N_SAMPLES)
amp = raw[:, 3].reshape(N_TRACES, N_SAMPLES)
hght = raw[:, 4].reshape(N_TRACES, N_SAMPLES)

# 一致性检查：每道内 lon/lat/elev/hght 应恒定
for name, a in (("lon", lon), ("lat", lat), ("elev", elev), ("hght", hght)):
    dev = float(np.abs(a - a[:, :1]).max())
    print(f"{name}: 道内最大偏差={dev:.3e}")

traces = amp  # (2378, 501) [trace, sample]
t_ns = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
dist_m = np.arange(N_TRACES) * 0.09093

lon1, lat1 = lon[:, 0], lat[:, 0]
elev1, hght1 = elev[:, 0], hght[:, 0]
print(f"经度范围 {lon1.min():.7f}~{lon1.max():.7f}, 纬度 {lat1.min():.7f}~{lat1.max():.7f}")
print(f"地表高程 {elev1.min():.2f}~{elev1.max():.2f} m, 飞行高度 {hght1.min():.3f}~{hght1.max():.3f} m")
print(f"振幅范围 [{traces.min():.4f}, {traces.max():.4f}], 全程RMS={np.sqrt((traces**2).mean()):.5f}")

# 包络平均剖面（能量随时间分布，定位直达波/地表反射）
from scipy.signal import hilbert
env_mean = np.abs(hilbert(traces, axis=1)).mean(axis=0)
pk = np.argsort(env_mean)[::-1][:6]
print("平均包络最强时刻(ns):", sorted(round(float(t_ns[i]), 1) for i in pk))

# ---------- 绘图 ----------
fig = plt.figure(figsize=(20, 13))
gs = fig.add_gridspec(3, 2, height_ratios=[3, 3, 1.2], hspace=0.32, wspace=0.18)

extent = [0, dist_m[-1], t_ns[-1], t_ns[0]]

# ① 原始符号振幅灰度（剔除 t<25 ns 后按 99 分位限幅）
w = t_ns >= 25.0
alim = float(np.percentile(np.abs(traces[:, w]), 99.0))
ax = fig.add_subplot(gs[0, :])
ax.imshow(traces.T, aspect="auto", cmap="gray", vmin=-alim, vmax=alim,
          extent=extent)
ax.set_title(f"① 9号测线原始 B-scan（符号振幅灰度，限幅=99分位={alim:.4f}，"
             f"{N_TRACES}道×{N_SAMPLES}采样，道距0.091 m，时窗700 ns）")
ax.set_ylabel("时间 (ns)")
ax.set_xlabel("测线距离 (m)")

# ② 增益恢复显示：逐采样 AGC（滑动 RMS 归一）看深部弱信号
ax = fig.add_subplot(gs[1, :])
env = np.abs(hilbert(traces, axis=1))
k = 25  # ~35 ns 滑动窗
kernel = np.ones(k) / k
rms = np.sqrt(np.apply_along_axis(
    lambda r: np.convolve(r, kernel, mode="same"), 1, env ** 2))
agc = traces / (rms + 1e-6)
alim2 = float(np.percentile(np.abs(agc[:, w]), 98.0))
ax.imshow(agc.T, aspect="auto", cmap="gray", vmin=-alim2, vmax=alim2,
          extent=extent)
ax.set_title("② AGC 增益恢复显示（滑动RMS归一，凸显深部弱反射）")
ax.set_ylabel("时间 (ns)")
ax.set_xlabel("测线距离 (m)")

# ③ 高程/飞行高度 + 平均包络
ax = fig.add_subplot(gs[2, 0])
ax.plot(dist_m, elev1, label="地表高程", lw=0.8)
ax.plot(dist_m, elev1 + hght1, label="天线高度(高程+飞行高度)", lw=0.8)
ax.set_xlabel("测线距离 (m)")
ax.set_ylabel("高程 (m)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
ax.set_title("③ 测线地形与飞行高度")

ax = fig.add_subplot(gs[2, 1])
ax.plot(t_ns, env_mean / env_mean.max(), lw=0.9)
ax.set_xlabel("时间 (ns)")
ax.set_ylabel("归一化平均包络")
ax.set_xlim(0, 300)
ax.grid(alpha=0.3)
ax.set_title("④ 全测线平均包络（0–300 ns）")

fig.savefig(OUT, dpi=150)
print(f"saved: {OUT}")
