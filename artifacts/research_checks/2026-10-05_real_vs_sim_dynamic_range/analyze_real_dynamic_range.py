"""实测营山 CSV 早/晚期动态范围体检：回答'实测为何没看到仿真的淹没/畸变'。

只读原始 CSV；输出到 artifacts/research_checks/2026-10-05_real_vs_sim_dynamic_range/。
- 早/晚期能量比（dB），与仿真'界面回波/直达 ≈ -72 dB'对比
- 首达时间 -> 推断航高/时零
- 幅值随时间衰减形态 -> 判断是否已含增益/AGC
- 灰度 B-scan + 平均道 dB 曲线
"""
import json
import os
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "real_data_yingshan/营山测线数据/Line9origin(36).csv"
OUT = "artifacts/research_checks/2026-10-05_real_vs_sim_dynamic_range"
os.makedirs(OUT, exist_ok=True)

# ---- 头部 ----
with open(CSV, encoding="utf-8") as f:
    hdr = [f.readline() for _ in range(4)]
n_samples = int(hdr[0].split("=")[1].split(",")[0])
t_window_ns = float(hdr[1].split("=")[1].split(",")[0])
n_traces = int(hdr[2].split("=")[1].split(",")[0])
dx = float(hdr[3].split("=")[1].split(",")[0])

data = np.loadtxt(CSV, delimiter=",", skiprows=4)
assert data.shape[0] == n_traces * n_samples, data.shape
amp = data[:, 3].reshape(n_traces, n_samples).T  # [sample, trace]
alt = data[4::n_samples, 4] if False else data[:, 4].reshape(n_traces, n_samples)[:, 0]

t_ns = np.linspace(0.0, t_window_ns, n_samples)

# ---- 包络（Hilbert 不行就用滑动 RMS） ----
try:
    from scipy.signal import hilbert
    env = np.abs(hilbert(amp, axis=0))
    env_method = "hilbert"
except Exception:
    k = max(3, int(round(10.0 / (t_ns[1] - t_ns[0]))))
    ker = np.ones(k) / k
    env = np.sqrt(np.apply_along_axis(lambda x: np.convolve(x * x, ker, "same"), 0, amp))
    env_method = f"moving_rms_k{k}"

# ---- 指标 ----
def band_rms(e, t0, t1):
    m = (t_ns >= t0) & (t_ns < t1)
    return np.sqrt(np.mean(e[m] ** 2, axis=0))  # per trace

rms_early = band_rms(env, 0, 100)
rms_mid = band_rms(env, 200, 400)
rms_late = band_rms(env, 400, 700)
eps = np.finfo(float).tiny

ratio_mid_db = 20 * np.log10(np.maximum(rms_mid, eps) / np.maximum(rms_early, eps))
ratio_late_db = 20 * np.log10(np.maximum(rms_late, eps) / np.maximum(rms_early, eps))

peak_t = t_ns[np.argmax(np.max(env, axis=1))]
mean_env = np.mean(env, axis=1)
mean_db = 20 * np.log10(np.maximum(mean_env, mean_env.max() * 1e-12) / mean_env.max())

# 晚期相对峰值的余量：晚期最强反射 vs 早期峰值
late_peak = np.max(env[t_ns >= 300], axis=0)
early_peak = np.max(env[t_ns < 150], axis=0)
ratio_peak_db = 20 * np.log10(np.maximum(late_peak, eps) / np.maximum(early_peak, eps))

summary = {
    "file": CSV,
    "n_traces": n_traces, "n_samples": n_samples,
    "t_window_ns": t_window_ns, "dx_m": dx,
    "envelope": env_method,
    "altitude_m_median": float(np.median(alt)),
    "altitude_m_minmax": [float(alt.min()), float(alt.max())],
    "peak_time_ns": float(peak_t),
    "rms_200_400_vs_0_100_dB": {
        "median": float(np.median(ratio_mid_db)),
        "p05": float(np.percentile(ratio_mid_db, 5)),
        "p95": float(np.percentile(ratio_mid_db, 95)),
    },
    "rms_400_700_vs_0_100_dB": {
        "median": float(np.median(ratio_late_db)),
        "p05": float(np.percentile(ratio_late_db, 5)),
        "p95": float(np.percentile(ratio_late_db, 95)),
    },
    "late_peak_vs_early_peak_dB": {
        "median": float(np.median(ratio_peak_db)),
        "p05": float(np.percentile(ratio_peak_db, 5)),
        "p95": float(np.percentile(ratio_peak_db, 95)),
    },
    "amp_abs_max": float(np.abs(amp).max()),
    "amp_abs_median": float(np.median(np.abs(amp))),
}
with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

# ---- 图 1：灰度 B-scan（包络 dB，各自归一） ----
env_db = 20 * np.log10(np.maximum(env, env.max() * 1e-10) / env.max())
fig, ax = plt.subplots(figsize=(14, 5))
ax.imshow(env_db, aspect="auto", cmap="gray", vmin=-70, vmax=0,
          extent=[0, n_traces * dx, t_window_ns, 0])
ax.set_xlabel("distance (m)")
ax.set_ylabel("time (ns)")
ax.set_title("Line9(36) envelope dB (self-normalized, [-70,0])")
fig.savefig(os.path.join(OUT, "bscan_gray_db.png"), dpi=120)
plt.close(fig)

# ---- 图 2：平均道 dB 曲线 + 原始首尾道波形 ----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
axes[0].plot(mean_db, t_ns)
axes[0].invert_yaxis()
axes[0].set_xlabel("mean envelope dB")
axes[0].set_ylabel("time (ns)")
axes[0].set_title("mean trace decay")
axes[0].grid(alpha=0.3)
axes[1].plot(amp[:, 0], t_ns, lw=0.6, label="trace0")
axes[1].plot(amp[:, n_traces // 2], t_ns, lw=0.6, label=f"trace{n_traces//2}")
axes[1].invert_yaxis()
axes[1].set_xlabel("amplitude (signed)")
axes[1].legend()
axes[1].set_title("raw A-scans")
axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "trace_decay.png"), dpi=120)
plt.close(fig)

print(json.dumps(summary, ensure_ascii=False, indent=2))
