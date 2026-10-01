# -*- coding: utf-8 -*-
"""振铃幅度比 alpha 的实测标定。

对营山 6 条测线逐道拟合加性振铃模板：
  ring(t') = A * exp(-t'/tau) * sin(2*pi*f0*t' + phi),  t' = t - t_direct
  - f0 = 107 MHz、tau = 73 ns 固定（前期已标定）
  - A、phi 通过线性最小二乘（sin/cos 双基）解出
  - alpha = A / 直达波峰值幅度
拟合窗：直达峰后 25~220 ns（避开主瓣与 >250 ns 深部反射带）。
输出：终端统计 + fig_ringing_alpha_calibration.png（alpha 直方图 + 单道拟合示例）
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

DIRS = [Path(r"E:\automation_djh\real_data_yingshan\营山测线数据"),
        Path(r"E:\automation_djh\temporary product")]
N_SAMPLES = 501
T_WINDOW_NS = 700.0
F0 = 0.107   # cycles / ns
TAU = 73.0   # ns
FIT_LO, FIT_HI = 25.0, 220.0   # 相对直达峰的拟合窗 (ns)

LINES = ["Line3", "Line6", "Line7", "Line9", "LineL1", "LineX1"]


def load(csv_path):
    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        hdr = [f.readline() for _ in range(4)]
    n_traces = int(hdr[2].split("=")[1].strip().rstrip(","))
    raw = np.loadtxt(csv_path, delimiter=",", skiprows=4, dtype=np.float64)
    return raw[:, 3].reshape(n_traces, N_SAMPLES)


def find_csv(name):
    for d in DIRS:
        p = d / f"{name}origin(36).csv"
        if p.exists():
            return p
    raise FileNotFoundError(name)


def fit_trace(tr, t_ns):
    """返回 (alpha, A, peak, t_direct, R2) 或 None。"""
    early = t_ns <= 60
    i0 = int(np.argmax(np.abs(tr[early])))
    t0 = t_ns[i0]
    peak = float(abs(tr[i0]))
    if peak <= 0:
        return None
    tp = t_ns - t0
    w = (tp >= FIT_LO) & (tp <= FIT_HI)
    x = tp[w]
    y = tr[w]
    env = np.exp(-x / TAU)
    b1 = env * np.sin(2 * np.pi * F0 * x)
    b2 = env * np.cos(2 * np.pi * F0 * x)
    B = np.stack([b1, b2], axis=1)
    coef, *_ = np.linalg.lstsq(B, y, rcond=None)
    A = float(np.hypot(*coef))
    yhat = B @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) + 1e-30
    return A / peak, A, peak, t0, 1 - ss_res / ss_tot


def main():
    t_ns = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    all_alpha, per_line = [], {}
    example = None
    for name in LINES:
        amp = load(find_csv(name))
        alphas, r2s = [], []
        for i in range(amp.shape[0]):
            r = fit_trace(amp[i], t_ns)
            if r is None:
                continue
            alpha, A, peak, t0, r2 = r
            alphas.append(alpha)
            r2s.append(r2)
            if name == "Line9" and example is None and i == amp.shape[0] // 2:
                example = (amp[i], t0, A, peak, r2)
        alphas = np.array(alphas)
        per_line[name] = alphas
        all_alpha.append(alphas)
        print(f"{name:7s} n={len(alphas):4d}  alpha 中位={np.median(alphas):.3f} "
              f"p25={np.percentile(alphas,25):.3f} p75={np.percentile(alphas,75):.3f} "
              f"拟合R2中位={np.median(r2s):.2f}")
    pooled = np.concatenate(all_alpha)
    print(f"\n汇总 n={len(pooled)}  alpha 中位={np.median(pooled):.3f} "
          f"均值={pooled.mean():.3f} p10={np.percentile(pooled,10):.3f} "
          f"p90={np.percentile(pooled,90):.3f}")

    # ---------- 图 ----------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
    ax1.hist(pooled, bins=40, color="#6BAED6", edgecolor="white")
    med = np.median(pooled)
    ax1.axvline(med, color="#C0392B", lw=1.5, label=f"中位 α = {med:.3f}")
    ax1.axvline(0.25, color="#888888", ls="--", lw=1.2, label="演示占位值 0.25")
    ax1.set_xlabel("α（振铃幅度 / 直达波峰值）")
    ax1.set_ylabel("道数")
    ax1.set_title(f"α 分布（6 条测线共 {len(pooled)} 道）")
    ax1.legend()

    tr, t0, A, peak, r2 = example
    tp = t_ns - t0
    w = (tp >= -10) & (tp <= 250)
    x = tp[w]
    env = np.exp(-np.clip(x, 0, None) / TAU)
    fit = np.where(x >= FIT_LO, 0.0, 0.0)  # 占位
    # 重建完整模板曲线用于显示
    b1 = env * np.sin(2 * np.pi * F0 * x)
    b2 = env * np.cos(2 * np.pi * F0 * x)
    # 用同一道重拟合拿系数
    wf = (tp >= FIT_LO) & (tp <= FIT_HI)
    xf = tp[wf]
    B = np.stack([np.exp(-xf / TAU) * np.sin(2 * np.pi * F0 * xf),
                  np.exp(-xf / TAU) * np.cos(2 * np.pi * F0 * xf)], axis=1)
    coef, *_ = np.linalg.lstsq(B, tr[wf], rcond=None)
    fit = np.where(x >= 0, A * np.exp(-np.clip(x, 0, None) / TAU) *
                   np.sin(2 * np.pi * F0 * x + np.arctan2(coef[1], coef[0])), 0.0)
    ax2.plot(x, tr[w], color="#9AA5B1", lw=1.0, label="实测道（Line9 中点）")
    ax2.plot(x, fit, color="#C0392B", lw=1.2,
             label=f"拟合振铃  A/峰={A/peak:.3f}, R²={r2:.2f}")
    ax2.axvspan(FIT_LO, FIT_HI, color="#1F6FEB", alpha=0.06, label="拟合窗 25–220 ns")
    ax2.set_xlabel("相对直达峰时间 (ns)")
    ax2.set_ylabel("振幅")
    ax2.set_title("单道拟合示例")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig.suptitle("振铃幅度比 α 实测标定（f0=107 MHz, τ=73 ns 固定）",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(r"E:\automation_djh\fig_ringing_alpha_calibration.png", dpi=150)
    print("fig saved")


if __name__ == "__main__":
    main()
