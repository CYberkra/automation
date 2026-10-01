# -*- coding: utf-8 -*-
"""直达波保真度量化分析：v0.3 卷积振铃是否扭曲了直达波？

对比三组：
  A 仿真原始 s0
  B 仿真 v0.3 增广（s0 + beta*(s0⊛g)）
  C 实测 Line9 全线

指标（逐道）：
  - 直达峰到时、峰值幅度
  - 直达窗(峰-5~峰+25ns)能量 E_d、振铃窗(峰+40~峰+100ns)能量 E_r，比值 E_r/E_d
  - 直达窗内过零次数（波形宽窄/振荡程度）
  - A→B：峰值变化率、到时漂移、直达窗波形相关系数（保真度）

输出：终端统计 + fig_direct_fidelity.png
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
from augment_sim2real_v0_3 import op1_ringing_conv, F0, TAU  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

LINE9 = Path(r"E:\automation_djh\temporary product\Line9origin(36).csv")
N_SAMPLES, T_WINDOW_NS = 501, 700.0


def stats_block(name, sigs, t_ns):
    """逐道直达波统计，返回 dict。"""
    n = sigs.shape[0]
    w0 = (t_ns >= 0) & (t_ns <= 60)
    idx = np.where(w0)[0]
    t0s, peaks, zcr, eratio = [], [], [], []
    for i in range(n):
        tr = sigs[i]
        i0 = idx[int(np.argmax(np.abs(tr[w0])))]
        t0 = t_ns[i0]
        t0s.append(t0)
        peaks.append(float(np.abs(tr[i0])))
        wd = (t_ns >= t0 - 5) & (t_ns <= t0 + 25)
        wr = (t_ns >= t0 + 40) & (t_ns <= t0 + 100)
        s_d = tr[wd]
        zcr.append(int(np.sum(np.diff(np.sign(s_d)) != 0)))
        eratio.append(float(np.sum(tr[wr] ** 2) / (np.sum(s_d ** 2) + 1e-30)))
    return {"name": name, "t0": np.array(t0s), "peak": np.array(peaks),
            "zcr": np.array(zcr), "eratio": np.array(eratio)}


def show(st):
    print(f"[{st['name']}] n={len(st['t0'])}")
    print(f"  直达峰到时: 中位 {np.median(st['t0']):.1f} ns (p10-p90 {np.percentile(st['t0'],10):.1f}-{np.percentile(st['t0'],90):.1f})")
    print(f"  峰值幅度:   中位 {np.median(st['peak']):.4g} 相对离散 ±{50*(np.percentile(st['peak'],90)/max(np.median(st['peak']),1e-30)-50):.0f}%")
    print(f"  直达窗过零次数: 中位 {np.median(st['zcr']):.0f} (p25-p75 {np.percentile(st['zcr'],25):.0f}-{np.percentile(st['zcr'],75):.0f})")
    print(f"  E_ring/E_direct: 中位 {np.median(st['eratio']):.3f} (p25-p75 {np.percentile(st['eratio'],25):.3f}-{np.percentile(st['eratio'],75):.3f})")


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    rng = np.random.default_rng(20261001)
    s1 = op1_ringing_conv(sigs, t_ns, rng)   # 只含振铃（不含 P2-2），隔离分析

    with open(LINE9, "r", encoding="utf-8", errors="replace") as f:
        hdr = [f.readline() for _ in range(4)]
    n_traces = int(hdr[2].split("=")[1].strip().rstrip(","))
    raw = np.loadtxt(LINE9, delimiter=",", skiprows=4, dtype=np.float64)
    real = raw[:, 3].reshape(n_traces, N_SAMPLES)
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)

    A = stats_block("仿真原始", sigs, t_ns)
    B = stats_block("仿真 v0.3 振铃", s1, t_ns)
    C = stats_block("实测 Line9", real, t_real)
    for st in (A, B, C):
        show(st)

    # A→B 保真度
    w0 = (t_ns >= 0) & (t_ns <= 40)
    ccs, dpeak, dt0 = [], [], []
    for i in range(sigs.shape[0]):
        a, b = sigs[i][w0], s1[i][w0]
        ccs.append(float(np.corrcoef(a, b)[0, 1]))
        dpeak.append(abs(np.abs(b).max() - np.abs(a).max()) / (np.abs(a).max() + 1e-30))
        dt0.append(abs(A["t0"][i] - B["t0"][i]))
    print("\nA→B 直达窗保真度：")
    print(f"  波形相关: 中位 {np.median(ccs):.3f} (min {np.min(ccs):.3f})")
    print(f"  峰值变化率: 中位 {100*np.median(dpeak):.1f}% (p90 {100*np.percentile(dpeak,90):.1f}%)")
    print(f"  到时漂移: 最大 {np.max(dt0):.2f} ns")

    # 图：直达窗波形叠加（道17）+ eratio 分布
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6))
    wz = t_ns <= 40
    i = 16
    ax1.plot(t_ns[wz], sigs[i][wz] / np.abs(sigs[i]).max(), color="#9AA5B1", lw=1.2, label="仿真原始")
    ax1.plot(t_ns[wz], s1[i][wz] / np.abs(s1[i]).max(), color="#C0392B", lw=1.2, label="仿真 v0.3 振铃")
    j = real.shape[0] // 2
    wrz = (t_real >= 0) & (t_real <= 45)
    ax1.plot(t_real[wrz] - 13.0, real[j][wrz] / np.abs(real[j]).max(), color="#1F6FEB", lw=1.0,
             alpha=0.85, label="实测 Line9 中点道（左移 13 ns 对齐）")
    ax1.set_xlabel("时间 (ns)")
    ax1.set_ylabel("归一化振幅")
    ax1.set_title("直达波窗波形对比")
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)

    bins = np.linspace(0, 0.6, 40)
    ax2.hist(C["eratio"], bins=bins, alpha=0.6, color="#1F6FEB", label=f"实测 Line9（中位 {np.median(C['eratio']):.2f}）")
    ax2.hist(B["eratio"], bins=bins, alpha=0.6, color="#C0392B", label=f"仿真 v0.3（中位 {np.median(B['eratio']):.2f}）")
    ax2.hist(A["eratio"], bins=bins, alpha=0.6, color="#9AA5B1", label=f"仿真原始（中位 {np.median(A['eratio']):.2f}）")
    ax2.set_xlabel("E_ring / E_direct（振铃窗/直达窗能量比）")
    ax2.set_ylabel("道数")
    ax2.set_title("振铃/直达能量比分布")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(r"E:\automation_djh\fig_direct_fidelity.png", dpi=150)
    print("fig saved")


if __name__ == "__main__":
    main()
