# -*- coding: utf-8 -*-
"""诊断2：t0 相对窗的振铃包络衰减剖面对比（v0.4 vs Line9）。

诊断1发现绝对时间窗受 time-zero 偏移污染（Line9 t0≈14ns, 仿真 t0≈1.7ns）。
本轮全部改为相对各道自身直达峰 t0 的窗，回答两个问题：
  A. 振铃幅值衰减剖面形状是否一致（τ=73 假设在实测里是否成立）；
  B. 带通(90-125MHz)后的振铃包络，v0.4 与实测各时段幅度比是多少。
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(r"E:\automation_djh\automation_repo")
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, r"E:\automation_djh")

from run_reward_protocol_b2_pilot_v0_1 import load_bscan  # noqa: E402
from augment_sim2real_v0_3 import (  # noqa: E402
    T_CROP_NS, N_SAMPLES, T_WINDOW_NS, op2_surface_edit, load_line9,
)
from augment_sim2real_v0_4 import op1_ringing_matched, trnorm, SEED  # noqa: E402


def t0_of(x, tt):
    w = (tt >= 0) & (tt <= 60)
    i0 = np.where(w)[0][int(np.argmax(np.abs(x[w])))]
    return tt[i0]


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0, t = sigs[:, crop], t_ns[crop]
    dt = float(t[1] - t[0])

    rng = np.random.default_rng(SEED)
    s_v04 = op2_surface_edit(op1_ringing_matched(sigs, t_ns, rng), t_ns, rng)[:, crop]

    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    rcrop = t_real <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real[rcrop]
    dtr = float(tr[1] - tr[0])

    v4n, r33n = trnorm(s_v04), trnorm(r33)

    from scipy.signal import hilbert, butter, sosfiltfilt
    sos = sos = butter(4, [0.090, 0.125], btype="band", fs=1.0 / dt, output="sos")
    sosr = butter(4, [0.090, 0.125], btype="band", fs=1.0 / dtr, output="sos")

    print("== t0 相对窗 |x| 的 p99（归一化幅值，全道中位±p25/p75） ==")
    rel_wins = [(10, 25), (25, 45), (45, 70), (70, 100), (100, 140), (140, 200)]
    for name, x, tt in [("v0.4", v4n, t), ("Line9", r33n, tr)]:
        rows = []
        for i in range(x.shape[0]):
            t0 = t0_of(x[i], tt)
            row = []
            for lo, hi in rel_wins:
                m = (tt >= t0 + lo) & (tt < t0 + hi)
                row.append(float(np.percentile(np.abs(x[i][m]), 99)) if m.sum() else np.nan)
            rows.append(row)
        rows = np.array(rows)
        med = np.nanmedian(rows, axis=0)
        p25 = np.nanpercentile(rows, 25, axis=0)
        p75 = np.nanpercentile(rows, 75, axis=0)
        print(f"  -- {name} --")
        for (lo, hi), a, b, c in zip(rel_wins, med, p25, p75):
            print(f"    t0+{lo:>3}~{hi:<3}ns  p99_med={a:.3f}  (p25={b:.3f}, p75={c:.3f})")

    print("== 带通 90-125MHz 包络：t0 相对窗包络中位（归一化） ==")
    for name, x, tt, s in [("v0.4", v4n, t, sos), ("Line9", r33n, tr, sosr)]:
        rows = []
        for i in range(x.shape[0]):
            t0 = t0_of(x[i], tt)
            xf = sosfiltfilt(s, x[i])
            env = np.abs(hilbert(xf))
            row = []
            for lo, hi in rel_wins:
                m = (tt >= t0 + lo) & (tt < t0 + hi)
                row.append(float(np.median(env[m])) if m.sum() else np.nan)
            rows.append(row)
        rows = np.array(rows)
        med = np.nanmedian(rows, axis=0)
        print(f"  -- {name} --  (直达峰=1 的归一化口径)")
        for (lo, hi), a in zip(rel_wins, med):
            print(f"    t0+{lo:>3}~{hi:<3}ns  env_med={a:.4f}")
        # 由 t0+45~70 与 t0+100~140 推等效 tau
        e1, e2 = med[2], med[4]
        if e1 > 0 and e2 > 0:
            tau_eff = (140 - 57.5) / np.log(e1 / e2)
            print(f"    -> 等效 tau ≈ {tau_eff:.0f} ns（由两段包络比推算）")


if __name__ == "__main__":
    main()
