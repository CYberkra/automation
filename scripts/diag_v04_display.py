# -*- coding: utf-8 -*-
"""诊断 v0.4 B-scan 图"振铃带偏淡"：显示口径 vs 数据本身。

1. 各面板 clim（t>=30ns 的 99 分位）实际取值，以及 v0.4 面板 t>30ns 最强事件是谁；
2. v0.4 vs Line9 逐道 E_ring/E_direct 分布（各自直达窗）；
3. 振铃带横向纹理（包络沿道起伏）对比。
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


def clim_val(sigs, t, q=99.0):
    return float(np.percentile(np.abs(sigs[:, t >= 30.0]), q))


def env_ratio(sigs, t, ring_win, direct_win):
    rm = (t >= ring_win[0]) & (t <= ring_win[1])
    dm = (t >= direct_win[0]) & (t <= direct_win[1])
    e_r = (sigs[:, rm] ** 2).sum(axis=1)
    e_d = (sigs[:, dm] ** 2).sum(axis=1)
    return e_r / np.maximum(e_d, 1e-30)


def main():
    sigs, t_ns = load_bscan("B2D-C1mX-BG", "2026-10-01")
    crop = t_ns <= T_CROP_NS
    s0, t = sigs[:, crop], t_ns[crop]

    rng = np.random.default_rng(SEED)
    s_v04 = op2_surface_edit(op1_ringing_matched(sigs, t_ns, rng), t_ns, rng)[:, crop]

    real = load_line9()
    t_real = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    mid = real.shape[0] // 2
    rcrop = t_real <= T_CROP_NS
    r33 = real[mid - 16: mid + 17][:, rcrop]
    tr = t_real[rcrop]

    s0n, v4n, r33n = trnorm(s0), trnorm(s_v04), trnorm(r33)

    print("== 1) clim (|x|, t>=30ns, q99) ==")
    print(f"  仿真原始 : {clim_val(s0n, t):.4f}")
    print(f"  v0.4     : {clim_val(v4n, t):.4f}")
    print(f"  Line9    : {clim_val(r33n, tr):.4f}")

    print("== 2) 分时段 |x| 统计（归一化后，全道） ==")
    for name, x, tt in [("v0.4", v4n, t), ("Line9", r33n, tr)]:
        a = np.abs(x)
        print(f"  -- {name} --")
        for lo, hi in [(30, 60), (60, 100), (100, 140), (140, 200), (200, 400)]:
            m = (tt >= lo) & (tt < hi)
            if m.sum() == 0:
                continue
            seg = a[:, m]
            print(f"    {lo:>3}-{hi:<3}ns  max={seg.max():.3f}  p99={np.percentile(seg, 99):.3f}  p50={np.median(seg):.4f}")

    print("== 3) 逐道 E_ring/E_direct ==")
    # v0.4：直达峰 ~1.7ns，窗 (-5,+25)/(+40,+100) 相对各自直达峰
    r_v4 = []
    for i in range(v4n.shape[0]):
        w60 = (t >= 0) & (t <= 60)
        i0 = np.where(w60)[0][int(np.argmax(np.abs(v4n[i][w60])))]
        t0 = t[i0]
        r_v4.append(env_ratio(v4n[i:i + 1], t, (t0 + 40, t0 + 100), (t0 - 5, t0 + 25))[0])
    r_v4 = np.array(r_v4)
    # Line9：直达峰 ~14ns
    r_l9 = []
    for i in range(r33n.shape[0]):
        w60 = (tr >= 0) & (tr <= 60)
        i0 = np.where(w60)[0][int(np.argmax(np.abs(r33n[i][w60])))]
        t0 = tr[i0]
        r_l9.append(env_ratio(r33n[i:i + 1], tr, (t0 + 40, t0 + 100), (t0 - 5, t0 + 25))[0])
    r_l9 = np.array(r_l9)
    for nm, r in [("v0.4", r_v4), ("Line9", r_l9)]:
        print(f"  {nm:>5}: med={np.median(r):.4f}  p25-75=({np.percentile(r,25):.4f},{np.percentile(r,75):.4f})  p10-90=({np.percentile(r,10):.4f},{np.percentile(r,90):.4f})")

    print("== 4) 振铃带横向纹理（包络沿道起伏） ==")
    from scipy.signal import hilbert

    def lat_env(x, tt, lo, hi, t0s):
        # 每道相对各自直达峰的振铃窗包络均值
        envs = []
        for i in range(x.shape[0]):
            m = (tt >= t0s[i] + lo) & (tt <= t0s[i] + hi)
            envs.append(np.abs(hilbert(x[i][m])).mean())
        e = np.array(envs)
        return e / e.mean()

    def t0_list(x, tt):
        out = []
        w60 = (tt >= 0) & (tt <= 60)
        for i in range(x.shape[0]):
            i0 = np.where(w60)[0][int(np.argmax(np.abs(x[i][w60])))]
            out.append(tt[i0])
        return out

    e_v4 = lat_env(v4n, t, 40, 100, t0_list(v4n, t))
    e_l9 = lat_env(r33n, tr, 40, 100, t0_list(r33n, tr))
    for nm, e in [("v0.4", e_v4), ("Line9", e_l9)]:
        d1 = np.diff(e)
        print(f"  {nm:>5}: std={e.std():.3f}  diff_rms={np.sqrt((d1 ** 2).mean()):.3f}  range=({e.min():.2f},{e.max():.2f})")


if __name__ == "__main__":
    main()
