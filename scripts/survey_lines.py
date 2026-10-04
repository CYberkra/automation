# -*- coding: utf-8 -*-
"""营山 6 条测线 CSV 统一勘察（只读分析）。

对每条测线报告：道数/采样/时窗、道距与线长、高程与真高范围、
直达波峰值到时、深部强反射带位置（>250 ns 最大包络到时）。
"""
from pathlib import Path

import numpy as np

DIR = Path(r"E:\automation_djh\real_data_yingshan\营山测线数据")
N_SAMPLES = 501
T_WINDOW_NS = 700.0


def survey(csv_path):
    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        hdr = [f.readline() for _ in range(4)]
    n_traces = int(hdr[2].split("=")[1].strip().rstrip(","))
    raw = np.loadtxt(csv_path, delimiter=",", skiprows=4, dtype=np.float64)
    assert raw.shape[0] == n_traces * N_SAMPLES, (csv_path, raw.shape, n_traces)
    lon = raw[:, 0].reshape(n_traces, N_SAMPLES)[:, 0]
    lat = raw[:, 1].reshape(n_traces, N_SAMPLES)[:, 0]
    elev = raw[:, 2].reshape(n_traces, N_SAMPLES)[:, 0]
    amp = raw[:, 3].reshape(n_traces, N_SAMPLES)
    hght = raw[:, 4].reshape(n_traces, N_SAMPLES)[:, 0]
    # 道间距：经纬度 -> 近似米（当地纬度 ~31.26N）
    dlat = np.diff(lat) * 111320.0
    dlon = np.diff(lon) * 111320.0 * np.cos(np.deg2rad(np.mean(lat)))
    step = np.hypot(dlon, dlat)
    t_ns = np.linspace(0, T_WINDOW_NS, N_SAMPLES)
    # 直达波峰值到时（每道 |amp| 在 0-60 ns 内的峰值）
    early = amp[:, t_ns <= 60]
    t_direct = t_ns[: early.shape[1]][np.argmax(np.abs(early), axis=1)]
    # 深部最强反射到时（>250 ns 包络峰值）
    deep = amp[:, t_ns >= 250]
    t_deep = t_ns[t_ns >= 250][np.argmax(np.abs(deep), axis=1)]
    # 深部/浅部功率比（>250 ns vs 30-120 ns）
    p_deep = float(np.mean(deep ** 2))
    p_shallow = float(np.mean(amp[:, (t_ns >= 30) & (t_ns <= 120)] ** 2))
    return {
        "n_traces": n_traces,
        "step_m": (float(np.median(step)), float(step.min()), float(step.max())),
        "length_m": float(step.sum()),
        "elev_m": (float(elev.min()), float(elev.max())),
        "alt_m": (float(hght.min()), float(np.median(hght)), float(hght.max())),
        "t_direct_ns": (float(np.median(t_direct)), float(t_direct.min()), float(t_direct.max())),
        "t_deep_ns": (float(np.median(t_deep)), float(np.percentile(t_deep, 10)), float(np.percentile(t_deep, 90))),
        "deep_shallow_db": float(10 * np.log10(p_deep / p_shallow)),
        "amp_absmax": float(np.abs(amp).max()),
    }


print(f"{'line':7s} {'traces':>6s} {'step_m(med/min/max)':>26s} {'len_m':>7s} {'elev_m':>16s} {'alt_m(min/med/max)':>24s} {'t_direct(med/min/max)':>24s} {'t_deep(med/p10/p90)':>24s} {'deep/sh dB':>10s} {'ampmax':>9s}")
for name in ["Line3", "Line6", "Line7", "Line9", "LineL1", "LineX1"]:
    p = DIR / f"{name}origin(36).csv"
    if not p.exists():
        p = Path(r"E:\automation_djh\temporary product") / f"{name}origin(36).csv"
    r = survey(p)
    st = r["step_m"]; al = r["alt_m"]; td = r["t_direct_ns"]; tk = r["t_deep_ns"]
    print(f"{name:7s} {r['n_traces']:6d} {st[0]:9.4f} {st[1]:7.4f} {st[2]:7.4f} {r['length_m']:7.1f} "
          f"{r['elev_m'][0]:7.2f}-{r['elev_m'][1]:7.2f} {al[0]:7.2f} {al[1]:7.2f} {al[2]:7.2f} "
          f"{td[0]:10.1f} {td[1]:6.1f} {td[2]:6.1f} {tk[0]:10.1f} {tk[1]:6.1f} {tk[2]:6.1f} "
          f"{r['deep_shallow_db']:10.1f} {r['amp_absmax']:9.4f}")
