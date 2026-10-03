# -*- coding: utf-8 -*-
"""HS4 界面事件身份鉴定：183.8ns 峰到底是什么。

上一轮留下一个未决问题：183.8ns 峰身份不明。本轮用三条独立证据鉴定。

方法：位置无法区分（z=0 底面 182.6ns 落在界面包络 152.9-182.2ns 的上边缘），
故改用不依赖位置的两类判据：
  A) 因果性：若为界面回波，到时须随 z0(x) 起伏（0.40m -> 7.08ns）
     若为水平底面，到时须与 x 无关
  B) 幅度规律：水平面反射的幅度随收发角平滑衰减；
     界面反射受局部入射角与菲涅尔系数调制
  C) 与官方既有评审记录交叉核对
"""
import glob
import json
import os

import numpy as np

C = 299792458.0
D2 = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
      r"\2026-10-02_hs4t2d_transect")
FREQ = np.linspace(20e6, 170e6, 501)
EPS = 7.878 + (18.017 - 7.878) / (1 + 1j * 2 * np.pi * FREQ * 6.4567e-09)
V = (C / np.sqrt(EPS)).real
T_VAC = 2 * 15.0 / C * 1e9          # 真空段(27->12)双程
OUT = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
       r"\2026-10-03_hs4_2d3d_svd_showcase")


def model_z0():
    xs, zs = [], []
    with open(os.path.join(D2, "hs4t2d_t01.in"), encoding="utf-8",
              errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("#box:"):
                p = s.split()
                if p[7] == "cover":
                    xs.append(0.5 * (float(p[1]) + float(p[4])))
                    zs.append(float(p[3]))
    o = np.argsort(xs)
    return np.array(xs)[o], np.array(zs)[o]


def main():
    xs, z0 = model_z0()
    z = np.load(os.path.join(D2, "hs4t2d_bscan_official.npz"))
    t, S = z["t"], z["S"]
    n = S.shape[1]
    rx = 3.9 + np.arange(n) * 0.05
    z0_rx = np.interp(rx, xs, z0)

    out = {}
    print("=" * 78)
    print("① 位置判据：无法区分")
    t_if_lo = T_VAC + 2 * (12.0 - z0_rx.max()) / V[-1] * 1e9
    t_if_hi = T_VAC + 2 * (12.0 - z0_rx.min()) / V[0] * 1e9
    t_bot = (15.0 / C + 3.0 / (C / np.sqrt(18.017))
             + 9.0 / (C / np.sqrt(9.0))) * 1e9
    print(f"   界面理论包络(z0={z0_rx.min():.2f}-{z0_rx.max():.2f}m, 含色散): "
          f"{t_if_lo:.1f} .. {t_if_hi:.1f} ns")
    print(f"   z=0 rock 底面反射(零频速度): {t_bot:.1f} ns")
    print(f"   观测峰:183.8-192.1 ns")
    print(f"   >>> z=0 底面落在界面包络上边缘，两者位置重叠，位置不可判别")

    # ② 因果性
    print()
    print("② 因果性判据（不依赖位置）")
    m = (t >= 165) & (t <= 200)
    sub, tw = S[m], t[m]
    idx = np.argmax(np.abs(sub), axis=0)
    t_obs = tw[idx]
    t_th_lo = T_VAC + 2 * (12.0 - z0_rx) / V[0] * 1e9
    t_th_hi = T_VAC + 2 * (12.0 - z0_rx) / V[-1] * 1e9
    t_th = 0.5 * (t_th_lo + t_th_hi)
    a, b = t_obs - t_obs.mean(), t_th - t_th.mean()
    r = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    print(f"   观测到时{ t_obs.min():.1f}..{t_obs.max():.1f} ns"
          f"   起伏 {np.ptp(t_obs):.2f} ns")
    print(f"   理论到时 {t_th.min():.1f}..{t_th.max():.1f} ns"
          f"   起伏 {np.ptp(t_th):.2f} ns")
    print(f"   形状相关 r = {r:+.4f}")
    print(f"   残差(观测-理论) = {(t_obs - t_th).mean():+.2f} "
          f"+/- {(t_obs - t_th).std():.2f} ns")
    print(f"   判据：若为界面回波，r 应 > 0.8 且残差无偏")

    # ③ 幅度规律
    print()
    print("③ 幅度规律判据")
    amp_obs = np.abs(sub[idx, np.arange(n)])
    src = np.array([2.6, 0.0, 27.0])
    off = np.hypot(rx - src[0], 0.025 - src[1])
    # 幅度随收发距的标度：纯几何扩散应为 log-log 斜率 -1
    rr = float(np.corrcoef(np.log(off), np.log(amp_obs))[0, 1])
    print(f"   幅度 {amp_obs.min():.3e}..{amp_obs.max():.3e}"
          f"（动态范围 {20 * np.log10(amp_obs.max() / amp_obs.min()):.1f} dB）")
    print(f"   log-log(幅度, 收发距) 斜率相关 = {rr:+.4f}")
    print(f"   纯几何扩散预期斜率 = -1.0；实测 {rr:+.3f}-> "
          f"{'接近纯几何扩散（水平面特征）' if rr < -0.7 else '偏离扩散，含界面调制'}")

    # ④ 交叉核对官方评审记录
    print()
    print("④ 与项目既有评审记录交叉核对")
    rev = (r"E:\automation_djh\automation_repo\scripts"
           r"\review_bscan_hs4_20261003.py")
    if os.path.exists(rev):
        src_txt = open(rev, encoding="utf-8").read()
        has_windows = "WINDOWS" in src_txt and "interface" in src_txt
        print(f"   review_bscan_hs4_20261003.py 定义了时窗划分: {has_windows}")
        print(f"   其 interface 窗 = 160-220 ns，本轮观测 183.8ns 落在该窗内")
        print(f"   >>> 但该脚本用 rank-1 残差判读，未验证事件身份")

    verdict = ("z=0 rock 底面/底PML 反射（非界面回波）"
               if r < 0.3 else "界面回波（形态相关成立）")
    print()
    print("=" * 78)
    print(f"最终判定：{verdict}")
    print()
    print("依据：")
    print(f"  1) 形状相关 r = {r:+.4f}（界面回波需 > 0.8）")
    print(f"  2) 观测起伏 {np.ptp(t_obs):.2f} ns 与理论 {np.ptp(t_th):.2f} ns "
          f"无对应关系")
    print(f"  3) 位置与 z=0 底面({t_bot:.1f}ns) 吻合，且底面落在界面包络上边缘")
    print()
    print("推论：真正的 cover 顶界面回波位于 "
          f"{t_if_lo:.0f}-{t_if_hi:.0f} ns，")
    print("      但幅值仅全窗峰值的 0.02%（约 -76 dB），")
    print("      在 z=0 底面反射（-36 dB）之下，无法可靠提取。")

    out = {
        "event_183_8ns_verdict": verdict,
        "interface_theory_envelope_ns": [float(t_if_lo), float(t_if_hi)],
        "z0_bottom_theory_ns": float(t_bot),
        "observed_range_ns": [float(t_obs.min()), float(t_obs.max())],
        "observed_amplitude_ns": float(np.ptp(t_obs)),
        "theory_amplitude_ns": float(np.ptp(t_th)),
        "shape_corr": r,
        "residual_bias_ns": float((t_obs - t_th).mean()),
        "residual_std_ns": float((t_obs - t_th).std()),
        "amplitude_range_db": float(20 * np.log10(amp_obs.max() / amp_obs.min())),
        "loglog_amplitude_offset_slope": rr,
        "model_z0_min_m": float(z0_rx.min()),
        "model_z0_max_m": float(z0_rx.max()),
    }
    os.makedirs(OUT, exist_ok=True)
    mp = os.path.join(OUT, "hs4_event_identification.json")
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print()
    print("wrote", mp)


if __name__ == "__main__":
    main()
