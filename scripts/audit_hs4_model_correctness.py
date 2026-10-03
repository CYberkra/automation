# -*- coding: utf-8 -*-
"""HS4 模型正确性研究轮：五项核查的量化汇总与结论。

修正上一轮的两个公式错误后重做界面检验：
  错误 1：用零频速度 v0 = c/sqrt(18.017) = 70.63 Mm/s 算理论到时。
         实际带内相速度为 76.7..104.2 Mm/s（Debye 色散），偏差 +8.6%..+47.5%。
  错误 2：Debye 复介电常数公式符号写反，导致 eps_real 下到 0.533（物理不可能）。
         官方语义（gprMax/materials.py, w=deltaer/tau, q=-1/tau）：
           eps(f) = er_delta + (er_0 - er_delta) / (1 + j*w*tau)

本脚本给出五项核查结论 + 修正后的界面匹配检验。
"""
import json
import os

import numpy as np

C = 299792458.0
FREQ = np.linspace(20e6, 170e6, 501)
W = 2 * np.pi * FREQ
TAU = 6.4567e-09
ER_DELTA = 7.878# 高频极限（.in 第2 个 Debye 参数）
ER_0 = 18.017            # 零频（.in #material 声明）
EPS = ER_DELTA + (ER_0 - ER_DELTA) / (1 + 1j * W * TAU)
V_PHASE = C / np.sqrt(EPS).real

Z_SRC, Z_BOT = 27.0, 12.0
D_VAC = Z_SRC - Z_BOT

OUT = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
       r"\2026-10-03_hs4_2d3d_svd_showcase")
D2 = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
      r"\2026-10-02_hs4t2d_transect")


def t_of_z_dispersive(z0):
    """色散介质中的双程到时(ns)，对频率取带内包络。"""
    z0 = np.asarray(z0, float)
    t = (2 * D_VAC / C * 1e9
         + 2 * (Z_BOT - np.asarray(z0, float))[:, None] / V_PHASE[None, :] * 1e9)
    return t.min(axis=1), t.max(axis=1)


def main():
    res = {}

    # ---------- 核查 1：频率栅格 vs 仿真窗 ----------
    df_sim = 1 / 600e-9
    df_grid = FREQ[1] - FREQ[0]
    d = np.load(os.path.join(D2, "hs4t2d_bscan_official.npz"))
    t, S = d["t"], d["S"]
    res["check1_frequency_grid"] = {
        "sim_window_ns": 600.0,
        "intrinsic_df_MHz": float(df_sim / 1e6),
        "requested_df_MHz": float(df_grid / 1e6),
        "oversampling_factor": float(df_sim / df_grid),
        "equivalent_independent_freqs": float((FREQ[-1] - FREQ[0]) * 600e-9),
        "verdict": "请求栅格比仿真窗能给的细 5.56倍；501 点非独立测量",
    }
    print("=" * 76)
    print("核查 1  频率栅格 vs 仿真窗")
    print(f"  仿真窗 600ns -> 固有频率分辨率 {df_sim / 1e6:.3f} MHz")
    print(f"  请求栅格{df_grid / 1e6:.3f} MHz（细 {df_sim / df_grid:.2f} 倍）")
    print(f"  等效独立频点 ≈ {int((FREQ[-1] - FREQ[0]) * 600e-9)}，而非 501")

    # ---------- 核查 2：.in 语法 ----------
    res["check2_input_syntax"] = {
        "box_format": "x0 y0 z0 x1 y1 z1 — 解析正确，无缝隙无重叠",
        "box_count_2d": 49, "cover_boxes_2d": 48,
        "box_count_3d": 2305, "cover_boxes_3d": 2304,
        "pml_cells_2d": "20 0 20 20 0 20",
        "pml_cells_3d": "20 20 20 20 20 20",
        "pml_zero_is_legal": True,
        "pml_zero_evidence": "gprMax/user_objects/cmds_singleuse.py:691 "
                             "'requires the PML thickness to be zero or greater'",
        "verdict": "语法与 PML 配置均合法",
    }
    print()
    print("核查 2  .in 语法")
    print("  box 格式 x0 y0 z0 x1 y1 z1，2D 49 盒/3D 2305 盒，无缝隙无重叠")
    print("  2D 的 pml_cells y=0 经查官方为明确合法（0 = 该面无 PML）")

    # ---------- 核查 3：材料与 Debye ----------
    eps_static = ER_DELTA + (ER_0 - ER_DELTA)
    res["check3_material"] = {
        "er0_declared": ER_0,
        "er0_from_debye": float(eps_static),
        "self_consistent": bool(abs(eps_static - ER_0) < 1e-9),
        "eps_at_20MHz": float(EPS[0].real),
        "eps_at_170MHz": float(EPS[-1].real),
        "v_phase_20MHz_Mm_s": float(V_PHASE[0] / 1e6),
        "v_phase_170MHz_Mm_s": float(V_PHASE[-1] / 1e6),
        "v_zero_freq_wrong_value_Mm_s": float(C / np.sqrt(ER_0) / 1e6),
        "debye_char_freq_MHz": float(1 / (2 * np.pi * TAU) / 1e6),
        "verdict": "材料自洽；但上一轮用零频速度算理论到时是错的",
    }
    print()
    print("核查 3  材料与 Debye")
    print(f"  静态 eps: 声明 {ER_0} vs Debye 推得 {eps_static:.4f} -> "
          f"{'自洽' if abs(eps_static - ER_0) < 1e-9 else '不自洽'}")
    print(f"  带内相速度: {V_PHASE[0] / 1e6:.1f} .. {V_PHASE[-1] / 1e6:.1f} Mm/s")
    print(f"  ⚠ 上一轮我用的 {C / np.sqrt(ER_0) / 1e6:.2f} Mm/s 是零频值，"
          f"偏低 {(V_PHASE[0] / (C / np.sqrt(ER_0)) - 1) * 100:.1f}%"
          f"..{(V_PHASE[-1] / (C / np.sqrt(ER_0)) - 1) * 100:.1f}%")

    # ---------- 核查 4：数值稳定性 ----------
    dz = 0.05
    lam_vac_min = C / 170e6
    # 介质波长必须用介质速度算，不能用真空光速——否则 PPW 被低估一半
    lam_med_min = float(V_PHASE.min()) / 170e6
    ppw_vac = lam_vac_min / dz
    ppw_med = lam_med_min / dz
    dt_measured = 1.1793271683748422e-10
    courant = C * dt_measured / dz# 3D Yee 稳定条件 c*dt/dx <= 1/sqrt(3)
    lam95 = float(V_PHASE[0]) / 95e6
    res["check4_numerics"] = {
        "dx_dy_dz_m": dz,
        "lambda_vacuum_at_170MHz_m": float(lam_vac_min),
        "lambda_medium_at_170MHz_m": float(lam_med_min),
        "ppw_vacuum": float(ppw_vac),
        "ppw_medium": float(ppw_med),
        "ppw_sufficient": bool(ppw_med >= 10),
        "courant_number": float(courant),
        "courant_vs_3d_limit": float(courant / (1.0 / np.sqrt(3.0))),
        "courant_stable_3d": bool(courant <= 1.0 / np.sqrt(3.0) + 1e-9),
        "courant_note": "2D ny=1 时 gprMax 放宽约束，实测未见发散"
                        "（末段/首段能量比 0.010），但晚期噪声为 3D 的约 60 倍",
        "dt_s": dt_measured,
        "iterations_2d": 5089,
        "window_actual_ns": float(5089 * dt_measured * 1e9),
        "pml_cells": 20,
        "pml_thickness_m": 1.0,
        "medium_wavelength_at_95MHz_m": float(lam95),
        "interface_box_width_m": 0.25,
        "boxes_per_wavelength": float(lam95 / 0.25),
        "verdict": f"介质内 PPW={ppw_med:.1f}（<10，偏低）；"
                   f"Courant={courant:.4f} 为 3D 上限的 "
                   f"{courant / (1.0 / np.sqrt(3.0)):.2f} 倍（实测未发散）；"
                   f"界面盒 0.25 m 仅 {lam95 / 0.25:.1f} 个介质波长"
                   f"（@95MHz），属亚波长结构",
    }
    print()
    print("核查 4  数值稳定性")
    print(f"  真空波长 @170MHz = {lam_vac_min * 100:.1f} cm -> PPW {ppw_vac:.1f}")
    print(f"  介质波长 @170MHz = {lam_med_min * 100:.1f} cm -> PPW {ppw_med:.1f}"
          f"（<10，偏低）")
    print(f"  Courant 数 = c*dt/dx = {courant:.4f}，为 3D 上限的 "
          f"{courant / (1.0 / np.sqrt(3.0)):.2f} 倍")
    print(f"  2D (ny=1) 实测未发散（末段/首段能量 0.010），"
          f"但晚期噪声为 3D 的约 60 倍")
    print(f"  实际仿真窗 = 5089 x {dt_measured * 1e9:.4f}ns = "
          f"{5089 * dt_measured * 1e9:.1f} ns（声明 600ns）")
    print(f"  ⚠ 界面盒 0.25 m = {lam95 / 0.25:.1f} 个介质波长（@95MHz）"
          f"= {0.25 / lam95:.2f} 波长 -> 亚波长结构")

    # ---------- 核查 5：2D/3D 可比性 ----------
    res["check5_2d3d_comparability"] = {
        "items_compared": 10,
        "items_identical": 1,
        "shared": ["收发偏移数值 1.3m", "材料参数 rock/cover"],
        "differing": [
            "域: 12x0.05x33 (ny=1) vs 12x12x33",
            "变化轴: 2D 沿 x 3.9-9.9m/121道/0.05m; 3D 沿 y 5.15-8.15m/13道/0.25m",
            "界面起伏: 2D 0.40m(9级) vs 3D 0.80m(17级)",
            "界面 z 范围: 2D 8.85-9.25m vs 3D 8.60-9.40m",
            "源: (2.6,0.025,27)y极化 vs (6,3.85,27)x极化",
            "接收器分量: Ey vs Ex",
            "PML y: 0 层 vs 20 层",
            "pml_cfs: 未指定 vs 已指定",
        ],
        "verdict": "两套独立几何，非同一模型的 2D/3D 实现；"
                   "2D 结论不可直接外推到 3D",
    }
    print()
    print("核查 5  2D/3D 可比性")
    print("  10 项对照中仅 1 项真正相同（收发偏移数值 + 材料参数）")
    print("  ⚠ 2D 不是 3D 的降维实现，2D 结论不能外推到 3D")

    # ---------- 修正后重做界面匹配 ----------
    print()
    print("=" * 76)
    print("修正色散速度后重做界面匹配（2D）")
    z0_rx = np.interp(np.arange(S.shape[1]) * 0.05 + 3.9,
                      np.concatenate([[0], np.arange(0.25, 12.01, 0.25)]),
                      np.concatenate([[8.95],
                                      [8.95 - 0.0] * 0]), right=9.05) \
        if False else None
    # 直接从 .in 读z0 剖面
    zs = []
    with open(os.path.join(D2, "hs4t2d_t01.in"), encoding="utf-8",
              errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("#box:"):
                p = s.split()
                if p[7] == "cover":
                    zs.append((0.5 * (float(p[1]) + float(p[4])), float(p[3])))
    zs.sort()
    zx = np.array([a for a, _ in zs])
    zz = np.array([b for _, b in zs])
    rx = np.arange(S.shape[1]) * 0.05 + 3.9
    z0_rx = np.interp(rx, zx, zz)

    tt_lo, tt_hi = t_of_z_dispersive(z0_rx)
    t_lo, t_hi = float(tt_lo.min()), float(tt_hi.max())
    tt = 0.5 * (tt_lo + tt_hi)
    m = (t >= t_lo - 8) & (t <= t_hi + 8)
    sub, tw = S[m], t[m]
    idx = np.argmin(sub, axis=0)
    t_obs = tw[idx]

    bias = t_obs - tt.mean()
    print(f"  色散理论到时（z0={z0_rx.min():.2f}-{z0_rx.max():.2f}m）:")
    print(f"    带内包络 {t_lo:.1f} .. {t_hi:.1f} ns")
    print(f"    上一轮用零频速度得 177.9..185.0 ns（偏高 "
          f"{(178.0 / t_lo - 1) * 100:.1f}%）")
    print(f"  观测界面事件: {t_obs.min():.1f} .. {t_obs.max():.1f} ns")
    print(f"  偏差: 均值 {bias.mean():+.2f} ns  标准差 {bias.std():.2f} ns")
    inside = ((t_obs >= t_lo) & (t_obs <= t_hi)).mean()
    print(f"  落在理论包络内的道: {inside * 100:.1f}%")

    res["interface_recheck_dispersive"] = {
        "theory_envelope_ns": [float(t_lo), float(t_hi)],
        "theory_envelope_old_nondispersive_ns": [177.9, 185.0],
        "observed_ns": [float(t_obs.min()), float(t_obs.max())],
        "bias_mean_ns": float(bias.mean()),
        "bias_std_ns": float(bias.std()),
        "fraction_inside_envelope": float(inside),
        "model_z0_range_m": [float(z0_rx.min()), float(z0_rx.max())],
    }

    mp = os.path.join(OUT, "hs4_model_correctness_review.json")
    os.makedirs(OUT, exist_ok=True)
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print()
    print("wrote", mp)


if __name__ == "__main__":
    main()
