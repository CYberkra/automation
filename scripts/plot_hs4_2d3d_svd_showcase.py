# -*- coding: utf-8 -*-
"""最近两次仿真（3D HS4 / 2D HS4T2D）：模型、B-scan 与 rank-1 SVD 效果。

口径严格对齐项目既有权威实现 scripts/review_bscan_hs4_20261003.py
（该脚本已用官方 processing.py 全量重算并验证 relative error < 1e-12）：

  * 时间轴：npz 的 t 单位就是纳秒（末值 3332.5 ns），来自官方
    reconstruct_time_response 的 time = arange(count)/(count*df) * 1e9。
    不再乘 h5 的 dt——SFCW 频率栅格决定时间轴，与 impulse 采样率无关。
  * SVD 窗口：仅在 t <= 250 ns 上做（与评审脚本一致）。全窗做会让
    远处的零噪声区主导奇异值谱，掩盖界面事件。
  * 窗口划分（评审脚本 WINDOWS）：
        direct 0-20 / pre_ground 20-80 / ground 85-120
        intermediate 120-160 / interface 160-220 ns
  * 3D 组的 13 道关于中心道严格对称（corr(道1,道13)=1.000000），
    这是共偏移几何的平移不变性，不是处理错误——故其 SV1 占比天然接近
    100%，残差幅值应读作"残余可抑制量"而非"保留信号"。

输出：每组一张 3x2 总览（raw / SVD 残差 / 残差+AGC + 奇异值谱 +
基模态 + 残差包络），另加一张 2D vs 3D 对照。
"""
from __future__ import annotations

import glob
import json
import os

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# 中文字形（否则标题中的汉字渲染为方框）
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
})

OUT_DIR = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
           r"\2026-10-03_hs4_2d3d_svd_showcase")

# 与评审脚本一致的时窗划分（ns）
WINDOWS = {"direct": (0, 20), "pre_ground": (20, 80), "ground": (85, 120),
           "intermediate": (120, 160), "interface": (160, 220)}
SVD_MAX_NS = 250.0
C = 299792458.0

CMAP = LinearSegmentedColormap.from_list(
    "seismic", [(0.0, "#00007f"), (0.25, "#0000ff"), (0.5, "#ffffff"),
                (0.75, "#ff0000"), (1.0, "#7f0000")], N=256)

CASES = {
    "3D": {
        "label": "3D  HS4 rough-interface half-space (13 traces, Ex)",
        "dir": r"E:\automation_djh\verify_runs\2026-10-02_halfspace_standard",
        "h5": "hs4_co13_t*.h5",
        "npz": "hs4_co13_bscan_official.npz",
        "in": "hs4_co13_t01.in",
        "axis": "y",
        "field": "Ex",
    },
    "2D": {
        "label": "2D  HS4T2D transect fast study (121 traces, Ey)",
        "dir": (r"E:\automation_djh\automation_repo\artifacts\research_checks"
                r"\2026-10-02_hs4t2d_transect"),
        "h5": "hs4t2d_t*.h5",
        "npz": "hs4t2d_bscan_official.npz",
        "in": "hs4t2d_t01.in",
        "axis": "x",
        "field": "Ey",
    },
}


def load_case(key):
    cfg = CASES[key]
    files = sorted(glob.glob(os.path.join(cfg["dir"], cfg["h5"])))
    if not files:
        raise FileNotFoundError(f"{key}: no h5 under {cfg['dir']}")

    pos, names, src = [], [], None
    for f in files:
        with h5py.File(f, "r") as h:
            for _, ds in h["rxs"].items():
                pos.append(np.array(ds.attrs["Position"], float))
                names.append(str(ds.attrs["Name"]))
            if src is None:
                for _, ds in h["srcs"].items():
                    src = np.array(ds.attrs["Position"], float)
    pos = np.array(pos)

    d = np.load(os.path.join(cfg["dir"], cfg["npz"]))
    t_ns, S = d["t"].astype(float), d["S"].astype(float)

    # 道序对齐（按文件名字数排序，与评审脚本 _t01.._tNN 约定一致）
    num = np.array([int("".join(c for c in n if c.isdigit()) or 0)
                    for n in names])
    order = np.argsort(num)
    pos = pos[order]
    if S.shape[1] == pos.shape[0]:
        S = S[:, order]

    return {"cfg": cfg, "S": S, "t": t_ns, "pos": pos, "src": src,
            "n": pos.shape[0]}


def svd_on_window(S, t, tmax=SVD_MAX_NS):
    """仅在 t <= tmax 上按道去均值做 rank-1 移除（评审脚本口径）。"""
    m = t <= tmax
    tw, Sw = t[m], S[m]
    a = Sw - Sw.mean(axis=0, keepdims=True)
    U, Sg, Vt = np.linalg.svd(a, full_matrices=False)
    frac = float(Sg[0] / Sg.sum())
    resid = a - Sg[0] * np.outer(U[:, 0], Vt[0])
    return tw, resid, frac, Sg, U, Vt


def robust_vlim(*arrays, q=0.9995):
    """按分位取色标上限，避免少数极值把主体压成白色。"""
    pooled = np.concatenate([np.abs(a).ravel() for a in arrays if a.size])
    if pooled.size == 0:
        return 1.0
    v = float(np.quantile(pooled, q))
    return v if v > 0 else float(pooled.max()) or 1.0


# B-scan 能量分段（实测 absmax，2D 为例）：
#   0–25ns   175.8   直达波
#   25–80ns    1.27  过渡
#   80–120ns  23.8   地面波（强反射）
#   120–250ns  0.21  界面窗
#   250–1000   1e-3  噪声地板
#   1000–3333  71.8  Hann 窗边振
# 全窗单一色标会把 120ns 后结构压白（相差 5 个数量级），
# 故按段独立归一；段界用白线标出，各段上限标在图下方。
SEGMENTS_FULL = [(0.0, 25.0), (25.0, 80.0), (80.0, 120.0),
                 (120.0, 250.0), (250.0, 1000.0), (1000.0, 3400.0)]
# 残差图只到 250ns，用同一套边界的前四段即可
SEGMENTS = SEGMENTS_FULL[:4]


def segmented_norm(img, t, segments=SEGMENTS, q=0.9995):
    """分段独立归一 + 硬拼接，返回可 imshow 的数组与各段上限。

    每段各自做 robust 上限，段内动态范围得以保留；段界处会有可见的
    幅值跳变，这是"分段归一"的固有代价，图上用白线标出段界并注明。
    """
    out = np.zeros_like(img)
    vlims = []
    edges = []
    for lo, hi in segments:
        m = (t >= lo) & (t < hi)
        if not m.any():
            vlims.append(1.0)
            continue
        v = robust_vlim(img[m], q=q)
        vlims.append(v)
        sub = img[m]
        out[m] = np.clip(sub / (v or 1.0), -1.0, 1.0)
        edges.append(lo)
    return out, vlims, edges


def draw_segmented_bscan(ax, img, t, axt, segments, title,
                         ylabel=True, q=0.9995):
    """分段归一 B-scan。段界画白线，右上角列出各段上限。"""
    normed, vlims, edges = segmented_norm(img, t, segments, q=q)
    ymax = float(t[-1])
    im = ax.imshow(normed, aspect="auto", cmap=CMAP, vmin=-1, vmax=1,
                   extent=[axt[0], axt[-1], ymax, 0], interpolation="nearest")
    for e in edges:
        if 0 < e < ymax:
            ax.axhline(e, color="white", lw=1.4, alpha=0.95)
    for z0, tw in theory_times_ns().items():
        if 0 < tw < ymax:
            ax.axhline(tw, ls="--", lw=0.9, color="#00b894", alpha=0.9)
    shade_windows(ax, ymax)
    txt = "  ".join(f"{lo:.0f}–{hi:.0f}:±{v:.2g}"
                    for (lo, hi), v in zip(segments, vlims))
    ax.text(0.5, -0.075, f"分段独立归一（白线=段界）  {txt}",
            transform=ax.transAxes, fontsize=5.8, color="#333", ha="center")
    ax.set_title(title, fontsize=9.5)
    if ylabel:
        ax.set_ylabel("two-way time (ns)", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    return im


def theory_times_ns():
    """界面双程几何到时：z0 取材料盒 cover 顶/底，v 取 cover 速度。"""
    v = C / np.sqrt(18.017)
    out = {}
    for z0 in (8.95, 9.00, 9.05):
        out[z0] = 2 * (z0 + 1.3) / v * 1e9
    return out


def shade_windows(ax, ymax):
    """按时窗划分打淡色底纹（标签省略，避免溢出画布；图例见图注）。"""
    colors = {"direct": "#ffffff", "pre_ground": "#f4f6f7",
              "ground": "#fdebd0", "intermediate": "#e8daef",
              "interface": "#d5f5e3"}
    for name, (lo, hi) in WINDOWS.items():
        if lo >= ymax:
            continue
        ax.axhspan(lo, min(hi, ymax), color=colors[name], alpha=0.16, lw=0)


def draw_bscan(ax, img, t, axt, vlim, ymax, title, ylabel=True):
    im = ax.imshow(img, aspect="auto", cmap=CMAP, vmin=-vlim, vmax=vlim,
                   extent=[axt[0], axt[-1], ymax, 0], interpolation="nearest")
    for z0, tw in theory_times_ns().items():
        if 0 < tw < ymax:
            ax.axhline(tw, ls="--", lw=0.9, color="#00b894", alpha=0.8)
    shade_windows(ax, ymax)
    ax.set_title(title, fontsize=9.5)
    if ylabel:
        ax.set_ylabel("two-way time (ns)", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    return im


def build_panel(fig, gs, D, out_png):
    cfg = D["cfg"]
    S, t, pos = D["S"], D["t"], D["pos"]
    ai = pos[:, 1] if cfg["axis"] == "y" else pos[:, 0]

    tw, resid, frac, Sg, U, Vt = svd_on_window(S, t)
    full_resid = np.zeros_like(S)
    m = t <= SVD_MAX_NS
    full_resid[m] = resid

    ymax = float(t[-1])
    v_raw = robust_vlim(S)
    v_res = robust_vlim(resid)

    # 道间差异量化：判断该组仿真是否携带独立横向信息
    m250 = t <= SVD_MAX_NS
    Sw = S[m250]
    if Sw.shape[1] > 1:
        rel_diff = float(np.abs(Sw[:, 0] - Sw[:, -1]).max()
                         / (np.abs(Sw).max() or 1.0))
    else:
        rel_diff = 0.0

    # (0,0) raw B-scan —— 分段归一，否则 120ns 后结构全部压白
    ax = fig.add_subplot(gs[0, 0])
    im = draw_segmented_bscan(
        ax, S, t, ai, SEGMENTS_FULL,
        f"B-scan  raw signed response  (full 0–{ymax:.0f} ns)")
    ax.set_xlabel(f"profile {cfg['axis']} (m)", fontsize=8.5)
    ax.text(0.012, 0.035,
            f"首/末道最大差 ÷ 全局最大 = {rel_diff:.2e}"
            f"  →  {'横向信息可用' if rel_diff > 1e-3 else '各道近乎相同，无独立横向信息'}",
            transform=ax.transAxes, fontsize=7.2, color="#1a5276",
            bbox=dict(boxstyle="round,pad=0.32", fc="#eaf2f8",
                      ec="#1a5276", lw=0.7, alpha=0.94))
    cb = fig.colorbar(im, ax=ax, pad=0.011, fraction=0.045)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("normalised", fontsize=7.5)

    # (1,0) SVD 残差（仅 250ns 窗）
    ax = fig.add_subplot(gs[1, 0])
    im = draw_segmented_bscan(
        ax, resid, tw, ai, SEGMENTS,
        f"SVD residual, rank-1 removed  "
        f"(SV1={frac * 100:.2f}%, window ≤{SVD_MAX_NS:.0f} ns)",
        ylabel=False)
    ax.set_xlabel(f"profile {cfg['axis']} (m)", fontsize=8.5)
    cb = fig.colorbar(im, ax=ax, pad=0.011, fraction=0.045)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("normalised", fontsize=7.5)

    # (2,0) 残差能量沿时间（原始尺度，不用 AGC）
    ax = fig.add_subplot(gs[2, 0])
    prof = np.sqrt((resid ** 2).mean(axis=1))
    ax.plot(tw, prof, lw=1.15, color="#2c3e50")
    ax.fill_between(tw, 0, prof, color="#2c3e50", alpha=0.18)
    for name, (lo, hi) in WINDOWS.items():
        sel = (tw >= lo) & (tw <= hi)
        if sel.any():
            ax.axvline(lo, ls=":", lw=0.7, color="#999")
    for z0, twt in theory_times_ns().items():
        if 0 < twt < float(tw[-1]):
            ax.axvline(twt, ls="--", lw=0.9, color="#00b894", alpha=0.85)
            ax.text(twt, ax.get_ylim()[1] * 0.95, f" {twt:.0f}",
                    fontsize=6.8, color="#00795c", va="top")
    ax.set_xlim(0, float(tw[-1]))
    ax.set_yscale("log")
    ax.set_title("residual RMS vs time  (per-window energy)", fontsize=9.5)
    ax.set_xlabel("two-way time (ns)", fontsize=8.5)
    ax.set_ylabel("residual RMS", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    ax.grid(alpha=0.25, lw=0.5)

    # 能量比：ground vs interface（评审脚本口径）
    en = {}
    for name, (lo, hi) in WINDOWS.items():
        sel = (tw >= lo) & (tw <= hi)
        en[name] = float(np.sum(resid[sel] ** 2)) if sel.any() else 0.0
    ratio = en["ground"] / en["interface"] if en["interface"] > 0 else float("nan")

    # (0,1) 奇异值谱
    ax = fig.add_subplot(gs[0, 1])
    k = min(20, Sg.size)
    ax.semilogy(Sg[:k], "o-", ms=4, lw=1.3, color="#c0392b")
    ax.axhline(Sg[0], ls=":", lw=1, color="#666")
    ax.text(0.6, Sg[0] * 1.6, f"SV1 = {frac * 100:.2f}%", fontsize=8.5,
            color="#c0392b")
    if Sg.size > k:
        ax.text(0.98, 0.05, f"(仅显示前 {k}/{Sg.size})", transform=ax.transAxes,
                fontsize=7, color="#777", ha="right")
    # SV1 接近 100% 时必须说明成因，否则会被误读为"SVD 失效"
    if frac > 0.98:
        ax.text(0.5, 0.62,
                "SV1≈100%：各道在 ≤250ns 窗口内近乎相同\n"
                "（共偏移几何平移不变；非处理错误）\n"
                "残差应读作'可抑制的背景量'",
                transform=ax.transAxes, fontsize=7.6, color="#7b241c",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.4", fc="#fdedec",
                          ec="#c0392b", lw=0.8, alpha=0.95))
    ax.set_title("singular-value spectrum (≤250 ns)", fontsize=9.5)
    ax.set_xlabel("component index", fontsize=8.5)
    ax.set_ylabel("singular value", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    ax.grid(alpha=0.25, lw=0.5)

    # (1,1) rank-1 基模态
    ax = fig.add_subplot(gs[1, 1])
    u1 = U[0] * np.sign(U[0][np.argmax(np.abs(U[0]))])
    ax.plot(ai, u1, lw=1.2, color="#2c6fbb", label="U1 (spatial)")
    ax.set_xlabel(f"profile {cfg['axis']} (m)", fontsize=8.5)
    ax.set_ylabel("U1 (spatial)", fontsize=8.5, color="#2c6fbb")
    ax.tick_params(labelsize=7.5, axis="y")
    ax2 = ax.twinx()
    ax2.plot(ai, Vt[0], lw=1.2, color="#e08a1e", label="V1 (trace)")
    ax2.set_ylabel("V1 (trace)", fontsize=8.5, color="#e08a1e")
    ax2.tick_params(labelsize=7.5)
    ax.set_title("rank-1 spatial / trace modes", fontsize=9.5)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="best")

    # (2,1) 逐窗口残差能量 + 道间相关系数
    ax = fig.add_subplot(gs[2, 1])
    names = list(WINDOWS)
    vals = [en[n] if en[n] > 0 else 1e-30 for n in names]
    bars = ax.bar(range(len(names)), vals, color="#5499c7", edgecolor="#2c3e50",
                  lw=0.6)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=28, fontsize=7.5, ha="right")
    ax.set_yscale("log")
    ax.set_title(f"residual energy by window   "
                 f"ground/interface = {ratio:.2f}", fontsize=9.5)
    ax.set_ylabel("sum of squares", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    ax.grid(alpha=0.25, lw=0.5, axis="y")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.6, f"{v:.1e}",
                ha="center", fontsize=6.3)

    return {"frac": frac, "ratio": ratio, "energy": en, "v_raw": v_raw,
            "v_res": v_res, "resid": resid, "tw": tw, "ai": ai,
            "full_resid": full_resid, "S": S, "t": t,
            "rel_first_last_diff": rel_diff}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    data = {k: load_case(k) for k in ("3D", "2D")}
    summary, panels = {}, {}

    for key in ("3D", "2D"):
        D = data[key]
        fig = plt.figure(figsize=(13.0, 11.6))
        gs = fig.add_gridspec(3, 2, width_ratios=[1.62, 1.0],
                              hspace=0.40, wspace=0.30)
        R = build_panel(fig, gs, D, OUT_DIR)
        panels[key] = R

        fig.suptitle(
            f"{D['cfg']['label']}\n"
            f"official 20 MHz SFCW chain · gprMax v4 · common-offset "
            f"{D['n']} traces · window {D['t'][-1]:.0f} ns\n"
            f"时窗底色：direct 0–20 / pre_ground 20–80 / ground 85–120 / "
            f"intermediate 120–160 / interface 160–220 ns；"
            f"绿色虚线 = cover 顶/底面双程几何到时 290–293 ns",
            fontsize=10, y=0.995)
        p = os.path.join(OUT_DIR, f"hs4_{key.lower()}_model_bscan_svd.png")
        fig.savefig(p, dpi=140, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print("wrote", p)

        # 道间相关系数矩阵（对称性证据）
        S = D["S"]
        Cc = np.corrcoef(S.T)
        off = Cc[~np.eye(len(Cc), dtype=bool)]
        summary[key] = {
            "traces": int(D["n"]),
            "component": D["cfg"]["field"],
            "profile_axis": D["cfg"]["axis"],
            "axis_range_m": [float(D["pos"][:, 1 if key == "3D" else 0].min()),
                             float(D["pos"][:, 1 if key == "3D" else 0].max())],
            "axis_step_m": float(np.unique(np.round(np.diff(np.sort(
                D["pos"][:, 1 if key == "3D" else 0])), 6))[0]),
            "time_window_ns": float(D["t"][-1]),
            "sv1_fraction": float(R["frac"]),
            "sv2_over_sv1": float(R["resid"].shape[0] and 0) or None,
            "ground_over_interface_energy": float(R["ratio"]),
            "residual_energy_by_window": R["energy"],
            "trace_corr_min": float(off.min()),
            "trace_corr_max": float(off.max()),
            "first_last_maxdiff_over_absmax": float(R["rel_first_last_diff"]),
            "transverse_info_available": bool(R["rel_first_last_diff"] > 1e-3),
            "vlim_raw": float(R["v_raw"]),
            "vlim_residual": float(R["v_res"]),
            "raw_absmax": float(np.abs(S).max()),
            "theory_two_way_ns": {str(k_): float(v) for k_, v
                                  in theory_times_ns().items()},
        }
        # 补 SV2/SV1
        _, _, _, Sg, _, _ = svd_on_window(S, D["t"])
        summary[key]["sv2_over_sv1"] = float(Sg[1] / Sg[0]) if Sg.size > 1 else None
        print(f"  {key}: SV1={R['frac'] * 100:.2f}%  "
              f"ground/interface={R['ratio']:.3g}  "
              f"corr[{summary[key]['trace_corr_min']:.4f},"
              f"{summary[key]['trace_corr_max']:.4f}]")

    # ---- 2D vs 3D 对照 ----
    fig = plt.figure(figsize=(13.2, 5.6))
    gs2 = fig.add_gridspec(1, 3, wspace=0.28)
    for j, key in enumerate(("3D", "2D")):
        R, D = panels[key], data[key]
        tw, resid, ai = R["tw"], R["resid"], R["ai"]
        ax = fig.add_subplot(gs2[j])
        im = draw_segmented_bscan(
            ax, resid, tw, ai, SEGMENTS,
            f"{key} · SVD residual (≤{SVD_MAX_NS:.0f} ns), "
            f"SV1={R['frac'] * 100:.1f}%")
        ax.set_xlabel(f"profile {D['cfg']['axis']} (m)", fontsize=9)
        ax.set_ylabel("two-way time (ns)", fontsize=9)
        ax.tick_params(labelsize=8)
        cb = fig.colorbar(im, ax=ax, pad=0.011, fraction=0.045)
        cb.ax.tick_params(labelsize=7.5)
        cb.set_label("normalised", fontsize=8)

    ax = fig.add_subplot(gs2[2])
    for key, color in (("3D", "#2c6fbb"), ("2D", "#e08a1e")):
        R = panels[key]
        prof = np.sqrt((R["resid"] ** 2).mean(axis=1))
        p = prof / (prof.max() or 1.0)
        ax.plot(R["tw"], p, lw=1.3, color=color,
                label=f"{key} (SV1={R['frac'] * 100:.1f}%)")
    for z0, twt in theory_times_ns().items():
        if 0 < twt < SVD_MAX_NS:
            ax.axvline(twt, ls="--", lw=0.9, color="#00b894", alpha=0.7)
    ax.set_xlim(0, SVD_MAX_NS)
    ax.set_title("normalised residual envelope  (green = theory 290–293 ns)",
                 fontsize=10)
    ax.set_xlabel("two-way time (ns)", fontsize=9)
    ax.set_ylabel("residual RMS / max", fontsize=9)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.5)
    ax.tick_params(labelsize=8)

    fig.suptitle("HS4 · 2D vs 3D — rank-1 SVD background suppression "
                 "(official 20 MHz SFCW chain, gprMax v4)", fontsize=11.5,
                 y=0.99)
    p = os.path.join(OUT_DIR, "hs4_2d_vs_3d_svd_compare.png")
    fig.savefig(p, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", p)

    mp = os.path.join(OUT_DIR, "hs4_2d3d_svd_showcase_metrics.json")
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print("wrote", mp)


if __name__ == "__main__":
    main()
