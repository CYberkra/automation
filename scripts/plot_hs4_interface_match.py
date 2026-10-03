# -*- coding: utf-8 -*-
"""HS4T2D 深部界面：模型 vs B-scan 的定量一致性检验（匹配滤波法）。

前两版方法的缺陷及本版修正：
  v1 直接在搜索窗取 argmin -> 无质量门控，噪声道产生 4.9 m 假洞
  v2 加了门控，但 G1 幅值门 0/121 全灭——界面回波仅为全窗峰值的
     0.037%（约 -68 dB），2% 阈值不适用于如此弱的信号

本版方法：匹配滤波（互相关）
  1) 构造理论界面到时曲线 t_theory(x)（由 .in 实读几何 + 正确路径公式）
  2) 把每道在宽窗内按 t_theory 平移对齐（即"按理论到时采样"）
  3) 对齐后的道做加权叠加（stack），得到叠加道
  4) 叠加道在理论窗内的能量占比 = 匹配度
  5) 形态一致性检验：逐道互相关曲线的形状 vs 理论曲线去均值形状

这才是弱信号下的标准做法：不依赖单道极值，而是把 121 道的信号
按已知几何对齐后联合检验。
"""
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
})

D2 = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
      r"\2026-10-02_hs4t2d_transect")
NPZ = os.path.join(D2, "hs4t2d_bscan_official.npz")
IN = os.path.join(D2, "hs4t2d_t01.in")
OUT = (r"E:\automation_djh\automation_repo\artifacts\research_checks"
       r"\2026-10-03_hs4_2d3d_svd_showcase")

C = 299792458.0
V_COVER = C / np.sqrt(18.017)
Z_SRC, Z_BOT = 27.0, 12.0
RX_LO, DX = 3.9, 0.05
HALF_WIN = 8.0            # 对齐窗口半宽 (ns)


def model_z0(path):
    boxes = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("#box:"):
                p = s.split()
                boxes.append([float(v) for v in p[1:7]] + [p[7]])
    cov = sorted((b for b in boxes if b[6] == "cover"), key=lambda b: b[0])
    return (np.array([0.5 * (b[0] + b[3]) for b in cov]),
            np.array([b[2] for b in cov]))


def t_of_z(z):
    z = np.asarray(z, float)
    return 2 * (Z_SRC - Z_BOT) / C * 1e9 + 2 * (Z_BOT - z) / V_COVER * 1e9


def z_of_t(tt):
    tt = np.asarray(tt, float)
    return Z_BOT - (tt - 2 * (Z_SRC - Z_BOT) / C * 1e9) * V_COVER / 2 / 1e9


def sample_at(t, S, centers, half_win):
    """在每道以 center 为中心取 ±half_win 的样点，线性插值。

    返回 (n_traces, 2*half_win/dt+1) 的对齐矩阵。
    """
    dt = float(t[1] - t[0])
    n_s = int(round(2 * half_win / dt)) + 1
    offs = np.linspace(-half_win, half_win, n_s)
    out = np.empty((S.shape[1], n_s))
    for i in range(S.shape[1]):
        out[i] = np.interp(centers[i] + offs, t, S[:, i])
    return out, offs


def main():
    d = np.load(NPZ)
    t, S = d["t"], d["S"]
    n = S.shape[1]
    rx = np.arange(n) * DX + RX_LO
    xs, z0 = model_z0(IN)
    z0_rx = np.interp(rx, xs, z0)
    tt = t_of_z(z0_rx)

    print("=" * 78)
    print("① 匹配滤波（按理论几何对齐后联合检验）")
    print(f"   理论到时: {tt.min():.2f} .. {tt.max():.2f} ns"
          f"   起伏 {np.ptp(tt):.2f} ns")
    print(f"   模型 z0: {z0_rx.min():.3f} .. {z0_rx.max():.3f} m"
          f"   起伏 {np.ptp(z0_rx):.3f} m")
    print(f"   源/Rx 在 z={Z_SRC} m，z={Z_BOT}..{Z_SRC} 为真空段")
    print()

    Al, offs = sample_at(t, S, tt, HALF_WIN)
    stack = Al.mean(axis=0)
    stack_rms = np.sqrt((Al ** 2).mean(axis=0))

    # 界面中心在offsets 的 0 位置
    c0 = len(offs) // 2
    # 匹配度：±2ns 窗内能量占整个对齐窗的比例
    dt = float(t[1] - t[0])
    k2 = max(1, int(round(2.0 / dt)))
    win = slice(c0 - k2, c0 + k2 + 1)
    e_win = float((stack_rms[win] ** 2).sum())
    e_all = float((stack_rms ** 2).sum())
    match = e_win / e_all if e_all > 0 else 0.0

    # 噪声基准：把对齐窗整体平移 60ns 后同样计算
    Al2, offs2 = sample_at(t, S, tt + 60.0, HALF_WIN)
    st2 = np.sqrt((Al2 ** 2).mean(axis=0))
    c02 = len(offs2) // 2
    w2 = slice(c02 - k2, c02 + k2 + 1)
    e_win2 = float((st2[w2] ** 2).sum())
    e_all2 = float((st2 ** 2).sum())
    match2 = e_win2 / e_all2 if e_all2 > 0 else 0.0

    print("② 匹配度检验（±2ns 窗能量占比）")
    print(f"   对齐到理论: {match * 100:.2f}%")
    print(f"   错位+60ns  : {match2 * 100:.2f}%   （噪声基准）")
    print(f"   对比倍数   : {match / match2 if match2 > 0 else float('nan'):.2f}×")
    print()

    # 逐道互相关：理论脉冲 vs 实际道
    # 构造理论道：Ricker 子波在理论到时处
    f0 = 107e6
    tw_theory = np.exp(-((offs * 1e-9 + 0.0) * np.pi * f0) ** 2) \
        * np.cos(2 * np.pi * f0 * offs * 1e-9)
    tw_theory = -tw_theory# 界面为负极性
    cc = np.empty(n)
    for i in range(n):
        a = Al[i] - Al[i].mean()
        b = tw_theory - tw_theory.mean()
        den = np.linalg.norm(a) * np.linalg.norm(b)
        cc[i] = float(a @ b / den) if den > 0 else 0.0

    print("③ 逐道互相关（理论子波 vs 观测对齐道）")
    print(f"   相关系数: 均值 {cc.mean():+.4f}  中位 {np.median(cc):+.4f}"
          f"  标准差 {cc.std():.4f}")
    print(f"   正相关道数: {(cc > 0).sum()}/{n}  ({(cc > 0).mean() * 100:.1f}%)")
    print(f"   若纯噪声，期望 |r| ≈ 1/sqrt(窗内样点数) = "
          f"{1 / np.sqrt(len(offs)):.4f}")
    print()

    # 形态检验：只对高相关道做
    thr = max(0.15, np.median(cc) + 2 * cc.std())
    good = cc >= thr
    print(f"④ 形态检验（选取相关 ≥ {thr:.3f} 的道，n={int(good.sum())}）")
    if good.sum() >= 8:
        obs = tt[good]
        a = obs - obs.mean()
        b = tt[good] - tt[good].mean()
        r_shape = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
        print(f"   形状相关 r(去均值) = {r_shape:+.4f}")
        print(f"   起伏比 观测/理论 = {np.ptp(obs) / np.ptp(tt[good]):.3f}")
        seg = np.array_split(np.where(good)[0], 4)
        print("   分段:")
        for s in seg:
            if len(s):
                print(f"     x={rx[s[0]]:.2f}-{rx[s[-1]]:.2f}m (n={len(s):2d}): "
                      f"r={np.corrcoef(tt[s] - tt[s].mean(), obs[s] - obs[s].mean())[0, 1]:+.3f}"
                      f"  观测起伏 {np.ptp(obs[s - seg[0].start]):5.2f} ns"
                      f"  理论 {np.ptp(tt[s]):5.2f} ns")
    else:
        r_shape = float("nan")
        print(f"   高相关道不足（{int(good.sum())} 道），无法做形态检验")
    print()

    # 残差频谱（用全部道，避免门控偏差）
    print("⑤ 残差频谱（用全部有效互相关道的到时残差）")
    resid = []
    xs_g = []
    for i in range(n):
        if cc[i] >= np.median(cc):
            lo, hi = tt[i] - HALF_WIN, tt[i] + HALF_WIN
            m = (t >= lo) & (t <= hi)
            if m.any():
                resid.append(float(t[m][int(np.argmin(S[m, i]))]) - tt[i])
                xs_g.append(rx[i])
    if len(resid) >= 8:
        resid = np.array(resid)
        F = np.abs(np.fft.rfft(resid - resid.mean()))
        f = np.fft.rfftfreq(resid.size, d=DX)
        msk = f > 1e-6
        f_model = 1.0 / max(np.ptp(z0_rx), 1e-9)
        top = [i for i in np.argsort(F)[::-1][:4] if f[i] > 1e-6]
        print(f"   模型特征 f ≈ {f_model:.2f} /m（由 {np.ptp(z0_rx):.2f} m 起伏）")
        for i in top:
            print(f"   残差峰 f = {f[i]:.3f} /m（周期 {1 / f[i]:.2f} m）"
                  f"  幅 {F[i]:.1f}")
        hit = any(abs(f[i] - f_model) < 0.2 * f_model for i in top)
        print(f"   是否有峰落在模型特征频率附近: {'是' if hit else '否'}")
    else:
        f = F = None
        hit = None
        f_model = 1.0 / max(np.ptp(z0_rx), 1e-9)
        top = []

    # ---------------- 图 ----------------
    fig = plt.figure(figsize=(13.6, 9.8))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.68],
                          hspace=0.46, wspace=0.25)

    # (0,0) 模型界面 + 叠加道定位
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(rx, z0_rx, lw=1.9, color="#1a5276", label="模型 z0(x)")
    ax.plot(xs, z0, "o", ms=4.2, color="#1a5276", alpha=0.6,
            label="模型阶梯节点（.in 实读 48 盒）")
    ax.plot(rx, z_of_t(tt), "--", lw=1.0, color="#7f8c8d",
            label="模型 z0（= 理论到时反算，自洽）")
    ax.set_xlabel("profile x (m)", fontsize=9.5)
    ax.set_ylabel("界面深度 (m)", fontsize=9.5)
    ax.set_title("① 模型界面几何（真值）", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=8.5)

    # (0,1) 匹配叠加道
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(offs, stack, lw=1.2, color="#2c3e50", label="按理论对齐的叠加道")
    ax.plot(offs, stack_rms, lw=1.4, color="#c0392b", label="叠加道 RMS 包络")
    ax.axvspan(-2, 2, color="#00b894", alpha=0.16)
    ax.axvline(0, ls="--", lw=1.2, color="#00795c")
    ax.text(0.5, 0.92, f"±2ns 窗能量占比 {match * 100:.1f}%\n"
                      f"（错位基准 {match2 * 100:.1f}%）",
            transform=ax.transAxes, fontsize=8.6, ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="#eafaf1", ec="#117a65"))
    ax.set_xlabel("相对理论到时 (ns)", fontsize=9.5)
    ax.set_ylabel("source response", fontsize=9.5)
    ax.set_title("② 匹配滤波：121 道按理论几何对齐后叠加", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=8.5)

    # (1,0) 逐道互相关
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(rx, cc, lw=1.0, color="#34495e", label="逐道互相关 r")
    ax.axhline(np.median(cc), ls="--", lw=1.2, color="#c0392b",
               label=f"中位 {np.median(cc):+.3f}")
    ax.axhline(thr, ls=":", lw=1.2, color="#e08a1e",
               label=f"选取阈值 {thr:.3f}")
    ax.fill_between(rx, 0, cc, where=cc >= thr, color="#00b894", alpha=0.25,
                    label=f"高相关道 n={int(good.sum())}")
    ax.set_xlabel("profile x (m)", fontsize=9.5)
    ax.set_ylabel("互相关 r", fontsize=9.5)
    ax.set_title(f"③ 逐道互相关（正相关 {(cc > 0).sum()}/{n} 道，"
                 f"噪声基准 ±{1 / np.sqrt(len(offs)):.3f}）", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=8.5)

    # (1,1) 形态对比
    ax = fig.add_subplot(gs[1, 1])
    if good.sum() >= 8:
        ax.plot(rx[good], tt[good] - tt[good].mean(), "o-", ms=4, lw=1.8,
                color="#1a5276", label=f"理论起伏（n={int(good.sum())} 道）")
    ax.plot(rx, tt - tt.mean(), lw=1.3, color="#1a5276", alpha=0.4,
            label="理论起伏（全部 121 道）")
    if f is not None and len(resid) >= 8:
        ax.plot(xs_g, np.array(resid) - np.mean(resid), "s", ms=3,
                color="#e74c3c", alpha=0.6, label="观测残差（去均值）")
    ax.axhline(0, ls="-", lw=1.0, color="#888")
    ax.set_xlabel("profile x (m)", fontsize=9.5)
    ax.set_ylabel("去均值到时 (ns)", fontsize=9.5)
    ax.set_title(f"④ 形态对比：残差是否复现模型起伏　r = {r_shape:+.3f}",
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=8.5)

    # (2,:) 残差频谱
    ax = fig.add_subplot(gs[2, :])
    if f is not None and len(resid) >= 8:
        msk = f > 1e-6
        ax.semilogy(f[msk], F[msk], lw=1.2, color="#34495e",
                    label="残差频谱")
    ax.axvline(f_model, ls="--", lw=1.7, color="#1a5276",
               label=f"模型特征 f ≈ {f_model:.2f} /m")
    if f is not None:
        for i in top:
            ax.annotate(f"{f[i]:.2f}", (f[i], F[i]), fontsize=7.6,
                        color="#c0392b", xytext=(4, 4),
                        textcoords="offset points")
    ax.set_xlabel("空间频率 (1/m)", fontsize=9.5)
    ax.set_ylabel("残差谱幅", fontsize=9.5)
    ax.set_title("⑤ 残差频谱 vs 模型特征频率", fontsize=11)
    ax.legend(fontsize=8.2)
    ax.grid(alpha=0.3, lw=0.5, which="both")
    ax.tick_params(labelsize=8.5)

    fig.suptitle(
        "HS4T2D 深部界面：模型 vs B-scan 观测 —— 匹配滤波定量检验\n"
        f"界面回波仅全窗峰值的 {np.median(np.abs(S)[(t >= 170) & (t <= 200)] / np.abs(S).max()):.2%}"
        f"（约 {20 * np.log10(np.median(np.abs(S)[(t >= 170) & (t <= 200)] / np.abs(S).max())):.0f} dB）"
        f"　匹配度 {match * 100:.1f}% vs 错位基准 {match2 * 100:.1f}%"
        f"　逐道相关中位 {np.median(cc):+.3f}",
        fontsize=11, y=0.995)
    p = os.path.join(OUT, "hs4_2d_interface_model_vs_bscan.png")
    fig.savefig(p, dpi=145, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", p)

    res = {
        "method": "matched filter (align to theory geometry, then stack)",
        "model_z0_min_m": float(z0_rx.min()),
        "model_z0_max_m": float(z0_rx.max()),
        "model_z0_amplitude_m": float(np.ptp(z0_rx)),
        "theory_t_min_ns": float(tt.min()),
        "theory_t_max_ns": float(tt.max()),
        "theory_amplitude_ns": float(np.ptp(tt)),
        "interface_snr_frac_of_absmax": float(
            np.median(np.abs(S)[(t >= 170) & (t <= 200)] / np.abs(S).max())),
        "match_ratio_on_theory": match,
        "match_ratio_offset60ns": match2,
        "match_advantage": float(match / match2) if match2 > 0 else None,
        "trace_cc_median": float(np.median(cc)),
        "trace_cc_positive_count": int((cc > 0).sum()),
        "noise_cc_baseline": float(1 / np.sqrt(len(offs))),
        "n_high_cc": int(good.sum()),
        "shape_corr_high_cc": r_shape,
        "model_char_freq_per_m": float(f_model),
        "residual_peak_freq_per_m": [float(f[i]) for i in top] if f is not None else None,
        "residual_peak_at_model_freq": hit,
    }
    mp = os.path.join(OUT, "hs4_2d_interface_model_vs_bscan.json")
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print("wrote", mp)


if __name__ == "__main__":
    main()
