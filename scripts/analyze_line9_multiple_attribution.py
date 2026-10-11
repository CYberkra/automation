"""A1+A2 多次波/侧偏归因分析（纯分析，不调用求解器）。

对 2026-10-11 夜间 41 道批（AGL 8m，契约 ed9791f1...）的官方 direct/Hann B-scan：
A1 逐道叠加几何预测的首次覆盖层返回、垂直二次往返（多次波）与底砂到时曲线，
  并拾取强斜带脊线，给出到时残差与路径比（垂直多次波路径比=2.0）。
A2 邻道包络互相关测强斜带斜率，与几何预测的首次返回斜率及其 2 倍（多次波）对照。

所有拾取窗口、足迹、显示裁剪均为下方预声明常量；结果只表征当前合成模型。
"""
import argparse, json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

EPS0 = 8.8541878128e-12
C0 = 299792458.0
F_REF = 95e6

# ---- 预声明常量 ----
PRIMARY_WINDOW_NS = (100.0, 280.0)      # 首次覆盖层返回拾取窗（覆盖170/190m已知的148/196ns）
BAND_LEAD_NS = 30.0                     # 强带拾取窗起点 = 首次拾取 + 30ns
BAND_SPAN_NS = 160.0                    # 强带拾取窗宽度
INTERFACE_FOOTPRINT_M = 0.5             # 每站界面厚度取 rx 列两侧各 0.5m 平均
DISPLAY_TMIN_NS = 80.0                  # 显示下限，避开直达/地表首峰大动态
DISPLAY_TMAX_NS = 600.0
DISPLAY_CLIP_PERCENTILE = 98.0
SLOPE_HALF_WIDTH = 2                    # 斜率中心差分半宽（站）
DX = 0.025


def group_index(material, f=F_REF):
    base = material['base']; sigma = base['electric_conductivity_s_per_m']
    er = base['relative_permittivity'] - 1j * sigma / (2 * np.pi * f * EPS0)
    derivative = 1j * sigma / (2 * np.pi * EPS0 * f * f)
    for pole in material.get('poles', []):
        dd = pole['relative_permittivity_difference']; tau = pole['relaxation_time_s']
        den = 1 + 2j * np.pi * f * tau
        er += dd / den; derivative -= 2j * np.pi * tau * dd / den ** 2
    n = np.sqrt(er)
    return float((n + f * derivative / (2 * n)).real)


def interfaces(grid):
    """每列表面 z、覆盖层/泥岩界面 z、泥岩/砂岩界面 z（取每种材料最顶部）。"""
    nx, ny = grid.shape
    surf = np.full(nx, np.nan); i1 = np.full(nx, np.nan); i2 = np.full(nx, np.nan)
    for c in range(nx):
        col = grid[c, :]
        top = {m: (np.flatnonzero(col == m).max() if np.any(col == m) else -1) for m in (1, 2, 3)}
        if top[1] >= 0:
            surf[c] = (top[1] + 1) * DX
        if 0 <= top[2] < top[1]:
            i1[c] = (top[2] + 1) * DX
        if 0 <= top[3] < top[2]:
            i2[c] = (top[3] + 1) * DX
    return surf, i1, i2


def pick_parabola(t, env, lo, hi):
    ids = np.flatnonzero((t >= lo) & (t <= hi))
    if len(ids) < 3:
        return None, None
    j = ids[np.argmax(env[ids])]
    if j <= ids[0] or j >= ids[-1]:
        return float(t[j]), float(env[j])
    y0, y1, y2 = env[j - 1], env[j], env[j + 1]
    denom = (y0 - 2 * y1 + y2)
    frac = 0.5 * (y0 - y2) / denom if denom > 0 else 0.0
    dt = t[1] - t[0]
    return float(t[j] + frac * dt), float(y1 - 0.25 * (y0 - y2) * frac)


def xcorr_lag(e0, e1, dt, max_lag_s):
    n = min(len(e0), len(e1))
    a = e0[:n] - e0[:n].mean(); b = e1[:n] - e1[:n].mean()
    cc = np.correlate(b, a, mode='full')
    lags = np.arange(-n + 1, n) * dt
    m = np.abs(lags) <= max_lag_s
    cc = cc[m]; lags = lags[m]
    j = int(np.argmax(cc))
    if 0 < j < len(cc) - 1:
        y0, y1, y2 = cc[j - 1], cc[j], cc[j + 1]
        denom = (y0 - 2 * y1 + y2)
        frac = 0.5 * (y0 - y2) / denom if denom > 0 else 0.0
    else:
        frac = 0.0
    return float(lags[j] + frac * dt), float(cc[j])


def main(a):
    out = a.out; out.mkdir(parents=True, exist_ok=False)
    mats = json.loads(Path(a.materials).read_text('utf-8'))['materials']
    ng_cover = group_index(mats['material_001_cover'])
    ng_mud = group_index(mats['material_002_mudstone'])

    with h5py.File(a.geometry, 'r') as f:
        grid = np.asarray(f['data'][:, :, 0])
    surf, i1, i2 = interfaces(grid)

    contract = json.loads(Path(a.contract).read_text('utf-8'))
    rx_by_x = {g['profile_x_m']: g['rx_m'] for g in contract['groups'] if g['role'] == 'H1'}

    with h5py.File(a.sfcw, 'r') as f:
        t = np.asarray(f['time_s']) * 1e9
        env = np.abs(np.asarray(f['hann/H1_complex_bandpass']))
        real = np.asarray(f['hann/H1_real_bandpass'])
        xs = np.asarray(f['profile_x_m'])

    fp = int(round(INTERFACE_FOOTPRINT_M / DX))
    rows = []
    for k, px in enumerate(xs):
        rx = rx_by_x.get(float(px))
        col = int(round((rx[0] if rx else px - 19.35) / DX))
        lo, hi = max(0, col - fp), min(grid.shape[0], col + fp + 1)
        sz = np.nanmean(surf[lo:hi]); z1 = np.nanmean(i1[lo:hi]); z2 = np.nanmean(i2[lo:hi])
        d1 = sz - z1; d2 = z1 - z2
        agl = (rx[1] - sz) if rx else 8.0
        sep = 1.3
        t_air = 2 * np.sqrt(agl ** 2 + (sep / 2) ** 2) / C0 * 1e9
        t_prim = t_air + 2 * d1 * ng_cover / C0 * 1e9
        t_mult = t_air + 4 * d1 * ng_cover / C0 * 1e9
        t_basal = t_prim + 2 * d2 * ng_mud / C0 * 1e9
        tp, ap = pick_parabola(t, env[:, k], *PRIMARY_WINDOW_NS)
        tb, ab = pick_parabola(t, env[:, k], tp + BAND_LEAD_NS, tp + BAND_LEAD_NS + BAND_SPAN_NS)
        rows.append(dict(profile_x_m=float(px), surface_z_m=float(sz), cover_m=float(d1),
                         mudstone_m=float(d2), agl_m=float(agl), t_air_ns=float(t_air),
                         t_prim_geom_ns=float(t_prim), t_mult_vertical_geom_ns=float(t_mult),
                         t_basal_geom_ns=float(t_basal),
                         t_prim_pick_ns=tp, prim_pick_abs=ap, t_band_pick_ns=tb, band_pick_abs=ab,
                         band_minus_vertical_multiple_ns=(tb - t_mult) if tb else None,
                         band_minus_calibrated_multiple_ns=(tb - (2 * tp - t_air)) if (tb and tp) else None,
                         band_path_ratio=((tb - t_air) / (tp - t_air)) if (tb and tp) else None))

    # A2 斜率：几何预测（中心差分）与互相关实测
    px_a = np.array([r['profile_x_m'] for r in rows])
    tg = np.array([r['t_prim_geom_ns'] for r in rows])
    tm = np.array([r['t_mult_vertical_geom_ns'] for r in rows])
    tb = np.array([r['t_band_pick_ns'] for r in rows])
    dt_s = (t[1] - t[0]) * 1e-9
    slope_rows = []
    for k in range(SLOPE_HALF_WIDTH, len(rows) - SLOPE_HALF_WIDTH):
        w = slice(k - SLOPE_HALF_WIDTH, k + SLOPE_HALF_WIDTH + 1)
        span = (px_a[k + SLOPE_HALF_WIDTH] - px_a[k - SLOPE_HALF_WIDTH])
        s_prim_geom = (tg[w][-1] - tg[w][0]) / span
        s_mult_geom = (tm[w][-1] - tm[w][0]) / span
        s_band_pick = (tb[w][-1] - tb[w][0]) / span
        lo_t = tb[k] + BAND_LEAD_NS - 30.0
        ids = np.flatnonzero((t >= lo_t) & (t <= tb[k] + 60.0))
        lag, cc = xcorr_lag(env[ids, k], env[ids, k + 1], dt_s, 20e-9)
        slope_rows.append(dict(profile_x_m=float(px_a[k]),
                               slope_primary_geom_ns_per_m=float(s_prim_geom),
                               slope_multiple_geom_ns_per_m=float(s_mult_geom),
                               slope_band_pick_ns_per_m=float(s_band_pick),
                               slope_band_xcorr_ns_per_m=float(lag * 1e9 / (px_a[k + 1] - px_a[k])),
                               xcorr_peak=float(cc)))

    # ---- 图 1：B-scan 叠加（上行全动态，下行弱带增强裁剪）----
    sel = (t >= DISPLAY_TMIN_NS) & (t <= DISPLAY_TMAX_NS)
    clip = np.percentile(np.abs(real[sel, :]), DISPLAY_CLIP_PERCENTILE)
    sel2 = (t >= 220.0) & (t <= 450.0)  # 预声明：弱带增强窗
    clip2 = np.percentile(np.abs(real[sel2, :]), DISPLAY_CLIP_PERCENTILE)
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    for ax, cl, ttl in ((axes[0], clip, '全动态裁剪（98百分位，80–600ns）'),
                        (axes[1], clip2, '弱带增强裁剪（98百分位，仅220–450ns窗统计）')):
        ax.imshow(real[sel, :], aspect='auto', cmap='gray', vmin=-cl, vmax=cl,
                  extent=[px_a[0], px_a[-1], t[sel][-1], t[sel][0]])
        ax.plot(px_a, tg, 'g-', lw=1.5, label='几何预测 首次覆盖层返回')
        ax.plot(px_a, tm, 'c--', lw=1.5, label='几何预测 垂直二次往返（多次波）')
        ax.plot(px_a, [r['t_basal_geom_ns'] for r in rows], 'm-', lw=1.5, label='几何预测 底砂顶')
        ax.plot(px_a, [r['t_prim_pick_ns'] for r in rows], 'y.', ms=4, label='实测拾取 首次返回')
        ax.plot(px_a, [r['t_band_pick_ns'] for r in rows], 'r.', ms=4, label='实测拾取 强斜带')
        ax.set_xlabel('剖面 X / m'); ax.set_ylabel('官方 direct SFCW 时间 / ns')
        ax.set_title(f'41道 H1 Hann B-scan 与几何预测到时叠加（AGL 8m 夜间批）；{ttl}')
        ax.legend(loc='lower right', fontsize=8)
    fig.tight_layout(); fig.savefig(out / 'bscan_multiple_overlay.png', dpi=150); plt.close(fig)

    # ---- 图 2：残差与路径比 ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(px_a, [r['band_minus_vertical_multiple_ns'] for r in rows], 'b.-')
    axes[0].axhline(0, color='k', lw=0.5)
    axes[0].set_ylabel('强带到时 − 垂直多次波预测 / ns')
    axes[0].set_title('强斜带相对垂直二次往返预测的到时残差（负=提前到达）')
    axes[1].plot(px_a, [r['band_path_ratio'] for r in rows], 'r.-')
    axes[1].axhline(2.0, color='k', ls='--', lw=0.8, label='垂直多次波=2.0')
    axes[1].set_ylabel('路径比 (t带−t空)/(t首−t空)'); axes[1].set_xlabel('剖面 X / m')
    axes[1].legend()
    fig.tight_layout(); fig.savefig(out / 'multiple_residual_ratio.png', dpi=150); plt.close(fig)

    # ---- 图 3：斜率对照 ----
    sx = [r['profile_x_m'] for r in slope_rows]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(sx, [r['slope_primary_geom_ns_per_m'] for r in slope_rows], 'g-', label='几何 首次返回斜率')
    ax.plot(sx, [r['slope_multiple_geom_ns_per_m'] for r in slope_rows], 'c--', label='几何 多次波斜率(=2×首次)')
    ax.plot(sx, [r['slope_band_pick_ns_per_m'] for r in slope_rows], 'r.-', ms=4, label='实测 强带斜率(脊线差分)')
    ax.set_xlabel('剖面 X / m'); ax.set_ylabel('斜率 / (ns/m)')
    ax.set_title('A2 强斜带斜率对照（邻道互相关因包络过宽不分辨，未采用；数值保留在 JSON）')
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / 'slope_comparison.png', dpi=150); plt.close(fig)

    resid = np.array([r['band_minus_vertical_multiple_ns'] for r in rows])
    resc = np.array([r['band_minus_calibrated_multiple_ns'] for r in rows])
    ratio = np.array([r['band_path_ratio'] for r in rows])
    result = dict(
        status='COMPLETE', schema='line9-multiple-attribution/1',
        predeclared=dict(primary_window_ns=PRIMARY_WINDOW_NS, band_lead_ns=BAND_LEAD_NS,
                         band_span_ns=BAND_SPAN_NS, interface_footprint_m=INTERFACE_FOOTPRINT_M,
                         display_clip_percentile=DISPLAY_CLIP_PERCENTILE,
                         display_window_ns=[DISPLAY_TMIN_NS, DISPLAY_TMAX_NS],
                         slope_half_width_stations=SLOPE_HALF_WIDTH, reference_frequency_Hz=F_REF),
        cover_group_index_95MHz=ng_cover, mudstone_group_index_95MHz=ng_mud,
        station_count=len(rows),
        band_minus_vertical_multiple_ns=dict(mean=float(resid.mean()), min=float(resid.min()),
                                             max=float(resid.max()), std=float(resid.std())),
        band_minus_calibrated_multiple_ns=dict(
            note='校准版垂直多次波预测 t=2*t_prim_pick−t_air，剔除速度模型系统偏差后的残余提前量',
            mean=float(resc.mean()), min=float(resc.min()),
            max=float(resc.max()), std=float(resc.std())),
        band_path_ratio=dict(mean=float(ratio.mean()), min=float(ratio.min()),
                             max=float(ratio.max()), std=float(ratio.std())),
        slope_band_xcorr_vs_multiple_geom_ratio=dict(
            note='邻道互相关在0.5m道距/0.6ns采样/Hann包络宽度下不分辨（零滞后占优），不作证据；仅保留数值',
            mean=float(np.mean([r['slope_band_xcorr_ns_per_m'] / r['slope_multiple_geom_ns_per_m']
                                for r in slope_rows if abs(r['slope_multiple_geom_ns_per_m']) > 1e-6]))),
        slope_band_pick_vs_multiple_geom_ratio=dict(
            mean=float(np.mean([r['slope_band_pick_ns_per_m'] / r['slope_multiple_geom_ns_per_m']
                                for r in slope_rows if abs(r['slope_multiple_geom_ns_per_m']) > 1e-6])),
            std=float(np.std([r['slope_band_pick_ns_per_m'] / r['slope_multiple_geom_ns_per_m']
                              for r in slope_rows if abs(r['slope_multiple_geom_ns_per_m']) > 1e-6]))),
        per_station=rows, slope_rows=slope_rows,
        inputs=dict(sfcw=str(a.sfcw), geometry=str(a.geometry), materials=str(a.materials),
                    contract=str(a.contract)),
        limits='仅当前合成模型（AGL8m夜间批，官方impulse/direct/Hann）；几何预测用法向局部平层近似，'
               '未计侧偏/角度/展宽；拾取窗均为预声明；不是实测结论。',
    )
    (out / 'attribution.json').write_text(json.dumps(result, ensure_ascii=False, indent=1), 'utf-8')
    print('residual mean/min/max ns:', resid.mean(), resid.min(), resid.max())
    print('path ratio mean/min/max:', ratio.mean(), ratio.min(), ratio.max())
    print('wrote', out)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--sfcw', type=Path, required=True)
    p.add_argument('--geometry', type=Path, required=True)
    p.add_argument('--materials', type=Path, required=True)
    p.add_argument('--contract', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
