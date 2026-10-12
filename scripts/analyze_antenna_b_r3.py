"""Stage B r3 (3D bowtie line segment) analysis: raw impulse A-scans, official SFCW
reconstruction, and same-chain comparison against the 2D overnight anchors (H1).

Inputs (all already local, hash-audited upstream):
- artifacts/local_checks/2026-10-12_antenna_b_exec_r3/prepared/cases/ant3d_x*/profile.h5
  (3D bowtie pair, 5 m AGL, 900 ns, rx1/Ey at the passive Rx gap)
- artifacts/local_checks/2026-10-11_line9_overnight_anchors_r1/line_x*_H1.h5
  (2D line source, ~8 m AGL, 2400 ns, rx1/Ez)

Official chain identical for both: load_source/load_receiver ->
direct_frequency_response (501 tones, 20-170 MHz, step 0.3 MHz) ->
reconstruct_time_response (hann, zero_pad 8). No new physics; amplitude
normalised per trace (2D line source vs 3D bowtie strengths are not comparable).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import hilbert, find_peaks

from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
                                              direct_frequency_response,
                                              reconstruct_time_response)
from hs_capsule_identity import sha256 as sha

FREQ = 20e6 + 0.3e6*np.arange(501)
STATIONS = [170.0, 180.0, 190.0]
AGL_3D, AGL_2D = 5.0, 8.0
C0 = 299792458.0
# observed 2D overnight-line features (docs/research/2026-10-11_overnight_bscan_execution.md)
FEATURES_2D_NS = {'surface': 53.5, 'cover_return': 232.0,
                  'bottom_sand': {170.0: 355.955, 180.0: 380.073, 190.0: 399.202}}


def load_sfcw(path, component):
    src = load_source(path)
    rx = load_receiver(path, component=component)
    fr = direct_frequency_response(src, rx, FREQ)
    tr = reconstruct_time_response(fr, window='hann', zero_pad_factor=8)
    return src, rx, fr, tr


def envelope_db(x):
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()  # drop the near-DC charge-relaxation tail before the envelope
    n = len(x)
    z = hilbert(np.pad(x, (0, n)))[:n]  # zero-pad so the FFT-Hilbert wrap lands outside
    env = np.abs(z)
    return 20*np.log10(env/env.max()+1e-12), env


def sliding_env_db(x, dt_ns, win_ns=2.0):
    """Windowed max-|.| envelope on the raw (uncompressed) trace.

    The raw 900 ns record ends in a decaying tail that is still well above the
    Hilbert transform's numeric floor; an FFT-Hilbert envelope of such a record
    shows a spurious end-of-trace upturn even after zero padding. A sliding
    window max avoids the circular-convolution boundary entirely. Returns
    (t_centers_ns, db normalised to the record max)."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    w = max(1, int(round(win_ns/dt_ns)))
    n = len(x)//w*w
    env = np.abs(x[:n]).reshape(-1, w).max(axis=1)
    t = (np.arange(len(env))+0.5)*w*dt_ns
    return t, 20*np.log10(env/env.max()+1e-12)


def peaks(t, db, tmin, height, dist_ns, dt_ns):
    m = t >= tmin
    pk, _ = find_peaks(db[m], height=height, distance=max(1, int(dist_ns/dt_ns)))
    return [(float(t[m][i]), float(db[m][i])) for i in pk]


def main(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    summary = {'stations': {}, 'chain': 'official direct 501 tones 20-170MHz + hann zpad8',
               'inputs': {}, 'limits': '3 stations only; 3D uses cross-track invariant terrain; '
                                       '2D/3D amplitudes normalised per trace, not comparable in absolute terms. '
                                       '2D: line source ~8m AGL 2400ns Ez; 3D: bowtie pair 5m AGL 900ns Ey.'}
    shift = 2*(AGL_2D-AGL_3D)/C0*1e9   # ns; 2D features appear this much earlier in 3D
    surf_3d = 2*np.hypot(AGL_3D, 0.65)/C0*1e9
    direct_3d = 1.3/C0*1e9
    summary['expected_3d_ns'] = dict(direct_tx_rx=direct_3d, surface_bounce=surf_3d,
                                     agl_shift_2d_to_3d=shift)

    data = {}
    for sx in STATIONS:
        cid = f'ant3d_x{round(sx*100):05d}'
        p3 = a.exec3d/'prepared/cases'/cid/'profile.h5'
        if sx < 185:
            p2 = a.exec2d/f'line_x{round(sx*100):05d}_H1.h5'
        else:
            # 190 m was reused into the overnight line from the paper-mudstone batch
            p2 = Path('artifacts/local_checks/2026-10-11_line9_paper_mudstone_results_r2/high_x19000_H1.h5')
        src3, rx3, fr3, tr3 = load_sfcw(p3, 'Ey')
        src2, rx2, fr2, tr2 = load_sfcw(p2, 'Ez')
        data[sx] = dict(p3=p3, p2=p2, src3=src3, rx3=rx3, tr3=tr3, tr2=tr2, fr3=fr3, fr2=fr2)
        summary['inputs'][cid] = dict(sha256=sha(p3))
        summary['inputs'][f'line_x{round(sx*100):05d}_H1'] = dict(sha256=sha(p2))

    # ---- fig1: raw 3D impulse A-scan envelopes ----
    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, sx in zip(axs, STATIONS):
        d = data[sx]
        dt3 = d['rx3'].dt
        tsw, dbsw = sliding_env_db(d['rx3'].samples, dt3*1e9)
        ax.plot(tsw, dbsw, lw=.6, color='0.2')
        ax.set_ylim(-80, 2)
        ax.axvline(direct_3d, color='tab:blue', ls='--', lw=.8)
        ax.axvline(surf_3d, color='tab:green', ls='--', lw=.8)
        bs = FEATURES_2D_NS['bottom_sand'][sx]-shift
        ax.axvspan(FEATURES_2D_NS['cover_return']-shift-6, FEATURES_2D_NS['cover_return']-shift+6,
                   color='tab:orange', alpha=.25)
        ax.axvspan(bs-6, bs+6, color='tab:red', alpha=.25)
        ax.set_ylabel('包络 dB')
        ax.grid(alpha=.3)
        ax.set_title(f'站位 x={sx:.0f} m', fontsize=10)
    axs[0].annotate('蓝虚线=Tx→Rx直达(4.3ns)；绿虚线=地表反射(33.6ns)\n'
                    '橙带=2D覆盖层返回平移后期望(212ns)；红带=2D底砂平移后期望',
                    xy=(0.01, 0.55), xycoords='axes fraction', fontsize=9)
    axs[-1].set_xlabel('时间 ns')
    fig.suptitle('B-r3 三维领结天线测线段：原始冲激 A-scan 滑动窗包络（Rx 间隙 Ey，5 m 离地，900 ns 窗）\n'
                 '激励=官方 impulse（t=0 单格冲激），未做 SFCW 频带压缩；包络=2ns 窗 max|.|（避免 Hilbert 边界假象）')
    fig.tight_layout()
    fig.savefig(a.out/'fig1_raw_ascan.png', dpi=150)
    plt.close(fig)

    # ---- fig2: SFCW reconstruction 3D vs 2D ----
    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, sx in zip(axs, STATIONS):
        d = data[sx]
        t3 = d['tr3'].time*1e9
        t2 = d['tr2'].time*1e9
        db3, _ = envelope_db(d['tr3'].real_bandpass)
        db2, _ = envelope_db(d['tr2'].real_bandpass)
        m3, m2 = t3 <= 900, t2 <= 900
        ax.plot(t3[m3], db3[m3], lw=.7, color='tab:red', label='3D 领结对 5m AGL（本批）')
        ax.plot(t2[m2], db2[m2], lw=.7, color='tab:blue', alpha=.75, label='2D 线源 8m AGL（夜间测线）')
        for name, tt in [('cover_return', FEATURES_2D_NS['cover_return']),
                         ('bottom_sand', FEATURES_2D_NS['bottom_sand'][sx])]:
            ax.axvline(tt, color='tab:blue', ls=':', lw=.8)
            ax.axvline(tt-shift, color='tab:red', ls=':', lw=.8)
        ax.set_ylim(-70, 2)
        ax.set_ylabel('包络 dB')
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_title(f'站位 x={sx:.0f} m（蓝点线=2D 实测特征到时；红点线=按航高差平移 20ns 后的 3D 期望）', fontsize=10)
    axs[-1].set_xlabel('时间 ns')
    fig.suptitle('同一官方 SFCW 链重建（direct 501 音调 20–170 MHz + Hann）：3D 领结 vs 2D 线源\n'
                 '2D 特征：覆盖层返回 232ns、底砂差场峰 356/380/399ns（170/180/190m）')
    fig.tight_layout()
    fig.savefig(a.out/'fig2_sfcw_compare.png', dpi=150)
    plt.close(fig)

    # ---- fig3: 3-station pseudo B-scan (3D SFCW envelope) ----
    fig, ax = plt.subplots(figsize=(7.5, 7))
    envs, taxis = [], None
    for sx in STATIONS:
        d = data[sx]
        t3 = d['tr3'].time*1e9
        m = t3 <= 900
        db, _ = envelope_db(d['tr3'].real_bandpass)
        envs.append(db[m])
        taxis = t3[m]
    img = np.array(envs).T
    ax.imshow(img, aspect='auto', cmap='gray', vmin=-55, vmax=0,
              extent=[0, 3, taxis[-1], 0])
    ax.set_xticks([0.5, 1.5, 2.5], [f'x={sx:.0f} m' for sx in STATIONS])
    ax.set_ylabel('时间 ns')
    ax.set_title('B-r3 三维领结测线段伪 B-scan（仅 3 站位，未插值；SFCW Hann 包络，dB）\n'
                 '横轨不变地形近似；不能据此判读连续同相轴')
    fig.tight_layout()
    fig.savefig(a.out/'fig3_pseudo_bscan.png', dpi=150)
    plt.close(fig)

    # ---- fig4: zoom 100-500 ns, feature visibility ----
    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, sx in zip(axs, STATIONS):
        d = data[sx]
        t3 = d['tr3'].time*1e9
        t2 = d['tr2'].time*1e9
        db3, _ = envelope_db(d['tr3'].real_bandpass)
        db2, _ = envelope_db(d['tr2'].real_bandpass)
        m3, m2 = (t3 >= 100) & (t3 <= 500), (t2 >= 100) & (t2 <= 500)
        ax.plot(t3[m3], db3[m3], lw=.8, color='tab:red', label='3D 领结对 5m AGL（本批）')
        ax.plot(t2[m2], db2[m2], lw=.8, color='tab:blue', alpha=.75, label='2D 线源 8m AGL（夜间测线）')
        cv, bs = FEATURES_2D_NS['cover_return'], FEATURES_2D_NS['bottom_sand'][sx]
        ax.axvspan(cv-6, cv+6, color='tab:blue', alpha=.12)
        ax.axvspan(bs-6, bs+6, color='tab:blue', alpha=.2)
        ax.axvspan(cv-shift-6, cv-shift+6, color='tab:red', alpha=.12)
        ax.axvspan(bs-shift-6, bs-shift+6, color='tab:red', alpha=.2)
        ax.set_ylim(-70, -25)
        ax.set_ylabel('包络 dB')
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_title(f'站位 x={sx:.0f} m（蓝底=2D 覆盖层返回 232ns 与底砂 {bs:.0f}ns；红底=3D 平移期望）', fontsize=10)
    axs[-1].set_xlabel('时间 ns')
    fig.suptitle('晚期窗 100–500 ns 放大：3D 领结在期望到时处是否出现可分辨回波？\n'
                 '（包络已去均值并零延拓，Hilbert 边界假象已消除）')
    fig.tight_layout()
    fig.savefig(a.out/'fig4_late_window_zoom.png', dpi=150)
    plt.close(fig)

    # ---- numbers ----
    for sx in STATIONS:
        d = data[sx]
        dt3 = d['rx3'].dt
        t3 = np.arange(len(d['rx3'].samples))*dt3*1e9
        db, _ = envelope_db(d['rx3'].samples)
        tsw, dbsw = sliding_env_db(d['rx3'].samples, dt3*1e9)
        row = {'raw_peaks_over_-40dB_first400ns': peaks(t3, db, 0, -40, 10, dt3*1e9)[:15],
               'raw_tail_sliding_env_600_900ns_db': float(dbsw[tsw >= 600].max())}
        t3r = d['tr3'].time*1e9
        db3r, env3r = envelope_db(d['tr3'].real_bandpass)
        m = t3r <= 900
        row['sfcw_peaks_over_-45dB'] = peaks(t3r[m], db3r[m], 0, -45, 8, t3r[1]-t3r[0])[:20]
        for name, tt in [('cover_return', FEATURES_2D_NS['cover_return']-shift),
                         ('bottom_sand', FEATURES_2D_NS['bottom_sand'][sx]-shift)]:
            w = (t3r >= tt-8) & (t3r <= tt+8)
            row[f'sfcw_env_max_near_{name}_{tt:.0f}ns_db'] = float(db3r[w].max())
        wfloor = (t3r >= 600) & (t3r <= 900)
        row['sfcw_floor_600_900ns_db'] = float(db3r[wfloor].mean())
        # 2D same-chain reference levels (observed 2D feature times, no shift)
        t2r = d['tr2'].time*1e9
        db2r, _ = envelope_db(d['tr2'].real_bandpass)
        for name, tt in [('cover_return', FEATURES_2D_NS['cover_return']),
                         ('bottom_sand', FEATURES_2D_NS['bottom_sand'][sx])]:
            w = (t2r >= tt-8) & (t2r <= tt+8)
            row[f'2d_env_max_at_{name}_{tt:.0f}ns_db'] = float(db2r[w].max())
        wfloor2 = (t2r >= 600) & (t2r <= 900)
        row['2d_floor_600_900ns_db'] = float(db2r[wfloor2].mean())
        # port ringing
        with h5py.File(d['p3']) as h:
            vt = h['ports/feed/Vtotal'][:]
            pt = h['ports/feed/time'][:]*1e9
        vdb, _ = envelope_db(vt)
        # residual port level AFTER each limit (window excludes the t=0 peak,
        # so this measures how far the feed ringing has decayed by then)
        row['port_Vtotal_residual_db'] = {
            f'after_{lim}ns': float(vdb[(pt >= lim) & (pt <= lim+100)].max())
            for lim in [50, 100, 200, 400]}
        summary['stations'][f'x{sx:.0f}'] = row
    (a.out/'analysis.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False), 'utf-8')
    print(json.dumps({'out': str(a.out)}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exec3d', type=Path,
                   default=Path('artifacts/local_checks/2026-10-12_antenna_b_exec_r3'))
    p.add_argument('--exec2d', type=Path,
                   default=Path('artifacts/local_checks/2026-10-11_line9_overnight_anchors_r1'))
    p.add_argument('--out', type=Path,
                   default=Path('artifacts/research_checks/2026-10-12_antenna_b_r3'))
    main(p.parse_args())
