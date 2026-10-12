"""B-r3 follow-up: background / coherent-clutter suppression on the 3-station
3D bowtie SFCW traces, to test whether the bottom-sand return becomes
resolvable once the system drag tail and antenna ringing comb are removed.

Background estimators (all within the existing scope of 'background and
coherent clutter suppression'):
- xmean  : per-time mean across the 3 stations (rank-1 cross-trace estimate)
- xmedian: per-time median across the 3 stations (robust variant)
- tavg   : along-trace moving average (80 ns window), removes the smooth
           exponential drag tail but keeps fast structure

Caveat: with only 3 stations the cross-trace estimators are weak; whatever is
common to all 3 stations (system tail, ringing, flat-lying geology) is treated
as background by construction. 2D same-station traces processed identically
for reference. Same official SFCW chain as analyze_antenna_b_r3.py.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
                                              direct_frequency_response,
                                              reconstruct_time_response)
from analyze_antenna_b_r3 import (FREQ, STATIONS, FEATURES_2D_NS,
                                  envelope_db)

TAVG_WIN_NS = 80.0


def load_sfcw(path, component):
    src = load_source(path)
    rx = load_receiver(path, component=component)
    fr = direct_frequency_response(src, rx, FREQ)
    tr = reconstruct_time_response(fr, window='hann', zero_pad_factor=8)
    return tr


def moving_average(x, dt_ns, win_ns):
    w = max(3, int(round(win_ns/dt_ns)) | 1)  # odd window
    k = np.ones(w)/w
    return np.convolve(x, k, mode='same')


def evaluate(t, sig, shift):
    """Peak dB of the envelope near expected cover/bottom windows and the
    late floor, plus protrusion over the floor."""
    db, _ = envelope_db(sig)
    out = {}
    for name, tt in [('cover', FEATURES_2D_NS['cover_return']-shift),
                     ('bottom', None)]:
        if tt is None:
            continue
        w = (t >= tt-8) & (t <= tt+8)
        out[f'{name}_peak_db'] = float(db[w].max())
    return db, out


def main(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    shift = 2*(8.0-5.0)/299792458.0*1e9  # 2D features appear this much earlier in 3D

    # load all traces first (need the full set before cross-trace backgrounds)
    tr3, tr2 = {}, {}
    for sx in STATIONS:
        cid = f'ant3d_x{round(sx*100):05d}'
        tr3[sx] = load_sfcw(a.exec3d/'prepared/cases'/cid/'profile.h5', 'Ey')
        if sx < 185:
            p2 = a.exec2d/f'line_x{round(sx*100):05d}_H1.h5'
        else:
            p2 = Path('artifacts/local_checks/2026-10-11_line9_paper_mudstone_results_r2/high_x19000_H1.h5')
        tr2[sx] = load_sfcw(p2, 'Ez')

    def prep(tr):
        t = tr.time*1e9
        m = t <= 900
        s = np.asarray(tr.real_bandpass, dtype=np.float64)[m]
        return t[m], s - s.mean()

    T3 = {sx: prep(tr3[sx]) for sx in STATIONS}
    T2 = {sx: prep(tr2[sx]) for sx in STATIONS}
    t3 = T3[STATIONS[0]][0]
    dt3 = t3[1]-t3[0]

    stack3 = np.array([T3[sx][1] for sx in STATIONS])
    stack2 = np.array([T2[sx][1] for sx in STATIONS])
    bg = {
        'xmean': (stack3.mean(axis=0), stack2.mean(axis=0)),
        'xmedian': (np.median(stack3, axis=0), np.median(stack2, axis=0)),
    }

    summary = {'estimators': ['xmean', 'xmedian', 'tavg80ns'],
               'shift_2d_to_3d_ns': shift,
               'note': 'cross-trace backgrounds from only 3 stations; common-mode '
                       'geology is removed by construction. dB normalised per trace '
                       'AFTER suppression; compare protrusion over own floor.',
               'stations': {}}

    def metrics(t, sig, sx, is3d):
        db, _ = envelope_db(sig)
        row = {}
        bs = FEATURES_2D_NS['bottom_sand'][sx] - (shift if is3d else 0.0)
        cv = FEATURES_2D_NS['cover_return'] - (shift if is3d else 0.0)
        for name, tt in [('cover', cv), ('bottom', bs)]:
            w = (t >= tt-8) & (t <= tt+8)
            idx = np.flatnonzero(w)[np.argmax(db[w])]
            row[f'{name}_peak_db'] = float(db[idx])
            row[f'{name}_peak_time_ns'] = float(t[idx])
        floor = float(db[(t >= 600) & (t <= 900)].mean())
        row['floor_600_900_db'] = floor
        row['bottom_protrusion_db'] = row['bottom_peak_db'] - floor
        row['cover_protrusion_db'] = row['cover_peak_db'] - floor
        return db, row

    # ---- per-estimator figures: 3D residual envelopes ----
    for est in ['xmean', 'xmedian']:
        fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        b3, b2 = bg[est]
        for ax, sx in zip(axs, STATIONS):
            r3 = T3[sx][1] - b3
            t2 = T2[sx][0]
            r2 = T2[sx][1] - b2
            db3, row3 = metrics(t3, r3, sx, True)
            db2, row2 = metrics(t2, r2, sx, False)
            summary['stations'].setdefault(f'x{sx:.0f}', {})[f'3d_{est}'] = row3
            summary['stations'][f'x{sx:.0f}'][f'2d_{est}'] = row2
            ax.plot(t3, db3, lw=.7, color='tab:red', label='3D 领结（去背景后）')
            ax.plot(t2, db2, lw=.7, color='tab:blue', alpha=.75, label='2D 线源（同法去背景）')
            cv, bs = FEATURES_2D_NS['cover_return'], FEATURES_2D_NS['bottom_sand'][sx]
            ax.axvspan(cv-6, cv+6, color='tab:blue', alpha=.12)
            ax.axvspan(bs-6, bs+6, color='tab:blue', alpha=.2)
            ax.axvspan(cv-shift-6, cv-shift+6, color='tab:red', alpha=.12)
            ax.axvspan(bs-shift-6, bs-shift+6, color='tab:red', alpha=.2)
            ax.set_xlim(60, 500)
            ax.set_ylim(-62, 2)
            ax.set_ylabel('包络 dB')
            ax.grid(alpha=.3)
            ax.legend(fontsize=8, loc='upper right')
            # mark the detected in-band bottom peak so weak residuals are visible
            bsr = FEATURES_2D_NS['bottom_sand'][sx]-shift
            wb = (t3 >= bsr-8) & (t3 <= bsr+8)
            i_pk = np.flatnonzero(wb)[np.argmax(db3[wb])]
            ax.plot(t3[i_pk], db3[i_pk], 'v', color='tab:red', ms=7, zorder=5)
            ax.annotate(f'{db3[i_pk]:.1f} dB @ {t3[i_pk]:.0f} ns',
                        xy=(t3[i_pk], db3[i_pk]), xytext=(t3[i_pk]-55, db3[i_pk]+7),
                        fontsize=8, color='tab:red',
                        arrowprops=dict(arrowstyle='-', color='tab:red', lw=.6))
            ax.set_title(f'站位 x={sx:.0f} m：底砂突出度 3D {row3["bottom_protrusion_db"]:+.1f} dB / '
                         f'2D {row2["bottom_protrusion_db"]:+.1f} dB（相对各自 600–900ns 本底；红三角=带内 3D 峰）',
                         fontsize=10)
        axs[-1].set_xlabel('时间 ns')
        name = {'xmean': '跨道均值（3 站）', 'xmedian': '跨道中位数（3 站）'}[est]
        fig.suptitle(f'背景抑制后包络：背景={name}\n'
                     '红=3D 本批，蓝=2D 同站同法；红底=3D 期望到时（2D 实测减 20ns），蓝底=2D 实测到时\n'
                     '注意：仅 3 站位，凡三站共有成分（系统拖尾/振铃/平缓地层）均被当作背景去除')
        fig.tight_layout()
        fig.savefig(a.out/f'fig5_bgsuppress_{est}.png', dpi=150)
        plt.close(fig)

    # ---- along-trace moving-average suppression (no cross-trace info) ----
    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, sx in zip(axs, STATIONS):
        t2 = T2[sx][0]
        dt2 = t2[1]-t2[0]
        r3 = T3[sx][1] - moving_average(T3[sx][1], dt3, TAVG_WIN_NS)
        r2 = T2[sx][1] - moving_average(T2[sx][1], dt2, TAVG_WIN_NS)
        db3, row3 = metrics(t3, r3, sx, True)
        db2, row2 = metrics(t2, r2, sx, False)
        summary['stations'][f'x{sx:.0f}']['3d_tavg80'] = row3
        summary['stations'][f'x{sx:.0f}']['2d_tavg80'] = row2
        ax.plot(t3, db3, lw=.7, color='tab:red', label='3D 领结（去背景后）')
        ax.plot(t2, db2, lw=.7, color='tab:blue', alpha=.75, label='2D 线源（同法去背景）')
        cv, bs = FEATURES_2D_NS['cover_return'], FEATURES_2D_NS['bottom_sand'][sx]
        ax.axvspan(cv-6, cv+6, color='tab:blue', alpha=.12)
        ax.axvspan(bs-6, bs+6, color='tab:blue', alpha=.2)
        ax.axvspan(cv-shift-6, cv-shift+6, color='tab:red', alpha=.12)
        ax.axvspan(bs-shift-6, bs-shift+6, color='tab:red', alpha=.2)
        ax.set_xlim(60, 500)
        ax.set_ylim(-45, 2)
        ax.set_ylabel('包络 dB')
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_title(f'站位 x={sx:.0f} m：底砂突出度 3D {row3["bottom_protrusion_db"]:+.1f} dB / '
                     f'2D {row2["bottom_protrusion_db"]:+.1f} dB', fontsize=10)
    axs[-1].set_xlabel('时间 ns')
    fig.suptitle(f'背景抑制后包络：背景=沿道 {TAVG_WIN_NS:.0f}ns 滑动平均（不用跨道信息，只去缓变拖尾）\n'
                 '红=3D 本批，蓝=2D 同站同法；红底=3D 期望到时，蓝底=2D 实测到时')
    fig.tight_layout()
    fig.savefig(a.out/'fig5_bgsuppress_tavg80.png', dpi=150)
    plt.close(fig)

    # ---- baseline (no suppression) metrics for the same windows ----
    for sx in STATIONS:
        t2 = T2[sx][0]
        db3, row3 = metrics(t3, T3[sx][1], sx, True)
        db2, row2 = metrics(t2, T2[sx][1], sx, False)
        summary['stations'][f'x{sx:.0f}']['3d_none'] = row3
        summary['stations'][f'x{sx:.0f}']['2d_none'] = row2

    (a.out/'bgsuppress_analysis.json').write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), 'utf-8')
    print(json.dumps({'out': str(a.out)}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exec3d', type=Path,
                   default=Path('artifacts/local_checks/2026-10-12_antenna_b_exec_r3'))
    p.add_argument('--exec2d', type=Path,
                   default=Path('artifacts/local_checks/2026-10-11_line9_overnight_anchors_r1'))
    p.add_argument('--out', type=Path,
                   default=Path('artifacts/research_checks/2026-10-12_antenna_b_r3_bgsuppress'))
    main(p.parse_args())
