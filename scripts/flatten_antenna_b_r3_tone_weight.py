"""Tone flattening for B-r3: emulate the real system's flat per-tone transmit
power (measured 39.3-40.6 dBm across 20-170 MHz) on the existing 3D bowtie
SFCW data, and re-evaluate the bottom-sand / direct-wave amplitude ratio
against the real-data reference of about -40 dB (raw data, user-provided).

Weight construction (empirical, no new solves):
- A2_ground_pair (stage A r3): Tx + passive Rx bowties, 1.3 m x-baseline,
  arms along y, 5 m AGL over a flat cover half-space. Its rxgap 501-tone
  response H_pair(f) contains the full two-way antenna system weighting
  (mismatch + patterns + effective lengths) plus direct/ground-bounce path
  interference. The path ripple is NOT system property, so |H_pair| is
  smoothed over ~15 MHz before inversion: w(f) = |Hsm(95MHz)| / |Hsm(f)|.
- Cross-check weight from A1_free_y port S11: w_s11(f) = round-trip mismatch
  factor ratio (1-|S11(fref)|^2)/(1-|S11(f)|^2).

Application: fr_flat = fr * w, then the unchanged official reconstruction
(hann, zero_pad 8). Direct wave peak = global envelope maximum; bottom-sand
and cover-return windows as in analyze_antenna_b_r3.py (20 ns AGL shift).

Limits: A2 is over flat cover half-space (no mudstone/sandstone); its record
is ~600 ns; smoothing cannot perfectly separate ripple from mismatch trend;
the real Rx antenna gain(f) is unknown beyond -2 dBi @100 MHz, so the
flattening approximates 'flat per-tone system response', not a measured
device transfer function.
"""
import argparse
import dataclasses
import json
from pathlib import Path

import h5py
import numpy as np

from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
                                              direct_frequency_response,
                                              reconstruct_time_response)
from analyze_antenna_b_r3 import (FREQ, STATIONS, FEATURES_2D_NS,
                                  envelope_db)

A2 = Path('artifacts/local_checks/2026-10-12_antenna_a_exec_r3/prepared/cases/A2_ground_pair/profile.h5')
A1 = Path('artifacts/local_checks/2026-10-12_antenna_a_exec_r3/prepared/cases/A1_free_y/profile.h5')
SMOOTH_MHZ = 15.0
FREF_MHZ = 95.0


def smooth(x, win):
    win = int(win) | 1
    return np.convolve(x, np.ones(win)/win, mode='same')


def main(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    # ---- weight from A2 pair link ----
    srcA = load_source(A2)
    rxA = load_receiver(A2, component='Ey', receiver_path='/rxs/rx1')
    frA = direct_frequency_response(srcA, rxA, FREQ)
    mag = np.abs(frA.response)
    wsmooth = smooth(mag, SMOOTH_MHZ/0.3)
    iref = int(round((FREF_MHZ-20.0)/0.3))
    w = wsmooth[iref]/wsmooth            # amplitude weight, 1.0 at 95 MHz

    # ---- cross-check weight from A1 free-space port S11 ----
    with h5py.File(A1) as h:
        fport = h['ports/feed/frequency'][:]
        s11 = h['ports/feed/S11'][:]
    mm = np.interp(FREQ, fport, np.abs(s11)**2)   # |S11|^2 at the 501 tones
    w_s11 = (1-mm[iref])/(1-mm)                    # round-trip mismatch ratio

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(FREQ/1e6, 20*np.log10(np.abs(w)), label='w(f)：A2 收发对链路（~15MHz 平滑）', color='tab:red')
    ax.plot(FREQ/1e6, 20*np.log10(w_s11), label='w_S11(f)：A1 端口失配双程比（交叉核对）', color='tab:blue', alpha=.8)
    ax.axhline(0, color='k', lw=.5)
    ax.set_xlabel('频率 MHz')
    ax.set_ylabel('幅度权重 dB（95 MHz 归一）')
    ax.set_title('等功率音调加权权重：把模型逐音调响应拉到"发射谱平坦"（对应实测发射功率 39.3–40.6 dBm 全带平坦）\n'
                 '低频端为正权重 = 补偿领结低频辐射塌陷；两曲线趋势一致则权重可信')
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(a.out/'fig1_tone_weight.png', dpi=150)
    plt.close(fig)

    # ---- apply to B-r3 stations ----
    shift = 2*(8.0-5.0)/299792458.0*1e9
    weights = {'s11': w_s11,   # mismatch flattening = flat accepted-power emulation (primary)
               'linkA2': w}    # full empirical link flattening (ablation: also removes channel lambda^2 physics)
    summary = {'weight': {'primary': 'w_s11: round-trip mismatch factor ratio from A1_free_y port S11, '
                                     'emulating the measured flat per-tone transmit power',
                          'ablation': 'w_linkA2: 1/|A2 pair link| smoothed 15 MHz (removes channel physics too; '
                                      'near-zero effect, kept as evidence)',
                          'w_db_at_20MHz_linkA2': float(20*np.log10(abs(w[0]))),
                          'w_s11_db_at_20MHz': float(20*np.log10(w_s11[0])),
                          'w_s11_db_at_170MHz': float(20*np.log10(w_s11[-1]))},
               'reference': 'real-data bottom-sand/direct about -40 dB (raw, user 2026-10-12)',
               'stations': {}}

    fig, axs = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for ax, sx in zip(axs, STATIONS):
        cid = f'ant3d_x{round(sx*100):05d}'
        p3 = a.exec3d/'prepared/cases'/cid/'profile.h5'
        src = load_source(p3)
        rx = load_receiver(p3, component='Ey')
        fr = direct_frequency_response(src, rx, FREQ)
        tr0 = reconstruct_time_response(fr, window='hann', zero_pad_factor=8)
        t = tr0.time*1e9
        m = t <= 900
        db0, _ = envelope_db(tr0.real_bandpass)
        dbs = {'raw': db0}
        for tag, wt in weights.items():
            frw = dataclasses.replace(fr, response=fr.response*wt)
            trw = reconstruct_time_response(frw, window='hann', zero_pad_factor=8)
            dbs[tag], _ = envelope_db(trw.real_bandpass)
        row = {}
        for tag, db in dbs.items():
            direct = float(db[:int(60/(t[1]-t[0]))].max())   # direct+surface region is the 0 dB anchor
            row[f'{tag}_direct_db'] = direct
            for name, tt in [('cover', FEATURES_2D_NS['cover_return']-shift),
                             ('bottom', FEATURES_2D_NS['bottom_sand'][sx]-shift)]:
                wnd = (t >= tt-8) & (t <= tt+8)
                row[f'{name}_{tag}_db'] = float(db[wnd].max()) - direct
            floor = float(db[(t >= 600) & (t <= 900)].mean()) - direct
            row[f'floor_{tag}_db'] = floor
            row[f'bottom_over_floor_{tag}_db'] = row['bottom_'+tag+'_db'] - floor
        summary['stations'][f'x{sx:.0f}'] = row
        ax.plot(t[m], db0[m], lw=.6, color='0.55', label='原样（未加权）')
        ax.plot(t[m], dbs['s11'][m], lw=.8, color='tab:red', label='等功率加权（失配补偿，主方案）')
        ax.plot(t[m], dbs['linkA2'][m], lw=.6, color='tab:blue', alpha=.6, label='全链路拉平（消融，含信道物理）')
        cv, bs = FEATURES_2D_NS['cover_return']-shift, FEATURES_2D_NS['bottom_sand'][sx]-shift
        ax.axvspan(cv-8, cv+8, color='tab:orange', alpha=.2)
        ax.axvspan(bs-8, bs+8, color='tab:red', alpha=.18)
        ax.axhline(-40, color='tab:green', ls='--', lw=.8)
        ax.set_xlim(0, 500)
        ax.set_ylim(-75, 3)
        ax.set_ylabel('包络 dB（直达=0）')
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_title(f'x={sx:.0f} m：底砂/直达 原样 {row["bottom_raw_db"]:.1f} → 失配补偿 {row["bottom_s11_db"]:.1f} dB'
                     f'（绿虚线=实测参照 −40 dB；橙带=覆盖层期望，红带=底砂期望）', fontsize=10)
    axs[-1].set_xlabel('时间 ns')
    fig.suptitle('等功率音调加权对 B-r3 三维领结 A-scan 的作用（同一官方 SFCW 链，仅插入逐音调权重）\n'
                 '目标：不差分条件下底砂回波达到实测量级（基覆界面约 −40 dB 相对直达波）')
    fig.tight_layout()
    fig.savefig(a.out/'fig2_flattened_ascan.png', dpi=150)
    plt.close(fig)

    (a.out/'toneflatten_analysis.json').write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), 'utf-8')
    print(json.dumps({'out': str(a.out)}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exec3d', type=Path,
                   default=Path('artifacts/local_checks/2026-10-12_antenna_b_exec_r3'))
    p.add_argument('--out', type=Path,
                   default=Path('artifacts/research_checks/2026-10-12_antenna_b_r3_toneflatten'))
    main(p.parse_args())
