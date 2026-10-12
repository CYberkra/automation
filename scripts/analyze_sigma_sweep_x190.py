"""Sigma sweep analysis (x190): bottom-sand / direct-wave ratio vs cover and
mudstone conductivity, same official SFCW chain as the B-r3 analysis.

Inputs:
- base point (mud 0.006, cover 0.001): B-r3 ant3d_x19000 (reused, hash-pinned)
- 3 new variants: artifacts/local_checks/2026-10-12_sigma_exec_r1/prepared/cases/sig_x19000_*/

Reference: real-data bottom-sand/direct about -40 dB (raw, user 2026-10-12).
Expected 3D feature times at x190 (5 m AGL): cover return 212 ns, bottom sand
379 ns (2D observed 399.202 ns minus 20 ns AGL shift).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
                                              direct_frequency_response,
                                              reconstruct_time_response)
from analyze_antenna_b_r3 import FREQ, FEATURES_2D_NS, envelope_db
from hs_capsule_identity import sha256 as sha

SHIFT = 2*(8.0-5.0)/299792458.0*1e9
BASE = ('artifacts/local_checks/2026-10-12_antenna_b_exec_r3/prepared/cases/ant3d_x19000/profile.h5',
        0.006, 0.001)
VARIANT_DIR = 'artifacts/local_checks/2026-10-12_sigma_exec_r1/prepared/cases'
VARIANTS = [('sig_x19000_m002_c001', 0.002, 0.001),
            ('sig_x19000_m002_c0005', 0.002, 0.0005),
            ('sig_x19000_m006_c0005', 0.006, 0.0005)]
COVER_T = FEATURES_2D_NS['cover_return']-SHIFT
BOTTOM_T = FEATURES_2D_NS['bottom_sand'][190.0]-SHIFT


def measure(path):
    src = load_source(path)
    rx = load_receiver(path, component='Ey')
    fr = direct_frequency_response(src, rx, FREQ)
    tr = reconstruct_time_response(fr, window='hann', zero_pad_factor=8)
    t = tr.time*1e9
    db, _ = envelope_db(tr.real_bandpass)
    out = {}
    for name, tt in [('surface', 33.6), ('cover', COVER_T), ('bottom', BOTTOM_T)]:
        w = (t >= tt-8) & (t <= tt+8)
        out[f'{name}_db'] = float(db[w].max())
        out[f'{name}_time_ns'] = float(t[w][np.argmax(db[w])])
    out['floor_db'] = float(db[(t >= 600) & (t <= 900)].mean())
    out['bottom_over_floor_db'] = out['bottom_db']-out['floor_db']
    return t, db, out


def main(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    rows = []
    curves = []
    for label, mud, cov, path in [('基准 m006_c001（B-r3 复用）', 0.006, 0.001, Path(BASE[0]))]+[
            (cid, mud, cov, Path(VARIANT_DIR)/cid/'profile.h5') for cid, mud, cov in VARIANTS]:
        t, db, m = measure(path)
        rows.append(dict(case=label, mudstone_sigma=mud, cover_sigma=cov,
                         input_sha256=sha(path), **m))
        curves.append((label, mud, cov, t, db))

    # ---- fig1: envelope overlay ----
    fig, ax = plt.subplots(figsize=(11, 6.5))
    colors = {'0.006': 'tab:red', '0.002': 'tab:blue'}
    for label, mud, cov, t, db in curves:
        ls = '-' if cov == 0.001 else '--'
        ax.plot(t[t <= 500], db[t <= 500], lw=.9, color=colors[str(mud)], ls=ls,
                label=f'泥岩σ={mud} 覆盖层σ={cov}')
    for tt, c, txt in [(33.6, 'tab:green', '地表'), (COVER_T, 'tab:orange', '覆盖层返回'),
                       (BOTTOM_T, 'tab:red', '底砂')]:
        ax.axvspan(tt-8, tt+8, color=c, alpha=.15)
    ax.axhline(-40, color='k', ls=':', lw=1)
    ax.text(405, -39, '实测参照：基覆界面 ≈ −40 dB（相对直达波）', fontsize=9)
    ax.set_xlim(0, 500)
    ax.set_ylim(-75, 3)
    ax.set_xlabel('时间 ns')
    ax.set_ylabel('包络 dB（直达=0）')
    ax.set_title('x190 σ 敏感性：四种（泥岩σ × 覆盖层σ）组合的 SFCW 包络（同一官方链，3D 领结，5 m AGL）\n'
                 '实线=覆盖层σ 0.001，虚线=0.0005；红=泥岩σ 0.006（基准），蓝=0.002')
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out/'fig1_sigma_envelopes.png', dpi=150)
    plt.close(fig)

    # ---- fig2: ratio vs sigma ----
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for cov, mk, ls in [(0.001, 'o', '-'), (0.0005, 's', '--')]:
        sel = sorted([r for r in rows if r['cover_sigma'] == cov], key=lambda r: r['mudstone_sigma'])
        ax.plot([r['mudstone_sigma']*1e3 for r in sel], [r['bottom_db'] for r in sel],
                mk+ls, color='tab:red', label=f'底砂/直达，覆盖层σ={cov}')
        ax.plot([r['mudstone_sigma']*1e3 for r in sel], [r['cover_db'] for r in sel],
                mk+ls, color='tab:orange', label=f'覆盖层返回/直达，覆盖层σ={cov}')
    ax.axhline(-40, color='k', ls=':', lw=1)
    ax.text(2.1, -41.5, '实测基覆界面 ≈ −40 dB', fontsize=9)
    ax.set_xlabel('泥岩 σ (mS/m)')
    ax.set_ylabel('包络峰 dB（直达=0）')
    ax.set_title('底砂与覆盖层返回电平 vs 电导率（x190，3D 领结，不差分）\n'
                 '判断：哪组 σ 使底砂自然落入实测量级')
    ax.grid(alpha=.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out/'fig2_ratio_vs_sigma.png', dpi=150)
    plt.close(fig)

    (a.out/'sigma_analysis.json').write_text(json.dumps(
        {'reference_bottom_over_direct_db': -40.0,
         'expected_times_ns': {'surface': 33.6, 'cover': COVER_T, 'bottom': BOTTOM_T},
         'rows': rows}, indent=1, ensure_ascii=False), 'utf-8')
    print(json.dumps({'out': str(a.out)}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path,
                   default=Path('artifacts/research_checks/2026-10-12_sigma_sweep_x190'))
    main(p.parse_args())
