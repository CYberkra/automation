"""Sparse three-anchor B-scan view for the 108m vs 36m large-domain control.

Only the three computed anchor stations are drawn (1.3 m column width, the
actual Tx/Rx offset); all other positions are masked grey -- no interpolation.
Total response and interface contrast (rough minus halfspace) are reconstructed
from the archived comparison arrays with the same source-normalized SFCW chain
used by the analysis script. Shared fixed symlog scales across all panels.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from analyze_hs4_cross_pc import response, reconstruct  # noqa: E402

RESULTS = ROOT / 'artifacts/research_checks/2026-10-04_hs4_large_domain_target_results_r1'
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-04_hs4_large_domain_target_r1'
STATIONS_36 = {'left': 14.6, 'centre': 17.6, 'right': 20.6}
SHIFT = 36.0
COLUMN_HALF_WIDTH_M = 0.65
T_LO, T_HI = 80.0, 240.0
LINTHRESH = 1e-4


def main():
    data = np.load(RESULTS / 'comparison_arrays.npz')
    reference = response(CAPSULE / 'centre_rough' / 'profile.h5', 'Ey')
    check = np.linalg.norm(reference.response - data['108m_centre_rough_frequency_response']) \
        / np.linalg.norm(data['108m_centre_rough_frequency_response'])
    if check > 1e-9:
        raise ValueError('archived arrays do not match capsule raw: %g' % check)
    traces = {}
    for label in ('36m', '108m'):
        for station in STATIONS_36:
            total = reconstruct(reference, data['%s_%s_rough_frequency_response' % (label, station)])
            traces[label, station, 'total'] = total.real_bandpass
            traces[label, station, 'contrast'] = data['%s_%s_signed_waveform' % (label, station)]
    time = total.time * 1e9
    mask = (time >= T_LO) & (time <= T_HI)
    t_edge = np.concatenate(([time[mask][0] - (time[1] - time[0]) / 2],
                             (time[mask][:-1] + time[mask][1:]) / 2,
                             [time[mask][-1] + (time[1] - time[0]) / 2]))
    vmax = {'total': max(np.max(np.abs(v[mask])) for v in
                         (traces[k] for k in traces if k[2] == 'total')),
            'contrast': max(np.max(np.abs(v[mask])) for v in
                            (traces[k] for k in traces if k[2] == 'contrast'))}
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharey=True, layout='constrained')
    for row, kind in enumerate(('total', 'contrast')):
        for col, label in enumerate(('36m', '108m')):
            ax = axes[row, col]
            stations = STATIONS_36 if label == '36m' else {s: x + SHIFT for s, x in STATIONS_36.items()}
            xs = sorted(stations.values())
            x0, x1 = (0.0, 36.0) if label == '36m' else (0.0, 108.0)
            ax.axvspan(x0, x1, color='0.92', zorder=0)
            for station, x in stations.items():
                strip = traces[label, station, kind][mask]
                ax.pcolormesh([x - COLUMN_HALF_WIDTH_M, x + COLUMN_HALF_WIDTH_M], t_edge,
                              strip[:, None], cmap='RdBu_r', shading='flat',
                              norm=matplotlib.colors.SymLogNorm(linthresh=LINTHRESH * vmax[kind],
                                                                vmin=-vmax[kind], vmax=vmax[kind]))
                ax.text(x, T_LO + 3, station, ha='center', fontsize=8, rotation=0)
            ax.set_xlim(x0, x1); ax.set_ylim(T_HI, T_LO)
            ax.set_title('%s  %s (fixed symlog |max|=%.2e)' % (label, kind, vmax[kind]), fontsize=10)
            ax.set_xlabel('x (m)')
            if col == 0:
                ax.set_ylabel('actual receiver time (ns)')
    fig.suptitle('Sparse 3-anchor view only; grey = not computed, no interpolation; '
                 'V4.0.0 double, source-normalized 20-170 MHz')
    out = RESULTS / 'large_domain_sparse_bscan.png'
    fig.savefig(out, dpi=150)
    print('saved', out, 'vmax', vmax)


if __name__ == '__main__':
    main()
