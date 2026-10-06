"""Validate and render the actual full-domain pulse wavefield; never run FDTD."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json
from line9_large_domain_snapshot import audit, SPACING
import analyze_line9_2d_sfcw as sfcw


def main(study, out):
    if out.exists():
        raise ValueError('Fresh wavefield report directory required')
    current = audit(study/'execution_contract.json', True)
    if current != json.loads((study/'completed_verification.json').read_text('utf-8')):
        raise ValueError('Completed wavefield evidence changed')
    c = json.loads((study/'execution_contract.json').read_text('utf-8'))
    files = current['snapshots']
    peak = 0.
    for row in files:
        with h5py.File(row['file']) as h:
            peak = max(peak, float(np.max(np.abs(h['Ez'][:]))))
    if not np.isfinite(peak) or peak <= 0:
        raise ValueError('Invalid all-frame field scale')
    out.mkdir(parents=True)
    # Delivers a correctly labelled single-station B-scan plus complete SFCW arrays.
    sfcw.main(study, out/'sfcw')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from PIL import Image
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    contacts = json.loads((study/'contacts_m.json').read_text('utf-8'))['points_x_elevation_m']
    selected = list(range(0, len(files), 3))
    if selected[-1] != len(files)-1:
        selected.append(len(files)-1)
    frames = []
    norm = SymLogNorm(linthresh=peak*1e-5, vmin=-peak, vmax=peak)
    colors = {'surface': 'black', 'clay_base': '#ef8a00', 'upper_mudstone_base': '#f600de',
              'lower_mudstone_top': '#00ad77', 'lower_mudstone_base': '#5066ff'}
    names = {'surface': '地表', 'clay_base': '粉质黏土底', 'upper_mudstone_base': '上段泥岩底',
             'lower_mudstone_top': '下段泥岩顶', 'lower_mudstone_base': '下段泥岩底'}
    for index in selected:
        with h5py.File(files[index]['file']) as h:
            field = h['Ez'][:, :, 0].T
            time_ns = float(h.attrs['time'])*1e9
        fig, ax = plt.subplots(figsize=(15, 5), layout='constrained')
        picture = ax.imshow(field, origin='lower', extent=[-50, 350, 405, 480],
                            aspect='equal', cmap='RdBu_r', norm=norm, interpolation='nearest')
        for name, points in contacts.items():
            p = np.asarray(points)
            ax.plot(p[:, 0], p[:, 1], color=colors[name], label=names[name], linewidth=.9)
        g = c['groups'][0]
        for key, marker, label in [('tx_m', '^', '发射点'), ('rx_m', 'v', '接收点')]:
            ax.plot(g[key][0]-50, g[key][1]+405, marker=marker, color='#fff700',
                    markeredgecolor='black', markersize=5, label=label)
        ax.set(xlim=(-50, 350), ylim=(405, 480), xlabel='原剖面横坐标 X / m', ylabel='高程 / m',
            title=f'大域400×75m；X{g["profile_x_m"]}；离地15m；Ricker脉冲总场 Ez；t={time_ns:.2f}ns\n'
                  '求解网格2.5cm；快照20cm/约2ns；固定SymLog色标；动画隔3帧播放，非SFCW波场')
        ax.legend(ncol=4, loc='upper right', fontsize=7)
        fig.colorbar(picture, ax=ax, label='瞬时 Ez / (V/m)，各帧同一色标')
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba()).copy()
        frames.append(Image.fromarray(rgba[:, :, :3]).convert('P', palette=Image.Palette.ADAPTIVE))
        if index == selected[len(selected)//3]:
            fig.savefig(out/'large_domain_wavefield_example.png', dpi=120)
        plt.close(fig)
    gif = out/'line9_large_domain_wavefield.gif'
    frames[0].save(gif, save_all=True, append_images=frames[1:], duration=65, loop=0)
    for frame in frames:
        frame.close()
    save_json(out/'wavefield_report.json', dict(status='PASS_NATIVE_PASSIVITY_AND_TRANSFORMS_NOT_CAUSAL_ATTRIBUTION',
        completed_verification_sha256=digest(study/'completed_verification.json'),
        calls_solver=False, full_snapshot_count=len(files), GIF_frame_count=len(selected),
        GIF_bytes=gif.stat().st_size, GIF_sha256=digest(gif), field_shared_peak_V_m=peak,
        passive_receiver_bit_identical=True, snapshot_grid_m=SPACING,
        native_grid_m=[.025]*3, snapshot_interval_s=c['dt_s']*34,
        SFCW_report_sha256=digest(out/'sfcw/analysis_report.json'),
        limitations='Raw pulse wavefront diagnostic; no per-frame normalization. Snapshot sampling does not certify all source frequencies or detailed subsurface phase. GIF decimates time further. Observed returns require path/control evidence before attribution to PML or geology.'))
    print(f'Wavefield report complete: {len(files)} native frames; {len(selected)} GIF frames')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    main(a.study.resolve(), a.out.resolve())
