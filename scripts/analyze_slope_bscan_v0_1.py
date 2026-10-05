"""Gentle-slope B-scan visualization: model view, total field, background removal, t^2 gain, GIF.

Panels (2x3): model cross-section; total field (shared fixed scale); ideal
background-removed (slope_rough minus slope_fullcover, declared pair);
catalogue mean-trace removal; catalogue SVD rank-2 removal; ideal + declared
t^2 analysis gain (non-amplitude product). Truth overlays: constant 8 m AGL
surface echo and per-station interface two-way prediction (95 MHz refractive
index of the frozen cover recipe). Independent 7-point DFT check on all 26
traces. GIF: 200 dense snapshots of slope_rough_s06 with relief overlays.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

import hs4_slope_bscan_v0_1 as sim
from analyze_uav_local_free_space import FREQ, transfer
from hs_capsule_identity import sha256
from research_operator_contract import apply_configuration

C_SI = 299792458.
E0_SI = 8.854187817e-12
EPS_INF, DEPS, TAU, SIGMA = 11., 0.5, 6.4567e-9, 0.001
AIR_NS = 2*sim.AGL/C_SI*1e9


def n_cover(freq):
    w = 2*np.pi*freq
    eps = EPS_INF + DEPS/(1+1j*w*TAU) - 1j*SIGMA/(w*E0_SI)
    return float(np.sqrt(eps.real))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('fresh analysis output required')
    a.out.mkdir(parents=True)
    c = json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    audit = json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha256(a.study/'execution_contract.json'):
        raise ValueError('completed raw identity required')
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}
    stations = c['geometry']['stations_tx_x_m']
    mats, times = {}, None
    dft_max = 0.
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        value, source, receiver = transfer(raw, 1, 'Ey', g['spacing_m'])
        take = np.array([0, 83, 167, 250, 333, 417, 500])
        sy = source.samples.copy(); ry = receiver.samples.copy()
        n = round(20e-9/receiver.dt)
        ry[-n:] *= .5*(1+np.cos(np.linspace(0, np.pi, n)))
        sft = np.exp(-2j*np.pi*FREQ[take, None]*(source.time_offset+source.dt*np.arange(len(sy))))@sy
        rft = np.exp(-2j*np.pi*FREQ[take, None]*(receiver.time_offset+receiver.dt*np.arange(len(ry))))@ry
        independent = rft/sft/g['spacing_m']
        dft_max = max(dft_max, float(np.linalg.norm(value[take]-independent)/np.linalg.norm(independent)))
        times = receiver.dt*np.arange(len(receiver.samples))*1e9
        mats.setdefault(g['case'], np.zeros((len(times), len(stations))))[:, g['station_index']] = receiver.samples
    if dft_max > 1e-9:
        raise ValueError('independent direct DFT mismatch')
    total = mats['slope_rough']
    ideal = total - mats['slope_fullcover']
    n95 = n_cover(95e6)
    tx_arr = np.array(stations)
    thick = (np.rint(sim.surface_z(tx_arr)/sim.DX) - np.rint(sim.interface_z(tx_arr)/sim.DX))*sim.DX
    t_iface = AIR_NS + 2*thick*n95/C_SI*1e9
    ops = {}
    op_diag = {}
    for cid, label in (('B3_G1_BG', 'mean_full'), ('B5_G1_BG', 'svd_rank2')):
        try:
            r = apply_configuration(total, cid)
            ops[label] = r['output']
            op_diag[label] = {'config_id': cid, 'diagnostics': [s['diagnostics'] for s in r['steps']]}
        except Exception as exc:
            ops[label] = None
            op_diag[label] = {'config_id': cid, 'unavailable': str(exc)}
    gain = (times/1.0)**2  # declared t^2 analysis gain, t in ns
    ideal_g = ideal*gain[:, None]
    # Metrics in the interface window
    win = np.zeros((len(times), len(stations)), dtype=bool)
    for j, ti in enumerate(t_iface):
        win[:, j] = (times >= ti-8) & (times <= ti+12)
    def wrms(m):
        return np.sqrt((m*m*win).sum(axis=0)/win.sum(axis=0))
    rows = {'predicted_interface_ns_per_station': t_iface.tolist(),
            'cover_thickness_m_per_station': thick.tolist(),
            'ideal_window_rms_per_station': wrms(ideal).tolist(),
            'total_window_rms_over_ideal_median': float(np.median(wrms(total)/wrms(ideal)))}
    for label, m in ops.items():
        if m is None:
            rows[f'{label}_window_rms_over_ideal_median'] = None
            continue
        rows[f'{label}_window_rms_over_ideal_median'] = float(np.median(wrms(m)/wrms(ideal)))
        num = (m*ideal*win).sum(axis=0)
        den = np.sqrt((m*m*win).sum(axis=0)*(ideal*ideal*win).sum(axis=0))
        rows[f'{label}_window_corr_with_ideal_median'] = float(np.median(num/den))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), layout='constrained')
    # Panel 1: model cross-section
    ax = axes[0, 0]
    ks, ki = sim.column_grids()
    xs_c = sim.DX*(np.arange(len(ks))+0.5)
    ax.fill_between(xs_c, 0, ki*sim.DX, color='#c2a36b', step='mid', label='基岩 εr9')
    ax.fill_between(xs_c, ki*sim.DX, ks*sim.DX, color='#8fb996', step='mid', label='覆盖层 εr11+Debye')
    for j, x in enumerate(stations):
        z = c['geometry']['antenna_z_m'][j]
        ax.plot([x, x+sim.BASELINE], [z, z], 'k^-', ms=3, lw=.5)
    ax.set_xlim(8, 28); ax.set_ylim(14, 31.5)
    ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
    ax.set_title('模型剖面（13站位地形跟随，离地8m）', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    extent = [stations[0], stations[-1], times[-1], 0]
    v_total = abs(total).max()*0.005
    v_iso = abs(ideal).max()
    def bscan(ax, m, v, title, cmap='RdBu_r', unit='Ey(V/m)'):
        im = ax.imshow(m, aspect='auto', cmap=cmap, vmin=-v, vmax=v, extent=extent)
        ax.plot(tx_arr, t_iface, 'k--', lw=1, label='界面双程预测(95MHz)')
        ax.axhline(AIR_NS, color='k', ls=':', lw=.8, label='地表回波预测(8m)')
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=6, loc='lower right')
        fig.colorbar(im, ax=ax, shrink=.8, label=unit)
        ax.set_ylim(200, 0)
    bscan(axes[0, 1], total, v_total, '总场（共享0.5%标尺）')
    bscan(axes[0, 2], ideal, v_iso, '理想去背景（坡模型−全覆盖背景）')
    if ops['mean_full'] is not None:
        bscan(axes[1, 0], ops['mean_full'], abs(ops['mean_full']).max(), '算子去背景：均值道(B3_G1_BG)')
    else:
        axes[1, 0].text(.5, .5, '均值道不可用', ha='center'); axes[1, 0].set_title('算子去背景：均值道')
    if ops['svd_rank2'] is not None:
        bscan(axes[1, 1], ops['svd_rank2'], abs(ops['svd_rank2']).max(), '算子去背景：SVD rank-2(B5_G1_BG)')
    else:
        axes[1, 1].text(.5, .5, 'SVD不可用', ha='center'); axes[1, 1].set_title('算子去背景：SVD rank-2')
    bscan(axes[1, 2], ideal_g, abs(ideal_g).max(), f'理想去背景+t²增益（t=200ns处×{gain[-1]:.3g}，分析处理非保幅）')
    for ax in (axes[1, 0], axes[1, 1], axes[1, 2]):
        ax.set_xlabel('发射站位 x(m)，接收=站位+1.3m')
    for ax in (axes[0, 1], axes[1, 0]):
        ax.set_ylabel('时间(ns)')
    fig.suptitle('2D缓坡地形B-scan：覆盖层ε∞11/σ0.001/Debye0.5，基岩εr9；地表12–24m降2m(~9.5°)，界面半坡+归档0.8m起伏；'
                 '增益曲线g(t)=(t/1ns)²已声明，属分析显示非保幅产品')
    fig.savefig(a.out/'slope_bscan_panels.png', dpi=140)
    plt.close(fig)
    # GIF from dense snapshots
    gif_bytes, n_frames = render_gif(a.study, a.out)
    summary = {'status': 'PASS', 'contract_sha256': sha256(a.study/'execution_contract.json'),
               'analysis_code_sha256': sha256(__file__), 'independent_DFT_relative_L2_max': dft_max,
               'metrics': rows, 'operator_diagnostics': op_diag,
               'gain_declaration': 'g(t)=(t/1ns)^2 applied to ideal background-removed B-scan; analysis display product, not amplitude fidelity',
               'isolation_definition': 'ideal = slope_rough minus slope_fullcover (declared target/background pair, same terrain)',
               'gif_frames': n_frames, 'gif_bytes': gif_bytes,
               'scope': 'single 2D synthetic gentle-slope mechanism model; terrain-following 8 m AGL; not site adaptation'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'dft_max': dft_max, 'metrics': rows,
                      'gif_frames': n_frames, 'gif_bytes': gif_bytes}, ensure_ascii=False, indent=2))


def render_gif(study, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from PIL import Image
    snap_dir = study/sim.SNAP_GROUP/'profile_snaps'
    snaps = sorted(snap_dir.glob('snap*.h5'))
    dt = sim.DX/(C_SI*np.sqrt(2))
    nx = round((sim.SNAP_REGION[3]-sim.SNAP_REGION[0])/sim.SNAP_STEP[0])
    nz = round((sim.SNAP_REGION[5]-sim.SNAP_REGION[2])/sim.SNAP_STEP[2])
    xs = sim.SNAP_REGION[0] + sim.SNAP_STEP[0]*(np.arange(nx)+0.5)
    zs = sim.SNAP_REGION[2] + sim.SNAP_STEP[2]*(np.arange(nz)+0.5)
    xq = np.linspace(sim.SNAP_REGION[0], sim.SNAP_REGION[3], 400)
    z_surf = np.rint(sim.surface_z(xq)/sim.DX)*sim.DX
    z_iface = np.rint(sim.interface_z(xq)/sim.DX)*sim.DX
    peak = 0.
    fields = []
    for f in snaps:
        with h5py.File(f) as h:
            ey = h['Ey'][:, 0, :]
            fields.append(ey)
            peak = max(peak, float(np.abs(ey).max()))
    pngs = []
    for k, ey in enumerate(fields):
        t = k*sim.SNAP_ITER_STEP*dt*1e9
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        norm = SymLogNorm(linthresh=peak*1e-3, vmin=-peak, vmax=peak)
        axes[0].pcolormesh(xs, zs, ey.T, cmap='RdBu_r', norm=norm, shading='auto')
        axes[0].set_title(f'Ey 总场 (SymLog 固定±{peak:.3g})', fontsize=9)
        env = np.log10(np.maximum(np.abs(ey), peak*1e-12)/peak)
        axes[1].pcolormesh(xs, zs, env.T, cmap='viridis', vmin=-8, vmax=0, shading='auto')
        axes[1].set_title('|Ey| log10 包络', fontsize=9)
        for ax in axes:
            ax.plot(xq, z_surf, 'w-', lw=1.)
            ax.plot(xq, z_iface, color='magenta', ls='--', lw=.8)
            ax.plot([18.0, 19.3], [28.5, 28.5], 'y^-', ms=4)
            ax.set_ylim(31.4, 14)
            ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
        fig.suptitle(f'2D缓坡波场快照（白=地表，品红虚线=基覆界面，黄=收发） t = {t:.1f} ns ({k+1}/{len(fields)})')
        fig.tight_layout()
        p = out/f'frame_{k:04d}.png'
        fig.savefig(p, dpi=80)
        plt.close(fig)
        pngs.append(p)
    images = [Image.open(p) for p in pngs]
    gif = out/'slope_wavefield.gif'
    images[0].save(gif, save_all=True, append_images=images[1:], duration=90, loop=0)
    for p in pngs:
        p.unlink()
    return gif.stat().st_size, len(pngs)


if __name__ == '__main__':
    main()
