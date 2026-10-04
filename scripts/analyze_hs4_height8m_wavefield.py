"""Analyse the dense ~1 ns 8 m height wavefields; receivers compared across heights.

Pure array analysis of completed, verified capsules; no solver. For the 8 m
pair: band-resolved rough-minus-halfspace Ey energy versus time (when and
where the interface-contrast field exists), footprint of that field on the
undulating interface, and a dense difference-field atlas across the echo
formation/return window. Across heights (2/8/15 m): native receiver traces
and their rough-minus-halfspace differences, since only the 8 m raw snapshot
frames survive locally (15 m/2 m frames were retained on the previous machine
only; their receiver H5 remain). Contract hashes recorded for provenance.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
CAP8 = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
CAPH = ROOT/'artifacts/research_checks/2026-10-04_hs4_height_wavefield_continuation'
OUT = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_analysis'
ORIGIN_Z = 8.0
SPACING = np.array([0.15, 0.05, 0.15])
XS = 13.5 + SPACING[0]*(np.arange(60)+0.5)
ZS = ORIGIN_Z + SPACING[2]*(np.arange(134)+0.5)
BANDS = {'interface': (8.5, 9.6), 'cover': (9.6, 12.0), 'air_low': (12.0, 19.0), 'antenna': (19.0, 21.0)}
C = 0.299792458  # m/ns


def interface_relief():
    xs, zs = [], []
    for line in (ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre/centre_rough/profile.in').read_text('utf-8').splitlines():
        if not line.startswith('#box:'):
            continue
        p = [float(v) for v in line.split()[1:7]]  # x0 y0 z0 x1 y1 z1
        if p[5] == 12.0 and p[1] == 0.0 and p[2] > 5.0:  # cover boxes only, not the rock base box
            xs.append(0.5*(p[0]+p[3])); zs.append(p[2])
    order = np.argsort(xs)
    return np.array(xs)[order], np.array(zs)[order]


IFACE_X, IFACE_Z = interface_relief()


def band_mask(z0, z1):
    return (ZS>=z0)&(ZS<z1)


def load_ey(path):
    with h5py.File(path) as h:
        return h['Ey'][:, 0, :]


def frames(group_dir):
    out = []
    for f in sorted(group_dir.glob('profile_snaps/*.h5')):
        with h5py.File(f) as h:
            out.append((int(h.attrs['iteration']), f))
    return out


def receiver_trace(group_dir):
    with h5py.File(group_dir/'profile.h5') as h:
        return h['rxs/rx1/Ey'][:], float(h.attrs['dt'])


def main():
    if OUT.exists():
        raise SystemExit('new analysis capsule required')
    c8 = json.loads((CAP8/'execution_contract.json').read_text('utf-8'))
    ch = json.loads((CAPH/'execution_contract.json').read_text('utf-8'))
    dt = c8['dt_s']
    fr = frames(CAP8/'mid8_rough')
    fh = frames(CAP8/'mid8_halfspace')
    assert [i for i,_ in fr]==[i for i,_ in fh]
    iters = np.array([i for i,_ in fr])
    t = iters*dt*1e9
    band_e = {k: [] for k in BANDS}
    footprint = np.zeros(60)
    global_peak = 0.0
    for (_, pr), (_, ph) in zip(fr, fh):
        diff = load_ey(pr) - load_ey(ph)
        global_peak = max(global_peak, float(np.abs(diff).max()))
        e = diff**2
        for k,(z0,z1) in BANDS.items():
            band_e[k].append(float(e[:, band_mask(z0,z1)].sum()))
        footprint = np.maximum(footprint, np.abs(diff[:, band_mask(*BANDS['interface'])]).max(axis=1))
    # Receiver traces across heights.
    rx = {}
    for h, cap, names in [(8, CAP8, ('mid8_rough','mid8_halfspace')), (15, CAPH, ('high_rough','high_halfspace')), (2, CAPH, ('low_rough','low_halfspace'))]:
        a, _ = receiver_trace(cap/names[0])
        b, dth = receiver_trace(cap/names[1])
        assert dth == dt
        rx[h] = {'rough': a, 'halfspace': b, 'diff': a-b}
    fp = footprint
    active = XS[fp >= 0.2*fp.max()]
    summary = {'status':'PASS','inputs':{
        '8m_contract_sha256': sha256(CAP8/'execution_contract.json'),
        'heights_contract_sha256': sha256(CAPH/'execution_contract.json'),
        'note':'15m/2m raw snapshot frames no longer on this machine (previous-machine retention); cross-height comparison uses surviving native receiver H5 only.'},
        'frame_interval_ns': float(np.median(np.diff(t))), 'frame_count': len(t)}
    summary['interface_footprint_8m'] = {'x_min_m': float(active.min()), 'x_max_m': float(active.max()),
        'width_m': float(active.max()-active.min())}
    summary['band_peaks_8m_ns'] = {k: float(t[int(np.argmax(v))]) for k, v in band_e.items()}
    # First time each band's diff energy rises above a numerical-noise-safe threshold.
    floor = (global_peak*1e-12)**2 * 60*134  # band-integrated floor proxy
    summary['global_diff_abs_max'] = global_peak
    summary['band_first_above_floor_ns'] = {}
    for k, v in band_e.items():
        v = np.array(v)
        above = np.nonzero(v > floor)[0]
        summary['band_first_above_floor_ns'][k] = float(t[above[0]]) if len(above) else None
    summary['predicted_vertical_interface_arrival_ns'] = 2*8/C + 2*3.0/(C/np.sqrt(18.017))
    for h in (2, 8, 15):
        d = rx[h]['diff']
        full_peak = float(np.abs(rx[h]['rough']).max())
        deep = np.abs(d)[int(100e-9/dt):int(220e-9/dt)].max()
        # Sliding 5 ns envelope peak after 60 ns: interface-echo arrival per height.
        k = int(2.5e-9/dt)
        env = np.array([np.abs(d)[max(0,i-k):i+k].max() for i in range(0, len(d), k)])
        te = np.arange(0, len(d), k)*dt*1e9
        late = te > 60
        j = int(np.argmax(env[late]))
        summary[f'height_{h}m_receiver'] = {
            'rough_full_record_peak': full_peak,
            'diff_peak_100_220ns': float(deep),
            'diff_over_full_peak_dB': float(20*np.log10(deep/full_peak)),
            'echo_envelope_peak_ns': float(te[late][j]),
            'echo_envelope_peak_amplitude': float(env[late][j])}
    summary['expected_two_way_shift_per_7m_ns'] = 2*7.0/C
    OUT.mkdir(parents=True)
    # Figure 1: 8 m band energy vs time, dB re per-band max, floored for readability.
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for ax, (k, v) in zip(axes, band_e.items()):
        v = np.array(v)
        db = 10*np.log10(np.maximum(v, v.max()*1e-30)/v.max())
        ax.plot(t, db)
        ax.set_ylabel(k); ax.grid(alpha=.3); ax.set_xlim(0, 300); ax.set_ylim(-140, 3)
        ax.axvline(2*8/C, color='r', ls=':', lw=.8)
        pk = summary['band_peaks_8m_ns'][k]
        ax.axvline(pk, color='g', ls=':', lw=.8)
        ax.annotate(f'peak {pk:.0f} ns', (pk, -10), fontsize=8, rotation=90)
    axes[0].set_title('8 m rough-minus-halfspace Ey band energy, dB re band max; red dotted = vertical ground two-way')
    axes[-1].set_xlabel('ns')
    fig.tight_layout()
    fig.savefig(OUT/'height8m_band_energy.png', dpi=110)
    plt.close(fig)
    # Figure 2: footprint with true relief.
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(XS, fp/fp.max())
    ax2 = ax.twinx(); ax2.plot(IFACE_X, IFACE_Z, 'k--', alpha=.5)
    ax2.set_ylabel('interface z (m)'); ax2.invert_yaxis()
    ax.set_xlabel('x (m)'); ax.set_ylabel('normalised peak |diff Ey| in interface band')
    ax.set_title('8 m interface-contrast footprint vs true relief')
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(OUT/'height8m_interface_footprint.png', dpi=110)
    plt.close(fig)
    # Figure 3: dense difference-field atlas, fixed global scale so decay is visible.
    picks = [30, 50, 53, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200]
    vmax = global_peak*0.2
    fig, axes = plt.subplots(4, 4, figsize=(17, 11))
    for ax, want in zip(axes.ravel(), picks):
        k = int(np.argmin(np.abs(t - want)))
        diff = load_ey(fr[k][1]) - load_ey(fh[k][1])
        ax.pcolormesh(XS, ZS, diff.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
        ax.plot(IFACE_X, IFACE_Z, 'k-', lw=.7)
        ax.axhline(12, color='k', ls=':', lw=.5)
        ax.plot([17.6, 18.9], [20, 20], 'k^', ms=3)
        ax.set_title(f'{t[k]:.1f} ns', fontsize=9)
        ax.set_ylim(22, 8)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle('8 m rough-minus-halfspace Ey, FIXED colour scale (±20%% of global diff peak); interface/surface/antenna marked')
    fig.tight_layout()
    fig.savefig(OUT/'height8m_diff_atlas.png', dpi=110)
    plt.close(fig)
    # Figure 3b: same frames, per-frame normalised log envelope to show weak propagation paths.
    fig, axes = plt.subplots(4, 4, figsize=(17, 11))
    for ax, want in zip(axes.ravel(), picks):
        k = int(np.argmin(np.abs(t - want)))
        diff = load_ey(fr[k][1]) - load_ey(fh[k][1])
        env = np.log10(np.maximum(np.abs(diff), global_peak*1e-12)/global_peak)
        ax.pcolormesh(XS, ZS, env.T, cmap='viridis', vmin=-8, vmax=0, shading='auto')
        ax.plot(IFACE_X, IFACE_Z, 'r-', lw=.7)
        ax.set_title(f'{t[k]:.1f} ns', fontsize=9)
        ax.set_ylim(22, 8)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle('8 m |rough-minus-halfspace Ey|, log10 re global diff peak, floor 1e-12')
    fig.tight_layout()
    fig.savefig(OUT/'height8m_diff_atlas_log.png', dpi=110)
    plt.close(fig)
    # Figure 4: receiver traces across heights, linear + log |diff|.
    fig, axes = plt.subplots(3, 2, figsize=(15, 9), sharex='col')
    tt = np.arange(len(rx[8]['rough']))*dt*1e9
    for row, h in enumerate((2, 8, 15)):
        ax = axes[row, 0]
        ax.plot(tt, rx[h]['rough'], lw=.6, label='rough')
        ax.plot(tt, rx[h]['diff'], lw=.6, label='rough-halfspace')
        ax.set_ylabel(f'h={h} m'); ax.set_xlim(0, 300); ax.grid(alpha=.3); ax.legend(fontsize=8)
        ax = axes[row, 1]
        d = np.abs(rx[h]['diff'])
        ax.semilogy(tt, np.maximum(d, d.max()*1e-14)/d.max(), lw=.8)
        ax.set_xlim(0, 300); ax.set_ylim(1e-8, 2); ax.grid(alpha=.3)
        ax.set_ylabel('|diff| re max')
    axes[-1, 0].set_xlabel('ns'); axes[-1, 1].set_xlabel('ns')
    fig.suptitle('Native receiver Ey (t01): linear traces (left), normalised |rough-halfspace| log (right)')
    fig.tight_layout()
    fig.savefig(OUT/'receiver_traces_by_height.png', dpi=110)
    plt.close(fig)
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=1)[:2500])


if __name__=='__main__':
    main()
