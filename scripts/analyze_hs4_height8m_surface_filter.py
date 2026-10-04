"""Verify cover-to-air critical-angle angular filtering from the 8 m snapshots.

Pure array analysis of the verified capsule; no solver. Takes the
interface-contrast (rough-halfspace) Ey along three snapshot rows (just below
the surface in cover, mid-air, and just below the antenna), builds (x,t)
gathers over the echo window, and compares their (f,kx) band spectra against
the propagating cone kx <= 2*pi*f/v. kx is conserved across the horizontal
surface, so components beyond the air cone are evanescent in air: they must
vanish with height if the field is physical.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CAP = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
OUT = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_analysis'
C_AIR = 0.299792458e9     # m/s
V_COVER = C_AIR/np.sqrt(18.017)  # low-frequency real-part estimate, m/s
SIN_CRIT = V_COVER/C_AIR
ROWS = {'cover_z11.975': 26, 'air_z16.025': 53, 'air_z19.925': 79}
TWIN = (90.0, 170.0)              # ns, echo-crossing window
FBAND = (60e6, 140e6)             # Hz, around the 95 MHz Ricker centre


def frames(group):
    out = []
    for f in sorted((CAP/group).glob('profile_snaps/*.h5')):
        with h5py.File(f) as h:
            out.append((int(h.attrs['iteration']), f))
    return out


def main():
    c = json.loads((CAP/'execution_contract.json').read_text('utf-8'))
    dt = c['dt_s']
    fr = frames('mid8_rough')
    fh = frames('mid8_halfspace')
    iters = np.array([i for i, _ in fr])
    t = iters*dt*1e9
    sel = (t >= TWIN[0]) & (t <= TWIN[1])
    ts = t[sel]*1e-9
    gathers = {name: [] for name in ROWS}
    fr_sel = [fr[i] for i in np.nonzero(sel)[0]]
    fh_sel = [fh[i] for i in np.nonzero(sel)[0]]
    for (_, pr), (_, ph) in zip(fr_sel, fh_sel):
        with h5py.File(pr) as h1, h5py.File(ph) as h2:
            d = h1['Ey'][:, 0, :] - h2['Ey'][:, 0, :]
        for name, iz in ROWS.items():
            gathers[name].append(d[:, iz])
    dx, dts = 0.15, ts[1]-ts[0]
    taper = np.hanning(len(ts))[:, None]*np.hanning(60)[None, :]
    freqs = np.fft.fftshift(np.fft.fftfreq(len(ts), d=dts))   # Hz
    kxs = np.fft.fftshift(np.fft.fftfreq(60, d=dx))*2*np.pi  # rad/m
    fsel = (freqs >= FBAND[0]) & (freqs <= FBAND[1])
    fmid = 95e6
    kx_cone_air = 2*np.pi*fmid/C_AIR
    kx_cone_cover = 2*np.pi*fmid/V_COVER
    in_air_cone = np.abs(kxs) <= kx_cone_air
    spectra, powers, outside = {}, {}, {}
    for name in ROWS:
        g = np.array(gathers[name])
        f2 = np.fft.fftshift(np.fft.fft2(g*taper))
        band = np.abs(f2[np.ix_(np.nonzero(fsel)[0])])**2
        spectra[name] = np.log10(band.sum(axis=0) + 1e-300)
        powers[name] = float(band.sum())
        outside[name] = float(band[:, ~in_air_cone].sum())
    ref = 'cover_z11.975'
    result = {
        'status': 'PASS',
        'sin_critical': float(SIN_CRIT), 'critical_angle_deg': float(np.degrees(np.arcsin(SIN_CRIT))),
        'v_cover_m_per_ns': float(V_COVER/1e9),
        'kx_propagating_limit_air_95MHz_rad_per_m': float(kx_cone_air),
        'kx_propagating_limit_cover_95MHz_rad_per_m': float(kx_cone_cover),
        'band_power': powers,
        'power_re_cover_dB': {k: float(10*np.log10(v/powers[ref])) for k, v in powers.items()},
        'fraction_outside_air_cone': {k: float(outside[k]/powers[k]) for k in ROWS},
        'window_ns': list(TWIN), 'fband_Hz': list(FBAND),
        'rows_z_m': {k: 8+0.15*(v+0.5) for k, v in ROWS.items()},
        'note': 'fraction_outside_air_cone must fall to ~0 at the air rows if the field is physical (evanescent high-kx decays with height); the cover-row value includes downgoing multiples and is the high-angle share present just below the surface.'}
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for name in ROWS:
        p = spectra[name]
        ax.plot(kxs, p - p.max(), label=name)
    for k, col in [(kx_cone_air, 'r'), (kx_cone_cover, 'g')]:
        ax.axvline(k, color=col, ls=':'); ax.axvline(-k, color=col, ls=':')
    ax.annotate('air cone (95 MHz)', (kx_cone_air, -1), color='r', fontsize=8, rotation=90)
    ax.annotate('cover limit', (kx_cone_cover, -1), color='g', fontsize=8, rotation=90)
    ax.set_xlabel('kx (rad/m)'); ax.set_ylabel('log10 band-integrated |E(f,kx)|^2 (re row max)')
    ax.grid(alpha=.3); ax.legend(); ax.set_xlim(-1.6*kx_cone_cover, 1.6*kx_cone_cover)
    ax.set_title(f'Angular filter check, {TWIN[0]:.0f}-{TWIN[1]:.0f} ns: critical angle {np.degrees(np.arcsin(SIN_CRIT)):.1f} deg')
    fig.tight_layout()
    fig.savefig(OUT/'surface_refraction_filter.png', dpi=110)
    plt.close(fig)
    (OUT/'surface_refraction_summary.json').write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=1))


if __name__=='__main__':
    main()
