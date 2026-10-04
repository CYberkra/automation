"""Diagnose why the bedrock-cover interface is hard to see in the 3D B-scan.

Checks:
1. Expected arrival times from geometry (direct / surface / interface).
2. Spectrum: full-band impulse vs device band 20-170 MHz.
3. Band-limited (20-170 MHz) B-scan raw + mean-removed, with event lines.
4. Coherence/energy in the interface prediction window per trace.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

CHECK = Path(r'E:\automation_djh\artifacts_check\r1')
OUTDIR = Path(r'E:\automation_djh')
B = np.load(CHECK / 'bscan_3D_5cm.npy')          # 13 x 6233
t_ns = np.load(CHECK / 't_ns_3D_5cm.npy')
ntr, nsteps = B.shape
dt = (t_ns[1] - t_ns[0]) * 1e-9
fs = 1 / dt
traces = np.arange(1, ntr + 1)

c = 299792458.0
# geometry (from frozen .in): Tx (6,3.85,27), Rx (6,5.15,27); surface z=12; cover/rock z~9
h_air = 27 - 12.0          # 15 m aircraft height
off = 5.15 - 3.85          # 1.3 m common offset
z_cover = 12.0 - 9.0       # ~3 m cover
er_cover = 18.017
v_cover = c / np.sqrt(er_cover)
t_direct = off / c * 1e9
t_surf = 2 * np.sqrt(h_air**2 + (off/2)**2) / c * 1e9
t_iface = t_surf + 2 * z_cover / v_cover * 1e9
print(f'expected: direct {t_direct:.1f} ns | surface {t_surf:.1f} ns | interface ~{t_iface:.1f} ns')

# ---- spectrum of a representative trace ----
a = B[6]
f, Pxx = signal.welch(a, fs=fs, nperseg=2048)
band = (f >= 20e6) & (f <= 170e6)
frac_band = float(np.trapezoid(Pxx[band], f[band]) / np.trapezoid(Pxx, f))
print(f'energy fraction inside 20-170 MHz: {frac_band*100:.1f}%')

# ---- band-limit 20-170 MHz (zero-phase, 4th order butterworth) ----
sos = signal.butter(4, [20e6, 170e6], btype='band', fs=fs, output='sos')
Bb = signal.sosfiltfilt(sos, B, axis=1)
mean_b = Bb.mean(axis=0, keepdims=True)
Br = Bb - mean_b

sel = t_ns <= 300
fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.4), dpi=130, sharey=True)
for ax, (M, ttl) in zip(axes, [
        (Bb[:, sel], 'band-limited 20-170 MHz (raw)'),
        (Br[:, sel], 'band-limited, mean-trace removed'),
        (B[:, sel], 'full-band impulse (ref, clipped 3%)')]):
    ref = np.max(np.abs(M)) or 1.0
    Mn = M / ref if M is not B[:, sel] else np.clip(M / (ref * 0.03), -1, 1)
    pc = ax.pcolormesh(traces, t_ns[sel], Mn.T, vmin=-1, vmax=1, cmap='RdBu_r', shading='auto')
    for tt, lab, col in [(t_direct, 'direct', 'k'), (t_surf, 'surface', 'g'),
                         (t_iface, 'interface', 'm')]:
        ax.axhline(tt, color=col, lw=1.0, ls='--')
        ax.text(13.1, tt, lab, color=col, fontsize=8, va='center')
    ax.set_xlabel('trace #'); ax.set_title(ttl, fontsize=10)
    ax.set_xticks(traces)
    fig.colorbar(pc, ax=ax)
axes[0].set_ylabel('two-way travel time (ns)')
axes[0].set_ylim(300, 0)
fig.suptitle(f'3D group B-scan | device-band energy = {frac_band*100:.0f}% of full-band impulse')
fig.tight_layout()
fig.savefig(OUTDIR / 'fig_3d_bscan_bandlimited.png')
print('saved fig_3d_bscan_bandlimited.png')

# ---- interface-window coherence (band-limited) ----
w = (t_ns >= t_iface - 25) & (t_ns <= t_iface + 25)
seg = Bb[:, w]
seg = seg - seg.mean(axis=1, keepdims=True)
nrm = np.linalg.norm(seg, axis=1)
seg = seg / np.where(nrm > 0, nrm, 1)[:, None]
C = seg @ seg.T
iu = np.triu_indices(ntr, 1)
print(f'interface-window ({t_iface-25:.0f}-{t_iface+25:.0f} ns) corr: mean {C[iu].mean():.3f} min {C[iu].min():.3f}')
# noise floor just before window
wn = (t_ns >= t_iface - 60) & (t_ns <= t_iface - 30)
snr = float(np.sqrt(np.mean(Bb[:, w]**2)) / np.sqrt(np.mean(Bb[:, wn]**2)))
print(f'interface-window RMS / pre-window noise RMS (band-limited): {snr:.2f}')
# peak per trace in interface window vs bandlimited global
for i in range(ntr):
    pk = np.max(np.abs(Bb[i, w]))
    print(f'  t{i+1:02d}: peak {pk:.3e} ({20*np.log10(pk/np.max(np.abs(Bb))+1e-12):.1f} dB re direct)')
