"""Render device-band (SFCW) B-scan comparison + interface-window quantification.

dB reference convention ('dB re direct'): per-trace max|complex envelope| in
the 0-50 ns direct-arrival window. The colorbar and the printed interface-window
metric share this single reference.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CHECK = Path(r'E:\automation_djh\artifacts_check\sfcw_r1')
OUT = Path(r'E:\automation_djh\fig_sfcw_bscan_compare.png')
GROUPS = [('3D_5cm', '3D 5cm'), ('2D_5cm', '2D 5cm'), ('2D_2p5cm', '2D 2.5cm')]
c = 299792458.0
t_surf = 100.2
t_iface = 185.1

d = np.load(CHECK / 'sfcw_responses.npz')
fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), dpi=130, sharey=True)
iface_peaks = {}
for ax, (tag, gname) in zip(axes, GROUPS):
    mats = []
    for k in range(1, 14):
        t = f't{k:02d}'
        env = d[f'{tag}_{t}_no_taper_hann_complex_envelope']
        tm = d[f'{tag}_{t}_no_taper_hann_envelope_time_s']
        mats.append(env)
    tm = tm * 1e9
    sel = tm <= 400
    B = np.stack([m[sel] for m in mats])
    direct = B[:, (tm[sel] >= 0) & (tm[sel] < 50)]
    ref = np.max(np.abs(direct), axis=1, keepdims=True)
    ref = np.where(ref > 0, ref, 1.0)
    im = 20 * np.log10(np.abs(B) / ref + 1e-12)
    pc = ax.pcolormesh(np.arange(1, 14), tm[sel], im.T, vmin=-60, vmax=0, cmap='viridis', shading='auto')
    for tt, lab, col in [(t_surf, 'surf', 'w'), (t_iface, 'iface', 'r')]:
        ax.axhline(tt, color=col, lw=1, ls='--')
        ax.text(13.1, tt, lab, color=col, fontsize=8, va='center')
    ax.set_xlabel('trace #'); ax.set_title(f'{gname} (SFCW 20-170 MHz, Hann)', fontsize=10)
    ax.set_xticks(range(1, 14))
    fig.colorbar(pc, ax=ax, label='dB re direct')
    w = (tm[sel] >= t_iface - 25) & (tm[sel] <= t_iface + 25)
    peaks = np.max(np.abs(B[:, w]), axis=1)
    iface_peaks[tag] = 20 * np.log10(peaks / ref[:, 0] + 1e-12)
axes[0].set_ylabel('two-way travel time (ns)')
axes[0].set_ylim(400, 0)
fig.suptitle('benchmark3d_r2_co device-band B-scans (official SFCW chain)')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
for tag, _ in GROUPS:
    p = iface_peaks[tag]
    print(f'{tag}: interface-window peak {p.mean():.1f} dB re direct (range {p.min():.1f}..{p.max():.1f})')
