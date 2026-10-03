"""Render benchmark B-scans, per-trace normalized to the direct arrival.

House convention: no absolute amplitudes across 2D/3D; each trace is
normalized by its direct-wave window peak so event structure is visible.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CHECK = Path(r'E:\automation_djh\artifacts_check\r1')
OUT = Path(r'E:\automation_djh\fig_3d_bscan_review.png')

panels = []
for tag, gname in [('3D_5cm', '3D 5cm (38M cells)'),
                   ('2D_5cm', '2D 5cm same section')]:
    p = CHECK / f'bscan_{tag}.npy'
    if p.exists():
        B = np.load(p)
        y = np.load(CHECK / f'rx_y_{tag}.npy')
        t_ns = np.load(CHECK / f't_ns_{tag}.npy')
        panels.append((B, y, t_ns, gname))

fig, axes = plt.subplots(1, len(panels), figsize=(6.4 * len(panels), 4.8), dpi=130)
if len(panels) == 1:
    axes = [axes]
for ax, (B, y, t_ns, gname) in zip(axes, panels):
    direct = B[:, (t_ns >= 0) & (t_ns < 50)]
    ref = np.max(np.abs(direct), axis=1, keepdims=True)
    ref = np.where(ref > 0, ref, 1.0)
    im = 20 * np.log10(np.abs(B) / ref + 1e-12)
    pc = ax.pcolormesh(t_ns, y, im, vmin=-60, vmax=0, cmap='viridis', shading='auto')
    ax.set_xlabel('t (ns)')
    ax.set_title(gname + ' (per-trace norm to direct)')
    ax.set_xlim(0, 400)
    ax.axvline(100, color='white', lw=0.6, ls='--', alpha=0.6)
    ax.axvline(185, color='red', lw=0.6, ls='--', alpha=0.6)
    fig.colorbar(pc, ax=ax, label='dB re direct peak')
axes[0].set_ylabel('rx y (m)')
axes[0].text(102, axes[0].get_ylim()[1], 'surf~100ns', color='white', fontsize=7, va='top')
axes[0].text(187, axes[0].get_ylim()[1], 'iface~185ns', color='red', fontsize=7, va='top')
fig.suptitle('benchmark3d_r2_co B-scans (13 common-offset traces)')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
