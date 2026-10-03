"""Render benchmark3d_r2_co B-scan figures from assembled .npy products."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CHECK = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r'E:\automation_djh\artifacts_check\r1')
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(r'E:\automation_djh\fig_benchmark3d_r2_bscans.png')

panels = []
for tag, gname in [('3D_5cm', '3D 5cm (38M cells)'),
                   ('2D_5cm', '2D 5cm same section'),
                   ('2D_2p5cm', '2D 2.5cm production grid')]:
    p = CHECK / f'bscan_{tag}.npy'
    if p.exists():
        B = np.load(p)
        ypath = CHECK / f'rx_y_{tag}.npy'
        tpath = CHECK / f't_ns_{tag}.npy'
        y = np.load(ypath) if ypath.exists() else np.arange(B.shape[0])
        t_ns = np.load(tpath) if tpath.exists() else np.arange(B.shape[1])
        panels.append((B, y, t_ns, gname))

if not panels:
    print('no bscan arrays found'); sys.exit(1)

fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.6), dpi=130, sharey=True)
if len(panels) == 1:
    axes = [axes]
for ax, (B, y, t_ns, gname) in zip(axes, panels):
    ref = np.percentile(np.abs(B), 99.5) or 1.0
    im = 20 * np.log10(np.abs(B) / ref + 1e-12)
    pc = ax.pcolormesh(t_ns, y, im, vmin=-60, vmax=0, cmap='viridis', shading='auto')
    ax.set_xlabel('t (ns)')
    ax.set_title(gname)
    ax.set_xlim(0, 400)
    fig.colorbar(pc, ax=ax, label='dB (ref p99.5)')
axes[0].set_ylabel('rx y (m)')
fig.suptitle('benchmark3d_r2_co B-scans (13 common-offset traces)')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
