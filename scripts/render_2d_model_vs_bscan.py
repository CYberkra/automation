"""2x2 panel: per 2D group - model cross-section (left) vs SFCW device-band B-scan (right).

Geometry parsed from each group's t01 .in (the two 2D groups share identical
geometry; only the Yee grid differs: 0.05 m vs 0.025 m).
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

CFG = Path(r'E:\automation_djh\automation_repo\configs\research\benchmark3d_r2_co')
SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_2d_model_vs_bscan.png')

GROUPS = [
    ('B2D5CM-C3mR2-BG-CO13', '2D_5cm', '2D 5cm (grid 0.05 m)'),
    ('B2D-C3mR2-BG-CO13', '2D_2p5cm', '2D 2.5cm (grid 0.025 m)'),
]
t_surf, t_iface = 100.2, 185.1


def parse_in(prefix):
    boxes = []
    src = rx = None
    dx = None
    for line in (CFG / f'{prefix}-t01.in').read_text().splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == '#box:':
            nums = []
            for v in p[1:]:
                if v == 'inf':
                    continue
                try:
                    nums.append(float(v))
                except ValueError:
                    pass
            if len(nums) == 4:
                boxes.append((nums[0], nums[1], nums[2], nums[3], p[-1]))
        elif p[0] == '#hertzian_dipole:':
            src = (float(p[3]), float(p[4]))
        elif p[0] == '#rx:':
            rx = (float(p[2]), float(p[3]))
        elif p[0] == '#dx_dy_dz:':
            dx = float(p[1])
    return boxes, src, rx, dx


d = np.load(SFCW)
fig, axes = plt.subplots(2, 2, figsize=(17, 10.5), dpi=130)

for row, (prefix, tag, gname) in enumerate(GROUPS):
    boxes, src, rx0, dx = parse_in(prefix)

    # ---- left: model cross-section ----
    ax = axes[row, 0]
    ax.set_facecolor('white')
    # air (nothing to draw); subsurface materials
    for y0, z0, y1, z1, mat in boxes:
        color = '#b0b0b0' if mat == 'rock' else '#c8863c'
        ax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0, facecolor=color,
                               edgecolor='none', zorder=1))
    # surface line
    ax.plot([0, 12], [12, 12], 'k-', lw=1.4, zorder=3)
    # PML frame (dashed), domain y 0..12, z 0..33 (top beyond view)
    ax.plot([0, 0, 12, 12], [0, 28, 28, 0], ls=(0, (5, 3)), color='purple',
            lw=1.2, zorder=2)
    ax.text(0.15, 4.4, 'PML', color='purple', fontsize=8, va='top')
    # material labels
    ax.text(7.5, 5.0, 'rock  er=9', ha='center', va='center', fontsize=11, color='k')
    ax.text(7.5, 10.9, 'cover  er=18 (Debye)', ha='center', va='center',
            fontsize=9, color='k')
    ax.text(11.8, 21.5, 'air  er=1', ha='right', va='center', fontsize=11, color='k')
    ax.annotate('rough interface\n(2.95-3.3 m deep)',
                xy=(9.3, 8.88), xytext=(11.4, 7.2), fontsize=8.5, color='darkred',
                ha='right', va='center',
                arrowprops=dict(arrowstyle='->', color='darkred',
                                connectionstyle='arc3,rad=-0.25'))
    # source/receiver of t01 + CO profile line
    ax.plot([src[0]], [src[1]], marker='*', ms=16, color='red', zorder=5)
    ax.plot([rx0[0]], [rx0[1]], marker='v', ms=10, color='blue', zorder=5)
    ys = np.linspace(3.85, 6.85, 13)
    ax.plot(ys, np.full(13, 27.0), '-', color='gray', lw=0.8, zorder=4)
    ax.plot(ys, np.full(13, 27.0), 'r*', ms=5, zorder=5)
    ax.plot(ys + 1.3, np.full(13, 27.0), 'bv', ms=4, zorder=5)
    ax.annotate('Tx (t01)', xy=src, xytext=(0.3, 24.2), fontsize=9, color='red',
                arrowprops=dict(arrowstyle='->', color='red'))
    ax.annotate('Rx (t01)', xy=rx0, xytext=(9.0, 24.2), fontsize=9, color='blue',
                arrowprops=dict(arrowstyle='->', color='blue'))
    ax.text(6.0, 26.1, 'CO profile: 13 traces, offset 1.3 m\n'
            'dy=0.25 m, height 15 m', fontsize=8, ha='center', va='top',
            color='dimgray')
    # axes cosmetics: z up in model coords -> invert so surface on top
    ax.set_xlim(0, 12)
    ax.set_ylim(28, 4)
    ax.set_xlabel('y (m)')
    ax.set_ylabel('z (m, model coords; surface z=12)')
    ax.set_title(f'{gname} | model section (parsed from {prefix}-t01.in)',
                 fontsize=10)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(alpha=0.15)
    # inset: Yee grid zoom at interface (upper-left, empty rock region)
    iax = ax.inset_axes([0.04, 0.60, 0.42, 0.34])
    iax.set_facecolor('white')
    for y0, z0, y1, z1, mat in boxes:
        if y1 < 4.0 or y0 > 5.2 or z1 < 8.2 or z0 > 9.4:
            continue
        color = '#b0b0b0' if mat == 'rock' else '#c8863c'
        iax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0, facecolor=color,
                                edgecolor='none'))
    gx = np.arange(4.0, 5.2 + 1e-9, dx)
    gz = np.arange(8.2, 9.4 + 1e-9, dx)
    for g in gx:
        iax.axvline(g, color='k', lw=0.2, alpha=0.6)
    for g in gz:
        iax.axhline(g, color='k', lw=0.2, alpha=0.6)
    iax.set_xlim(4.0, 5.2)
    iax.set_ylim(9.4, 8.2)
    iax.set_title(f'Yee grid dx={dx} m', fontsize=8)
    iax.set_xticks([]); iax.set_yticks([])

    # ---- right: SFCW device-band B-scan ----
    ax = axes[row, 1]
    mats = []
    tm = None
    for k in range(1, 14):
        t = f't{k:02d}'
        mats.append(d[f'{tag}_{t}_no_taper_hann_complex_envelope'])
        tm = d[f'{tag}_{t}_no_taper_hann_envelope_time_s']
    tm = tm * 1e9
    sel = tm <= 400
    B = np.stack([m[sel] for m in mats])
    direct = B[:, (tm[sel] >= 0) & (tm[sel] < 50)]
    ref = np.max(np.abs(direct), axis=1, keepdims=True)
    ref = np.where(ref > 0, ref, 1.0)
    im = 20 * np.log10(np.abs(B) / ref + 1e-12)
    pc = ax.pcolormesh(np.arange(1, 14), tm[sel], im.T, vmin=-60, vmax=0,
                       cmap='viridis', shading='auto')
    for tt, lab, col in [(t_surf, 'surf', 'w'), (t_iface, 'iface', 'r')]:
        ax.axhline(tt, color=col, lw=1, ls='--')
        ax.text(13.1, tt, lab, color=col, fontsize=8, va='center')
    ax.set_xlabel('trace #')
    ax.set_ylim(400, 0)
    ax.set_title(f'{gname} | simulated B-scan (SFCW 20-170 MHz, Hann,'
                 ' dB re direct)', fontsize=10)
    ax.set_xticks(range(1, 14))
    fig.colorbar(pc, ax=ax, label='dB re direct')

fig.suptitle('benchmark3d_r2_co: 2D model section vs device-band B-scan'
             ' (identical geometry, two grids)')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
