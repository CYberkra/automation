"""3D group: model views (y-z section at x=6.0 + plan-view interface inset) vs SFCW B-scan."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

CFG = Path(r'E:\automation_djh\automation_repo\configs\research\benchmark3d_r2_co')
SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_3d_model_vs_bscan.png')

PREFIX = 'B3D5CM-C3mR2-BG-CO13'
TAG = '3D_5cm'
SLICE_X = 6.0
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
            if len(nums) == 6:
                boxes.append(tuple(nums) + (p[-1],))
        elif p[0] == '#hertzian_dipole:':
            src = (float(p[2]), float(p[3]), float(p[4]))
        elif p[0] == '#rx:':
            rx = (float(p[1]), float(p[2]), float(p[3]))
        elif p[0] == '#dx_dy_dz:':
            dx = float(p[1])
    return boxes, src, rx, dx


boxes, src, rx0, dx = parse_in(PREFIX)
slice_boxes = [b for b in boxes if b[0] <= SLICE_X < b[3]]

# plan-view interface depth from pillar tops
tops = {}
for x0, y0, z0, x1, y1, z1, mat in boxes:
    if mat == 'cover':
        tops[(x0, y0)] = z0
xs = sorted({k[0] for k in tops})
ys = sorted({k[1] for k in tops})
Z = np.full((len(ys), len(xs)), np.nan)
for j, y0 in enumerate(ys):
    for i, x0 in enumerate(xs):
        Z[j, i] = tops.get((x0, y0), np.nan)
depth = 12.0 - Z

d = np.load(SFCW)
fig, axes = plt.subplots(1, 2, figsize=(17, 6.2), dpi=130)

# ---- left: y-z section at x=6 ----
ax = axes[0]
ax.set_facecolor('white')
for x0, y0, z0, x1, y1, z1, mat in slice_boxes:
    color = '#b0b0b0' if mat == 'rock' else '#c8863c'
    ax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0, facecolor=color,
                           edgecolor='none', zorder=1))
ax.plot([0, 12], [12, 12], 'k-', lw=1.4, zorder=3)
ax.plot([0, 0, 12, 12], [0, 28, 28, 0], ls=(0, (5, 3)), color='purple',
        lw=1.2, zorder=2)
ax.text(0.15, 4.4, 'PML', color='purple', fontsize=8, va='top')
ax.text(3.0, 5.2, 'rock  er=9', ha='center', va='center', fontsize=11, color='k')
ax.text(3.0, 10.9, 'cover  er=18 (Debye)', ha='center', va='center',
        fontsize=9, color='k')
ax.text(11.8, 21.5, 'air  er=1', ha='right', va='center', fontsize=11, color='k')
ax.annotate('rough interface (2.7-3.4 m deep)\nRMS 0.16 m, CL 2.5 m',
            xy=(9.3, 8.88), xytext=(11.4, 7.2), fontsize=8.5, color='darkred',
            ha='right', va='center',
            arrowprops=dict(arrowstyle='->', color='darkred',
                            connectionstyle='arc3,rad=-0.25'))
ax.plot([src[1]], [src[2]], marker='*', ms=16, color='red', zorder=5)
ax.plot([rx0[1]], [rx0[2]], marker='v', ms=10, color='blue', zorder=5)
yy = np.linspace(3.85, 6.85, 13)
ax.plot(yy, np.full(13, 27.0), '-', color='gray', lw=0.8, zorder=4)
ax.plot(yy, np.full(13, 27.0), 'r*', ms=5, zorder=5)
ax.plot(yy + 1.3, np.full(13, 27.0), 'bv', ms=4, zorder=5)
ax.annotate('Tx (t01)', xy=(src[1], src[2]), xytext=(0.3, 24.2), fontsize=9,
            color='red', arrowprops=dict(arrowstyle='->', color='red'))
ax.annotate('Rx (t01)', xy=(rx0[1], rx0[2]), xytext=(9.0, 24.2), fontsize=9,
            color='blue', arrowprops=dict(arrowstyle='->', color='blue'))
ax.text(6.0, 25.6, 'CO profile at x=6.0 m: 13 traces, offset 1.3 m\n'
        'dy=0.25 m, height 15 m', fontsize=8, ha='center', va='top',
        color='dimgray')
ax.set_xlim(0, 12)
ax.set_ylim(28, 4)
ax.set_xlabel('y (m)')
ax.set_ylabel('z (m, model coords; surface z=12)')
ax.set_title(f'3D 5cm | y-z section at x={SLICE_X} m (Tx/Rx plane),'
             f' parsed from {PREFIX}-t01.in', fontsize=10)
ax.set_aspect('equal', adjustable='box')
ax.grid(alpha=0.15)
# inset: plan view of interface depth (placed in empty air column, lower-left)
iax = ax.inset_axes([0.05, 0.30, 0.40, 0.28])
iax.set_facecolor('white')
pc = iax.pcolormesh(np.array(xs) + 0.125, np.array(ys) + 0.125, depth,
                    cmap='terrain', vmin=2.5, vmax=3.6, shading='nearest')
iax.plot([SLICE_X, SLICE_X], [1, 11], 'w--', lw=1.2)
iax.text(SLICE_X + 0.15, 10.5, 'section\nplane', color='w', fontsize=7)
iax.set_xlim(1, 11); iax.set_ylim(1, 11)
iax.set_title('interface depth (m)', fontsize=8)
iax.text(0.03, 0.03, '2.6 (blue) - 3.4 (brown)', transform=iax.transAxes,
         fontsize=6.5, color='k',
         bbox=dict(facecolor='white', alpha=0.75, pad=1, edgecolor='none'))
iax.set_xticks([]); iax.set_yticks([])

# ---- right: SFCW device-band B-scan ----
ax = axes[1]
mats = []
tm = None
for k in range(1, 14):
    t = f't{k:02d}'
    mats.append(d[f'{TAG}_{t}_no_taper_hann_complex_envelope'])
    tm = d[f'{TAG}_{t}_no_taper_hann_envelope_time_s']
tm = tm * 1e9
sel = tm <= 400
B = np.stack([m[sel] for m in mats])
direct = B[:, (tm[sel] >= 0) & (tm[sel] < 50)]
ref = np.max(np.abs(direct), axis=1, keepdims=True)
ref = np.where(ref > 0, ref, 1.0)
im = 20 * np.log10(np.abs(B) / ref + 1e-12)
pc2 = ax.pcolormesh(np.arange(1, 14), tm[sel], im.T, vmin=-60, vmax=0,
                    cmap='viridis', shading='auto')
for tt, lab, col in [(t_surf, 'surf', 'w'), (t_iface, 'iface', 'r')]:
    ax.axhline(tt, color=col, lw=1, ls='--')
    ax.text(13.1, tt, lab, color=col, fontsize=8, va='center')
ax.set_xlabel('trace #')
ax.set_ylim(400, 0)
ax.set_ylabel('two-way time (ns)')
ax.set_title('3D 5cm | simulated B-scan (SFCW 20-170 MHz, Hann, dB re direct)',
             fontsize=10)
ax.set_xticks(range(1, 14))
fig.colorbar(pc2, ax=ax, label='dB re direct')

fig.suptitle('benchmark3d_r2_co: 3D model section + interface plan view vs'
             ' device-band B-scan')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
print('boxes:', len(boxes), 'slice boxes:', len(slice_boxes),
      'depth range:', np.nanmin(depth), np.nanmax(depth))
