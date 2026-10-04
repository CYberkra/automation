# -*- coding: utf-8 -*-
"""Schematic of the BF 2D scan model (h=7 m, measured-wavelet excitation,
sigma=2) + the three candidate excitation wavelets."""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from gen_bf2d_scan import (build_geometry, cfg_seed, ER_COVER, SIGMA_COVER,
                           Y0, YBIN, ZI_R2, Z_SURF, TX_Y0, RX_Y0, PML_CFS)

H, SIGMA, SEED = 7.0, 2.0, cfg_seed('bf2d-phase2-h7-meas-s2')
ZTOP = Z_SURF + H + 6.0
bins = build_geometry(H, SIGMA, SEED)

fig = plt.figure(figsize=(15, 7.2), dpi=130)

# ---- panel A: model cross-section ---------------------------------------
ax = fig.add_subplot(1, 2, 1)
# bedrock
ax.add_patch(Rectangle((0, 1), 12, Z_SURF - 1, fc='#d9c9a8', ec='none',
                       zorder=0))
ers = [b[4] for b in bins]
vmin, vmax = ER_COVER - 3 * SIGMA, ER_COVER + 3 * SIGMA
cmap = plt.get_cmap('RdYlBu_r')
for (y0, y1, z0, z1, er) in bins:
    ax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0,
                           fc=cmap((er - vmin) / (vmax - vmin)),
                           ec='none', zorder=1))
# interface outline
yi = [Y0 + i * YBIN for i in range(len(ZI_R2) + 1)]
ax.step(yi, ZI_R2 + [ZI_R2[-1]], where='post', color='k', lw=1.2, zorder=3)
# surface
ax.axhline(Z_SURF, color='k', lw=1.5, zorder=3)
# PML frame
for x0, w in ((0, 1), (11, 1)):
    ax.add_patch(Rectangle((x0, 1), w, ZTOP - 1, fill=False, hatch='///',
                           ec='gray', lw=0.4, zorder=4))
ax.add_patch(Rectangle((0, 1), 12, 1, fill=False, hatch='///',
                       ec='gray', lw=0.4, zorder=4))
ax.add_patch(Rectangle((0, ZTOP - 1), 12, 1, fill=False, hatch='///',
                       ec='gray', lw=0.4, zorder=4))
# Tx / Rx
ztx = Z_SURF + H
ax.plot([TX_Y0], [ztx], marker='*', ms=16, color='crimson', zorder=5)
ax.plot([RX_Y0], [ztx], marker='^', ms=11, color='navy', zorder=5)
ax.annotate('Tx', (TX_Y0, ztx), xytext=(TX_Y0 - 0.55, ztx + 0.9),
            fontsize=10, color='crimson')
ax.annotate('Rx', (RX_Y0, ztx), xytext=(RX_Y0 + 0.15, ztx + 0.9),
            fontsize=10, color='navy')
# rays: direct / surface / interface
ymid_s, ymid_i = 6.0, 6.0
ax.plot([TX_Y0, RX_Y0], [ztx, ztx], color='crimson', lw=1.0, ls='-')
ax.plot([TX_Y0, ymid_s, RX_Y0], [ztx, Z_SURF, ztx], color='seagreen',
        lw=1.0, ls='--')
zi_mid = float(np.interp(ymid_i, [b[0] for b in bins], [b[2] for b in bins]))
ax.plot([TX_Y0, ymid_i, RX_Y0], [ztx, zi_mid, ztx], color='darkorange',
        lw=1.0, ls='-.')
ax.annotate('surface\n~46 ns rel t0', (7.9, Z_SURF + 1.6), fontsize=8,
            color='seagreen')
ax.annotate('bedrock interface\n(~86 ns rel t0)', (7.4, zi_mid - 1.7),
            fontsize=8, color='darkorange')
# geometry labels
ax.annotate('', (11.7, Z_SURF), (11.7, ztx),
            arrowprops=dict(arrowstyle='<->', color='k', lw=0.9))
ax.annotate(f'h = {H:g} m', (10.75, Z_SURF + H / 2), fontsize=9,
            ha='right')
ax.annotate('', (TX_Y0, ztx - 1.4), (RX_Y0, ztx - 1.4),
            arrowprops=dict(arrowstyle='<->', color='k', lw=0.9))
ax.annotate('offset 1.3 m', (5.5, ztx - 2.3), fontsize=8, ha='center')
ax.text(2.7, Z_SURF + 2.4, 'air', ha='center', fontsize=10)
ax.annotate('direct', (5.05, ztx - 0.9), fontsize=8, color='crimson')
ax.text(3.1, 10.6, f'cover ~3 m\n$\\epsilon_r$ = 18.0$\\pm${SIGMA:g}'
        f'\n$\\sigma$ = {SIGMA_COVER} S/m', fontsize=9, ha='center')
ax.text(6.0, 4.5, 'bedrock  $\\epsilon_r$ = 9', fontsize=10, ha='center')
sm = plt.cm.ScalarMappable(cmap=cmap,
                           norm=plt.Normalize(vmin, vmax))
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, fraction=0.038, pad=0.02)
cb.set_label('cover $\\epsilon_r$')
ax.set_xlim(0, 12)
ax.set_ylim(ZTOP, 0)
ax.set_xlabel('y (m)')
ax.set_ylabel('z (m, depth downward from surface z=12)')
ax.set_title(f'BF 2D model | h={H:g} m | measured-wavelet excitation | '
             f'seed {SEED}', fontsize=10)
ax.set_aspect('equal')

# ---- panel B: excitation candidates --------------------------------------
ax2 = fig.add_subplot(1, 2, 2)
t = np.arange(0, 120, 0.02)
# impulse (as it would appear after the 20-170 MHz SFCW band: ~7 ns sinc)
from scipy.signal import butter, sosfiltfilt
sos = butter(4, [20e6, 170e6], btype='band', fs=1 / 0.02e-9, output='sos')
imp = np.zeros_like(t)
imp[np.argmin(np.abs(t - 10))] = 1.0
imp_bp = sosfiltfilt(sos, imp)
ax2.plot(t, imp_bp / np.abs(imp_bp).max() * 0.5 + 2.0, lw=1.0,
         color='gray', label='impulse (r2 baseline, band-limited)')
# ricker 80 MHz
f0 = 80e6
r = (1 - 2 * (np.pi * f0 * (t - 20) * 1e-9) ** 2) * \
    np.exp(-(np.pi * f0 * (t - 20) * 1e-9) ** 2)
ax2.plot(t, r / np.abs(r).max() * 0.5 + 1.0, lw=1.0, color='steelblue',
         label='ricker 80 MHz')
# measured wavelet
d = np.load(r'E:\automation_djh\artifacts_check\ringing_field'
            r'\line9_wavelet.npz')
wav, tw = d['wavelet'], d['t_ns'] - d['t_ns'][0]
wf = interp1d(tw, wav / np.abs(wav).max(), kind='linear', fill_value=0.0,
              bounds_error=False)(t)
from scipy.signal import hilbert
env = np.abs(hilbert(wf))
ax2.plot(t, wf, lw=1.0, color='crimson', label='measured Line9 wavelet')
ax2.plot(t, env, lw=0.8, color='crimson', ls='--', alpha=0.6)
ax2.plot(t, -env, lw=0.8, color='crimson', ls='--', alpha=0.6)
ax2.axhline(2.0, color='gray', lw=0.3)
ax2.axhline(1.0, color='gray', lw=0.3)
ax2.axhline(0.0, color='gray', lw=0.3)
ax2.set_yticks([0.0, 1.0, 2.0],
               ['meas wavelet', 'ricker 80', 'impulse (r2)'])
ax2.set_xlabel('time (ns)')
ax2.set_xlim(0, 92)
ax2.set_title('excitation candidates (v3 uses measured wavelet)',
              fontsize=10)
ax2.legend(loc='upper right', fontsize=8)
ax2.annotate('envelope tail $\\geq$ 70 ns:\nthis is the missing 25-100 ns '
             'energy', (72, 0.45), fontsize=8, color='crimson')

fig.suptitle('benchmark v3 (BF) 2D scan model schematic | '
             'r2 skeleton + h=7 m + clutter sigma=2', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = r'E:\automation_djh\fig_bf2d_model_schematic.png'
fig.savefig(out)
print('saved:', out)
