"""Add the house V6 ringing operator (augment_sim2real_v0_6.op1_ringing_v06)
to the 3D_5cm SFCW B-scan and display before/after/component.

Ringing params (house v0.6): F0=107 MHz, TAU=73 ns, onset gate t0+12~22 ns,
beta solved per trace to hit Line9-calibrated q (90-125 MHz envelope median
in t0+45~70 ns / direct peak, p25..p75), phase jitter +/-3.5 ns, amplitude
jitter sigma=0.16, seed 20261001. Surface-edit op is NOT applied.
"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from augment_sim2real_v0_6 import calibrate_q, op1_ringing_v06, SEED  # noqa: E402
from augment_sim2real_v0_3 import load_line9, T_WINDOW_NS, N_SAMPLES  # noqa: E402

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_3d_v6_ringing.png')

TAG = '3D_5cm'
t_surf, t_iface = 100.2, 185.1
T_REF_NS, GMAX = 10.0, 30.0

# ---- load SFCW real bandpass (signed, full window) ----
d = np.load(SFCW)
mats, tm = [], None
for k in range(1, 14):
    tt = f't{k:02d}'
    mats.append(d[f'{TAG}_{tt}_tail_200ns_hann_real_bandpass'])
    tm = d[f'{TAG}_{tt}_tail_200ns_hann_envelope_time_s']
B0 = np.stack(mats)
tfull = tm * 1e9

# ---- calibrate q from measured Line9, apply V6 ringing ----
real = load_line9()
calib = calibrate_q(real, np.linspace(0, T_WINDOW_NS, N_SAMPLES))
rng = np.random.default_rng(SEED)
B1, n_clip = op1_ringing_v06(B0, tfull, rng, calib['p25'], calib['p75'])
C = B1 - B0  # ringing component

# ---- crop to 400 ns ----
sel = tfull <= 400
B0, B1, C, t = B0[:, sel], B1[:, sel], C[:, sel], tfull[sel]

# ---- t-gain with cap (display only) ----
g = np.minimum(t / T_REF_NS, GMAX)
clip = np.nanpercentile(np.abs(B0 * g[None, :]), 99.5)
p0 = np.clip(B0 * g[None, :] / clip, -1, 1)
p1 = np.clip(B1 * g[None, :] / clip, -1, 1)  # same normalization as baseline
clipc = np.nanpercentile(np.abs(C * g[None, :]), 99.5)
pc = np.clip(C * g[None, :] / clipc, -1, 1)

# ---- diagnostics ----
pk0 = np.max(np.abs(B0), axis=1)
m_ring = (t >= 45) & (t <= 70)
ring_rms = np.sqrt((C[:, m_ring] ** 2).mean(axis=1))
i0 = np.argmax(np.abs(B0), axis=1)
cor = []
for i in range(13):
    m_dir = np.zeros_like(t, dtype=bool)
    m_dir[max(0, i0[i] - 3):i0[i] + 4] = True
    cor.append(np.corrcoef(B0[i][m_dir], B1[i][m_dir])[0, 1])
print(f'calib q: ' + ', '.join(f'{k}={v:.4f}' for k, v in calib.items()))
print(f'beta clipped traces: {n_clip}/13')
print(f'direct-window corr med/min: {np.median(cor):.4f}/{np.min(cor):.4f}')
print('ringing rms (t0+45~70) re direct peak per trace (%):',
      np.round(100 * ring_rms / pk0, 2).tolist())
print(f'ringing period check: F0=107 MHz -> 9.35 ns')

fig, axes = plt.subplots(1, 3, figsize=(19, 6.2), dpi=130, sharey=True)
panels = [
    (p0, f'baseline: t-gain g=min(t/{T_REF_NS:.0f}ns, x{GMAX:.0f}),'
         ' 99.5-pct norm\n(no ringing)'),
    (p1, f'+ V6 ringing (F0=107 MHz, tau=73 ns, onset t0+12 ns,\n'
         f'q in [{calib["p25"]:.3f}, {calib["p75"]:.3f}], seed {SEED})\n'
         'same normalization as baseline panel'),
    (pc, 'ringing component only (after - before)\nown 99.5-pct norm,'
         ' same t-gain'),
]
for ax, (M, ttl) in zip(axes, panels):
    pcm = ax.pcolormesh(np.arange(1, 14), t, M.T, vmin=-1, vmax=1,
                        cmap='gray', shading='auto')
    for tt, lab in [(t_surf, 'surf'), (t_iface, 'iface')]:
        ax.axhline(tt, color='red', lw=0.8, ls='--')
        ax.text(13.15, tt, lab, color='red', fontsize=8, va='center')
    ax.set_xlabel('trace #')
    ax.set_title(ttl, fontsize=9.5)
    ax.set_xticks(range(1, 14))
axes[0].set_ylabel('two-way time (ns)')
axes[0].set_ylim(400, 0)
cax = fig.add_axes([0.92, 0.18, 0.010, 0.6])
fig.colorbar(pcm, cax=cax, label='normalized (white=+, black=-)')
fig.suptitle(f'benchmark3d_r2_co {TAG}: V6 ringing injection on SFCW'
             ' device-band B-scan (bipolar, capped t-gain)')
fig.tight_layout(rect=[0, 0, 0.91, 1])
fig.savefig(OUT)
print('saved:', OUT)
