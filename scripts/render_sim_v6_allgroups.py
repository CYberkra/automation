"""Inject house V6 ringing into all three benchmark3d_r2_co SFCW B-scans.

Same operator/parameters as the field-calibrated analysis: op1_ringing_v06,
q targets from Line9 p25..p75, seed 20261001. Surface-edit op NOT applied.
Display: bipolar grayscale, capped t-gain g=min(t/10ns, x30), per-row shared
99.5-pct normalization so before/after are directly comparable.
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
OUT = Path(r'E:\automation_djh\fig_sim_v6_ringing_all.png')

GROUPS = [('3D_5cm', '3D 5cm'), ('2D_5cm', '2D 5cm'),
          ('2D_2p5cm', '2D 2.5cm')]
t_surf, t_iface = 100.2, 185.1
win = (160.0, 210.0)
T_REF_NS, GMAX = 10.0, 30.0

real = load_line9()
calib = calibrate_q(real, np.linspace(0, T_WINDOW_NS, N_SAMPLES))
print(f'q target band: p25={calib["p25"]:.3f} p75={calib["p75"]:.3f}')

d = np.load(SFCW)
fig, axes = plt.subplots(3, 2, figsize=(16, 13.5), dpi=130, sharey=True)
gains = {}

for row, (tag, gname) in enumerate(GROUPS):
    mats, tm = [], None
    for k in range(1, 14):
        t = f't{k:02d}'
        mats.append(d[f'{tag}_{t}_tail_200ns_hann_real_bandpass'])
        tm = d[f'{tag}_{t}_tail_200ns_hann_envelope_time_s']
    B0 = np.stack(mats)
    tfull = tm * 1e9
    rng = np.random.default_rng(SEED)
    B1, n_clip = op1_ringing_v06(B0, tfull, rng, calib['p25'], calib['p75'])

    sel = tfull <= 400
    B0, B1, t = B0[:, sel], B1[:, sel], tfull[sel]
    g = np.minimum(t / T_REF_NS, GMAX)
    clip = np.nanpercentile(np.abs(B1 * g[None, :]), 99.5)  # shared per row
    p0 = np.clip(B0 * g[None, :] / clip, -1, 1)
    p1 = np.clip(B1 * g[None, :] / clip, -1, 1)

    pk = np.max(np.abs(B0), axis=1)
    m_r = (t >= 45) & (t <= 70)
    rr = np.sqrt(((B1 - B0)[:, m_r] ** 2).mean(axis=1)) / pk * 100
    print(f'{tag}: beta clipped {n_clip}/13, ringing rms % per trace:',
          np.round(rr, 1).tolist())

    for col, (M, ttl) in enumerate([
            (p0, f'{gname} | sim raw (no ringing)'),
            (p1, f'{gname} | + V6 ringing (107 MHz, tau=73 ns)')]):
        ax = axes[row, col]
        pc = ax.pcolormesh(np.arange(1, 14), t, M.T, vmin=-1, vmax=1,
                           cmap='gray', shading='auto')
        for tt, lab in [(t_surf, 'surf'), (t_iface, 'iface')]:
            ax.axhline(tt, color='red', lw=0.8, ls='--')
            ax.text(13.15, tt, lab, color='red', fontsize=8, va='center')
        ax.axhspan(win[0], win[1], color='blue', alpha=0.05, lw=0)
        ax.set_title(ttl, fontsize=10)
        ax.set_xticks(range(1, 14))
        if row == 2:
            ax.set_xlabel('trace #')
        if col == 0:
            ax.set_ylabel('two-way time (ns)')
            ax.set_ylim(400, 0)
axes[0, 0].set_ylim(400, 0)
cax = fig.add_axes([0.905, 0.12, 0.010, 0.76])
fig.colorbar(pc, cax=cax, label='normalized (white=+, black=-)')
fig.suptitle('benchmark3d_r2_co: V6 field ringing injected into simulated'
             f' device-band B-scans (q in [{calib["p25"]:.3f},'
             f' {calib["p75"]:.3f}], seed {SEED}; bipolar, capped t-gain)')
fig.tight_layout(rect=[0, 0, 0.90, 0.97])
fig.savefig(OUT)
print('saved:', OUT)
