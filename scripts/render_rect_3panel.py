# -*- coding: utf-8 -*-
"""Rectangular-window SFCW B-scans: raw / mean-removed / +revised ringing.

Chain per group (3D_5cm, 2D_5cm, 2D_2p5cm), 13 traces:
  col1: rectangular-window device-band B-scan, raw
  col2: col1 minus mean trace (background subtraction, no ringing)
  col3: revised ringing (tau=58 ns, jitter=0.20, pooled 6-line q band)
        injected into raw, then mean-removed -> full processing chain
Rectangular window preserves the 185 ns interface (-36 dB re direct) that the
Hann tail taper was suppressing by ~44 dB.
Display: bipolar gray, capped t-gain g=min(t/10ns, x30), per-row shared 99.5-pct
normalization across all three panels.
"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from v6r_op import op1_ringing_v07, TAU_REV, JITTER_REV, pooled_q_band  # noqa: E402
from augment_sim2real_v0_6 import SEED  # noqa: E402

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_rect_3panel_all.png')

GROUPS = [('3D_5cm', '3D 5cm'), ('2D_5cm', '2D 5cm'),
          ('2D_2p5cm', '2D 2.5cm')]
t_surf, t_iface = 100.2, 185.1
win = (160.0, 210.0)
T_REF_NS, GMAX = 10.0, 30.0

q_lo, q_hi, nq = pooled_q_band()
print(f'pooled q target band: n={nq} p25={q_lo:.4f} p75={q_hi:.4f}')

d = np.load(SFCW)
fig, axes = plt.subplots(3, 3, figsize=(19, 13.5), dpi=130, sharey=True)

for row, (tag, gname) in enumerate(GROUPS):
    mats, tm = [], None
    for k in range(1, 14):
        t = f't{k:02d}'
        mats.append(d[f'{tag}_{t}_tail_200ns_rectangular_real_bandpass'])
        tm = d[f'{tag}_{t}_tail_200ns_rectangular_envelope_time_s']
    B0 = np.stack(mats)
    tfull = tm * 1e9

    rng = np.random.default_rng(SEED)
    Br, n_clip = op1_ringing_v07(B0, tfull, rng, q_lo, q_hi,
                                 TAU_REV, JITTER_REV)
    mean0 = B0.mean(axis=0, keepdims=True)
    meanr = Br.mean(axis=0, keepdims=True)
    panels = [B0, B0 - mean0, Br - meanr]

    sel = tfull <= 400
    panels = [M[:, sel] for M in panels]
    t = tfull[sel]
    g = np.minimum(t / T_REF_NS, GMAX)
    clip = np.nanpercentile(np.abs(np.concatenate(panels) * g[None, :]), 99.5)
    pn = [np.clip(M * g[None, :] / clip, -1, 1) for M in panels]

    pk = np.max(np.abs(B0), axis=1)
    m_r = (tfull >= 45) & (tfull <= 70)
    rr = np.sqrt(((Br - B0)[:, m_r] ** 2).mean(axis=1)) / pk * 100
    print(f'{tag}: beta clipped {n_clip}/13, ringing rms % per trace:',
          np.round(rr, 1).tolist())

    titles = [f'{gname} | rectangular raw',
              f'{gname} | minus mean trace',
              f'{gname} | + revised ringing, minus mean']
    for col, (M, ttl) in enumerate(zip(pn, titles)):
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
cax = fig.add_axes([0.925, 0.12, 0.010, 0.76])
fig.colorbar(pc, cax=cax, label='normalized (white=+, black=-)')
fig.suptitle('benchmark3d_r2_co: rectangular-window device-band B-scans,'
             f' background subtraction and field-revised ringing'
             f' (tau={TAU_REV:.0f} ns, jitter={JITTER_REV:.2f},'
             f' pooled q [{q_lo:.4f}, {q_hi:.4f}], seed {SEED};'
             ' bipolar, capped t-gain)')
fig.tight_layout(rect=[0, 0, 0.92, 0.97])
fig.savefig(OUT)
print('saved:', OUT)
