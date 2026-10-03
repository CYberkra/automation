# -*- coding: utf-8 -*-
"""Inject FIELD-REVISED ringing into all three benchmark3d_r2_co SFCW B-scans.

Revision vs the stock V6 operator (augment_sim2real_v0_6):
  - tau  73 -> 58 ns   (field envelope-decay slopes: 49.7/50.4/50.8/59.4/62.1/62.4)
  - amplitude white jitter sigma 0.16 -> 0.20 (field lateral diff_rms 0.21-0.33)
  - q target band: pooled 6 lines (10394 traces) p25..p75 = 0.0456..0.0899
    (stock: Line9 only, 0.0530..0.0910)
Unchanged: F0=107 MHz, onset gate t0+12~22 ns, delta ±3.5 ns, q-match window
t0+45~70 ns, 90-125 MHz analysis band, seed 20261001. Surface-edit op NOT applied.
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
from augment_sim2real_v0_6 import (  # noqa: E402
    SEED, _sos, smooth_series, t0_index, win_mask,
    BP_LO, BP_HI, QW_LO, QW_HI, ONSET0, ONSET1, DELTA_NS,
)

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
QDIR = Path(r'E:\automation_djh\artifacts_check\ringing_field')
OUT = Path(r'E:\automation_djh\fig_sim_v6r_ringing_all.png')

TAU_REV = 58.0       # ns, field-revised ring-down
JITTER_REV = 0.20    # amplitude white jitter sigma, field-revised

GROUPS = [('3D_5cm', '3D 5cm'), ('2D_5cm', '2D 5cm'),
          ('2D_2p5cm', '2D 2.5cm')]
t_surf, t_iface = 100.2, 185.1
win = (160.0, 210.0)
T_REF_NS, GMAX = 10.0, 30.0


def onset_gate(t_ns, t0):
    g = np.ones_like(t_ns)
    trel = t_ns - t0
    g[trel < ONSET0] = 0.0
    ramp = (trel >= ONSET0) & (trel < ONSET1)
    u = (t_ns[ramp] - t0 - ONSET0) / (ONSET1 - ONSET0)
    g[ramp] = 0.5 - 0.5 * np.cos(np.pi * u)
    return g


def op1_ringing_v07(sigs, t_ns, rng, q_lo, q_hi, tau, jitter_std):
    """v06 body with parameterized tau / jitter_std (field-revised)."""
    from scipy.signal import hilbert, sosfiltfilt
    n, nt = sigs.shape
    dt = float(t_ns[1] - t_ns[0])
    tg = np.arange(0, 5 * tau, dt)
    g = np.exp(-tg / tau) * np.sin(2 * np.pi * 0.107 * tg)
    sos = _sos(dt)
    qs = smooth_series(rng, n, q_lo, q_hi, 5.0)
    deltas = smooth_series(rng, n, -DELTA_NS, DELTA_NS, 3.0)
    jitter = np.clip(1.0 + jitter_std * rng.standard_normal(n), 0.4, None)
    beta_grid = np.concatenate([[0.0], np.logspace(-3, np.log10(3.0), 96)])
    out = sigs.copy()
    n_clip = 0
    for i in range(n):
        i0 = t0_index(sigs[i], t_ns)
        t0 = t_ns[i0]
        peak = abs(sigs[i][i0]) + 1e-30
        c = np.convolve(sigs[i], g)[:nt] * dt
        c = np.interp(t_ns - deltas[i], t_ns, c, left=0.0, right=0.0)
        c = c * onset_gate(t_ns, t0)
        Ya = hilbert(sosfiltfilt(sos, sigs[i]))
        Ca = hilbert(sosfiltfilt(sos, c))
        m = win_mask(t_ns, t0, QW_LO, QW_HI)
        env_med = np.median(np.abs(Ya[m][None, :] + beta_grid[:, None] * Ca[m][None, :]), axis=1)
        q_curve = env_med / peak
        target = qs[i]
        hit = np.where(q_curve >= target)[0]
        if len(hit) == 0:
            beta = beta_grid[-1]
            n_clip += 1
        else:
            j = hit[0]
            if j == 0:
                beta = 0.0
            else:
                b0, b1 = beta_grid[j - 1], beta_grid[j]
                q0, q1 = q_curve[j - 1], q_curve[j]
                beta = b0 + (b1 - b0) * (target - q0) / max(q1 - q0, 1e-30)
        out[i] = sigs[i] + beta * jitter[i] * c
    return out, n_clip


# --- pooled q target band from the 6 field lines -------------------------
q_all = np.concatenate([np.load(p) for p in sorted(QDIR.glob('*_q.npy'))])
q_lo, q_hi = (float(np.percentile(q_all, 25)), float(np.percentile(q_all, 75)))
print(f'pooled q target band: n={len(q_all)} p25={q_lo:.4f} p75={q_hi:.4f}')

d = np.load(SFCW)
fig, axes = plt.subplots(3, 2, figsize=(16, 13.5), dpi=130, sharey=True)

for row, (tag, gname) in enumerate(GROUPS):
    mats, tm = [], None
    for k in range(1, 14):
        t = f't{k:02d}'
        mats.append(d[f'{tag}_{t}_tail_200ns_hann_real_bandpass'])
        tm = d[f'{tag}_{t}_tail_200ns_hann_envelope_time_s']
    B0 = np.stack(mats)
    tfull = tm * 1e9
    rng = np.random.default_rng(SEED)
    B1, n_clip = op1_ringing_v07(B0, tfull, rng, q_lo, q_hi,
                                 TAU_REV, JITTER_REV)

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
            (p1, f'{gname} | + revised ringing (107 MHz, tau=58 ns, jit=0.20)')]):
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
fig.suptitle('benchmark3d_r2_co: FIELD-REVISED ringing injected into simulated'
             f' device-band B-scans (pooled 6-line q in [{q_lo:.4f},'
             f' {q_hi:.4f}], tau=58 ns, jitter=0.20, seed {SEED};'
             ' bipolar, capped t-gain)')
fig.tight_layout(rect=[0, 0, 0.90, 0.97])
fig.savefig(OUT)
print('saved:', OUT)
