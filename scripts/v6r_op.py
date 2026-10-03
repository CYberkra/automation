# -*- coding: utf-8 -*-
"""Shared field-revised ringing operator (single source of truth).

Revision vs stock V6 (automation_repo/scripts/augment_sim2real_v0_6.py):
  - tau  73 -> 58 ns   (field envelope-decay slopes: 49.7/50.4/50.8/59.4/62.1/62.4)
  - amplitude white jitter sigma 0.16 -> 0.20 (field lateral diff_rms 0.21-0.33)
  - q target band: pooled 6 field lines (10394 traces) p25..p75 = 0.0456..0.0899
    (stock: Line9 only, 0.0530..0.0910)
Unchanged: F0=107 MHz, onset gate t0+12~22 ns, delta ±3.5 ns, q-match window
t0+45~70 ns, 90-125 MHz analysis band, seed 20261001.
"""
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from augment_sim2real_v0_6 import (  # noqa: E402
    SEED, _sos, smooth_series, t0_index, win_mask,
    BP_LO, BP_HI, QW_LO, QW_HI, ONSET0, ONSET1, DELTA_NS,
)

QDIR = Path(r'E:\automation_djh\artifacts_check\ringing_field')

TAU_REV = 58.0       # ns, field-revised ring-down
JITTER_REV = 0.20    # amplitude white jitter sigma, field-revised


def pooled_q_band():
    """q target band p25..p75 pooled across the 6 field lines."""
    q_all = np.concatenate([np.load(p) for p in sorted(QDIR.glob('*_q.npy'))])
    return (float(np.percentile(q_all, 25)), float(np.percentile(q_all, 75)),
            len(q_all))


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
