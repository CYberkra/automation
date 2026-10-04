# -*- coding: utf-8 -*-
"""Debug: where does the huge early-time energy in S1 come from?"""
from pathlib import Path
import sys

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from v6r_op import op1_ringing_v07, pooled_q_band, TAU_REV, JITTER_REV  # noqa: E402
from augment_sim2real_v0_6 import SEED  # noqa: E402

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
d = np.load(SFCW)
mats, tm = [], None
for k in range(1, 14):
    t_ = f't{k:02d}'
    mats.append(d[f'3D_5cm_{t_}_tail_200ns_rectangular_real_bandpass'])
    tm = d[f'3D_5cm_{t_}_tail_200ns_rectangular_envelope_time_s']
B0 = np.stack(mats)
tfull = tm * 1e9
print('sim dt_ns =', tfull[1] - tfull[0], 'n =', len(tfull))

q_lo, q_hi, _ = pooled_q_band()
B1, n_clip = op1_ringing_v07(B0, tfull, np.random.default_rng(SEED),
                             q_lo, q_hi, TAU_REV, JITTER_REV)

# reproduce internal beta by re-running one pass with instrumentation
from v6r_op import onset_gate  # noqa: E402
from augment_sim2real_v0_6 import _sos, smooth_series, t0_index, win_mask, QW_LO, QW_HI, DELTA_NS  # noqa: E402

rng = np.random.default_rng(SEED)
n, nt = B0.shape
dt = float(tfull[1] - tfull[0])
tg = np.arange(0, 5 * TAU_REV, dt)
g = np.exp(-tg / TAU_REV) * np.sin(2 * np.pi * 0.107 * tg)
sos = _sos(dt)
qs = smooth_series(rng, n, q_lo, q_hi, 5.0)
deltas = smooth_series(rng, n, -DELTA_NS, DELTA_NS, 3.0)
jitter = np.clip(1.0 + JITTER_REV * rng.standard_normal(n), 0.4, None)

betas = []
for i in range(n):
    i0 = t0_index(B0[i], tfull)
    peak = abs(B0[i][i0]) + 1e-30
    c = np.convolve(B0[i], g)[:nt] * dt
    c = np.interp(tfull - deltas[i], tfull, c, left=0.0, right=0.0)
    c = c * onset_gate(tfull, tfull[i0])
    Ya = hilbert(sosfiltfilt(sos, B0[i]))
    Ca = hilbert(sosfiltfilt(sos, c))
    m = win_mask(tfull, tfull[i0], QW_LO, QW_HI)
    bg = np.concatenate([[0.0], np.logspace(-3, np.log10(3.0), 96)])
    env_med = np.median(np.abs(Ya[m][None, :] + bg[:, None] * Ca[m][None, :]), axis=1)
    q_curve = env_med / peak
    hit = np.where(q_curve >= qs[i])[0]
    j = hit[0]
    if j == 0:
        beta = 0.0
    else:
        b0, b1 = bg[j - 1], bg[j]
        q0, q1 = q_curve[j - 1], q_curve[j]
        beta = b0 + (b1 - b0) * (qs[i] - q0) / max(q1 - q0, 1e-30)
    betas.append(beta)

betas = np.array(betas)
print('beta per trace:', np.round(betas, 3).tolist())
print('jitter:', np.round(jitter, 2).tolist())
print('conv template peak |c| rel direct peak:',
      np.round([np.abs(np.convolve(B0[i], g)[:nt] * dt).max() / np.abs(B0[i]).max()
                for i in range(3)], 3))

R = B1 - B0
print('max|R| / direct pk per trace:',
      np.round(np.abs(R).max(axis=1) / np.abs(B0).max(axis=1), 3).tolist())
m25 = (tfull >= 25) & (tfull <= 45)
print('max|R| in 25-45abs / direct pk:',
      np.round(np.abs(R[:, m25]).max(axis=1) / np.abs(B0).max(axis=1), 3).tolist())

S1 = B1 - B1.mean(axis=0, keepdims=True)
print('max|S1| in 25-45abs / direct pk:',
      np.round(np.abs(S1[:, m25]).max(axis=1) / np.abs(B0).max(axis=1), 3).tolist())
