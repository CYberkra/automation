# -*- coding: utf-8 -*-
"""Estimate Line9 flight height and antenna wavelet from field data.

Height: surface echo delay rel t0 -> h = dt*c/2 (air, terrain-following).
Wavelet: median direct arrival (t0-5 .. t0+60 ns) as antenna IR estimate.
"""
from pathlib import Path
import sys

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from analyze_field_ringing import load_csv  # noqa: E402

DATA = Path(r'E:\automation_djh\real_data_yingshan\营山测线数据')
x, t, dt, hdr = load_csv(DATA / 'Line9origin(36).csv')
c = 0.299792458  # m/ns
print('header:', hdr)

i0 = np.argmax(np.abs(x[:, t <= 60]), axis=1)
pk = np.abs(x[np.arange(len(x)), i0])
valid = pk > 0
x, i0, pk = x[valid], i0[valid], pk[valid]

sos = butter(4, [20e6, 170e6], btype='band', fs=1 / (dt * 1e-9), output='sos')
env = np.abs(hilbert(sosfiltfilt(sos, x, axis=1), axis=1))
nrm_env = env / pk[:, None]

# surface echo: first peak with envelope > 30% of direct, at t > t0+15 ns
tsurf = np.full(len(x), np.nan)
for j in range(len(x)):
    seg_t = t - t[i0[j]]
    m = (seg_t > 15) & (seg_t < 160)
    idx = np.where(m & (nrm_env[j] > 0.30))[0]
    if len(idx):
        tsurf[j] = seg_t[idx[0]]

med = np.nanmedian(tsurf)
print('surface echo delay rel t0: median=%.1f ns  p25=%.1f  p75=%.1f  '
      'frac found=%.2f' % (med, np.nanpercentile(tsurf, 25),
                           np.nanpercentile(tsurf, 75),
                           np.mean(~np.isnan(tsurf))))
print('=> flight height: median=%.1f m  IQR=[%.1f, %.1f] m' %
      (med * c / 2, np.nanpercentile(tsurf, 25) * c / 2,
       np.nanpercentile(tsurf, 75) * c / 2))

# antenna wavelet estimate: median of aligned traces, t0-5..t0+70 ns
nrm = x / pk[:, None]
wlen = int(round(75 / dt))
W = np.stack([nrm[j][max(i0[j] - 4, 0):i0[j] + wlen] for j in range(len(x))])
W = W[:, :wlen]
wav = np.median(W, axis=0)
tw = (np.arange(wlen) - 4) * dt
e = np.abs(hilbert(wav))
print('wavelet peak at t=%.1f ns; envelope >10%% of peak until t=%.1f ns; '
      '>5%% until t=%.1f ns' % (tw[np.argmax(np.abs(wav))],
                               tw[np.where(e > 0.1 * e.max())[0][-1]],
                               tw[np.where(e > 0.05 * e.max())[0][-1]]))
np.savez(r'E:\automation_djh\artifacts_check\ringing_field\line9_wavelet.npz',
         wavelet=wav, t_ns=tw)
print('saved line9_wavelet.npz')
