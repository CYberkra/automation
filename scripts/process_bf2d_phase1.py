# -*- coding: utf-8 -*-
"""Process BF phase1 (h=5/7/10/15, t07) through the official SFCW chain.

Steps per H5: official API (20-170 MHz x 501, rectangular window, zero-pad 4,
no tail taper -- 1200 ns window per audit) -> real_bandpass trace; optional
device-domain convolution with the measured Line9 antenna wavelet (manifest
wavelet flag). Checks the official -60 dB tail guideline on the raw record,
measures surface-echo arrival (wavelet-off panels), and renders a
wavelet-off/on x height comparison plus 6-window profiles vs real Line9.
"""
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import hilbert

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from gprMax.toolboxes.SFCW.processing import (  # noqa: E402
    direct_frequency_response, load_receiver, load_source,
    reconstruct_time_response)
from analyze_field_ringing import load_csv, WINS  # noqa: E402

PDIR = Path(r'E:\automation_djh\verify_runs\bf2d_scan_20261002\phase1')
WAVZ = Path(r'E:\automation_djh\artifacts_check\ringing_field'
            r'\line9_wavelet.npz')
REAL = Path(r'E:\automation_djh\real_data_yingshan\营山测线数据'
            r'\Line9origin(36).csv')
OUTFIG = Path(r'E:\automation_djh\fig_bf2d_phase1.png')
FREQS = np.linspace(20e6, 170e6, 501)
HEIGHTS = (5, 7, 10, 15)
T_REF_NS, GMAX, T_SHOW = 10.0, 30.0, 300.0

# ---- official-chain per H5 ------------------------------------------------
proc = {}
for h in HEIGHTS:
    p = PDIR / f'BF2D-h{h}-s0-t07.in'
    h5 = p.with_suffix('.h5')
    src = load_source(h5)
    with h5py.File(h5, 'r') as f:
        rx_name = str(f['rxs']['rx1'].attrs['Name'])
    rx = load_receiver(h5, receiver_path=f'name:{rx_name}', component='Ex')
    raw = rx.samples
    pk_raw = np.abs(raw).max()
    tail_db = 20 * np.log10(np.abs(raw[int(0.95 * len(raw)):]).max()
                            / pk_raw + 1e-30)
    print(f'h{h}: record {len(raw)} samples ({len(raw)*rx.dt*1e9:.0f} ns), '
          f'last-5% peak {tail_db:.1f} dB (official guideline -60)')
    fr = direct_frequency_response(src, rx, FREQS, tail_taper_fraction=0.0)
    assert bool(np.all(fr.source_valid))
    tr = reconstruct_time_response(fr, window='rectangular',
                                   zero_pad_factor=4, time_shift=0.0)
    proc[h] = (tr.real_bandpass, tr.time * 1e9)

# ---- measured wavelet (device-domain injection) ---------------------------
d = np.load(WAVZ)
wav, tw = d['wavelet'], d['t_ns'] - d['t_ns'][0]
wav = wav / np.abs(wav).max()

# ---- real Line9 reference -------------------------------------------------
x, t, dt, hdr = load_csv(REAL)
i0 = np.argmax(np.abs(x[:, t <= 60]), axis=1)
pk = np.abs(x[np.arange(len(x)), i0])
m = pk > 0
x, i0, pk = x[m], i0[m], pk[m]
nrm = x / pk[:, None]
FULL = np.stack([np.roll(nrm[j], -i0[j]) for j in range(len(nrm))])
t_full = np.arange(FULL.shape[1]) * dt
sel_r = t_full <= T_SHOW


def band_env(M, tt):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [20e6, 170e6], btype='band', fs=1 / (tt[1] - tt[0]) * 1e9,
                 output='sos')
    return np.abs(hilbert(sosfiltfilt(sos, M, axis=1), axis=1))


def prof6(M, tt):
    env = band_env(M, tt)
    m_d = tt <= 12.0
    pk_t = env[:, m_d].max(axis=1, keepdims=True) + 1e-30
    env = env / pk_t
    out = []
    for lo, hi in WINS:
        w = (tt >= lo) & (tt <= hi)
        out.append(float(np.median(env[:, w].mean(axis=1))))
    return out


def surf_delay(M, tt):
    """Envelope argmax in (t0+15, t0+170) ns (wavelet-off panels)."""
    env = band_env(M, tt)
    w = (tt > 15) & (tt < 170)
    idx = np.where(w)[0]
    return float(tt[idx[np.argmax(env[0][idx])]])


fig = plt.figure(figsize=(18.5, 10.5), dpi=130)
gs = fig.add_gridspec(3, 4, height_ratios=[2.4, 2.4, 1.6])
centers = [(lo + hi) / 2 for lo, hi in WINS]
colors = {5: '#1f77b4', 7: '#d62728', 10: '#2ca02c', 15: '#9467bd'}

for row, use_wav in enumerate((False, True)):
    for col, h in enumerate(HEIGHTS):
        B0, tt = proc[h]
        B0 = B0[None, :] if B0.ndim == 1 else B0
        if use_wav:
            wf = interp1d(tw, wav, kind='linear', fill_value=0.0,
                          bounds_error=False)(tt)
            B1 = np.array([np.convolve(B0[i], wf)[:B0.shape[1]]
                           for i in range(B0.shape[0])])
        else:
            B1 = B0
        i0s = np.argmax(np.abs(B1), axis=1)
        B1 = np.stack([np.roll(B1[j], -i0s[j]) for j in range(B1.shape[0])])
        t_rel = tt - tt[0]
        s = t_rel <= T_SHOW
        B1, t_rel = B1[:, s], t_rel[s]
        g = np.minimum(t_rel / T_REF_NS, GMAX)
        m_d = t_rel <= 12.0
        pk_t = np.abs(B1[:, m_d]).max(axis=1, keepdims=True) + 1e-30
        P = np.clip(B1 * g[None, :] / pk_t, -1, 1)
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(P, aspect='auto', vmin=-1, vmax=1, cmap='gray',
                  extent=[0.5, 1.5, t_rel[-1], t_rel[0]])
        if row == 0:
            sd = surf_delay(np.stack([np.roll(proc[h][0], -np.argmax(
                np.abs(proc[h][0])))]), (proc[h][1] - proc[h][1][0]))
            ax.set_title(f'h={h} m | wavelet OFF\nsurface echo ~{sd:.0f} ns '
                         f'rel t0', fontsize=9)
        else:
            ax.set_title(f'h={h} m | wavelet ON', fontsize=9)
        ax.set_ylim(T_SHOW, 0)
        if col == 0:
            ax.set_ylabel('t rel t0 (ns)')
        if row == 1:
            ax.set_xlabel('trace')

# ---- profiles --------------------------------------------------------------
axp = fig.add_subplot(gs[2, :])
env_real = band_env(FULL, t_full)
pr = []
for lo, hi in WINS:
    w = (t_full >= lo) & (t_full <= hi)
    pr.append(float(np.median(env_real[:, w].mean(axis=1))))
axp.plot(centers, pr, 'o-', color='k', lw=1.5,
         label='Line9 real (full line)')
for h in HEIGHTS:
    B0, tt = proc[h]
    wf = interp1d(tw, wav, kind='linear', fill_value=0.0,
                  bounds_error=False)(tt)
    B1 = np.convolve(B0, wf)[:len(B0)][None, :]
    t_rel = tt - tt[0]
    s = t_rel <= T_SHOW
    ps = prof6(B1[:, s], t_rel[s])
    axp.plot(centers, ps, 's--', ms=4, color=colors[h],
             label=f'sim h={h} m + wavelet')
axp.set_yscale('log')
axp.set_xlabel('window center (ns rel t0)')
axp.set_ylabel('median 20-170 MHz envelope\n(/ direct peak)')
axp.grid(True, which='both', alpha=0.3)
axp.legend(fontsize=8, ncol=3)
axp.set_title('ring-down profile: real Line9 vs BF phase1 (wavelet on)',
              fontsize=10)

fig.suptitle('BF phase1 | official SFCW chain (rect window, zp4) | '
             'measured-antenna-wavelet injection post-chain', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig(OUTFIG)
print('saved:', OUTFIG)
