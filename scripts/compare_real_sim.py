# -*- coding: utf-8 -*-
"""Side-by-side: real Line9 field B-scan vs simulated 3D final chain.

Same display contract for both: per-trace normalization by direct-arrival
peak, trace alignment to t0 (roll), bipolar grayscale, capped t-gain
g=min(rel_t/10ns, x30), 99.5-pct normalization per panel.
Panels:
  1. Line9 traces 1000-1200 (raw device output, 36 dBm)
  2. sim 3D_5cm rectangular, minus mean trace (no ringing)
  3. sim 3D_5cm rectangular, + revised ringing, minus mean trace
  4. 6-window median bandpass-envelope profile (rel t0) :
     Line9 (full line) vs sim raw vs sim + ringing (log scale)
"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from analyze_field_ringing import load_csv, WINS  # noqa: E402
from v6r_op import op1_ringing_v07, pooled_q_band, TAU_REV, JITTER_REV  # noqa: E402
from augment_sim2real_v0_6 import SEED  # noqa: E402

DATA = Path(r'E:\automation_djh\real_data_yingshan\营山测线数据')
SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_real_vs_sim.png')

T_REF_NS, GMAX, T_SHOW = 10.0, 30.0, 400.0

# ---- real Line9 ----------------------------------------------------------
x, t, dt, hdr = load_csv(DATA / 'Line9origin(36).csv')
i0 = np.argmax(np.abs(x[:, t <= 60]), axis=1)
pk = np.abs(x[np.arange(len(x)), i0])
valid = pk > 0
x, i0, pk = x[valid], i0[valid], pk[valid]
nrm = x / pk[:, None]
seg = nrm[1000:1200]
i0s = i0[1000:1200]
SEG = np.stack([np.roll(seg[j], -i0s[j]) for j in range(len(seg))])
t_rel = np.arange(SEG.shape[1]) * dt
sel_r = t_rel <= T_SHOW
SEG = SEG[:, sel_r]
t_rel = t_rel[sel_r]

# ---- sim 3D final chain --------------------------------------------------
d = np.load(SFCW)
mats, tm = [], None
for k in range(1, 14):
    t_ = f't{k:02d}'
    mats.append(d[f'3D_5cm_{t_}_tail_200ns_rectangular_real_bandpass'])
    tm = d[f'3D_5cm_{t_}_tail_200ns_rectangular_envelope_time_s']
B0 = np.stack(mats)
tfull = tm * 1e9
q_lo, q_hi, nq = pooled_q_band()
B1, n_clip = op1_ringing_v07(B0, tfull, np.random.default_rng(SEED),
                             q_lo, q_hi, TAU_REV, JITTER_REV)
S0 = B0 - B0.mean(axis=0, keepdims=True)
S1 = B1 - B1.mean(axis=0, keepdims=True)
i0_sim = np.argmax(np.abs(B0), axis=1)
S0 = np.stack([np.roll(S0[j], -i0_sim[j]) for j in range(13)])
S1 = np.stack([np.stack([np.roll(S1[j], -i0_sim[j]) for j in range(13)])])[0]
t_sim = tfull - tfull[0]
sel_s = t_sim <= T_SHOW
S0, S1, t_sim = S0[:, sel_s], S1[:, sel_s], t_sim[sel_s]
# per-trace normalization by direct peak (same contract as real data)
pk_s = np.abs(B0).max(axis=1)
S0 = S0 / pk_s[:, None]
S1 = S1 / pk_s[:, None]

# ---- profiles ------------------------------------------------------------
def make_sos(tt):
    dt_ = float(tt[1] - tt[0])
    return butter(4, [0.090, 0.125], btype='band', fs=1 / dt_, output='sos')


def prof(matr, tt, sos):
    """6-window median envelope profile. matr MUST already be per-trace
    direct-peak normalized by the caller; never self-normalize here (after
    mean removal a trace's own max IS the ringing, which would divide
    ringing by itself)."""
    env = np.abs(hilbert(sosfiltfilt(sos, matr, axis=1), axis=1))
    out = []
    for lo, hi in WINS:
        m = (tt >= lo) & (tt <= hi)
        out.append(float(np.median(env[:, m].mean(axis=1))))
    return out


# real: roll full normalized line to t0
FULL = np.stack([np.roll(nrm[j], -i0[j]) for j in range(len(nrm))])
t_full = np.arange(FULL.shape[1]) * dt
pr = prof(FULL, t_full, make_sos(t_full))
ps0 = prof(S0, t_sim, make_sos(t_sim))
ps1 = prof(S1, t_sim, make_sos(t_sim))
print('6-window median envelope profile (rel t0):')
print('windows            :', WINS)
print('Line9 (real)       :', [round(v, 4) for v in pr])
print('sim raw (no ring)  :', [round(v, 4) for v in ps0])
print('sim + revised ring :', [round(v, 4) for v in ps1])
centers = [(lo + hi) / 2 for lo, hi in WINS]

# ---- figure --------------------------------------------------------------
fig = plt.figure(figsize=(19, 13), dpi=130)
gs = fig.add_gridspec(2, 3, height_ratios=[3, 2])


def show(ax, M, tt, pk, ntr, title):
    """Fixed absolute scale: per-trace direct-peak reference (NOT a per-panel
    percentile -- an empty static zone would let ringing set the scale)."""
    g = np.minimum(tt / T_REF_NS, GMAX)
    P = np.clip(M * g[None, :] / pk[:, None], -1, 1)
    pc = ax.pcolormesh(np.arange(1, ntr + 1), tt, P.T, vmin=-1, vmax=1,
                       cmap='gray', shading='auto')
    ax.set_title(title, fontsize=10)
    ax.set_ylim(T_SHOW, 0)
    return pc


ax1 = fig.add_subplot(gs[0, 0])
show(ax1, SEG, t_rel, pk[1000:1200], SEG.shape[0],
     f'Line9 real | traces 1001-1200 | dt={dt:.1f} ns')
ax1.set_ylabel('two-way time rel t0 (ns)')
ax2 = fig.add_subplot(gs[0, 1])
show(ax2, S0, t_sim, pk_s, 13, 'sim 3D 5cm | rectangular, minus mean (no ringing)')
ax3 = fig.add_subplot(gs[0, 2])
pc = show(ax3, S1, t_sim, pk_s, 13,
          'sim 3D 5cm | + revised ringing, minus mean')
for ax in (ax2, ax3):
    for tt_, lab in [(100.2, 'surf'), (185.1, 'iface')]:
        ax.axhline(tt_, color='red', lw=0.8, ls='--')
        ax.text(13.2, tt_, lab, color='red', fontsize=8, va='center')
cax = fig.add_axes([0.935, 0.56, 0.008, 0.36])
fig.colorbar(pc, cax=cax, label='normalized (white=+, black=-)')

ax4 = fig.add_subplot(gs[1, :])
ax4.plot(centers, pr, 'o-', label='Line9 real (full line, n=%d)' % len(nrm))
ax4.plot(centers, ps0, 's--', label='sim raw (no ringing)')
ax4.plot(centers, ps1, '^-', label='sim + revised ringing (tau=58 ns)')
for lo, hi in WINS:
    ax4.axvspan(lo, hi, alpha=0.06, color='gray')
ax4.set_yscale('log')
ax4.set_xlabel('window center (ns rel t0)')
ax4.set_ylabel('median bandpass envelope\n(per-trace / direct peak)')
ax4.set_title('ring-down profile: real vs sim', fontsize=10)
ax4.grid(True, which='both', alpha=0.3)
ax4.legend()

fig.suptitle('real Line9 vs simulated 3D final chain | identical display:'
             ' per-trace peak-norm, t0-aligned, bipolar gray, capped t-gain')
fig.tight_layout(rect=[0, 0, 0.93, 0.97])
fig.savefig(OUT)
print('saved:', OUT)
