"""One-off visual check figure for the batch2d_v1 A0 anchor case (B2D-C3m-D10m-W4m-T0.5m-E20-S0.02).

Read-only over archived h5; no solver, no new analysis claims. Panel 1: YZ geometry schematic of the
A0 model (parsed values, not a re-simulation). Panel 2: raw Ex A-scan with the frozen event windows
(event table v0.1: ti_i_mul_dtoffset0) and self-computed two-way arrival times overlaid.
Panel 3: paired TGT-BG difference with the same markers, so the target response is visible.
"""

from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'artifacts/research_checks'
TGT_DIR = RUNS / '2026-09-26_B2D-C3m-D10m-W4m-T0.5m-E20-S0.02'
BG_DIR = RUNS / '2026-09-26_B2D-C3m-BG'
OUT_DIR = RUNS / '2026-09-26_a0_visual_check'
OUT_DIR.mkdir(parents=True, exist_ok=True)

DT_S = 5.896635841874211e-11
N = 20351  # analysis window, floor(1200ns/dt)+1

# frozen event table v0.1 values (C3 family)
T_COV_NS = 180.124611          # self-computed cover-bottom two-way time
T_TGT_NS = 320.221531          # self-computed D10m target two-way time
COV_WIN = (100.124611, 260.124611)
TGT_WIN = (240.221531, 400.221531)


def load_ex(d: Path):
    h5 = next(d.glob('*.h5'))
    with h5py.File(h5, 'r') as h:
        ex = np.asarray(h['rxs/rx1/Ex'][:N], dtype=np.float64).ravel()
    assert ex.shape[0] == N
    return ex


tgt = load_ex(TGT_DIR)
bg = load_ex(BG_DIR)
diff = tgt - bg
t_ns = np.arange(N) * DT_S * 1e9
env = np.abs(hilbert(diff))

fig = plt.figure(figsize=(13, 11), dpi=130)

# ---- Panel 1: geometry schematic (YZ cross-section, A0 / C3 family) ----
ax1 = fig.add_subplot(3, 1, 1)
Z_TOP, Z_BOT = 46.0, 8.0   # display crop
Y_L, Y_R = 8.0, 25.0
# air
ax1.add_patch(plt.Rectangle((Y_L, 30.0), Y_R - Y_L, Z_TOP - 30.0, fc='#eef4fb', ec='none'))
# cover (silty clay, 3 m)
ax1.add_patch(plt.Rectangle((Y_L, 27.0), Y_R - Y_L, 3.0, fc='#d9c9a3', ec='none'))
# sandstone
ax1.add_patch(plt.Rectangle((Y_L, Z_BOT), Y_R - Y_L, 27.0 - Z_BOT, fc='#c8d6c0', ec='none'))
# target box y14-18, z19.75-20.25
ax1.add_patch(plt.Rectangle((14.0, 19.75), 4.0, 0.5, fc='#b04a3a', ec='k', lw=0.8))
# surface & interfaces
for z, ls in ((30.0, '-'), (27.0, '-')):
    ax1.plot([Y_L, Y_R], [z, z], 'k', lw=1.0, ls=ls)
ax1.annotate('surface z=30 m', xy=(Y_L + 0.3, 30.0), xytext=(Y_L + 0.3, 31.5),
             ha='left', fontsize=8)
ax1.annotate('cover bottom z=27 m (clay er16 / sandstone er9)', xy=(Y_L + 0.3, 27.0),
             xytext=(Y_L + 0.3, 24.0), ha='left', fontsize=8)
# Tx/Rx
ax1.plot(15.35, 45.0, '^', ms=10, c='#1f5fa8')
ax1.plot(16.65, 45.0, 's', ms=9, c='#d07b1f')
ax1.annotate('Tx y15.35 z45', xy=(15.35, 45.0), xytext=(10.0, 37.8), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.7))
ax1.annotate('Rx y16.65 z45\n(baseline 1.3 m)', xy=(16.65, 45.0), xytext=(19.6, 37.8), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.7))
ax1.annotate('target y14-18, z19.75-20.25\ner20 / s0.02, D10m W4m T0.5m', xy=(17.5, 20.0),
             xytext=(11.5, 12.0), ha='left', fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.7))
ax1.set_xlim(Y_L, Y_R)
ax1.set_ylim(Z_BOT, Z_TOP)
ax1.set_xlabel('y (m)')
ax1.set_ylabel('z (m)')
ax1.set_title('Model schematic — A0 anchor  B2D-C3m-D10m-W4m-T0.5m-E20-S0.02  (YZ 2D TM, cover C3 = 3 m)')
ax1.set_aspect('equal')
ax1.grid(alpha=0.25)

# ---- Panel 2: raw A-scan with frozen windows ----
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(t_ns, tgt, lw=0.7, c='#1f5fa8', label='raw Ex A-scan (TGT)')
ax2.axvspan(*COV_WIN, color='#d9c9a3', alpha=0.45, label='frozen COV window [100.12, 260.12] ns')
ax2.axvspan(*TGT_WIN, color='#b04a3a', alpha=0.20, label='frozen TGT window [240.22, 400.22] ns')
ax2.axvline(T_COV_NS, color='#8a6d1f', ls='--', lw=1.0)
ax2.axvline(T_TGT_NS, color='#b04a3a', ls='--', lw=1.0)
ax2.annotate(f'cover bottom (self-computed {T_COV_NS:.1f} ns)', xy=(T_COV_NS, 0),
             xytext=(T_COV_NS - 150, ax2.get_ylim()[1] * 0.6), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.7))
ax2.annotate(f'target D10m (self-computed {T_TGT_NS:.1f} ns)', xy=(T_TGT_NS, 0),
             xytext=(T_TGT_NS + 15, ax2.get_ylim()[1] * 0.8), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.7))
ax2.set_xlim(0, 1200)
ax2.set_xlabel('two-way time (ns)')
ax2.set_ylabel('Ex (V/m, 2D line-source normalisation)')
ax2.set_title('Raw A-scan — strong early direct-coupling/air wave expected; interface & target arrivals marked')
ax2.grid(alpha=0.25)
ax2.legend(fontsize=7, loc='upper right')

# ---- Panel 3: paired difference TGT - BG ----
ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(t_ns, diff, lw=0.7, c='#444444', label='paired difference TGT - BG')
ax3.plot(t_ns, env, lw=0.9, c='#b04a3a', label='|envelope|')
ax3.axvspan(*TGT_WIN, color='#b04a3a', alpha=0.20)
ax3.axvspan(*COV_WIN, color='#d9c9a3', alpha=0.45)
ax3.axvline(T_TGT_NS, color='#b04a3a', ls='--', lw=1.0)
ax3.set_xlim(0, 1200)
ax3.set_xlabel('two-way time (ns)')
ax3.set_ylabel('difference Ex (V/m)')
ax3.set_title('Paired difference — common cover interface cancels; target response should sit inside the TGT window')
ax3.grid(alpha=0.25)
ax3.legend(fontsize=7, loc='upper right')

fig.tight_layout()
out = OUT_DIR / 'a0_model_and_ascan.png'
fig.savefig(out)
print('wrote', out)
print('raw Ex range:', float(tgt.min()), float(tgt.max()))
print('diff peak |Ex|:', float(np.abs(diff).max()), 'at ns:', float(t_ns[np.argmax(np.abs(diff))]))
print('envelope peak ns:', float(t_ns[np.argmax(env)]))
