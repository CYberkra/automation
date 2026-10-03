"""Generate intuition figures for the status widget (reads archived CO11 h5)."""
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(r'E:\automation_djh\automation_repo')
OUT = Path(r'C:\Users\ROG\AppData\Roaming\kimi-desktop\daimon-share\daimon\agents\main\blueprint\widgets\widget_d07af666-1a90-44a8-96a9-8de58971ab05\assets')
OUT.mkdir(parents=True, exist_ok=True)

DT = 5.896635841874211e-11
N = 20352
t_ns = np.arange(N) * DT * 1e9

TGT = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02'
BG = 'B2D-C3m-BG'

def read_trace(run_id, tag):
    p = REPO / 'artifacts/simulations' / f'2026-09-26_{run_id}-CO11-{tag}' / f'{run_id}-CO11-{tag}.h5'
    with h5py.File(p, 'r') as f:
        return f['rxs/rx1/Ex'][:].astype(np.float64)

tags = [f't{k:02d}' for k in range(1, 12)]
tgt = np.stack([read_trace(TGT, t) for t in tags], axis=1)  # [samples, 11]
bg = np.stack([read_trace(BG, t) for t in tags], axis=1)
diff = tgt - bg
y_pos = 12.65 + np.arange(11) * 1.0  # Rx y, m

# ---------- fig 1: b-scan panel ----------
def gain(a, pct=99.5):
    c = np.percentile(np.abs(a), pct)
    return np.clip(a / c, -1, 1)

windows = [
    (100.1, 260.1, 'cover interface'),
    (240.2, 400.2, 'target D10m window'),
]
fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.1), sharey=True)
for ax, data, title in (
    (axes[0], gain(tgt), 'raw target gather (TGT)'),
    (axes[1], gain(diff), 'background removed: TGT - BG'),
):
    ax.imshow(data, aspect='auto', cmap='gray', extent=[y_pos[0], y_pos[-1], t_ns[-1], t_ns[0]],
              vmin=-1, vmax=1, interpolation='nearest')
    for lo, hi, name in windows:
        ax.axhspan(lo, hi, color='tab:orange', alpha=0.13, lw=0)
    ax.set_xlabel('Rx position y (m)')
    ax.set_title(title, fontsize=10)
    ax.set_xlim(y_pos[0], y_pos[-1])
    ax.set_ylim(1200, 0)
axes[0].set_ylabel('two-way time (ns)')
axes[1].text(12.8, 300, 'difference energy 100% in 250-450 ns', color='tab:red', fontsize=8)
axes[0].text(12.8, 700, 'target response buried\nin cover/background energy', color='tab:blue', fontsize=8)
fig.tight_layout()
fig.savefig(OUT / 'bscan_panel.png', dpi=150)
plt.close(fig)

# ---------- fig 2: a-scan t05 ----------
fig, axes = plt.subplots(2, 1, figsize=(11.0, 3.4), sharex=True)
ax = axes[0]
ax.plot(t_ns, bg[:, 4], color='0.45', lw=0.7, label='BG (no target)')
ax.plot(t_ns, tgt[:, 4], color='tab:blue', lw=0.7, alpha=0.85, label='TGT (target at 10 m)')
for lo, hi, name in windows:
    ax.axvspan(lo, hi, color='tab:orange', alpha=0.13, lw=0)
ax.set_ylabel('raw Ex (V/m)')
ax.legend(fontsize=8, loc='upper right')
ax.set_title('single trace t05 (y = 16.65 m): raw A-scans of BG vs TGT nearly overlap', fontsize=10)
ax = axes[1]
amp = 800.0
d = np.clip(diff[:, 4] * amp, -1, 1)
ax.plot(t_ns, d, color='tab:blue', lw=0.7)
for lo, hi, name in windows:
    ax.axvspan(lo, hi, color='tab:orange', alpha=0.13, lw=0)
ax.text(410, 0.72, 'simulated target response lands inside the D10m event window', color='tab:red', fontsize=8)
ax.set_ylabel('diff x800 (clipped)')
ax.set_xlabel('two-way time (ns)')
ax.set_ylim(-1.05, 1.05)
fig.tight_layout()
fig.savefig(OUT / 'ascan_t05.png', dpi=150)
plt.close(fig)

# ---------- fig 3: damage ladder retention ----------
ops = ['identity', 'mean lam0.25', 'mean lam1.0', 'svd k1']
dmg_types = ['amplitude\nscale', 'polarity\nflip', 'sample\nshift', 'trace\ndeletion']
D = np.array([
    [0.0, 0.227, 0.909, 0.916],
    [0.0, 0.189, 0.758, 0.839],
    [0.0, 0.227, 0.906, 0.916],
    [0.0, 0.233, 0.931, 0.936],
])
colors = ['0.55', 'tab:blue', 'tab:red', 'tab:purple']
x = np.arange(len(dmg_types))
w = 0.2
fig, ax = plt.subplots(figsize=(7.4, 3.0))
for i, (op, c) in enumerate(zip(ops, colors)):
    ax.bar(x + (i - 1.5) * w, D[:, i], width=w, label=op, color=c, alpha=0.9)
ax.set_xticks(x)
ax.set_xticklabels(dmg_types, fontsize=9)
ax.set_ylabel('median D (constructed ref)')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8, ncol=4, loc='upper left')
ax.set_title('damage ladder: stronger background removal erases more local damage', fontsize=10)
ax.text(0.99, 0.40, 'D ~ 0: damage preserved\nD -> 1: damage erased', transform=ax.transAxes,
        ha='right', fontsize=8, color='0.25')
fig.tight_layout()
fig.savefig(OUT / 'ladder_retention.png', dpi=150)
plt.close(fig)

for p in sorted(OUT.glob('*.png')):
    print(p.name, p.stat().st_size)
