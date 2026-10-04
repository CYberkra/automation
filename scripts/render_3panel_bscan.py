"""Three-panel B-scan per user spec: raw / mean-trace removed / SVD rank-1 removed.

Axes: x = trace number, y = two-way travel time (ns). 3D group, 13 traces.
Each panel normalized to its own max |value|; energy notes in titles.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CHECK = Path(r'E:\automation_djh\artifacts_check\r1')
OUT = Path(r'E:\automation_djh\fig_3d_bscan_3panel.png')

B = np.load(CHECK / 'bscan_3D_5cm.npy')          # 13 x nsteps
t_ns = np.load(CHECK / 't_ns_2D_5cm.npy') if (CHECK / 't_ns_2D_5cm.npy').exists() else np.load(CHECK / 't_ns_3D_5cm.npy')
t_ns = np.load(CHECK / 't_ns_3D_5cm.npy')
ntr, nsteps = B.shape
traces = np.arange(1, ntr + 1)

# panel 1: raw
raw = B

# panel 2: minus mean trace (background = ensemble mean)
mean_trace = B.mean(axis=0, keepdims=True)
bg_removed = B - mean_trace
removed_frac = float(np.sum((B - bg_removed) ** 2) / np.sum(B ** 2))

# panel 3: SVD, remove first singular component
U, S, Vt = np.linalg.svd(B, full_matrices=False)
e_total = float(np.sum(S ** 2))
rank1_frac = float(S[0] ** 2) / e_total
svd_filtered = (U[:, 1:] * S[1:]) @ Vt[1:, :]

sel = t_ns <= 400
CLIP = 0.03  # raw panel clip at 3% of max |amp| so late events are visible
raw_clipped = np.clip(raw[:, sel] / (np.max(np.abs(raw)) * CLIP), -1, 1)
panels = [
    (raw_clipped, f'raw B-scan (clipped at {CLIP*100:.0f}% of max |amp|)'),
    (bg_removed[:, sel], f'mean-trace removed (removed {removed_frac*100:.2f}% energy = common direct wave)'),
    (svd_filtered[:, sel], f'SVD rank-1 removed (1st SV = {rank1_frac*100:.2f}% energy)'),
]

fig, axes = plt.subplots(1, 3, figsize=(17, 6.4), dpi=130, sharey=True)
for ax, (M, title) in zip(axes, panels):
    ref = np.max(np.abs(M)) or 1.0
    Mn = M / ref
    pc = ax.pcolormesh(traces, t_ns[sel], Mn.T, vmin=-1, vmax=1, cmap='RdBu_r', shading='auto')
    ax.set_xlabel('trace #')
    ax.set_title(title, fontsize=10)
    ax.set_xticks(traces)
    fig.colorbar(pc, ax=ax, label='normalized amplitude')
axes[0].set_ylabel('two-way travel time (ns)')
axes[0].set_ylim(400, 0)
fig.suptitle('benchmark3d_r2_co 3D group (38M cells, 5 cm) — per-panel max|amp| normalization')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
print('rank-1 SV energy frac: %.4f' % rank1_frac)
print('top-5 SV energy fractions:', (S[:5] ** 2 / e_total).round(4).tolist())
