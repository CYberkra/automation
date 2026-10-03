"""2x3 panel: 2D-5cm and 2D-2.5cm groups, raw / mean-removed / SVD rank-1 removed."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CHECK = Path(r'E:\automation_djh\artifacts_check\r1')
OUT = Path(r'E:\automation_djh\fig_2d_bscan_3panel.png')

GROUPS = [('2D_5cm', '2D 5cm same section'), ('2D_2p5cm', '2D 2.5cm production grid')]
CLIP = 0.03

fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.2), dpi=130, sharey=True)
for row, (tag, gname) in enumerate(GROUPS):
    B = np.load(CHECK / f'bscan_{tag}.npy')
    t_ns = np.load(CHECK / f't_ns_{tag}.npy')
    ntr = B.shape[0]
    traces = np.arange(1, ntr + 1)

    mean_trace = B.mean(axis=0, keepdims=True)
    bg_removed = B - mean_trace
    removed_frac = float(np.sum((B - bg_removed) ** 2) / np.sum(B ** 2))
    U, S, Vt = np.linalg.svd(B, full_matrices=False)
    rank1_frac = float(S[0] ** 2) / float(np.sum(S ** 2))
    svd_filtered = (U[:, 1:] * S[1:]) @ Vt[1:, :]

    sel = t_ns <= 400
    raw_clipped = np.clip(B[:, sel] / (np.max(np.abs(B)) * CLIP), -1, 1)
    panels = [
        (raw_clipped, f'{gname} | raw (clip 3%)'),
        (bg_removed[:, sel], f'mean-trace removed (removed {removed_frac*100:.2f}%)'),
        (svd_filtered[:, sel], f'SVD rank-1 removed (1st SV {rank1_frac*100:.2f}%)'),
    ]
    for col, (M, ttl) in enumerate(panels):
        ax = axes[row, col]
        ref = np.max(np.abs(M)) or 1.0
        Mn = M / ref if M is not raw_clipped else M
        pc = ax.pcolormesh(traces, t_ns[sel], Mn.T, vmin=-1, vmax=1, cmap='RdBu_r', shading='auto')
        for tt, lab, colr in [(100.2, 'surf', 'g'), (185.1, 'iface', 'm')]:
            ax.axhline(tt, color=colr, lw=0.9, ls='--')
            ax.text(13.1, tt, lab, color=colr, fontsize=8, va='center')
        ax.set_title(ttl, fontsize=9)
        ax.set_xticks(traces)
        if row == 1:
            ax.set_xlabel('trace #')
        if col == 0:
            ax.set_ylabel('two-way time (ns)')
            ax.set_ylim(400, 0)
fig.suptitle('benchmark3d_r2_co 2D groups (full-band impulse H5, per-panel normalization)')
fig.tight_layout()
fig.savefig(OUT)
print('saved:', OUT)
