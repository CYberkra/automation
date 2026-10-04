"""Bipolar grayscale B-scans for 3D_5cm: t-gain (capped) raw vs interface-only.

Left: full B-scan with deterministic time gain g(t)=min(t/t_ref, GMAX)
(spreading-compensation style display gain, capped at +29.5 dB), normalized
to the 99.5th percentile. The ~100 MHz stripes are Hann-window spectral
leakage ringing of the strong events (raw device-band late-time floor is
-126 dB, i.e. numerically empty), not simulation content and not V6 ringing.

Right: 'interface only' = per-time-sample mean-trace removal (kills all flat
events: direct, surface, leakage stripes) + hard window 185+/-25 ns.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
OUT = Path(r'E:\automation_djh\fig_3d_bipolar_interface.png')

TAG = '3D_5cm'
t_surf = 100.2
t_iface = 185.1
win = (160.0, 210.0)
T_REF_NS = 10.0
GMAX = 30.0  # linear gain cap (+29.5 dB)

d = np.load(SFCW)
mats = []
tm = None
for k in range(1, 14):
    t = f't{k:02d}'
    mats.append(d[f'{TAG}_{t}_tail_200ns_hann_real_bandpass'])
    tm = d[f'{TAG}_{t}_tail_200ns_hann_envelope_time_s']
tm = tm * 1e9
sel = tm <= 400
B = np.stack([m[sel] for m in mats])
t = tm[sel]

ref = np.max(np.abs(B))

# deterministic t-gain with cap (display only)
g = np.minimum(t / T_REF_NS, GMAX)
Bg = B * g[None, :]
clip = np.nanpercentile(np.abs(Bg), 99.5)
raw_disp = np.clip(Bg / clip, -1, 1)

# interface isolation: mean-trace removal + time window
Bmr = B - B.mean(axis=0, keepdims=True)
w = (t >= win[0]) & (t <= win[1])
iface_only = np.where(w, Bmr, 0.0)
iface_max = np.max(np.abs(iface_only))

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), dpi=130, sharey=True)
panels = [
    (raw_disp,
     f'raw B-scan with t-gain g=min(t/{T_REF_NS:.0f}ns, x{GMAX:.0f})'
     f' (+{20*np.log10(GMAX):.1f} dB cap), 99.5-pct norm\n'
     '200 ns tail taper. white=+, black=-'),
    (iface_only / iface_max,
     f'interface only: mean-trace removed + window {win[0]:.0f}-{win[1]:.0f}'
     ' ns\n(flat events cancelled; sloping arrival = rough interface)\n'
     'normalized to window max'),
]
for ax, (M, ttl) in zip(axes, panels):
    pc = ax.pcolormesh(np.arange(1, 14), t, M.T, vmin=-1, vmax=1,
                       cmap='gray', shading='auto')
    for tt, lab in [(t_surf, 'surf'), (t_iface, 'iface')]:
        ax.axhline(tt, color='red', lw=0.8, ls='--')
        ax.text(13.15, tt, lab, color='red', fontsize=8, va='center')
    ax.axhspan(win[0], win[1], color='blue', alpha=0.05, lw=0)
    ax.set_xlabel('trace #')
    ax.set_title(ttl, fontsize=9.5)
    ax.set_xticks(range(1, 14))
axes[0].set_ylabel('two-way time (ns)')
axes[0].set_ylim(400, 0)
fig.suptitle(f'benchmark3d_r2_co {TAG}: bipolar grayscale, t-gain (capped)'
             ' raw vs bedrock-cover interface (mean-removed + windowed)')
fig.tight_layout()
fig.savefig(OUT)

idx = np.argmax(np.abs(iface_only)[:, w], axis=1)
tw = t[w]
print('saved:', OUT)
print(f'raw global max |A| = {ref:.4e}')
print(f'interface-window max |A| = {iface_max:.4e}'
      f' = {20*np.log10(iface_max/ref):.1f} dB re raw max')
print('gain at surf/iface/window-end:',
      round(np.interp(100.2, t, g), 1), round(np.interp(185.1, t, g), 1),
      round(np.interp(210, t, g), 1))
print('per-trace argmax time in window (ns):',
      np.round(tw[idx], 1).tolist())
