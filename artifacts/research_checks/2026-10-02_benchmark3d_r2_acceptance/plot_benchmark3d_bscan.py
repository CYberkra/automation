"""3D benchmark B-scan figure v2 (envelope display, no moire).

Raw impulse response is oscillatory at grid scale -> display Hilbert envelope.
Panels: (a) full 600 ns envelope, t^1.5 gain; (b) zoom 0-260 ns envelope,
global clip; (c) zoom, per-trace surface-normalized. Arrival markers on b/c.
"""
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

font_manager.fontManager.addfont(r'C:\Windows\Fonts\msyh.ttc')
plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

SIMS = Path(r'E:\automation_djh\automation_repo\artifacts\simulations')
P = 'B3D5CM-C3mR2-BG-CO13'
OUT = Path(r'E:\automation_djh\fig_benchmark3d_r2_bscan.png')

traces = []
dt = None
for k in range(1, 14):
    rid = f'{P}-t{k:02d}'
    with h5py.File(SIMS / f'2026-10-01_{rid}' / f'{rid}.h5', 'r') as f:
        traces.append(np.asarray(f['rxs/rx1/Ex']).ravel())
        dt = float(f.attrs['dt'])
B = np.stack(traces)
t_ns = np.arange(B.shape[1]) * dt * 1e9
mid_y = 4.5 + 0.25 * np.arange(13)


def hilbert_env(a):
    X = np.fft.fft(a, axis=1)
    h = np.zeros(a.shape[1])
    h[0] = 1
    h[1:a.shape[1] // 2] = 2
    if a.shape[1] % 2 == 0:
        h[a.shape[1] // 2] = 1
    return np.abs(np.fft.ifft(X * h, axis=1))


E = hilbert_env(B)                       # 13 x nsteps envelope
# light time smoothing (2 ns boxcar) to kill residual sample-scale texture
ns = max(1, int(round(2.0 / (dt * 1e9))))
ker = np.ones(ns) / ns
E = np.apply_along_axis(lambda r: np.convolve(r, ker, mode='same'), 1, E)

fig, axes = plt.subplots(1, 3, figsize=(19.5, 8.2), constrained_layout=True)
ext = [mid_y[0] - 0.125, mid_y[-1] + 0.125, None, 0]

# (a) full window, t^1.5 gain on envelope
g = (np.maximum(t_ns, 1.0) / 100.0) ** 1.5
Ea = E * g[None, :]
va = np.percentile(Ea, 99.0)
axes[0].imshow(Ea.T, aspect='auto', cmap='gray_r', vmin=0, vmax=va,
               extent=[ext[0], ext[1], t_ns[-1], 0])
axes[0].set_title('(a) 包络全窗口 0–600 ns（t$^{1.5}$ 增益）')
axes[0].set_xlabel('测线位置 y（m）')
axes[0].set_ylabel('双程走时（ns）')

# (b) zoom, global clip
mz = t_ns <= 260
Ez = E[:, mz]
tz = t_ns[mz]
axes[1].imshow(Ez.T, aspect='auto', cmap='gray_r', vmin=0, vmax=0.45,
               extent=[ext[0], ext[1], tz[-1], 0])
axes[1].set_title('(b) 包络 0–260 ns，统一截幅 0.45')
axes[1].set_xlabel('测线位置 y（m）')
axes[1].set_ylabel('双程走时（ns）')

# (c) zoom, per-trace surface-normalized
surf = np.max(E[:, (t_ns >= 85) & (t_ns <= 125)], axis=1)
Ec = Ez / surf[:, None]
axes[2].imshow(Ec.T, aspect='auto', cmap='gray_r', vmin=0, vmax=0.85,
               extent=[ext[0], ext[1], tz[-1], 0])
axes[2].set_title('(c) 包络逐道按地表波归一（看界面带纹理）')
axes[2].set_xlabel('测线位置 y（m）')
axes[2].set_ylabel('双程走时（ns）')

for ax in axes[1:]:
    for tt, lab in [(4.3, '直达'), (100, '地表'), (185, '界面理论')]:
        ax.axhline(tt, color='red', lw=0.8, ls='--', alpha=0.75)
        ax.text(4.32, tt + 4, f'{lab} {tt}ns', color='red', fontsize=9)
    ax.axhspan(174, 196, color='blue', alpha=0.10)

pk_t = [tz := None]
env_pk_t = []
for k in range(13):
    m = (t_ns >= 165) & (t_ns <= 215)
    env_pk_t.append(t_ns[m][int(np.argmax(E[k][m]))])
axes[2].plot(mid_y, env_pk_t, 'b.-', ms=5, lw=0.9, label='包络峰 165–215ns')
axes[2].legend(loc='lower right', fontsize=9)

fig.suptitle('benchmark3d_r2_co 3D 组 13 道 B-scan（冲激响应包络，灰度）— 2026-10-02 验收',
             fontsize=13)
fig.savefig(OUT, dpi=150)
print('wrote', OUT)
