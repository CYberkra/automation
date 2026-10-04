"""House FFT band-pass synthesis on the delivered 13-trace 3D impulse responses.

NOT the official gprMax SFCW toolbox chain (gprMax.toolboxes.SFCW.processing):
this script band-weights the raw Ex rfft with a Hann mask and irffts back --
no source-spectrum normalisation, no tail taper, Hann not mean-normalised,
and the dense rfft grid is not the 501-point 20-170 MHz instrument grid.
Absolute amplitudes therefore do not follow the official convention; use for
display/diagnosis only. Official-chain products live in
artifacts_check/sfcw_r1 (sfcw_adapter_benchmark3d.py).
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
OUT = Path(r'E:\automation_djh\fig_benchmark3d_sfcw_bscan.png')

F0, F1 = 20e6, 170e6          # instrument sweep range (user-provided settings)

traces = []
dt = None
for k in range(1, 14):
    rid = f'{P}-t{k:02d}'
    with h5py.File(SIMS / f'2026-10-01_{rid}' / f'{rid}.h5', 'r') as f:
        traces.append(np.asarray(f['rxs/rx1/Ex']).ravel())
        dt = float(f.attrs['dt'])
B = np.stack(traces)                       # 13 x nsteps raw impulse response
n = B.shape[1]
t_ns = np.arange(n) * dt * 1e9
mid_y = 4.5 + 0.25 * np.arange(13)

# ---- house synthesis: band-weight the raw response rfft, Hann mask ----
freqs = np.fft.rfftfreq(n, d=dt)
H = np.fft.rfft(B, axis=1)
W = np.zeros_like(freqs)
mb = (freqs >= F0) & (freqs <= F1)
# Hann across the sweep band (assumption label per instrument follow-up doc)
W[mb] = np.hanning(mb.sum())
S = np.fft.irfft(H * W[None, :], n=n, axis=1)


def env(a):
    X = np.fft.fft(a, axis=1)
    h = np.zeros(a.shape[1])
    h[0] = 1
    h[1:a.shape[1] // 2] = 2
    if a.shape[1] % 2 == 0:
        h[a.shape[1] // 2] = 1
    return np.abs(np.fft.ifft(X * h, axis=1))


ES = env(S)
ns_sm = max(1, int(round(1.5 / (dt * 1e9))))
ker = np.ones(ns_sm) / ns_sm
ES = np.apply_along_axis(lambda r: np.convolve(r, ker, mode='same'), 1, ES)

fig, axes = plt.subplots(1, 3, figsize=(19.5, 8.2), constrained_layout=True)
ext = [mid_y[0] - 0.125, mid_y[-1] + 0.125]

# (a) raw impulse-response B-scan envelope (the ringy one) zoom
mz = t_ns <= 300
axes[0].imshow(env(B)[:, mz].T, aspect='auto', cmap='gray_r', vmin=0, vmax=0.45,
               extent=[ext[0], ext[1], t_ns[mz][-1], 0])
axes[0].set_title('(a) 原始冲激响应（之前给你看的：全是振铃）')
axes[0].set_xlabel('测线位置 y（m）')
axes[0].set_ylabel('双程走时（ns）')

# (b) SFCW-synthesized B-scan envelope, global scale
axes[1].imshow(ES[:, mz].T, aspect='auto', cmap='gray_r', vmin=0, vmax=np.percentile(ES, 99.5),
               extent=[ext[0], ext[1], t_ns[mz][-1], 0])
axes[1].set_title('(b) house 带通合成后（20–170 MHz，Hann 掩码，非官方链）')
axes[1].set_xlabel('测线位置 y（m）')
axes[1].set_ylabel('双程走时（ns）')

# (c) synthesized, per-trace surface-normalized to hunt for interface texture
surf = ES[:, (t_ns >= 85) & (t_ns <= 140)].max(axis=1)
SN = ES / surf[:, None]
axes[2].imshow(SN[:, mz].T, aspect='auto', cmap='gray_r', vmin=0, vmax=0.6,
               extent=[ext[0], ext[1], t_ns[mz][-1], 0])
axes[2].set_title('(c) 合成后逐道按地表波归一（找界面纹理）')
axes[2].set_xlabel('测线位置 y（m）')
axes[2].set_ylabel('双程走时（ns）')

for ax in axes:
    for tt, lab in [(4.3, '直达'), (100, '地表'), (185, '界面理论')]:
        ax.axhline(tt, color='red', lw=0.8, ls='--', alpha=0.75)
        ax.text(4.28, tt + 5, f'{lab} {tt}ns', color='red', fontsize=9)
axes[2].axhspan(165, 215, color='blue', alpha=0.10)

fig.suptitle('house FFT 带通合成（非官方 SFCW 工具箱）：13 道 3D B-scan（冲激响应 → 20–170 MHz 合成）— 2026-10-02',
             fontsize=13)
fig.savefig(OUT, dpi=150)
print('wrote', OUT)

# quantitative check: coda level in-band
m170 = (t_ns >= 170) & (t_ns <= 220)
print('synth surface env mean: %.4e' % surf.mean())
print('synth env mean 170-220ns: %.4e (ratio to surface %.4f)' % (
    ES[:, m170].mean(), ES[:, m170].mean() / surf.mean()))
