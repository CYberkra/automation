"""Impulse vs Ricker-100MHz comparison at t07 (rough interface, 12x12 domain).

F = rough + Ricker100, G = rock-only + Ricker100, A = rough + impulse.
F-G isolates the cover-layer contribution under a band-limited source.
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

V = Path(r'E:\automation_djh\verify_runs\2026-10-02_iface_sensitivity')
OUT = Path(r'E:\automation_djh\fig_benchmark3d_ricker_vs_impulse.png')


def load(p):
    with h5py.File(p, 'r') as f:
        return np.asarray(f['rxs/rx1/Ex']).ravel(), float(f.attrs['dt'])


def env(a):
    X = np.fft.fft(a)
    h = np.zeros(len(a))
    h[0] = 1
    h[1:len(a) // 2] = 2
    if len(a) % 2 == 0:
        h[len(a) // 2] = 1
    return np.abs(np.fft.ifft(X * h))


A, dt = load(V / 'ctl3dA_t07_asis.h5')
F, _ = load(V / 'ctl3dF_t07_ricker100.h5')
G, _ = load(V / 'ctl3dG_t07_rockonly_ricker100.h5')
t = np.arange(len(A)) * dt * 1e9
m = t <= 320

fig, axes = plt.subplots(1, 3, figsize=(18, 6.4), constrained_layout=True)

axes[0].plot(t[m], A[m], lw=0.5, color='steelblue')
axes[0].set_title('(a) 冲激源（现标杆）：长振铃拖尾 0.1 量级')
axes[0].set_xlabel('双程走时（ns）')
axes[0].set_ylabel('Ex（V/m）')
axes[0].set_ylim(-0.6, 0.6)

axes[1].plot(t[m], F[m], lw=0.6, color='darkgreen')
axes[1].set_title('(b) Ricker 100MHz：拖尾降约 160 倍')
axes[1].set_xlabel('双程走时（ns）')
axes[1].set_ylim(-0.6, 0.6)

dFG = F - G
axes[2].plot(t[m], F[m] * 20, lw=0.6, color='darkgreen', label='F ×20（糙界面 Ricker）')
axes[2].plot(t[m], dFG[m], lw=0.8, color='purple', label='F−G = 覆盖层+界面贡献')
axes[2].axvspan(165, 215, color='blue', alpha=0.10)
axes[2].axhline(0, color='k', lw=0.5)
axes[2].set_title('(c) 放大 20 倍：界面窗内回波仅 ~1.2e-4（地表的 0.6%）')
axes[2].set_xlabel('双程走时（ns）')
axes[2].legend(fontsize=9)

for ax in axes:
    for tt in (4.3, 100, 185):
        ax.axvline(tt, color='red', ls='--', lw=0.8, alpha=0.6)

fig.suptitle('去掉冲激源振铃后的 3D 波形（t07 糙界面，12×12 域，GPU double）— 2026-10-02', fontsize=13)
fig.savefig(OUT, dpi=150)
print('wrote', OUT)
