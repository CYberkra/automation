"""Final verdict figure: boundary-contamination dominance + buried interface echo.

Controls (GPU double, t07 geometry):
  A  rough interface, 12x12 domain (bit-exact reproduction of delivered t07)
  B  flat interface z=9, 12x12 domain
  C2 rock only, 12x12 domain
  E  flat interface z=9, 12x16 domain (y walls +2 m each side)
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
OUT = Path(r'E:\automation_djh\fig_benchmark3d_iface_echo_anatomy.png')


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
B, _ = load(V / 'ctl3dB_t07_flat9p00.h5')
C, _ = load(V / 'ctl3dC2_t07_rockonly.h5')
E, _ = load(V / 'ctl3dE_t07_flat_bigy.h5')
t = np.arange(len(A)) * dt * 1e9

fig, axes = plt.subplots(2, 2, figsize=(15.5, 9.5), constrained_layout=True)

ax = axes[0][0]
m = t <= 280
ax.semilogy(t[m], env(A)[m], lw=0.7, label='A 糙界面 12x12（=交付t07）')
ax.semilogy(t[m], env(E)[m], lw=0.7, label='E 平界面 12x16（侧墙外移2m）')
ax.semilogy(t[m], env(C)[m], lw=0.7, alpha=0.7, label='C2 纯岩石 12x12')
ax.set_xlim(0, 280)
ax.set_ylim(1e-4, 30)
ax.set_title('(a) 包络对比：拖尾在纯岩石模型中同样存在 → 非地质信号')
ax.set_xlabel('双程走时（ns）')
ax.set_ylabel('包络（V/m）')
for tt in (4.3, 100, 185):
    ax.axvline(tt, color='red', ls='--', lw=0.8, alpha=0.6)
ax.legend(fontsize=9)

ax = axes[0][1]
dBE = B - E
ax.plot(t[m], dBE[m], lw=0.6, color='teal', label='B−E（仅差侧墙位置）')
ax.set_xlim(0, 280)
ax.set_title('(b) 移动侧墙即改变 0.19 量级 → 中后期响应主要是边界混响')
ax.set_xlabel('双程走时（ns）')
ax.axhline(0, color='k', lw=0.5)
ax.legend(fontsize=9)

ax = axes[1][0]
BC = B - C
ax.plot(t[m], BC[m], lw=0.7, color='darkgreen', label='B−C2 = 覆盖层全部贡献（边界相消）')
ax.plot(t[m], env(BC)[m], lw=1.0, color='orange', label='包络(B−C2)')
ax.axvspan(174, 196, color='blue', alpha=0.10)
ax.set_xlim(60, 280)
ax.set_title('(c) 覆盖层真实贡献：界面理论窗内包络峰仅 0.016 @174.6ns')
ax.set_xlabel('双程走时（ns）')
ax.legend(fontsize=9)

ax = axes[1][1]
AB = A - B
ax.plot(t[m], AB[m], lw=0.7, color='purple', label='A−B = ±0.4m 界面起伏的影响')
ax.set_xlim(60, 280)
ax.set_title('(d) 起伏影响 ≤5.4e-5（直达波的 5e-6）→ 起伏纹理不可测')
ax.set_xlabel('双程走时（ns）')
ax.axhline(0, color='k', lw=0.5)
ax.legend(fontsize=9)

fig.suptitle('3D 标杆验收决定性对照（t07，GPU double）— 2026-10-02', fontsize=13)
fig.savefig(OUT, dpi=150)
print('wrote', OUT)
