"""Render GIF animations from the verified 8 m dense wavefield capsule.

Pure visualisation of completed evidence; no solver. Three panels per frame:
rough-model Ey (SymLog, fixed scale), interface-contrast field
(rough-halfspace, SymLog fixed scale), and its log10 envelope. Every second
frame is used (2.0049 ns step, 300 frames). Physical axes: z=8 m (deep) at
top, z=22 m (air) at bottom; true relief, surface and antenna marked.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
CAP = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
OUT = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_visual'
XS = 13.5 + 0.15*(np.arange(60)+0.5)
ZS = 8.0 + 0.15*(np.arange(134)+0.5)
STEP = 2  # frames


def interface_relief():
    xs, zs = [], []
    for line in (ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre/centre_rough/profile.in').read_text('utf-8').splitlines():
        if not line.startswith('#box:'):
            continue
        p = [float(v) for v in line.split()[1:7]]
        if p[5] == 12.0 and p[1] == 0.0 and p[2] > 5.0:
            xs.append(0.5*(p[0]+p[3])); zs.append(p[2])
    order = np.argsort(xs)
    return np.array(xs)[order], np.array(zs)[order]


IFACE_X, IFACE_Z = interface_relief()


def frames(group):
    out = []
    for f in sorted((CAP/group).glob('profile_snaps/*.h5')):
        with h5py.File(f) as h:
            out.append((int(h.attrs['iteration']), f))
    return out


def load_ey(path):
    with h5py.File(path) as h:
        return h['Ey'][:, 0, :]


def main():
    if OUT.exists():
        raise SystemExit('new visual capsule required')
    c = json.loads((CAP/'execution_contract.json').read_text('utf-8'))
    dt = c['dt_s']
    fr = frames('mid8_rough')[::STEP]
    fh = frames('mid8_halfspace')[::STEP]
    assert [i for i,_ in fr] == [i for i,_ in fh]
    # Fixed scales from global extrema over the rendered sequence.
    rough_peak = 0.0
    diff_peak = 0.0
    for (ir, pr), (_, ph) in zip(fr, fh):
        a = load_ey(pr)
        rough_peak = max(rough_peak, float(np.abs(a).max()))
        diff_peak = max(diff_peak, float(np.abs(a - load_ey(ph)).max()))
    OUT.mkdir(parents=True)
    pngs = []
    for k, ((ir, pr), (_, ph)) in enumerate(zip(fr, fh)):
        rough = load_ey(pr)
        diff = rough - load_ey(ph)
        t = ir*dt*1e9
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        for ax, data, peak, title in [
                (axes[0], rough, rough_peak, f'rough Ey (SymLog, fixed ±{rough_peak:.3g})'),
                (axes[1], diff, diff_peak, f'rough-halfspace Ey (±{diff_peak:.3g})')]:
            norm = SymLogNorm(linthresh=peak*1e-3, vmin=-peak, vmax=peak)
            ax.pcolormesh(XS, ZS, data.T, cmap='RdBu_r', norm=norm, shading='auto')
            ax.set_title(title, fontsize=9)
        env = np.log10(np.maximum(np.abs(diff), diff_peak*1e-12)/diff_peak)
        axes[2].pcolormesh(XS, ZS, env.T, cmap='viridis', vmin=-8, vmax=0, shading='auto')
        axes[2].set_title('|rough-halfspace| log10 re peak', fontsize=9)
        for ax in axes:
            ax.plot(IFACE_X, IFACE_Z, 'k-', lw=.8)
            ax.axhline(12, color='k', ls=':', lw=.6)
            ax.plot([17.6, 18.9], [20, 20], 'k^', ms=4)
            ax.set_ylim(22, 8)
            ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
        fig.suptitle(f'8 m antenna height, t = {t:.1f} ns (frame {k+1}/{len(fr)})')
        fig.tight_layout()
        p = OUT/f'frame_{k:04d}.png'
        fig.savefig(p, dpi=85)
        plt.close(fig)
        pngs.append(p)
    images = [Image.open(p) for p in pngs]
    gif = OUT/'height8m_wavefield.gif'
    images[0].save(gif, save_all=True, append_images=images[1:], duration=100, loop=0)
    for p in pngs:
        p.unlink()
    print('frames', len(pngs), 'rough_peak', rough_peak, 'diff_peak', diff_peak,
          'gif_bytes', gif.stat().st_size)


if __name__=='__main__':
    main()
