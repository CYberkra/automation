"""Render per-case GIF animations from the verified 15 m per-permittivity wavefield capsule.

Pure visualisation of completed evidence; no solver. One GIF per non-dispersive
scan variant (eps6, eps9, eps18nd): three panels per frame — rough-model Ey
(SymLog), interface-contrast field (rough-halfspace, SymLog) and its log10
envelope. Scales are fixed across ALL cases (global extrema) so the GIFs are
directly comparable. Every second frame (2.0049 ns step, 300 frames). Physical
axes: z=8 m (deep) at top, z=22 m (air) at bottom; true relief, surface and
the 15 m antenna pair marked.
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
CAP = ROOT/'artifacts/research_checks/2026-10-04_hs4_permittivity_wavefield'
OUT = ROOT/'artifacts/research_checks/2026-10-04_hs4_permittivity_wavefield_visual'
XS = 13.5 + 0.15*(np.arange(60)+0.5)
ZS = 8.0 + 0.15*(np.arange(134)+0.5)
STEP = 2
CASES = ['eps6', 'eps9', 'eps18nd']


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
    import hashlib
    import sys
    sys.path.insert(0, str(ROOT/'scripts'))
    def sha256(p):
        h = hashlib.sha256()
        with open(p, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
        return h.hexdigest()
    c = json.loads((CAP/'execution_contract.json').read_text('utf-8'))
    v = json.loads((CAP/'completed_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(CAP/'execution_contract.json'):
        raise ValueError('completed verified wavefield capsule required')
    dt = c['dt_s']
    seq = {}
    for tag in CASES:
        fr = frames(f'w15_{tag}_rough')[::STEP]
        fh = frames(f'w15_{tag}_halfspace')[::STEP]
        assert [i for i, _ in fr] == [i for i, _ in fh]
        seq[tag] = (fr, fh)
    nframes = len(seq['eps6'][0])
    assert all(len(f) == nframes for f, _ in seq.values())
    # Fixed scales from global extrema over ALL rendered frames of ALL cases.
    rough_peak = diff_peak = 0.0
    for tag in CASES:
        for (ir, pr), (_, ph) in zip(*seq[tag]):
            a = load_ey(pr)
            rough_peak = max(rough_peak, float(np.abs(a).max()))
            diff_peak = max(diff_peak, float(np.abs(a - load_ey(ph)).max()))
    OUT.mkdir(parents=True)
    for tag in CASES:
        pngs = []
        for k, ((ir, pr), (_, ph)) in enumerate(zip(*seq[tag])):
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
                ax.plot([17.6, 18.9], [27, 27], 'k^', ms=4)
                ax.set_ylim(28.1, 8)
                ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)')
            fig.suptitle(f'15 m antenna height, cover {tag} (non-dispersive), t = {t:.1f} ns (frame {k+1}/{nframes})')
            fig.tight_layout()
            p = OUT/f'frame_{tag}_{k:04d}.png'
            fig.savefig(p, dpi=85)
            plt.close(fig)
            pngs.append(p)
        images = [Image.open(p) for p in pngs]
        gif = OUT/f'wavefield_15m_{tag}.gif'
        images[0].save(gif, save_all=True, append_images=images[1:], duration=100, loop=0)
        for p in pngs:
            p.unlink()
        print(tag, 'gif_bytes', gif.stat().st_size, flush=True)
    print('frames_per_gif', nframes, 'rough_peak', rough_peak, 'diff_peak', diff_peak)


if __name__ == '__main__':
    main()
