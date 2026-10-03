"""Show actual stepped profiles alongside shared-scale underground B-scans."""
import argparse
import json
from pathlib import Path

import numpy as np

from hs_capsule_identity import sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', type=Path, required=True)
    ap.add_argument('--profile', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new output directory required')
    records = json.loads((args.results / 'manifest.json').read_text('utf-8'))
    for row in records:
        if sha256(args.results / row['file']) != row['sha256']:
            raise ValueError('analysis capsule identity changed')
    profile_digest = sha256(args.profile)
    profile = np.genfromtxt(args.profile, delimiter=',', names=True)
    with np.load(args.results / 'comparison_arrays.npz') as data:
        t, x = data['time_ns'], data['midpoint_m']
        mask = (t >= 160) & (t <= 220)
        matrices = [data[name + '_raw_signed'][mask].copy()
                    for name in ('baseline', 'relief2p0')]
    limit = float(np.quantile(np.abs(np.stack(matrices)), .999)) or 1.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Microsoft YaHei'
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True,
                             gridspec_kw={'height_ratios': [1, 2]})
    edges = np.r_[profile['x0_m'], profile['x1_m'][-1]]
    for j, (z, title) in enumerate(zip(
            (profile['original_z_m'], profile['grid_z_m']), ('原版 0.4 m 起伏', '新版 0.8 m 起伏'))):
        ax = axes[0, j]
        ax.stairs(12-z, edges, baseline=None, color='tab:brown', linewidth=2)
        ax.axhline(2.95, color='gray', linestyle='--', label='平界面对照：2.95 m')
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(3.5, 2.4)
        ax.set_ylabel('距水平地表深度（m）')
        ax.set_title(title + '：输入台阶剖面')
        ax.legend(loc='lower left')
        ax = axes[1, j]
        im = ax.pcolormesh(x, t[mask], matrices[j], cmap='gray', shading='nearest',
                           vmin=-limit, vmax=limit, rasterized=True)
        ax.set_ylim(220, 160)
        ax.set_xlabel('实际收发中点 X（m）')
        ax.set_ylabel('重建时间（ns）')
        ax.set_title(title + '：带限响应，未去背景')
        fig.colorbar(im, ax=ax, label='场/源响应（左右共用色标）')
    fig.suptitle('gprMax V4.0.0 / CUDA FP64；20–170 MHz / 501频点 / Hann\n'
                 '起伏指全模型高差；图中截取收发中点范围；160–220 ns 局部共享色标\n'
                 '未去背景、未增益；上方深度不转换为下方时间')
    args.out.mkdir(parents=True)
    fig.savefig(args.out / 'geometry_and_bscan.png', dpi=145)
    plt.close(fig)
    if sha256(args.profile) != profile_digest:
        raise ValueError('profile changed during plotting')
    report = {'results_manifest_sha256': sha256(args.results / 'manifest.json'),
              'profile_sha256': profile_digest, 'code_sha256': sha256(__file__),
              'window_ns': [160, 220], 'symmetric_color_limit': limit,
              'role': 'window-specific shared-scale diagnostic; no time-depth calibration'}
    (args.out / 'provenance.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    (args.out / 'manifest.json').write_text(json.dumps([
        {'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha256(p)}
        for p in sorted(args.out.iterdir())], indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
