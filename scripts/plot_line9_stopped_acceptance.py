"""Audited partial-line 300--450 ns view, with explicit independent display scales."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read = lambda p: json.loads(p.read_text('utf-8'))
    p = read(a.public / 'analysis.json')
    audit = read(a.public / 'independent_audit.json')
    m = read(a.package / 'manifest.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(a.public / 'analysis.json')
    assert p['numerical_sha256'] == sha(a.numerical) and p['manifest_sha256'] == sha(a.package / 'manifest.json')
    target = a.public / 'acceptance_hann_300_450.png'
    assert not target.exists()
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]
        z = h['response'][:]
        ids = [x.decode() if isinstance(x, bytes) else x for x in h['ids'][:]]
    np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
    t = np.arange(4008) / (4008 * 300000.)
    k = (t * 1e9 >= 300) & (t * 1e9 <= 450)
    w = np.hanning(501)
    w /= w.mean()
    profiles = np.exp(2j * np.pi * t[k, None] * f) @ (z * w[:, None]) / 501
    lookup = {sid: profiles[:, j] for j, sid in enumerate(ids)}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    cmap = plt.get_cmap('gray').copy()
    cmap.set_bad('#d6e5ee')
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
    for j, (v, label) in enumerate([('high', '原高损耗泥岩'), ('low', '低损耗诊断泥岩')]):
        total = np.full((k.sum(), len(m['chainage_m'])), np.nan + 1j * np.nan)
        delta = np.full_like(total, np.nan + 1j * np.nan)
        paired = []
        for col, xpos in enumerate(m['chainage_m']):
            prefix = f'{v}_x{round(xpos*100):05d}_'
            if prefix + 'H1' in lookup:
                total[:, col] = lookup[prefix + 'H1']
                if prefix + 'H0' in lookup:
                    delta[:, col] = lookup[prefix + 'H1'] - lookup[prefix + 'H0']
                    paired.append(xpos)
        for row, data, title in [(0, total, f'H1总场（{sum(p["availability"][v])}站）'),
                                 (1, delta, f'H1−H0底砂差场（仅{len(paired)}站）')]:
            lim = float(np.nanmax(abs(data.real)))
            im = axes[row, j].imshow(np.ma.masked_invalid(data.real), cmap=cmap, vmin=-lim, vmax=lim,
                aspect='auto', interpolation='nearest', extent=[190.125, 184.875, t[k][-1]*1e9, t[k][0]*1e9])
            axes[row, j].set(title=label + ' / ' + title, xlabel='剖面里程 / m（190→185）', ylabel='SFCW时间 / ns')
            axes[row, j].plot(m['chainage_m'], [s['templates'][v]['peak_ns'] for s in m['stations']],
                color='#008a47', ls='--', lw=1, label='局部底砂一次模板（非真值）')
            axes[row, j].legend(fontsize=8)
            fig.colorbar(im, ax=axes[row, j], label='(V/m)/(A·m)')
    fig.suptitle('用户停止后的局部验收 / 4.0.1原生FP64 / AGL约8m / 501频点 / Hann\n'
        '各面板独立灰度，仅比较波形位置；幅度请读各自色标。蓝灰为未算/无配对，不插值、AGC或拟合')
    fig.savefig(target, dpi=130)
    plt.close(fig)
    metadata = dict(status='AUDITED_INPUT_PARTIAL_DISPLAY', script_sha256=sha(__file__),
        analysis_sha256=sha(a.public / 'analysis.json'), numerical_sha256=sha(a.numerical),
        figure_sha256=sha(target), displayed_gate_ns=[300, 450], window='hann',
        scales='Independent symmetric linear scales; no cross-panel brightness comparison.',
        limits='Plot of audited complex profiles; H1-H0 is complete material-contrast difference, not pure primary or field clean truth.')
    (a.public / 'acceptance_figure.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(metadata))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['package', 'public', 'numerical']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
