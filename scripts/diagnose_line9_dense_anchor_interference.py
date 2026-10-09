"""Paired anchor interference in predeclared union and fixed wide windows; no fit."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from analyze_line9_v401_version_controls import inverse, FREQ
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    p = read(a.public / 'analysis.json')
    audit = read(a.public / 'independent_audit.json')
    m = read(a.package / 'manifest.json')
    assert p['manifest_sha256'] == sha(a.package / 'manifest.json')
    assert audit['analysis_sha256'] == sha(a.public / 'analysis.json') and audit['status'].startswith('PASS')
    assert p['numerical_sha256'] == sha(a.numerical)
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:], FREQ)
        ids = h['ids'].asstr()[:].tolist()
        z = h['response'][:]
    anchors = []
    for chain in m['anchors_m']:
        names = [f'{v}_x{round(chain * 100):05d}_{role}' for v in ['high', 'low'] for role in ['H0', 'H1']]
        if all(sid in ids for sid in names):
            anchors.append(next(s for s in m['stations'] if s['chainage_m'] == chain))
    assert [s['chainage_m'] for s in anchors] == [190., 187.5]
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    metrics = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean()
        profiles, t = inverse(z, w)
        tn = t * 1e9
        show = (tn >= 300) & (tn <= 450)
        fig, axes = plt.subplots(2, 2, figsize=(15, 9), layout='constrained')
        rows = []
        shared = {v: max(float(abs(profiles[show, ids.index(f'{v}_x{round(s["chainage_m"] * 100):05d}_{role}')].real).max())
                         for s in anchors for role in ['H0', 'H1']) for v in ['high', 'low']}
        for j, station in enumerate(anchors):
            chain = station['chainage_m']
            union = [min(station['templates'][v]['basal_gate_ns'][0] for v in ['high', 'low']),
                     max(station['templates'][v]['basal_gate_ns'][1] for v in ['high', 'low'])]
            for i, variant in enumerate(['high', 'low']):
                prefix = f'{variant}_x{round(chain * 100):05d}_'
                b = profiles[:, ids.index(prefix + 'H0')]
                q = profiles[:, ids.index(prefix + 'H1')]
                d = q - b
                for name, bounds in [('basal_union', union), ('fixed_wide', [300., 450.])]:
                    keep = (tn >= bounds[0]) & (tn <= bounds[1])
                    bb, dd, qq = b[keep], d[keep], q[keep]
                    norm = np.linalg.norm
                    nb, nd, nq = (float(norm(x)) for x in [bb, dd, qq])
                    inner = np.vdot(dd, bb)
                    incoherent = nb * nb + nd * nd
                    closure = abs(nq * nq - incoherent - 2 * inner.real) / (nq * nq)
                    assert closure < 1e-12
                    rows.append(dict(chainage_m=chain, variant=variant, window=name, gate_ns=bounds,
                        H0_norm=nb, delta_norm=nd, H1_norm=nq, H0_over_delta=nb / nd,
                        H0_delta_inner_product_real=float(inner.real), H0_delta_inner_product_imag=float(inner.imag),
                        H0_vs_delta_inner_phase_deg=float(np.angle(inner) * 180 / np.pi),
                        complex_coherence_magnitude=float(abs(inner) / (nb * nd)),
                        cross_term_over_incoherent_squared_norm=float(2 * inner.real / incoherent),
                        squared_norm_identity_relative_error=float(closure),
                        H1_peak_ns=float(tn[keep][np.argmax(abs(qq))]), delta_peak_ns=float(tn[keep][np.argmax(abs(dd))])))
                ax = axes[i, j]
                for value, label, color in [(q, 'H1原总场', 'black'), (b, 'H0保留覆盖层、去底砂', 'orange'), (d, 'H1−H0底砂对比（含交互）', '#1674bc')]:
                    ax.plot(tn[show], value[show].real, color=color, lw=1, label=label)
                ax.axvspan(*union, color='green', alpha=.09, label='事前两配方底砂窗并集')
                ax.axvline(station['templates'][variant]['peak_ns'], color='green', ls='--', lw=.8, label='局部一次模板，非真值')
                ax.set(title=f'{chain:g}m / ' + ('原高损耗' if variant == 'high' else '低损耗诊断'),
                    xlabel='SFCW时间 / ns', ylabel='(V/m)/(A·m)', ylim=(-1.05 * shared[variant], 1.05 * shared[variant]))
                ax.legend(fontsize=8)
        metrics[window] = rows
        fig.suptitle(f'两配对锚点相干叠加 / {window} / 精确501点 / 原生FP64 / 无AGC或幅相拟合\n同行两站共享纵轴；高低损耗两行纵轴不同。H0含真实地层，差场不等于纯一次波；平方范数恒等式不是物理能量标定。')
        fig.savefig(a.out / f'anchor_interference_{window}.png', dpi=130)
        plt.close(fig)
    result = dict(status='TWO_ANCHOR_PAIRED_INTERFERENCE_DIAGNOSTIC_NOT_H0_PATH_IDENTIFICATION',
        script_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        independent_source_audit_sha256=sha(a.public / 'independent_audit.json'),
        numerical_sha256=sha(a.numerical), manifest_sha256=sha(a.package / 'manifest.json'),
        anchors_m=[s['chainage_m'] for s in anchors], metrics=metrics, calls_solver=False,
        limits='Stored complex-array identity within fixed declared gates,not unique multiple/PML/path attribution. Phase alone does not measure interference when coherence is small; report actual cross term. Not a certified evaluation direction cosine; FDTD error bound remains unresolved. H0 retains cover interfaces. Low mud is not field calibrated. No3D/field/training certification.')
    (a.out / 'analysis.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['public', 'package', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
