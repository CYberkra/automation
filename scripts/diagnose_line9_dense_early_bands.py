"""Refine exploratory temporal accounting into four native receiver-time bands."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from analyze_line9_v401_version_controls import inverse
from diagnose_line9_dense_time_support import early_weight
from hs_capsule_identity import sha256 as sha

CONFIGS = [[40., 130., 260.], [30., 130., 260.], [50., 130., 260.],
           [40., 110., 260.], [40., 150., 260.], [40., 130., 240.], [40., 130., 280.]]
LABELS = ['早期段（含空气直达主波包）', '中早段（含地表候选波包）',
          '中段（含覆盖层底候选波包）', '晚段（尚未认证路径）']


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    r = read(a.support / 'analysis.json'); audit = read(a.support / 'independent_audit.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(a.support / 'analysis.json')
    assert r['numerical_sha256'] == sha(a.support_numerical)
    with h5py.File(a.support_numerical) as h:
        f = h['frequency_Hz'][:]; ids = list(h['ids'].asstr()[:]); full = h['full'][:]
        old = h['partition_response'][:]
    centers = sorted({v for config in CONFIGS for v in config})
    cumulative = np.empty((501, 4, len(centers)), complex)
    hashes = {row['id']: row['native_sha256'] for row in r['native']}
    for j, sid in enumerate(ids):
        path = a.source / (sid + '.h5'); assert sha(path) == hashes[sid]
        with h5py.File(path) as h:
            x = h['rxs/rx1/Ez'][:]; s = h['srcs/src1/excitation/samples'][:]
            tr = np.arange(len(x)) * h.attrs['dt'] + h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
            ts = np.arange(len(s)) * h.attrs['dt'] + h['srcs/src1/excitation'].attrs['TimeSampleOffset']
        new = [c for c in centers if c not in r['split_centers_native_ns']]
        weighted = x[:, None] * np.column_stack([early_weight(tr * 1e9, c) for c in new])
        values = np.empty((501, len(new)), complex)
        for first in range(0, 501, 16):
            ff = f[first:first + 16, None]
            values[first:first + len(ff)] = (np.exp(-2j * np.pi * ff * tr) @ weighted) / (np.exp(-2j * np.pi * ff * ts) @ s)[:, None] / .025
        for c, center in enumerate(centers):
            cumulative[:, j, c] = old[:, j, r['split_centers_native_ns'].index(center), 0] if center in r['split_centers_native_ns'] else values[:, new.index(center)]
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    metrics = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean(); total, t = inverse(full, w)
        cum, _ = inverse(cumulative.reshape(501, -1), w); cum = cum.reshape(4008, 4, len(centers))
        lo, hi = r['basal_common_gate_ns']; gate = (t * 1e9 >= lo) & (t * 1e9 <= hi); rows = []
        fig, axes = plt.subplots(2, 2, figsize=(15, 9), layout='constrained')
        show = (t * 1e9 >= 300) & (t * 1e9 <= 410)
        for j, sid in enumerate(ids):
            for n, config in enumerate(CONFIGS):
                c0, c1, c2 = [cum[:, j, centers.index(c)] for c in config]
                bands = np.column_stack([c0, c1 - c0, c2 - c1, total[:, j] - c2])
                norm = np.linalg.norm; q = total[gate, j]
                ratios = [float(norm(bands[gate, b]) / norm(q)) for b in range(4)]
                closure = float(norm(bands[gate].sum(axis=1) - q) / norm(q)); assert closure < 1e-10
                rows.append(dict(id=sid, native_cutoff_centers_ns=config, band_L2_over_full=ratios,
                                 closure_over_full_L2=closure))
                if n == 0:
                    ax = axes.flat[j]; ax.plot(t[show] * 1e9, total[show, j].real, color='black', lw=1.7, label='完整场')
                    for b, label in enumerate(LABELS):
                        ax.plot(t[show] * 1e9, bands[show, b].real, lw=1, label=label)
                    ax.axvspan(lo, hi, color='green', alpha=.07)
                    ax.set(title=('高损耗' if sid.startswith('high') else '低损耗诊断') + (' / 无底砂H0' if sid.endswith('H0') else ' / 总场H1'), xlabel='SFCW时间 / ns', ylabel='实部 / (V/m)/(A·m)')
                    ax.legend(fontsize=7)
        metrics[window] = rows
        fig.suptitle(f'190m原生时间四段贡献 / {window} / 原生分界中心40、130、260ns，各20ns互补过渡\n同一完整源归一；各面板独立纵轴；含某波包不等于纯物理分量，范数比不是能量份额，无生产数据修正')
        fig.savefig(a.out / f'early_bands_{window}.png', dpi=130); plt.close(fig)
    with h5py.File(a.numerical, 'x') as h:
        h.create_dataset('frequency_Hz', data=f); h.create_dataset('full', data=full)
        h.create_dataset('cumulative_early_response', data=cumulative)
        h.create_dataset('ids', data=np.array(ids, dtype=h5py.string_dtype('utf-8')))
    result = dict(status='EXPLORATORY_FOUR_NATIVE_TIME_BANDS_NOT_PURE_COMPONENTS',
        script_sha256=sha(__file__), weight_helper_sha256=sha(Path(__file__).with_name('diagnose_line9_dense_time_support.py')),
        support_analysis_sha256=sha(a.support / 'analysis.json'), support_audit_sha256=sha(a.support / 'independent_audit.json'),
        support_numerical_sha256=sha(a.support_numerical), numerical_sha256=sha(a.numerical),
        cumulative_centers_native_ns=centers, transition_width_ns=20., configurations=CONFIGS,
        labels=LABELS, basal_common_gate_ns=r['basal_common_gate_ns'], native=r['native'], metrics=metrics,
        new_solves=0, production_data_changed=False,
        limits='Posthoc temporal decomposition only. Ideal 2D direct field has a causal tail; each time band can contain several paths, so labels are descriptive candidates. A strong early-band contribution at deep SFCW time is finite-band/source-normalized leakage from those native samples, not independently isolated antenna or surface reflection. No claim of energy percentage, noise identity,field match or full-line attribution.')
    (a.out / 'analysis.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], central_metrics={w: [q for q in rows if q['native_cutoff_centers_ns'] == CONFIGS[0]] for w, rows in metrics.items()})))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'support', 'support_numerical', 'out', 'numerical']:
        p.add_argument('--' + key.replace('_', '-'), type=Path, required=True)
    main(p.parse_args())
