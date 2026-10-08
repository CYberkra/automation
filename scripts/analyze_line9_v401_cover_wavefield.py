"""Exact501-tone single-station geometry controls; no per-column normalization or fit."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import spectral_window

from analyze_line9_v401_version_controls import response, inverse, FREQ
from hs_capsule_identity import sha256 as sha

IDS = ['base_snapshot_H0', 'flat_snapshot_H0', 'far_cover_removed_H0', 'near_cover_removed_H0']
LABELS = ['原几何H0＋快照', '全横向平层H0＋快照', '远处174–184m覆盖层底移除', '近处184–194m覆盖层底移除']
GATES = dict(early=[0., 120.], wide=[300., 450.], inherited_basal=[345.6180971390552, 374.6081170991351], late=[450., 1100.])


def metrics(profiles, time_ns):
    result = {}
    for name, bounds in GATES.items():
        keep = (time_ns >= bounds[0]) & (time_ns <= bounds[1])
        x = profiles[keep]
        ref = x[:, 0]
        norm = np.linalg.norm
        assert norm(ref) > 0
        rows = []
        for j, sid in enumerate(IDS):
            q = x[:, j]
            assert norm(q) > 0
            rows.append(dict(id=sid, norm=float(norm(q)), norm_over_base=float(norm(q) / norm(ref)),
                difference_over_base=float(norm(ref - q) / norm(ref)),
                correlation_with_base=float(abs(np.vdot(ref, q)) / (norm(ref) * norm(q))),
                peak_ns=float(time_ns[keep][np.argmax(abs(q))])))
        result[name] = dict(gate_ns=bounds, rows=rows)
    return result


def main(a):
    read = lambda p: json.loads(p.read_text('utf-8'))
    assert not a.out.exists() and not a.numerical.exists()
    c = read(a.source / 'execution_contract.json')
    v = read(a.source / 'completed_verification.json')
    m = c['study_manifest']
    assert v['completed'] and v['contract_sha256'] == sha(a.source / 'execution_contract.json')
    assert [g['id'] for g in c['groups']] == [r['id'] for r in v['groups']] == IDS
    assert v['groups'][0]['observer_receiver_and_source_bitwise_equal']
    assert [r['snapshot_count'] for r in v['groups']] == [250, 250, 0, 0]
    zs = []
    records = []
    reference_source = None
    for g, r in zip(c['groups'], v['groups']):
        path = a.source / (g['id'] + '.h5')
        assert sha(path) == r['native_sha256']
        with h5py.File(path) as h:
            assert h.attrs['gprMax'] == '4.0.1' and h.attrs['Iterations'] == 20352
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [8400, 1700, 1])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025] * 3)
            assert abs(h.attrs['dt'] / m['dt_s'] - 1) < 1e-14
            for node, role in [('srcs/src1', 'tx_m'), ('rxs/rx1', 'rx_m')]:
                np.testing.assert_allclose(h[node].attrs['Position'], g[role], rtol=0, atol=1e-12)
            x = h['srcs/src1/excitation/samples'][:]
            assert x.dtype == np.float64 and np.isfinite(x).all()
            if reference_source is None:
                reference_source = x
            else:
                assert x.tobytes() == reference_source.tobytes()
        z, error = response(path, .025)
        zs.append(z)
        records.append(dict(id=g['id'], native_sha256=r['native_sha256'], direct_DFT_relative_L2=error))
    z = np.column_stack(zs)
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    summaries = {}
    for window in ['hann', 'blackman']:
        profiles, t = inverse(z, spectral_window(window, 501))
        tn = t * 1e9
        summaries[window] = metrics(profiles, tn)
        for name, bounds in [('early', GATES['early']), ('main', [150., 450.]), ('late', GATES['late'])]:
            keep = (tn >= bounds[0]) & (tn <= bounds[1])
            x = profiles[keep]
            difference = x[:, :1] - x
            lim = float(abs(x.real).max())
            dlim = float(abs(difference.real).max())
            fig, axes = plt.subplots(2, 2, figsize=(14, 9), layout='constrained')
            for ax, data, title, limit in [(axes[0, 0], x.real, '四配置H0剩余场：共同绝对灰度', lim),
                                          (axes[0, 1], difference.real, '原H0−控制：共同差场灰度（另设范围）', dlim)]:
                im = ax.imshow(data, cmap='gray', vmin=-limit, vmax=limit, aspect='auto', interpolation='nearest',
                               extent=[-.5, 3.5, tn[keep][-1], tn[keep][0]])
                ax.set_xticks(range(4), LABELS, rotation=15, fontsize=8)
                ax.set(title=title, ylabel='SFCW时间 / ns')
                fig.colorbar(im, ax=ax, label='(V/m)/(A·m)')
            for j, label in enumerate(LABELS):
                axes[1, 0].plot(tn[keep], abs(x[:, j]), lw=1, label=label)
                if j:
                    axes[1, 1].plot(tn[keep], abs(difference[:, j]), lw=1, label='原H0−' + label)
            for ax, title in [(axes[1, 0], 'H0复包络幅值（线性物理量）'), (axes[1, 1], '先复数相减，再取幅值；非纯多次波')]:
                ax.set(xlabel='SFCW时间 / ns', ylabel='(V/m)/(A·m)', title=title)
                ax.legend(fontsize=8)
                if name == 'main':
                    ax.axvspan(*GATES['inherited_basal'], color='green', alpha=.08, label='原底砂到时窗')
            fig.suptitle(f'190m单站四配置 / 高损耗H0 / {window} / 4.0.1原生FP64 / 约8m航高\n所有配置均去底砂；横轴是配置，不是空间测线。区域替换改变边缘与全部交互，无AGC或幅相拟合。')
            fig.savefig(a.out / f'cover_controls_{window}_{name}.png', dpi=130)
            plt.close(fig)
    with h5py.File(a.numerical, 'x') as h:
        h['frequency_Hz'] = FREQ
        h['response'] = z
        h['ids'] = np.array(IDS, dtype=h5py.string_dtype('utf-8'))
    report = dict(status='COMPLETE_FOUR_SINGLE_STATION_GEOMETRY_CONTROLS_NOT_UNIQUE_PATH_OR_FIELD_VALIDATION',
        script_sha256=sha(__file__), contract_sha256=sha(a.source / 'execution_contract.json'),
        verification_sha256=sha(a.source / 'completed_verification.json'), numerical_sha256=sha(a.numerical),
        native=records, metrics=summaries,
        limits='H0 retains real cover interface. Regional edits add edges/change continuation/all interactions. Flat changes whole lateral geology. No AGC, phase/delay/amplitude fit, per-column normalization, independent real-data validation or unique multiple attribution.')
    (a.out / 'analysis.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=report['status'], metrics=summaries)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'out', 'numerical']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
