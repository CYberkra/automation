"""Exploratory native-time partition: early band tails versus late native support.

This is a linear accounting diagnostic, never a production cleaning operator or
an identification of a reflected path. All partitions use the same full source.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from analyze_line9_v401_version_controls import FREQ, inverse
from hs_capsule_identity import sha256 as sha

CENTERS_NS = [240., 260., 280.]
TRANSITION_NS = 20.


def early_weight(time_ns, center):
    u = np.clip((time_ns - center + TRANSITION_NS / 2) / TRANSITION_NS, 0, 1)
    return .5 * (1 + np.cos(np.pi * u))


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    p = read(a.public / 'analysis.json'); audit = read(a.public / 'independent_audit.json')
    m = read(a.package / 'manifest.json')
    assert audit['analysis_sha256'] == sha(a.public / 'analysis.json')
    assert audit['status'].startswith('PASS') and p['numerical_sha256'] == sha(a.response)
    identities = {r['id']: r['native_sha256'] for r in p['native']}
    ids = [f'{v}_x19000_{role}' for v in ['high', 'low'] for role in ['H0', 'H1']]
    assert all(sid in identities for sid in ids)
    with h5py.File(a.response) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:], FREQ)
        names = list(h['ids'].asstr()[:]); full = h['response'][:, [names.index(s) for s in ids]]
    partitions = np.empty((501, 4, 3, 2), complex); closure = []; raw_info = []
    for j, sid in enumerate(ids):
        path = a.source / (sid + '.h5'); assert sha(path) == identities[sid]
        with h5py.File(path) as h:
            x = h['rxs/rx1/Ez'][:]; s = h['srcs/src1/excitation/samples'][:]
            assert x.dtype == s.dtype == np.float64 and np.isfinite(x).all()
            dt = float(h.attrs['dt'])
            tr = np.arange(len(x)) * dt + h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
            ts = np.arange(len(s)) * dt + h['srcs/src1/excitation'].attrs['TimeSampleOffset']
        ew = np.column_stack([early_weight(tr * 1e9, c) for c in CENTERS_NS])
        split = x[:, None, None] * np.stack([ew, 1 - ew], axis=-1)
        raw_info.append(dict(id=sid, native_sha256=identities[sid],
                             source_abs_peak_ns=float(ts[np.argmax(abs(s))] * 1e9),
                             raw_partition_closure_max=float(abs(split.sum(axis=2) - x[:, None]).max())))
        for first in range(0, 501, 16):
            f = FREQ[first:first + 16, None]
            denominator = np.exp(-2j * np.pi * f * ts) @ s
            assert np.all(abs(denominator) > 0)
            parts = np.exp(-2j * np.pi * f * tr) @ split.reshape(len(x), 6)
            partitions[first:first + len(f), j] = (parts / denominator[:, None] / .025).reshape(len(f), 3, 2)
        err = np.linalg.norm(partitions[:, j].sum(axis=-1) - full[:, j, None]) / (np.linalg.norm(full[:, j]) * np.sqrt(3))
        assert err < 1e-9; closure.append(float(err))
    station = next(s for s in m['stations'] if s['chainage_m'] == 190.)
    gate = [min(station['templates'][v]['basal_gate_ns'][0] for v in ['high', 'low']),
            max(station['templates'][v]['basal_gate_ns'][1] for v in ['high', 'low'])]
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    metrics = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean()
        total, t = inverse(full, w)
        split, _ = inverse(partitions.reshape(501, -1), w)
        split = split.reshape(4008, 4, 3, 2)
        k = (t * 1e9 >= gate[0]) & (t * 1e9 <= gate[1]); rows = []
        for j, sid in enumerate(ids):
            norm = np.linalg.norm; q = total[k, j]
            for c, center in enumerate(CENTERS_NS):
                e, late = split[k, j, c, 0], split[k, j, c, 1]
                rows.append(dict(id=sid, split_center_native_ns=center,
                    early_over_total_L2=float(norm(e) / norm(q)),
                    late_over_total_L2=float(norm(late) / norm(q)),
                    early_late_normalized_real_inner_product=float(np.vdot(e, late).real / (norm(e) * norm(late))),
                    closure_over_total_L2=float(norm(e + late - q) / norm(q))))
        metrics[window] = rows
        fig, axes = plt.subplots(2, 2, figsize=(14, 9), layout='constrained')
        show = (t * 1e9 >= 300) & (t * 1e9 <= 410)
        for j, sid in enumerate(ids):
            ax = axes.flat[j]
            ax.plot(t[show] * 1e9, total[show, j].real, color='black', lw=1.8, label='原完整场')
            ax.plot(t[show] * 1e9, split[show, j, 1, 0].real, color='#d58000', label='早段贡献（原生时间≤250ns，250–270ns过渡）')
            ax.plot(t[show] * 1e9, split[show, j, 1, 1].real, color='#0079a5', label='晚段贡献（原生时间≥270ns，过渡互补）')
            ax.axvspan(*gate, alpha=.07, color='green')
            ax.set(title=('高损耗' if sid.startswith('high') else '低损耗诊断') + (' / 无底砂H0' if sid.endswith('H0') else ' / 原总场H1'),
                   xlabel='SFCW时间 / ns', ylabel='实部 / (V/m)/(A·m)')
            ax.legend(fontsize=7)
        fig.suptitle(f'190m站位原生早晚时间支持诊断 / {window} / 同一完整源归一 / 20–170MHz精确501点\n各面板独立纵轴；早+晚=完整复场，范数比不是能量份额；时间分段不认证多次波或散射路径')
        fig.savefig(a.out / f'native_time_support_{window}.png', dpi=130); plt.close(fig)
    with h5py.File(a.numerical, 'x') as h:
        h.create_dataset('frequency_Hz', data=FREQ); h.create_dataset('full', data=full)
        h.create_dataset('partition_response', data=partitions)
        h.create_dataset('ids', data=np.array(ids, dtype=h5py.string_dtype('utf-8')))
    result = dict(status='EXPLORATORY_LINEAR_NATIVE_TIME_SUPPORT_NOT_PATH_IDENTIFICATION',
        script_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        independent_analysis_audit_sha256=sha(a.public / 'independent_audit.json'),
        manifest_sha256=sha(a.package / 'manifest.json'), response_sha256=sha(a.response),
        numerical_sha256=sha(a.numerical), split_centers_native_ns=CENTERS_NS,
        transition_width_ns=TRANSITION_NS, basal_common_gate_ns=gate,
        source_time_offsets_retained=True, full_source_used_for_all_parts=True,
        native=raw_info, spectral_closure_relative_L2=closure, metrics=metrics,
        new_solves=0, production_data_changed=False,
        limits='Exploratory receiver-time partition, not physical component separation or a proposed cleaning operator. Late-native support can contain true geology,multiples,lateral paths,boundary effects or numerical errors. Early support can contain multiple events. Source deconvolution and finite-band tails retained. L2 ratios are not additive energy fractions; no field validation.')
    (a.out / 'analysis.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], gate_ns=gate, metrics=metrics)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'package', 'public', 'response', 'out', 'numerical']:
        parser.add_argument('--' + key, type=Path, required=True)
    main(parser.parse_args())
