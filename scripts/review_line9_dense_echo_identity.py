"""CPU-only review: distinguish wide-window H0 packets from basal response."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from analyze_line9_v401_version_controls import inverse
from audit_line9_v401_delivery import native_response, sha


def measurements(q, b, times):
    d = q - b
    norm = np.linalg.norm
    return dict(H1_peak_ns=float(times[abs(q).argmax()]), H0_peak_ns=float(times[abs(b).argmax()]),
        delta_peak_ns=float(times[abs(d).argmax()]), H1_norm=float(norm(q)), H0_norm=float(norm(b)),
        delta_norm=float(norm(d)), H1_H0_complex_correlation=float(abs(np.vdot(q, b))/(norm(q)*norm(b))),
        H0_peak_amplitude=float(max(abs(b))), H1_peak_amplitude=float(max(abs(q))))


def main(a):
    read = lambda p: json.loads(p.read_text('utf-8'))
    assert not a.out.exists()
    p = read(a.public / 'analysis.json')
    audit = read(a.public / 'independent_audit.json')
    m = read(a.package / 'manifest.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(a.public / 'analysis.json')
    assert p['numerical_sha256'] == sha(a.numerical) and p['manifest_sha256'] == sha(a.package / 'manifest.json')
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]
        z = h['response'][:]
        ids = [x.decode() if isinstance(x, bytes) else x for x in h['ids'][:]]
    np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
    receipts = {r['id']: r['native_sha256'] for r in p['native']}
    assert set(ids) == set(receipts)
    anchors = [190., 187.5]
    independent = {}
    dft_errors = {}
    for xpos in anchors:
        for variant in ['high', 'low']:
            for role in ['H0', 'H1']:
                sid = f'{variant}_x{round(xpos*100):05d}_{role}'
                independent[sid] = native_response(a.source / (sid + '.h5'), .025, receipts[sid])
                err = float(np.linalg.norm(independent[sid]-z[:, ids.index(sid)])/np.linalg.norm(independent[sid]))
                assert err < 1e-9
                dft_errors[sid] = err
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    windows = {}
    inverse_errors = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean()
        x, t = inverse(z, w)
        ns = t * 1e9
        k = (ns >= 300) & (ns <= 450)
        lookup = {sid: x[:, j] for j, sid in enumerate(ids)}
        direct = {sid: np.exp(2j*np.pi*t[k, None]*f) @ (q*w)/501 for sid, q in independent.items()}
        errors = {sid: float(np.linalg.norm(q-lookup[sid][k])/np.linalg.norm(q)) for sid, q in direct.items()}
        assert max(errors.values()) < 1e-9
        inverse_errors[window] = errors
        rows = []
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
        common = max(float(max(abs(lookup[f'{v}_x{round(s*100):05d}_H0'][k]))) for v in ['high', 'low'] for s in anchors)
        for j, xpos in enumerate(anchors):
            station = next(s for s in m['stations'] if s['chainage_m'] == xpos)
            lo = min(s['basal_gate_ns'][0] for s in station['templates'].values())
            hi = max(s['basal_gate_ns'][1] for s in station['templates'].values())
            kk = (ns >= lo) & (ns <= hi)
            row = dict(chainage_m=xpos, wide_gate_ns=[300, 450], basal_union_gate_ns=[lo, hi], wide={}, basal={})
            for v, label, color in [('high', '原高损耗', '#b3423c'), ('low', '低损耗诊断', '#24779c')]:
                q = lookup[f'{v}_x{round(xpos*100):05d}_H1']
                b = lookup[f'{v}_x{round(xpos*100):05d}_H0']
                row['wide'][v] = measurements(q[k], b[k], ns[k])
                check = measurements(direct[f'{v}_x{round(xpos*100):05d}_H1'], direct[f'{v}_x{round(xpos*100):05d}_H0'], ns[k])
                for name, value in check.items():
                    assert abs(value-row['wide'][v][name])/max(1, abs(value)) < 1e-8
                row['basal'][v] = measurements(q[kk], b[kk], ns[kk])
                axes[0, j].plot(ns[k], abs(b[k]), color=color, label=label+' H0')
            for gate in ['wide', 'basal']:
                row[gate]['low_over_high_H0_norm'] = row[gate]['low']['H0_norm']/row[gate]['high']['H0_norm']
                row[gate]['low_over_high_delta_norm'] = row[gate]['low']['delta_norm']/row[gate]['high']['delta_norm']
            q = lookup[f'high_x{round(xpos*100):05d}_H1'][k]
            b = lookup[f'high_x{round(xpos*100):05d}_H0'][k]
            for value, label, color, style in [(q, '高损耗H1总场', '#222222', '-'), (b, '高损耗H0（已移除底砂）', '#bd4b40', '--'), (q-b, '底砂对比差场H1−H0', '#16845a', '-')]:
                axes[1, j].plot(ns[k], abs(value), color=color, ls=style, label=label)
            for r in [0, 1]:
                axes[r, j].set(xlabel='SFCW时间 / ns', ylabel='复包络幅度 / (V/m)/(A·m)', xlim=[300, 450], ylim=[0, common*1.12],
                    title=f'{xpos:g}m / '+('H0高低配方绝对尺度对照' if r == 0 else '高损耗总场强峰是否随底砂移除'))
                axes[r, j].axvspan(lo, hi, color='#cde5d7', alpha=.25, label='事前底砂共同窗')
                axes[r, j].legend(fontsize=8)
                axes[r, j].grid(alpha=.15)
            rows.append(row)
        tracks = {}
        for v in ['high', 'low']:
            tracks[v] = []
            for station in m['stations']:
                sid = f'{v}_x{round(station["chainage_m"]*100):05d}_H1'
                if sid in lookup:
                    tracks[v].append(dict(chainage_m=station['chainage_m'], wide_peak_ns=float(ns[k][abs(lookup[sid][k]).argmax()]),
                        local_template_peak_ns=station['templates'][v]['peak_ns']))
        windows[window] = dict(anchors=rows, wide_total_peak_tracks=tracks)
        fig.suptitle(f'已有两锚点再审 / {window} / 4.0.1原生FP64 / 20–170MHz、501点\n'
            '全部同一绝对幅度尺度；H0含真实上覆地质，差场含交互；包络幅度不可相加当能量')
        fig.savefig(a.out / f'echo_identity_{window}.png', dpi=130)
        plt.close(fig)
    result = dict(status='PASS_EXISTING_DATA_REVIEW_CAUSAL_PATH_OPEN', script_sha256=sha(__file__),
        source_analysis_sha256=sha(a.public / 'analysis.json'), source_audit_sha256=sha(a.public / 'independent_audit.json'),
        numerical_sha256=sha(a.numerical), native_sha256={sid: receipts[sid] for sid in independent},
        independent_native_DFT_relative_L2=dft_errors, independent_direct_inverse_relative_L2=inverse_errors,
        windows=windows, limits='CPU reanalysis only. Same pre-existing windows, no fit or simulation. H0 not noise; removing basal contrast does not isolate a pure primary. Correlation is descriptive, not causal-path certification.')
    (a.out / 'review.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], max_DFT_error=max(dft_errors.values()), windows=windows), ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'package', 'numerical', 'out']:
        p.add_argument('--'+key, type=Path, required=True)
    main(p.parse_args())
