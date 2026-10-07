"""Independent raw-to-deliverable checks and alternate-window layer diagnostics."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from review_line9_result_packages import FREQ, C0, EPS0, sha, save, rel


def basal_primary(geometry, materials):
    # Explicit scalar Debye evaluation, independent of indices()/primary_paths().
    refraction = []
    refraction95 = []
    for m in materials.values():
        base = m['base']
        def dielectric(f):
            e = base['relative_permittivity'] - 1j*base['electric_conductivity_s_per_m']/(2*np.pi*f*EPS0)
            for p in m.get('poles', []):
                e = e + p['relative_permittivity_difference']/(1+2j*np.pi*f*p['relaxation_time_s'])
            return np.sqrt(e)
        refraction.append(dielectric(FREQ))
        refraction95.append(dielectric(95e6).real)
    n = np.array(refraction)
    lengths = np.zeros(len(materials))
    lengths[0] = geometry['midpoint_agl_m']
    trans = np.ones(len(FREQ), complex)
    previous = 0.
    for b in geometry['boundaries']:
        a, c = b['above'], b['below']
        lengths[a] += b['depth_m'] - previous
        previous = b['depth_m']
        if b == geometry['basal_sand']:
            reflection = (n[a]-n[c])/(n[a]+n[c])
            extra = b['time95_ns']*1e-9 - 2*np.dot(lengths, refraction95)/C0
            phase_distance = np.sum(n*lengths[:, None], axis=0)
            return trans*reflection*np.exp(-2j*np.pi*FREQ*(2*phase_distance/C0+extra))
        trans *= (2*n[a]/(n[a]+n[c]))*(2*n[c]/(n[a]+n[c]))
    raise ValueError('Basal contact absent')


def corr(a, b):
    return float(abs(np.sum(a.conj()*b))/np.sqrt(np.sum(abs(a)**2)*np.sum(abs(b)**2)))


def main(root, review, products_report, kinematics, out):
    if out.exists():
        raise ValueError('Fresh independent output directory required')
    product_report = json.loads(products_report.read_text('utf-8'))
    kin = json.loads(kinematics.read_text('utf-8'))
    kin_by_package = {p['package']: p for p in kin['packages']}
    results = []
    out.mkdir(parents=True)
    for p in product_report['packages']:
        package = p['package']
        ap = review/(package+'_audit.json')
        if sha(ap) != p['audit_sha256']:
            raise ValueError('Audit identity mismatch')
        audit = json.loads(ap.read_text('utf-8'))
        rows = audit['records']
        derived = Path(p['private_hdf5_path'])
        if sha(derived) != p['private_hdf5_sha256']:
            raise ValueError('Deliverable changed')
        checks = []
        with h5py.File(derived) as h:
            np.testing.assert_array_equal(h['frequency'][:], FREQ)
            if h.attrs['NativeDtype'] != 'float32':
                raise ValueError('Expected explicit native FP32 provenance')
            response = h['response'][:]
            time = h['time_response/time'][:]
            actual_hann = h['time_response/complex_bandpass'][:]
            blackman = h['diagnostic_blackman/complex_bandpass'][:]
            gaussian = h['diagnostic_gaussian/complex_bandpass'][:]
            assert response.shape == (501, len(rows))
            assert h['planned_completed_mask'][:].sum() == len(rows)
            np.testing.assert_array_equal(h['chainage_m'][:], [r['chainage_m'] for r in rows])
        # ALL 501 physical tones, first/middle/last station of each original package.
        for j in sorted(set([0, len(rows)//2, len(rows)-1])):
            raw = root/package/'cases'/rows[j]['id']/'profile.h5'
            if sha(raw) != rows[j]['native_sha256']:
                raise ValueError('Original changed')
            with h5py.File(raw) as h:
                y = h['rxs/rx1/Ez'][:]
                x = h['srcs/src1/excitation/samples'][:]
                dt = float(h.attrs['dt'])
                yt = np.arange(len(y))*dt+float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset'])
                xt = np.arange(len(x))*dt+float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
                scale = float(h['srcs/src1/excitation'].attrs['SpatialScale'])
            pieces = []
            for start in range(0, 501, 61):
                f = FREQ[start:start+61, None]
                pieces.append((np.exp(-2j*np.pi*f*yt)@y)/(np.exp(-2j*np.pi*f*xt)@x)/scale)
            error = rel(np.concatenate(pieces), response[:, j])
            assert error < 1e-9
            checks.append(dict(station=rows[j]['id'], independent_all501_DFT_relative_L2=error))
        mf = next((root/package/'geometries').glob('*.json'))
        if sha(mf) != audit['materials_sha256']:
            raise ValueError('Materials changed')
        materials = json.loads(mf.read_text('utf-8'))['materials']
        primary = np.column_stack([basal_primary(r['geometry'], materials) for r in rows])
        metrics = {}
        models = {}
        for window, w, actual in [('hann', np.hanning(501), actual_hann),
                                   ('blackman', np.blackman(501), blackman)]:
            w /= w.mean()
            # Explicit inverse sum in the whole region containing all basal gates.
            centers = np.array(kin_by_package[package]['interfaces']['basal_sand']['hann']['primary_peak_ns'])
            selected = (time*1e9 >= centers.min()-20)&(time*1e9 <= centers.max()+20)
            small_t = time[selected]
            exponent = np.exp(2j*np.pi*small_t[:, None]*FREQ)
            model = exponent@(w[:, None]*primary)/501
            exact_data = exponent@(w[:, None]*response)/501
            inverse_error = rel(exact_data, actual[selected])
            assert inverse_error < 1e-9
            peaks = small_t[np.argmax(abs(model), axis=0)]*1e9
            mask = abs(small_t[:, None]*1e9-peaks[None, :])<=12
            a, b = model[mask], actual[selected][mask]
            coefficient = np.sum(a.conj()*b)/np.sum(abs(a)**2)
            n = len(rows)//2
            first = np.arange(len(rows))[None, :]<n
            aa, bb = model[mask&first], actual[selected][mask&first]
            first_fit = np.sum(aa.conj()*bb)/np.sum(abs(aa)**2)
            other = mask&~first
            metrics[window] = dict(complex_correlation=corr(a, b),
                explicit_inverse_relative_L2=inverse_error,
                peak_ns_endpoints=peaks[[0,-1]].tolist(),
                first_half_fit_second_half_relative_L2=rel(first_fit*model[other], actual[selected][other]),
                second_half_complex_correlation=corr(model[other], actual[selected][other]))
            models[window] = (small_t, model, coefficient, peaks)
            if window == 'hann':
                reported = kin_by_package[package]['interfaces']['basal_sand']['hann']['complex_correlation_primary']
                assert abs(corr(a, b)-reported)<1e-10
        item = dict(package=package, deliverable_sha256=sha(derived), checks=checks, basal_diagnostics=metrics)
        if 'pkg3' in package:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            plt.rcParams['font.sans-serif']=['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
            plt.rcParams['axes.unicode_minus']=False
            fig, axes = plt.subplots(2, 3, figsize=(15, 8), layout='constrained')
            reference = float(np.max(abs(gaussian)))
            display_peak = 0.
            for window, actual in [('hann', actual_hann), ('blackman', blackman)]:
                tt, model, coefficient, peaks = models[window]
                use = (time>=tt[0])&(time<=tt[-1])
                for j in [0, len(rows)//2, len(rows)-1]:
                    local = abs(tt*1e9-peaks[j])<=20
                    display_peak = max(display_peak,
                        float(max(abs(actual[use,j][local]).max(), abs(coefficient*model[:,j][local]).max())/reference*1e4))
            for row, window, actual in [(0, 'hann', actual_hann), (1, 'blackman', blackman)]:
                tt, model, coefficient, peaks = models[window]
                use = (time>=tt[0])&(time<=tt[-1])
                for col, j in enumerate([0, len(rows)//2, len(rows)-1]):
                    ax = axes[row, col]
                    ax.plot(tt*1e9, abs(actual[use,j])/reference*1e4, label='实际总场幅度（无背景扣除）', lw=1.5)
                    ax.plot(tt*1e9, abs(coefficient*model[:,j])/reference*1e4, '--', label='局部一次反射近似×全局复比例', lw=1.2)
                    ax.axvline(peaks[j], color='#1ec766', lw=.8, label='材料谱预测峰位')
                    ax.set(xlim=(peaks[j]-20, peaks[j]+20), ylim=(0, display_peak*1.1), xlabel='时间 / ns',
                           ylabel='幅度 / 原Gaussian峰值 × 10⁴',
                           title=window+f'；里程{rows[j]["chainage_m"]:.2f}m')
                    ax.legend(fontsize=7);ax.grid(alpha=.15)
            fig.suptitle('1m包底砂参考窗：实际原始复谱与独立材料谱近似；两个频窗、三个固定站位\n近似模型每个窗仅拟合一个全局复比例，无逐道拟合/延时搜索；非独立消融认证')
            fig.savefig(out/'pkg3_basal_actual_and_primary.png', dpi=135)
            plt.close(fig)
        results.append(item)
        print(json.dumps(item, ensure_ascii=False), flush=True)
    save(out/'independent_checks.json', dict(calls_solver=False, calls_training=False,
        script_sha256=sha(__file__), products_report_sha256=sha(products_report),
        kinematics_sha256=sha(kinematics), packages=results,
        limits='Exact output consistency and geometry-informed waveform agreement only; not blind detection, numerical convergence, target ablation, true SNR or field validation. First-half fit validation holds out station coefficients, not geology/model development.'))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','review','products-report','kinematics','out']:
        p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    main(args.root,args.review,args.products_report,args.kinematics,args.out)
