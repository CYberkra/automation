"""Compare matched joint-grid transfer functions, preserving complex phase."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

from analyze_hs4t2d_cause_controls import load_response, weighted_comparison
from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import verify_official_runtime, OFFICIAL_PROCESSING_SHA256

ROOT = Path(__file__).resolve().parents[1]


def compare(a, b):
    norms = np.linalg.norm(a, axis=0)
    if np.any(norms == 0):
        raise ValueError('zero reference; relative comparison undefined')
    return {'relative_L2': float(np.linalg.norm(b-a)/np.linalg.norm(a)),
            'relative_L2_by_trace': (np.linalg.norm(b-a, axis=0)/norms).tolist(),
            'reference_L2': float(np.linalg.norm(a)), 'difference_L2': float(np.linalg.norm(b-a)),
            'reference_peak_abs': float(np.max(np.abs(a))),
            'difference_peak_abs': float(np.max(np.abs(b-a)))}


def completed_study(path, stage):
    c = json.loads((path/'execution_contract.json').read_text('utf-8'))
    check = json.loads((path/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (path/'execution.jsonl').read_text('utf-8').splitlines()]
    digest = sha256(path/'execution_contract.json')
    if c['stage'] != stage or check['status'] != 'PASS' or check['contract_sha256'] != digest or events[0]['contract_sha256'] != digest or events[-1].get('status') != 'COMPLETED' or events[-1].get('traces') != c['max_runs']:
        raise ValueError('complete independently verified frozen study required')
    for g in c['groups']:
        if sha256(path/g['id']/'profile.in') != g['sha256']:
            raise ValueError('input identity differs')
    return c


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--centre', type=Path, required=True)
    p.add_argument('--remaining', type=Path)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise ValueError('new result directory required')
    verify_official_runtime()
    pilot = completed_study(args.centre, 'centre')
    studies = [(args.centre, pilot)]
    if args.remaining:
        full = completed_study(args.remaining, 'remaining')
        if full['prerequisite']['contract_sha256'] != sha256(args.centre/'execution_contract.json'):
            raise ValueError('wrong pilot for remaining stations')
        studies.append((args.remaining, full))
    base = Path(pilot['baseline_directory'])
    if sha256(base/'execution_contract.json') != pilot['baseline_contract_sha256']:
        raise ValueError('coarse reference identity differs')
    responses, audit = {}, {}
    source_info = {}
    for role in ('rough', 'halfspace'):
        r, rows = load_response(base/('x40wide_'+role))
        if len(rows) != 13 or not r.source_valid.all():
            raise ValueError('complete matched coarse reference required')
        responses[role], audit[role] = r, rows
    coarse = responses['rough']
    fine = {'rough': {}, 'halfspace': {}}
    for directory, contract in studies:
        for g in contract['groups']:
            r, rows = load_response(directory/g['id'])
            if not r.source_valid.all() or not np.array_equal(r.frequency, coarse.frequency) or r.source.spatial_scale != coarse.source.spatial_scale or r.source.quantity != coarse.source.quantity or r.source.units != coarse.source.units:
                raise ValueError('incompatible source/grid normalization')
            role = 'halfspace' if g['halfspace'] else 'rough'
            source_info[g['id']] = {'source_quantity': r.source.quantity, 'source_units': r.source.units,
                'spatial_scale': r.source.spatial_scale, 'dt_s': r.source.dt,
                'source_all501bins_valid': True, 'audit_sha256': sha256(directory/g['id']/'audit.json')}
            for k, row in enumerate(rows):
                i = row['old_station_index']
                if i in fine[role]:
                    raise ValueError('duplicate station')
                old = audit[role][(i-1)//10]
                if not np.allclose(row['source_m'], old['source_m'], rtol=0, atol=1e-12) or not np.allclose(row['receiver_m'], old['receiver_m'], rtol=0, atol=1e-12):
                    raise ValueError('actual station differs')
                fine[role][i] = r.response[:, k]
    indices = sorted(fine['rough'])
    expected = list(range(1, 122, 10)) if args.remaining else [61]
    if indices != expected or sorted(fine['halfspace']) != expected:
        raise ValueError('incomplete paired station coverage')
    take = [(i-1)//10 for i in indices]
    spectra = {'coarse': (responses['rough'].response-responses['halfspace'].response)[:, take],
               'fine': np.stack([fine['rough'][i]-fine['halfspace'][i] for i in indices], axis=1)}
    products = {k: reconstruct_time_response(replace(coarse, response=v), window='hann', zero_pad_factor=8)
                for k,v in spectra.items()}
    time = products['coarse'].time*1e9
    if not np.array_equal(time, products['fine'].time*1e9):
        raise ValueError('reconstruction axes differ')
    windows = {'ground': [85,120], 'underground': [160,220]}
    metrics = {'full_hann_complex_spectrum': weighted_comparison(spectra['coarse'], spectra['fine'])}
    peaks = {}
    for label, (lo,hi) in windows.items():
        mask = (time>=lo)&(time<=hi)
        metrics[label] = {'signed_band': compare(products['coarse'].real_bandpass[mask], products['fine'].real_bandpass[mask]),
                         'phase_aware_complex_envelope': compare(products['coarse'].complex_envelope[mask], products['fine'].complex_envelope[mask]),
                         'envelope_magnitude': compare(2*np.abs(products['coarse'].complex_envelope[mask]), 2*np.abs(products['fine'].complex_envelope[mask]))}
        if label == 'underground':
            for name, product in products.items():
                times = time[mask][np.argmax(np.abs(product.complex_envelope[mask]), axis=0)]
                peaks[name] = {'times_ns': times.tolist(), 'span_ns': float(np.ptp(times)),
                               'scope': 'fixed-window envelope maxima may switch lobes; not certified interface arrival/path identity'}
    args.out.mkdir(parents=True)
    x = np.array([(audit['rough'][k]['source_m'][0]+audit['rough'][k]['receiver_m'][0])/2-12 for k in take])
    arrays = {'time_ns': time, 'midpoint_m': x, 'old_station_indices': indices, 'frequency_Hz': coarse.frequency}
    for name, product in products.items():
        arrays[name+'_frequency_complex'] = spectra[name]
        arrays[name+'_signed'] = product.real_bandpass
        arrays[name+'_complex_envelope'] = product.complex_envelope
    np.savez_compressed(args.out/'joint_grid_arrays.npz', **arrays)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Microsoft YaHei'
    mask = (time>=160)&(time<=220)
    centre_k = indices.index(61)
    fig, axes = plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for name, product in products.items():
        label = 'X/Z 5cm' if name=='coarse' else 'X/Z 2.5cm'
        axes[0].plot(time[mask],product.real_bandpass[mask,centre_k],label=label)
        axes[1].plot(time[mask],2*np.abs(product.complex_envelope[mask,centre_k]),label=label)
    for ax in axes:
        ax.legend(); ax.set_xlim(160,220); ax.set_xlabel('时间（ns）'); ax.set_ylabel('场/电流响应')
    axes[0].set_title('中心道带符号响应（无归一化、无时间平移）')
    axes[1].set_title('中心道复包络模值（不替代相位比较）')
    fig.suptitle('同一15m航高 / 36m域 / 0.8m全剖面起伏 / 相同物理PML厚度\n匹配全覆盖层复谱差分 / 20–170MHz / Hann / V4.0.0 CUDA double')
    fig.savefig(args.out/'joint_grid_center.png',dpi=145); plt.close(fig)
    if len(indices)>1:
        a,b = products['coarse'].real_bandpass[mask],products['fine'].real_bandpass[mask]
        limit = float(np.max(np.abs(np.stack([a,b]))))
        fig,axes = plt.subplots(2,3,figsize=(14,8),constrained_layout=True,gridspec_kw={'height_ratios':[1,2]})
        profile = np.genfromtxt(ROOT/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',delimiter=',',names=True)
        for j,(matrix,title) in enumerate(zip([a,b,b-a],['X/Z 5cm','X/Z 2.5cm','细网格减粗网格（同色标）'])):
            ax=axes[0,j]
            ax.stairs(12-profile['grid_z_m'],np.r_[profile['x0_m'],profile['x1_m'][-1]],baseline=None)
            ax.set_xlim(x[0],x[-1]); ax.set_ylim(3.4,2.4); ax.set_ylabel('距地表深度（m）'); ax.set_title('相同实际台阶界面')
            ax=axes[1,j]
            im=ax.pcolormesh(x,time[mask],matrix,cmap='gray',vmin=-limit,vmax=limit,shading='nearest')
            ax.set_ylim(220,160); ax.set_xlabel('原模型收发中点X（m）'); ax.set_ylabel('时间（ns）'); ax.set_title(title)
            fig.colorbar(im,ax=ax,label='场/电流响应，无归一化')
        fig.suptitle('15m航高 / 36m域 / 左右2m、上下1m PML / 13匹配站位\n全覆盖层复谱差分 / 无SVD、AGC、时间对齐或成像 / 20–170MHz Hann')
        fig.savefig(args.out/'joint_grid_bscan.png',dpi=145); plt.close(fig)
    report = {'status': 'COMPLETED', 'station_count': len(indices), 'new_traces': sum(c['max_runs'] for _,c in studies),
              'studies': [{'directory': str(d.resolve()),'contract_sha256': sha256(d/'execution_contract.json'),
                          'completed_verification_sha256': sha256(d/'completed_verification.json')} for d,c in studies],
              'official_processing_sha256': OFFICIAL_PROCESSING_SHA256, 'code_sha256': sha256(__file__),
              'baseline_contract_sha256': pilot['baseline_contract_sha256'], 'source_audit': source_info,
              'fixed_windows_ns': windows, 'comparisons': metrics, 'window_peak_diagnostic': peaks,
              'scope': 'Joint X/Z grid sensitivity at fixed15m geometry. Two resolutions do not prove complete convergence or a unique physical cause; maxima are not event labels.'}
    (args.out/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (args.out/'manifest.json').write_text(json.dumps([{'file':f.name,'bytes':f.stat().st_size,'sha256':sha256(f)}
        for f in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'station_count':len(indices),'comparisons':metrics,'window_peak_diagnostic':peaks},indent=2))


if __name__ == '__main__':
    main()
