"""Two bounded X220 domain-width controls; preserve full-domain reference and complex phase."""
import argparse
import copy
import json
import shutil
from pathlib import Path

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json
import run_line9_2d_first as base
from analyze_line9_2d_sfcw import process_trace, FREQUENCY, PROCESSING_SHA

ROOT = Path(__file__).resolve().parents[1]


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    base.check_files(c)
    if digest(Path(c['reference_raw'])) != c['reference_raw_sha256']:
        raise ValueError('Full-domain reference changed')
    rows = []
    with h5py.File(c['reference_raw']) as ref:
        for g in c['groups']:
            p = Path(g['input'])
            if digest(p) != g['input_sha256']:
                raise ValueError('Frozen crop input changed')
            row = dict(id=g['id'], width_m=g['width_m'], input_sha256=digest(p))
            if completed:
                raw = p.with_suffix('.h5')
                with h5py.File(raw) as h:
                    v = h['rxs/rx1/Ez'][:]
                    if str(h.attrs['gprMax']) != '4.0.0' or v.dtype != np.float64 or not np.isfinite(v).all():
                        raise ValueError('Native crop runtime/precision mismatch')
                    np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025]*3)
                    np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [round(g['width_m']/.025), 3000, 1])
                    for k in ['dt','Iterations']:
                        np.testing.assert_array_equal(h.attrs[k], ref.attrs[k])
                    if v.shape != ref['rxs/rx1/Ez'].shape:
                        raise ValueError('Native history length changed')
                    source = h['srcs/src1/excitation/samples']
                    if source.dtype != np.float64:
                        raise ValueError('Original source precision changed')
                    np.testing.assert_array_equal(source[:], ref['srcs/src1/excitation/samples'][:])
                    for k in ['SpatialScale','Polarisation','WaveformAmplitude','WaveformFrequency']:
                        np.testing.assert_array_equal(h['srcs/src1/excitation'].attrs[k], ref['srcs/src1/excitation'].attrs[k])
                    for key, position in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                        actual = h[key].attrs['GridPosition']
                        np.testing.assert_array_equal(actual[:2], np.rint(np.asarray(position[:2])/.025).astype(int))
                        if actual[2] not in (0,1):
                            raise ValueError('Invariant-axis placement changed')
                row.update(raw_path=str(raw.resolve()), raw_sha256=digest(raw), dtype='float64')
            rows.append(row)
    return dict(status='PASS_NATIVE_IDENTITIES_NOT_DOMAIN_ACCEPTANCE', completed=completed,
                contract_sha256=digest(path), groups=rows)


def freeze(package, out, pilots):
    if out.exists():
        raise ValueError('Fresh crop attempt required')
    old_path = pilots/'execution_contract.json'
    old = json.loads(old_path.read_text('utf-8'))
    base.check_files(old)
    group = next(g for g in old['groups'] if g['id']=='full2d_pilot_0')
    events = [json.loads(v) for v in (pilots/'execution.jsonl').read_text('utf-8').splitlines()]
    completion = [v for v in events if v.get('status')=='COMPLETED' and v.get('group')==group['id']]
    raw = Path(group['input']).with_suffix('.h5')
    if len(completion)!=1 or digest(raw)!=completion[0]['raw_sha256']:
        raise ValueError('Require recorded completed X220 full-domain reference')
    hardware = base.resources(json.loads(base.TASK.read_text('utf-8')))
    if hardware['gpu_uuid'] != old['hardware_at_freeze']['gpu_uuid'] or hardware['driver_version'] != old['hardware_at_freeze']['driver_version']:
        raise ValueError('Reference GPU/driver changed')
    m = json.loads((package/'manifest.json').read_text('utf-8'))
    cases = [next(r for r in m['cases'] if r['id']==f'w{w}_full2d_pilot_0') for w in [200,160]]
    db = package/'geometries/line9_research_materials_v1.json'
    if digest(db) != m['material_database_sha256']:
        raise ValueError('Crop material database changed')
    # Recheck every retained voxel against the original frozen full H5.
    with h5py.File(pilots/'geometries/full2d.h5') as full:
        for r in cases:
            if digest(package/r['input'])!=r['input_sha256'] or digest(package/r['geometry'])!=r['geometry_sha256']:
                raise ValueError('Prepared crop changed')
            if r['source_input_sha256'] != group['input_sha256']:
                raise ValueError('Wrong source station')
            def commands(file):
                return {s.split(':',1)[0]:s.split(':',1)[1].split() for s in file.read_text('utf-8').splitlines() if s.startswith('#')}
            before, after = commands(Path(group['input'])), commands(package/r['input'])
            allowed={'#domain','#hertzian_dipole','#rx','#geometry_objects_read'}
            if set(before)!=set(after) or any(before[k]!=after[k] for k in before.keys()-allowed):
                raise ValueError('Crop changes more than horizontal extent/translation')
            np.testing.assert_array_equal(np.asarray(after['#domain'],dtype=float),[r['width_m'],75,.025])
            for k,command,begin in [('tx_m','#hertzian_dipole',1),('rx_m','#rx',0)]:
                np.testing.assert_allclose(np.asarray(after[command][begin:begin+3],dtype=float),r[k],rtol=0,atol=1e-9)
                restored=np.asarray(r[k])+[r['profile_x_range_m'][0]+50,0,0]
                np.testing.assert_allclose(restored,group[k],rtol=0,atol=1e-9)
            with h5py.File(package/r['geometry']) as crop:
                first = round((r['profile_x_range_m'][0]+50)/.025)
                np.testing.assert_array_equal(crop['material_keys'][:],full['material_keys'][:])
                for j in range(0,crop['data'].shape[0],64):
                    stop=min(j+64,crop['data'].shape[0])
                    np.testing.assert_array_equal(crop['data'][j:stop], full['data'][first+j:first+stop])
                np.testing.assert_array_equal(crop.attrs['dx_dy_dz'],[.025]*3)
    out.mkdir(parents=True)
    (out/'geometries').mkdir()
    shutil.copyfile(db,out/'geometries'/db.name)
    groups=[]
    for r in cases:
        dst=out/r['input'];dst.parent.mkdir(parents=True)
        shutil.copyfile(package/r['input'],dst)
        shutil.copyfile(package/r['geometry'],out/r['geometry'])
        g=copy.deepcopy(group)
        g.update(r, input=str(dst.resolve()),input_sha256=digest(dst))
        groups.append(g)
    view=copy.deepcopy(old)
    view.update(status='READONLY_COMPLETED_GROUP_AUDIT',groups=[group],max_runs=1)
    save_json(out/'reference_audit_view.json',view)
    save_json(out/'reference_verification.json',base.audit(out/'reference_audit_view.json',True))
    c=copy.deepcopy(old)
    sources=['line9_crop_controls.py','run_line9_2d_first.py','hs4_station_grid_controls.py',
             'gprmax_cached_cuda_entry.py','hs_capsule_identity.py','build_pdf_profile_geometry.py',
             'analyze_line9_2d_sfcw.py','run_line9_2d_v4.cmd']
    c.update(stage='crop_controls',approval_basis='User requests domain acceleration preserving quality; prior autonomous simulation authorization; two X220 width controls only, no full line',
        groups=groups,max_runs=2,reuse={},hardware_at_freeze=hardware,
        reference_raw=str(raw.resolve()),reference_raw_sha256=digest(raw),
        reference_completion_event=completion[0],reference_contract_sha256=digest(old_path),
        candidate_manifest_sha256=digest(package/'manifest.json'),
        code_identities={str(ROOT/'scripts'/n):digest(ROOT/'scripts'/n) for n in sources},
        file_identities={str(p.resolve()):digest(p) for p in out.rglob('*') if p.is_file()},
        max_group_wall_s=1200,max_batch_wall_s=2400,
        quality_threshold=None,quality_acceptance='Not declared; report raw, windowed weak returns and501 complex tones, never auto-promote whole line')
    save_json(out/'execution_contract.json',c)
    save_json(out/'preflight_verification.json',audit(out/'execution_contract.json'))


def analyze(study, out):
    from dataclasses import replace
    from gprMax.toolboxes.SFCW import processing as sf
    if out.exists() or digest(Path(sf.__file__))!=PROCESSING_SHA:
        raise ValueError('Fresh report and reviewed processing required')
    current=audit(study/'execution_contract.json',True)
    if current!=json.loads((study/'completed_verification.json').read_text('utf-8')):
        raise ValueError('Completed crop verification changed')
    c=json.loads((study/'execution_contract.json').read_text('utf-8'))
    paths=[Path(c['reference_raw'])]+[Path(g['input']).with_suffix('.h5') for g in c['groups']]
    names=['完整域400m','精确切片200m','精确切片160m']
    products=[];histories=[];signed=[];records=[]
    for path in paths:
        p,r,record=process_trace(path)
        result=sf.reconstruct_time_response(replace(p,response=p.response),window='hann',zero_pad_factor=8,normalise_window=True,time_shift=0)
        products.append(p.response);histories.append(r.samples);signed.append(result.real_bandpass);records.append(record)
    windows=[]
    for lo,hi in [(0,100),(100,300),(300,600),(600,1200)]:
        row=dict(window_ns=[lo,hi],comparisons=[])
        for i in [1,2]:
            metrics={}
            for key,t,values in [('native',r.times,histories),('Hann_SFCW',result.time,signed)]:
                mask=(t>=lo*1e-9)&(t<hi*1e-9)
                ref=values[0][mask];delta=values[i][mask]-ref
                metrics[key]=dict(relative_L2=float(np.linalg.norm(delta)/np.linalg.norm(ref)) if np.linalg.norm(ref)>0 else None,
                    residual_peak=float(np.max(abs(delta))),reference_peak=float(np.max(abs(ref))),
                    residual_peak_over_local_reference_peak=float(np.max(abs(delta))/np.max(abs(ref))) if np.max(abs(ref))>0 else None)
            row['comparisons'].append(dict(width_m=c['groups'][i-1]['width_m'],**metrics))
        windows.append(row)
    out.mkdir(parents=True)
    np.savez_compressed(out/'crop_comparison_arrays.npz',native_time_s=r.times,native_Ez_V_m=np.stack(histories,axis=1),
        frequency_Hz=FREQUENCY,response_complex=np.stack(products,axis=1),sfcw_time_s=result.time,Hann_signed=np.stack(signed,axis=1))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(2,3,figsize=(13,9),layout='constrained')
    for row,(t,values,label) in enumerate([(r.times,histories,'原始Ricker总场Ez / V/m'),(result.time,signed,'SFCW Hann总场；源归一化')]):
        mask=t<=1.2e-6;scale=max(float(np.max(abs(v[mask]))) for v in values)
        for j in range(3):
            axes[row,j].imshow(values[j][mask,None],extent=[-.25,.25,t[mask][-1]*1e9,t[mask][0]*1e9],aspect='auto',cmap='gray',vmin=-scale,vmax=scale)
            axes[row,j].set(title=names[j]+'；'+label,ylabel='时间 / ns',xlabel='单站s=0（X220）；无相邻道',xticks=[0])
            axes[row,j].axhline(100.16,color='red',ls='--',lw=.7)
    fig.suptitle('同一X220单道B-scan对照；15m离地；2.5cm/FP64/材料/1200ns全保留\n每行共享色标；未去背景、未增益；仅显示已计算站位，不能代表整条测线等价')
    fig.savefig(out/'crop_single_station_bscan.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    for row,(t,values,label) in enumerate([(r.times,histories,'原始总场Ez / V/m'),(result.time,signed,'SFCW Hann带符号总场')]):
        for col,(lo,hi) in enumerate([(100,600),(600,1200)]):
            mask=(t>=lo*1e-9)&(t<hi*1e-9)
            for j in range(3):
                axes[row,col].plot(t[mask]*1e9,values[j][mask],label=names[j],lw=.8)
            axes[row,col].set(title=label+f'；{lo}–{hi}ns',xlabel='时间 / ns',ylabel='物理幅值（无逐道归一化）')
            axes[row,col].legend(fontsize=8)
    fig.suptitle('缩域弱回波检查；保留完整域实幅，不由直耦波全局L2决定质量')
    fig.savefig(out/'crop_weak_return_comparison.png',dpi=140);plt.close(fig)
    save_json(out/'comparison_report.json',dict(status='COMPLETED_X220_CONTROL_NOT_WHOLE_LINE_ACCEPTANCE',
        calls_solver=False,contract_sha256=digest(study/'execution_contract.json'),records=records,windows=windows,
        complex501_relative_L2=[float(np.linalg.norm(v-products[0])/np.linalg.norm(products[0])) for v in products[1:]],
        arrays_sha256=digest(out/'crop_comparison_arrays.npz'),
        limitations='One X220 station only. Threshold not declared; late/weak windows separately reported. No automatic 99/391-trace promotion. Shared scales and signed complex phase retained.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['freeze','run','analyze'])
    p.add_argument('--package',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--pilots',type=Path)
    p.add_argument('--report',type=Path)
    a=p.parse_args()
    if a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.pilots.resolve())
    elif a.action=='analyze':analyze(a.out.resolve(),a.report.resolve())
    else:
        import hs4_station_grid_controls as supervisor
        c=json.loads((a.out/'execution_contract.json').read_text('utf-8'))
        hw=base.resources(json.loads(base.TASK.read_text('utf-8')))
        if any(hw[k]!=c['hardware_at_freeze'][k] for k in ['gpu_uuid','driver_version']):
            raise ValueError('Frozen crop GPU/driver changed')
        supervisor.audit=audit
        supervisor.run(a.out.resolve()/'execution_contract.json')
