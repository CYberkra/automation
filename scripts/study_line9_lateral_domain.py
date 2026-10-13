"""Authorized up-to-six 1/30m lateral-domain solves with a174m stop gate."""
from __future__ import annotations
import argparse
import copy
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from study_line9_v31_smoothing_domain import DB, load, save, sha, write_geometry
from analyze_line9_v31_smoothing_domain import read_sf, independent_check, compare, coefficients, relative, plot_sparse_compare


def prepare(a):
    root=a.root.resolve();root.mkdir(parents=True,exist_ok=False);geo=root/'geometries';geo.mkdir()
    prior=load(a.baseline/'contract.json')
    assert prior['grid_fraction_m']==[1,30] and prior['record_ns']==1100
    contacts_path=a.package/'model/contacts_m.json';selection_path=a.package/'config/selection_pdf_reconstructed.json'
    cs=load(contacts_path);sel=load(selection_path);contacts=cs['contact_points_x_elevation_m']
    assert cs['selection_sha256']==sha(selection_path)
    assert cs['geometry_model_version']=='V3.1_VERTICAL_PDF_SCALE_FIXED'
    shutil.copy2(a.inputs/'geometries'/f'{DB}.json',geo/f'{DB}.json')
    assert sha(geo/f'{DB}.json')==prior['material_database_sha256']
    cases=[];references=[];checks=[];dl=1/30
    for station in [174.,150.,162.]:
        base=next(v for v in prior['cases'] if v['station_m']==station)
        source=a.inputs/base['geometry'];assert sha(source)==base['geometry_sha256']
        with h5py.File(source) as h:old=h['data'][:];attrs=dict(h.attrs)
        ox,oy=base['origin_global_x_elevation_m']
        xs=ox+(np.arange(-600,2100)+.5)*dl;ys=oy+(np.arange(1500)+.5)*dl
        lines={k:np.interp(xs,np.array(v)[:,0],np.array(v)[:,1]) for k,v in contacts.items()}
        surface,clay,upper,lower_top,lower_bottom=[lines[k][:,None] for k in ['surface','clay_base','upper_mudstone_base','lower_mudstone_top','lower_mudstone_base']]
        lo,hi=sel['upper_mudstone_extent_m'];upper_on=((xs>=lo)&(xs<hi))[:,None];y=ys[None,:]
        data=np.where(y>=surface,0,np.where(y>=clay,1,np.where(upper_on&(y>=upper)&(y<clay),2,
            np.where((y>=lower_bottom)&(y<np.minimum(lower_top,clay)),2,3)))).astype(np.int16)[:,:,None]
        assert data.shape==(2700,1500,1) and np.array_equal(data[600:2100],old)
        for k,v in contacts.items():
            active_x=xs[(xs>=lo)&(xs<hi)] if k=='upper_mudstone_base' else xs
            assert active_x.min()>=np.array(v)[:,0].min() and active_x.max()<=np.array(v)[:,0].max(),k
        refdir=root/'references'/'cases'/base['id'];refdir.mkdir(parents=True)
        narrow=a.narrow/'cases'/base['id'];rec=load(a.baseline/'cases'/base['id']/'receipt.json')
        assert sha(narrow/'impulse.h5')==rec['raw_sha256']
        for filename in ['impulse.h5','sfcw.h5']:
            shutil.copy2(narrow/filename,refdir/filename)
        shutil.copy2(a.baseline/'cases'/base['id']/'receipt.json',refdir/'original_receipt.json')
        for filename in ['surface_yee.h5','surface_yee_materials.json']:
            shutil.copy2(a.baseline/'cases'/base['id']/filename,refdir/filename)
        shutil.copy2(source,geo/source.name)
        references.append(dict(case=base,raw_sha256=sha(refdir/'impulse.h5'),sfcw_sha256=sha(refdir/'sfcw.h5')))
        for width,start,stop in [(70,300,2400),(90,0,2700)]:
            case=copy.deepcopy(base);cid=f'width{width}_x{station:.1f}';shift=(width-50)/2
            subset=data[start:stop];origin=ox-shift;path=geo/f'{cid}.h5'
            write_geometry(path,subset,attrs,dl,origin)
            case.update(id=cid,stage='A' if station==174 else 'B',family='lateral_domain',
                geometry=path.relative_to(root).as_posix(),geometry_sha256=sha(path),domain_m=[width,50.],
                origin_global_x_elevation_m=[origin,oy],baseline_case=base['id'])
            for key in ['tx_m','rx_m']:
                case[key][0]+=shift
                assert np.allclose(np.array(case[key][:2])+[origin,oy],case[key.replace('_m','_global_x_elevation_m')],atol=1e-10,rtol=0)
                assert np.allclose(np.array(case[key][:2])/dl,np.round(np.array(case[key][:2])/dl),atol=1e-8,rtol=0)
            for bounds in case['surface_coefficient_export_bounds']:bounds[0]+=shift
            offset=round(shift/dl);assert np.array_equal(subset[offset:offset+1500],old)
            checks.append(dict(case=cid,retained50m_voxels_identical=True,actual_global_pose_identical=True,grid_phase_shift_cells=offset))
            cases.append(case)
    contract=copy.deepcopy(prior)
    contract['official_sfcw'].update(f_stop_hz=110_000_000,steps=301)
    contract.update(id='line9_lateral_domain_20261013_r1',schema='line9-lateral-domain/1',
        authorization='User explicitly approved at most6 lateral-domain solves,174m70/90 gate then150/162m',
        cases=cases,references=references,reference_checks=checks,reuse_references=[],max_new_solves=6,max_wall_seconds=600,
        parent_contract_sha256=sha(a.baseline/'contract.json'),
        geometry_source_sha256=sha(contacts_path),selection_sha256=sha(selection_path),
        comparison_scope='Fixed exact1/30m grid, retained voxel geometry and actual global TX/RX; domain extension includes lateral geology and PML relocation',
        gate='174m width70 vs width90 must pass original frozen screen; otherwise stop without stageB',
        current_model_modified=False,automatic_domain_adoption=False)
    save(root/'contract.json',contract)
    for filename in ['study_line9_v31_smoothing_domain.py','analyze_line9_v31_smoothing_domain.py']:
        shutil.copy2(Path(__file__).with_name(filename),root/filename)
    shutil.copy2(Path(__file__).with_name('study_line9_v31_smoothing_domain.py'),root/'study.py')
    shutil.copy2(__file__,root/'lateral.py')
    save(root/'preparation_audit.json',dict(status='PASS',reference_checks=checks,
        contact_source_sha256=sha(contacts_path),contract_sha256=sha(root/'contract.json')))
    save(root/'input_manifest.json',dict(files=[dict(path=f.relative_to(root).as_posix(),sha256=sha(f),bytes=f.stat().st_size)
        for f in sorted(root.rglob('*')) if f.is_file()]))
    print('Prepared6 candidates,3 reused50m references; overlap voxels and actual pose identical')


def audit_case(root,case):
    dest=root/'cases'/case['id'];rec=load(dest/'receipt.json')
    assert rec['status']=='PASS'
    for f in rec['files']:
        p=root/f['path'];assert sha(p)==f['sha256'] and p.stat().st_size==f['bytes']
    with h5py.File(dest/'impulse.h5') as h:
        y=h['rxs/rx1/Ez'];x=h['srcs/src1/excitation']
        assert y.dtype==np.dtype('float32') and y.shape==(case['expected_samples'],)
        assert np.isfinite(y[:]).all() and np.max(np.abs(y[:]))>0
        assert np.count_nonzero(x['samples'][:])==1 and np.max(x['samples'][:])==1
        assert np.isclose(y.attrs['SampleInterval'],case['dt_s'],rtol=1e-12,atol=0)
        assert x.attrs['TimeSampleOffset']==float(h.attrs['dt'])/2 and y.attrs['TimeSampleOffset']==0
        for group,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:
            assert np.allclose(h[group].attrs['Position'],case[key],atol=1e-8,rtol=0)
    log=(dest/'solver.log').read_text(encoding='utf-8')
    assert 'CUDA backend using Device 0' in log and 'Interface averaging is enabled' in log
    sf=read_sf(dest/'sfcw.h5',frequency_spec=(20e6,110e6,301))
    check=independent_check(dest/'impulse.h5',sf)
    co=coefficients(root,case['id']);ref=coefficients(root/'references',case['baseline_case'])
    assert np.array_equal(co['epsilon'],ref['epsilon']) and np.array_equal(co['sigma'],ref['sigma'])
    assert co['mixed_ez_edge_count']>0
    return sf,rec,dict(**check,retained_surface_coefficients_identical=True)


def run(a):
    root=a.root.resolve();c=load(root/'contract.json');started=time.monotonic()
    for f in load(root/'input_manifest.json')['files']:
        p=root/f['path'];assert sha(p)==f['sha256'] and p.stat().st_size==f['bytes']
    try:
        subprocess.run([sys.executable,str(root/'study.py'),'run-stage','--root',str(root),'--stage','A'],check=True,timeout=600)
        cases={v['id']:v for v in c['cases']}
        sf70,_,check70=audit_case(root,cases['width70_x174.0'])
        sf90,_,check90=audit_case(root,cases['width90_x174.0'])
        base=next(v['case'] for v in c['references'] if v['case']['station_m']==174)
        stats=compare(sf70,sf90,base,c['domain_screen'])
        proceed=stats['passes_frozen_domain_screen']
        save(root/'gate_174.json',dict(proceed=proceed,comparison=stats,independent_checks=[check70,check90],
            contract_sha256=sha(root/'contract.json')))
        print('174m gate: '+('PASS; proceeding with4 remaining solves' if proceed else 'FAILED; stopping after2 solves'),flush=True)
        if proceed:
            remaining=600-(time.monotonic()-started)
            assert remaining>0
            subprocess.run([sys.executable,str(root/'study.py'),'run-stage','--root',str(root),'--stage','B'],check=True,timeout=remaining)
        save(root/'batch_complete.json',dict(status='COMPLETE' if proceed else 'STOPPED_AT_FROZEN_GATE',
            new_solves=6 if proceed else 2,gate_pass=proceed,elapsed_seconds=time.monotonic()-started,
            contract_sha256=sha(root/'contract.json'),retries=0))
    except Exception as e:
        save(root/'batch_failure.json',dict(error_type=type(e).__name__,error=str(e),elapsed_seconds=time.monotonic()-started))
        raise


def analyze(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=a.root;out=a.out;out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    c=load(root/'contract.json');completion=load(root/'batch_complete.json')
    assert sha(root/'contract.json')==completion['contract_sha256']==load(root/'execution_start.json')['contract_sha256']
    assert sha(root/'study.py')==load(root/'execution_start.json')['script_sha256']
    for f in load(root/'input_manifest.json')['files']:
        p=root/f['path'];assert sha(p)==f['sha256'] and p.stat().st_size==f['bytes']
    results={};receipts={};checks={};refs={};rows=[]
    for ref in c['references']:
        case=ref['case'];dest=root/'references'/'cases'/case['id']
        assert sha(dest/'impulse.h5')==ref['raw_sha256'] and sha(dest/'sfcw.h5')==ref['sfcw_sha256']
        refs[case['station_m']]=read_sf(dest/'sfcw.h5',frequency_spec=(20e6,110e6,301))
        checks[case['id']]=independent_check(dest/'impulse.h5',refs[case['station_m']])
    completed_cases=[case for case in c['cases'] if (root/'cases'/case['id']/'receipt.json').exists()]
    assert len(completed_cases)==completion['new_solves']
    if not completion['gate_pass']:
        assert len(completed_cases)==2
        assert not any((root/'cases'/v['id']/'attempt.json').exists() for v in c['cases'] if v['stage']=='B')
    for case in completed_cases:
        sf,rec,check=audit_case(root,case);results[case['id']]=sf;receipts[case['id']]=rec;checks[case['id']]=check
    stations=sorted({v['station_m'] for v in completed_cases})
    for s in stations:
        base=next(v['case'] for v in c['references'] if v['case']['station_m']==s)
        sf70,sf90=results[f'width70_x{s:.1f}'],results[f'width90_x{s:.1f}']
        for label,new,ref in [('70_vs90',sf70,sf90),('50_vs90',refs[s],sf90),('50_vs70',refs[s],sf70)]:
            stats=compare(new,ref,base,c['domain_screen'])
            mask=(ref['t']>=100)&(ref['t']<400)
            rows.append(dict(station_m=s,comparison=label,**stats,
                magnitude_100_400_relative_l2=relative(np.abs(new['env'][mask]),np.abs(ref['env'][mask]))))
            if s==174 and label=='70_vs90':
                assert stats==load(root/'gate_174.json')['comparison']
    save(out/'comparisons.json',dict(batch=completion,cases=rows,thresholds=c['domain_screen'],
        automatic_adoption=False,reference='Finite90m/1100ns reference; absolute convergence not certified'))
    save(out/'audit.json',dict(status='PASS',new_solves=len(completed_cases),reused_references=3,checks=checks,
        verified_receipt_files=sum(len(v['files']) for v in receipts.values()),input_manifest_verified=True,
        contract_sha256=sha(root/'contract.json'),analysis_script_sha256=sha(__file__),
        total_fdtd_seconds=sum(v['official_fdtd_seconds'] for v in receipts.values()),
        total_build_solver_seconds=sum(v['official_build_and_solver_seconds'] for v in receipts.values())))
    for width in [50,70]:
        reference=np.array([results[f'width90_x{s:.1f}']['real'] for s in stations]).T
        candidate=np.array([(refs[s] if width==50 else results[f'width70_x{s:.1f}'])['real'] for s in stations]).T
        plot_sparse_compare(out,f'width{width}_vs90_sparse_bscan.png',stations,reference,candidate,sf90['t'],
            f'90m基准 → {width}m候选｜精确1/30m／20–110MHz／同一实际收发坐标｜仅完成站位',
            reference_events={s:next(v['case']['vertical_reference_events'] for v in c['references'] if v['case']['station_m']==s) for s in stations})
    fig,axes=plt.subplots(2,len(stations),figsize=(max(8.5,5*len(stations)),7),squeeze=False,constrained_layout=True)
    for j,s in enumerate(stations):
        for sf,label in [(refs[s],'50m旧基准'),(results[f'width70_x{s:.1f}'],'70m'),(results[f'width90_x{s:.1f}'],'90m')]:
            mask=(sf['t']>=100)&(sf['t']<400)
            axes[0,j].plot(sf['t'][mask],sf['real'][mask],lw=1,label=label)
            axes[1,j].plot(sf['t'][mask],np.abs(sf['env'][mask]),lw=1,label=label)
        for ax in axes[:,j]:ax.set(title=f'{s:g}m｜同一1/30m网格',xlabel='时间（ns）');ax.grid(alpha=.2);ax.legend()
        axes[0,j].set_ylabel('带符号总场响应');axes[1,j].set_ylabel('复包络模')
    fig.suptitle('20–110MHz／301点｜原生平滑y\n无增益、背景抑制、时间对齐或幅度拟合')
    fig.savefig(out/'ascan_width_comparison.png',dpi=170);plt.close(fig)
    q=results['width70_x174.0'];r=results['width90_x174.0'];m=(r['t']>=100)&(r['t']<400)
    delta=np.abs(q['env']-r['env']);t=r['t'][m]
    save(out/'residual_summary.json',dict(station_m=174,
        difference='abs(complex_envelope70-complex_envelope90)',
        peak_time_within100_400_ns=float(t[np.argmax(delta[m])]),peak_amplitude=float(delta[m].max()),
        deep308_338_window_relative_l2=next(v['complex_relative_l2'] for v in rows[0]['windows'] if v['lo_ns']>300),
        attribution='Residual timing alone cannot separate PML error from omitted lateral geology'))
    fig,ax=plt.subplots(figsize=(11,4.5),constrained_layout=True)
    ax.plot(t,delta[m],label='|70m复包络－90m复包络|')
    ax.plot(t,np.abs(refs[174.]['env']-r['env'])[m],label='|50m复包络－90m复包络|')
    for width in [70,90]:
        ax.axvline((width-6)/299792458*1e9,ls='--',lw=.9,label=f'{width}m侧面空气往返几何尺度（非回波证据）')
    ax.set(xlabel='时间（ns）',ylabel='复包络差的模；独立纵轴、未加增益',
        title='174m｜20–110MHz／精确1/30m｜域差量的时间分布；不自动称为PML反射')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(out/'domain_residual_timing.png',dpi=170);plt.close(fig)
    print(completion)
    for row in rows:print(row['station_m'],row['comparison'],row['passes_frozen_domain_screen'],max(v['complex_relative_l2'] for v in row['windows']))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','run','analyze']);p.add_argument('--root',type=Path,required=True)
    for arg in ['baseline','inputs','narrow','package','out']:p.add_argument('--'+arg,type=Path)
    a=p.parse_args();{'prepare':prepare,'run':run,'analyze':analyze}[a.action](a)
