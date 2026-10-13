"""Reprocess existing FDTD only; research lateral domains in20–110MHz."""
from __future__ import annotations
import argparse
import copy
from pathlib import Path
import shutil
import sys
import time
import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_line9_v31_smoothing_domain import load, save, sha, process_sfcw


def prepare(a):
    a.root.mkdir(parents=True, exist_ok=False)
    specs = []
    for parent, fraction in [(a.domain_root, False), (a.fraction_root, True)]:
        contract = load(parent/'contract.json')
        for case in contract['cases']:
            keep = fraction or (case['averaging']=='y' and not case['geometry_only'] and (
                (case['id'].startswith('full_x') and case['domain_m']==[50.,50.]) or
                (case['family'] in ['expanded_domain','reduced_domain'] and case['domain_m'][1]==50.)))
            if not keep:
                continue
            src = parent/'cases'/case['id']
            old = load(src/'receipt.json')
            assert sha(src/'impulse.h5')==old['raw_sha256']
            assert sha(src/'sfcw.h5')==old['sfcw_sha256']
            dest = a.root/'cases'/case['id']
            dest.mkdir(parents=True)
            for filename, newname in [('impulse.h5','impulse.h5'),('sfcw.h5','sfcw_20_170.h5'),('receipt.json','original_receipt.json')]:
                shutil.copy2(src/filename,dest/newname)
            specs.append(dict(case=case, group='selected_grid_baseline' if fraction else 'historical_domain_screen',
                raw_sha256=old['raw_sha256'], wide_sfcw_sha256=old['sfcw_sha256'],
                source_contract_sha256=sha(parent/'contract.json')))
    assert len(specs)==20
    spec=copy.deepcopy(contract['official_sfcw']);spec.update(f_stop_hz=110_000_000,steps=301)
    save(a.root/'contract.json',dict(id='line9_lateral_domain_research_20261013_r1',
        authorization='User selected20–110MHz and1/30m, requested lateral-domain research; existing FDTD reused only',
        new_fdtd=0, official_sfcw=spec, cases=specs,
        domain_screen=contract['domain_screen'],
        processing_source_sha256='adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b',
        limits='2.5cm domain results inform planning only; do not certify1/30m domain convergence',
        max_processing_wall_seconds=300))
    shutil.copy2(__file__,a.root/'research.py')
    shutil.copy2(Path(__file__).with_name('study_line9_v31_smoothing_domain.py'),a.root/'study_line9_v31_smoothing_domain.py')
    save(a.root/'input_manifest.json',dict(files=[dict(path=f.relative_to(a.root).as_posix(),sha256=sha(f),bytes=f.stat().st_size)
        for f in sorted(a.root.rglob('*')) if f.is_file()]))
    print('Prepared20 existing records; no FDTD entrypoint in this research script')


def process(a):
    import gprMax
    from gprMax.toolboxes.SFCW import processing
    assert gprMax.__version__=='4.0.1'
    c=load(a.root/'contract.json')
    assert sha(processing.__file__)==c['processing_source_sha256']
    started=time.monotonic();receipts=[]
    for definition in c['cases']:
        assert time.monotonic()-started<c['max_processing_wall_seconds']
        cid=definition['case']['id'];dest=a.root/'cases'/cid
        assert sha(dest/'impulse.h5')==definition['raw_sha256']
        assert sha(dest/'sfcw_20_170.h5')==definition['wide_sfcw_sha256']
        process_sfcw(dest/'impulse.h5',dest/'sfcw.h5',dest/'sfcw.log',c['official_sfcw'])
        receipts.append(dict(case=cid,files=[dict(path=f.relative_to(a.root).as_posix(),sha256=sha(f),bytes=f.stat().st_size)
            for f in sorted(dest.glob('*')) if f.is_file()]))
        print(f'{len(receipts)}/20 official SFCW: {cid}',flush=True)
    save(a.root/'completion.json',dict(status='COMPLETE',new_fdtd=0,reprocessed_existing=len(receipts),
        elapsed_seconds=time.monotonic()-started,receipts=receipts,contract_sha256=sha(a.root/'contract.json'),
        processing_source_sha256=sha(processing.__file__),research_script_sha256=sha(__file__)))


def analyze(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze_line9_v31_smoothing_domain import read_sf, independent_check, compare, relative, plot_sparse_compare
    a.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    c=load(a.root/'contract.json');completion=load(a.root/'completion.json')
    assert completion['status']=='COMPLETE' and completion['new_fdtd']==0 and completion['reprocessed_existing']==20
    assert sha(a.root/'contract.json')==completion['contract_sha256']
    repair=load(a.root/'processing_repair.json')
    assert sha(a.root/'research_r2.py')==repair['corrected_worker_sha256']==completion['research_script_sha256']
    assert sha(a.root/'research.py')==repair['old_worker_sha256']
    hashed=0
    for r in completion['receipts']:
        for f in r['files']:
            p=a.root/f['path'];assert sha(p)==f['sha256'] and p.stat().st_size==f['bytes'];hashed+=1
    for f in load(a.root/'input_manifest.json')['files']:
        p=a.root/f['path'];assert sha(p)==f['sha256'] and p.stat().st_size==f['bytes']
    results={};checks={}
    for definition in c['cases']:
        case=definition['case'];cid=case['id'];dest=a.root/'cases'/cid
        with h5py.File(dest/'impulse.h5') as h:
            y=h['rxs/rx1/Ez'];x=h['srcs/src1/excitation']
            assert y.dtype==np.dtype('float32') and y.shape==(case['expected_samples'],)
            assert np.isfinite(y[:]).all() and np.count_nonzero(x['samples'][:])==1
            assert np.isclose(y.attrs['SampleInterval'],case['dt_s'],rtol=1e-12,atol=0)
            assert x.attrs['TimeSampleOffset']==float(h.attrs['dt'])/2
            for group,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:
                assert np.allclose(h[group].attrs['Position'],case[key],atol=1e-8,rtol=0)
        sf=read_sf(dest/'sfcw.h5',frequency_spec=(20e6,110e6,301));results[cid]=sf
        checks[cid]=independent_check(dest/'impulse.h5',sf)
        wide=read_sf(dest/'sfcw_20_170.h5')
        assert relative(sf['H'],wide['H'][:301])<1e-10
    comparisons=[]
    cases={d['case']['id']:d['case'] for d in c['cases']}
    path_rows=[]
    for case in cases.values():
        if not case['id'].startswith('grid1over30'):
            continue
        tx=np.array(case['tx_global_x_elevation_m']);rx=np.array(case['rx_global_x_elevation_m'])
        for width in [40,50,70,90]:
            origin=case['origin_global_x_elevation_m'][0]-(width-50)/2
            left,right=origin+3,origin+width-3
            lengths=[np.hypot(tx[0]+rx[0]-2*left,tx[1]-rx[1]),np.hypot(2*right-tx[0]-rx[0],tx[1]-rx[1])]
            path_rows.append(dict(station_m=case['station_m'],width_m=width,
                origin_x_m=origin,left_inner_pml_x_m=left,right_inner_pml_x_m=right,
                left_air_mirror_time_ns=float(lengths[0]/299792458*1e9),
                right_air_mirror_time_ns=float(lengths[1]/299792458*1e9)))
    save(a.out/'boundary_paths.json',dict(cases=path_rows,
        centered_approximation_formula='t=(W-6)/c; potential boundary interaction, not measured reflection',
        conservative_width_for_air_return_after400ns_m=299792458*400e-9+6,
        cells_per110MHz_wavelength_in_er18_at1over30=299792458/(110e6*np.sqrt(18)*(1/30)),
        proposed_new_fdtd=6,proposed_widths_m=[70,90],proposed_stations_m=[150,162,174],
        execution_authorized=False))
    for cid,case in cases.items():
        if not cid.startswith(('crop40x50','wide70')):
            continue
        station=case['station_m'];baseid=f'full_x{station:.1f}_y'
        base=cases[baseid]
        for key in ['tx_global_x_elevation_m','rx_global_x_elevation_m']:
            assert np.allclose(case[key],base[key],atol=1e-10,rtol=0)
        referenceid,candidateid=(cid,baseid) if cid.startswith('wide70') else (baseid,cid)
        ref,new=results[referenceid],results[candidateid]
        stats=compare(new,ref,base,c['domain_screen'])
        masks=[(lo,hi,(ref['t']>=lo)&(ref['t']<hi)) for lo,hi in c['domain_screen']['common_windows_ns']]
        comparisons.append(dict(station_m=station,reference=referenceid,candidate=candidateid,
            reference_width_m=cases[referenceid]['domain_m'][0],candidate_width_m=cases[candidateid]['domain_m'][0],
            actual_global_pose_equal=True,**stats,
            magnitude_windows=[dict(lo_ns=lo,hi_ns=hi,relative_l2=relative(np.abs(new['env'][m]),np.abs(ref['env'][m]))) for lo,hi,m in masks]))
    summaries=[]
    for candidate_width,reference_width in [(40,50),(50,70)]:
        rows=[r for r in comparisons if r['candidate_width_m']==candidate_width]
        summaries.append(dict(candidate_width_m=candidate_width,reference_width_m=reference_width,stations=len(rows),
            all_screen_pass=all(r['passes_frozen_domain_screen'] for r in rows),
            max_band_relative_l2=max(r['complex_band_relative_l2'] for r in rows),
            max_window_relative_l2=max(w['complex_relative_l2'] for r in rows for w in r['windows']),
            max_100_400_complex_relative_l2=max(r['windows'][0]['complex_relative_l2'] for r in rows),
            max_100_400_magnitude_relative_l2=max(r['magnitude_windows'][0]['relative_l2'] for r in rows)))
    save(a.out/'domain_narrow_results.json',dict(summary=summaries,comparisons=comparisons,
        new_fdtd=0,grid_m=.025,selected_grid_m=1/30,
        selected_grid_domain_converged=False,thresholds=c['domain_screen']))
    save(a.out/'audit.json',dict(status='PASS',new_fdtd=0,receipts_verified=hashed,independent_checks=checks,
        input_manifest_verified=True,contract_sha256=sha(a.root/'contract.json'),
        analysis_script_sha256=sha(__file__),official_processing_source_sha256=completion['processing_source_sha256']))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for rows,label in [([r for r in comparisons if r['candidate_width_m']==40],'40m对50m（7站）'),
                       ([r for r in comparisons if r['candidate_width_m']==50],'50m对70m（3站）')]:
        rows=sorted(rows,key=lambda r:r['station_m'])
        axes[0].plot([r['station_m'] for r in rows],[100*max(w['complex_relative_l2'] for w in r['windows']) for r in rows],'-o',label=label)
        axes[1].plot([r['station_m'] for r in rows],[100*r['magnitude_windows'][0]['relative_l2'] for r in rows],'-o',label=label)
    axes[0].axhline(1,ls='--',color='black',label='原项目窗复差容差1%')
    axes[0].set(ylabel='最差预定义窗复差（%）',title='全窗复响应筛查')
    axes[1].set(ylabel='包络模相对L2（%）',title='100–400ns地下时段；单列模差')
    for ax in axes:ax.set_xlabel('站位（m）');ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.suptitle('20–110MHz／301点｜复用旧2.5cm／高度50m／PML3m／平滑y；不认证1/30m域宽')
    fig.savefig(a.out/'lateral_domain_screen.png',dpi=170);plt.close(fig)
    ss=[150.,174.,198.];t=results['full_x150.0_y']['t']
    plot_sparse_compare(a.out,'width50_vs70_sparse_bscan.png',ss,
        np.array([results[f'wide70_x{s:.1f}']['real'] for s in ss]).T,
        np.array([results[f'full_x{s:.1f}_y']['real'] for s in ss]).T,t,
        '70m基准 → 50m候选｜旧2.5cm／20–110MHz／同一实际收发坐标',
        reference_events={s:cases[f'full_x{s:.1f}_y']['vertical_reference_events'] for s in ss})
    print(summaries)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','process','analyze'])
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--domain-root',type=Path);p.add_argument('--fraction-root',type=Path);p.add_argument('--out',type=Path)
    a=p.parse_args();{'prepare':prepare,'process':process,'analyze':analyze}[a.action](a)
