"""Audit the six-case authorized coarse-grid study, preserving official products."""
from __future__ import annotations
import argparse
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

from study_line9_v31_smoothing_domain import load,save,sha
from analyze_line9_v31_smoothing_domain import read_sf,independent_check,compare,coefficients,plot_sparse_compare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['root','inputs','reference','reference_inputs','out']:
        p.add_argument('--'+n.replace('_','-'),type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    contract=load(a.root/'contract.json');complete=load(a.root/'stage_A_complete.json')
    assert sha(a.root/'contract.json')==sha(a.inputs/'contract.json')==load(a.root/'execution_start.json')['contract_sha256']
    assert complete['status']=='COMPLETE' and complete['cases']==6
    prior=load(a.reference/'contract.json');assert sha(a.reference/'contract.json')==contract['parent_contract_sha256']
    rows=[];replays={};results={};refs={};hashed=0
    for case in contract['cases']:
        cid=case['id'];receipt=next(v for v in complete['receipts'] if v['case']==cid)
        assert receipt['status']=='PASS'
        for rec in receipt['files']:
            f=a.root/rec['path'];assert sha(f)==rec['sha256'] and f.stat().st_size==rec['bytes'];hashed+=1
        assert sha(a.inputs/case['geometry'])==case['geometry_sha256']
        with h5py.File(a.root/'cases'/cid/'impulse.h5') as h:
            ds=h['rxs/rx1/Ez'];ex=h['srcs/src1/excitation'];src=ex['samples'][:]
            assert ds.dtype==np.dtype('float32') and ds.shape==(case['expected_samples'],)
            assert np.isfinite(ds[:]).all() and np.count_nonzero(src)==1 and np.max(src)==1
            assert np.isclose(ds.attrs['SampleInterval'],case['dt_s'],rtol=1e-12,atol=0)
            assert ex.attrs['TimeSampleOffset']==float(h.attrs['dt'])/2 and ds.attrs['TimeSampleOffset']==0
            for group,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:
                pos=h[group].attrs['Position'];assert np.allclose(pos,case[key],atol=1e-8,rtol=0)
                assert np.allclose(pos[:2]+case['origin_global_x_elevation_m'],case[key.replace('_m','_global_x_elevation_m')],atol=1e-8,rtol=0)
        log=(a.root/'cases'/cid/'solver.log').read_text(encoding='utf-8')
        assert 'CUDA backend using Device 0' in log and 'Interface averaging is enabled' in log
        sf=read_sf(a.root/'cases'/cid/'sfcw.h5');results[cid]=sf
        replays[cid]=independent_check(a.root/'cases'/cid/'impulse.h5',sf)
        refid=case['baseline_case'];refcase=next(v for v in prior['cases'] if v['id']==refid)
        refrec=load(a.reference/'cases'/refid/'receipt.json')
        if refid not in refs:
            for rec in refrec['files']:
                f=a.reference/rec['path'];assert sha(f)==rec['sha256'] and f.stat().st_size==rec['bytes'];hashed+=1
            refs[refid]=read_sf(a.reference/'cases'/refid/'sfcw.h5')
            replays[refid]=independent_check(a.reference/'cases'/refid/'impulse.h5',refs[refid])
        co=coefficients(a.root,cid);assert co['mixed_ez_edge_count']>0
        stats=compare(sf,refs[refid],refcase,contract['domain_screen'])
        rows.append(dict(case=cid,station_m=case['station_m'],grid_m=case['grid_m'],**stats,
            tail_db=receipt['tail_db'],reference_tail_db=refrec['tail_db'],
            new_official_fdtd_seconds=receipt['official_fdtd_seconds'],
            new_build_solver_seconds=receipt['official_build_and_solver_seconds'],
            fdtd_speedup=refrec['official_fdtd_seconds']/receipt['official_fdtd_seconds'],
            build_solver_speedup=refrec['official_build_and_solver_seconds']/receipt['official_build_and_solver_seconds'],
            tx_coordinate_delta_m=case['tx_coordinate_delta_m'],rx_coordinate_delta_m=case['rx_coordinate_delta_m'],
            actual_horizontal_separation_m=case['actual_horizontal_separation_m'],domain_m=case['domain_m'],
            pml_physical_m=case['pml_physical_m'],mixed_ez_edge_count=co['mixed_ez_edge_count'],
            min_center_sandstone_interval_m=min(v['thickness_m'] for v in case['center_column_intervals'] if v['material_id']==3)))
    summary=[]
    for dl in [.0325,.0375]:
        rr=[v for v in rows if v['grid_m']==dl]
        summary.append(dict(grid_m=dl,all_three_stations_pass=all(v['passes_frozen_domain_screen'] for v in rr),
            max_spectrum_error=max(v['complex_band_relative_l2'] for v in rr),
            max_window_error=max(w['complex_relative_l2'] for v in rr for w in v['windows']),
            max_phase_rms_degrees=max(v['strong_bin_phase_rms_degrees'] for v in rr),
            median_fdtd_speedup=float(np.median([v['fdtd_speedup'] for v in rr])),
            median_build_solver_speedup=float(np.median([v['build_solver_speedup'] for v in rr])),
            timing_scope='Single-run current coarse cases vs reused earlier baseline; not repeated benchmark'))
        ordered=sorted(rr,key=lambda r:r['station_m']);stations=[v['station_m'] for v in ordered]
        base=np.array([refs[f'full_x{s:.1f}_y']['real'] for s in stations]).T
        new=np.array([results[v['case']]['real'] for v in ordered]).T
        plot_sparse_compare(a.out,f'grid_{dl*100:g}cm_sparse_bscan.png',stations,base,new,sf['t'],
            f'网格2.5→{dl*100:g} cm｜原生平滑y｜仅150/162/174 m三站，包含声明的位置及几何量化')
    save(a.out/'comparisons.json',dict(summary=summary,cases=rows,thresholds=contract['domain_screen'],
        reference_status='finite-grid finite-domain baseline, tail screen unresolved',
        interpretation=contract['comparison_scope'],current_model_modified=False))
    save(a.out/'audit.json',dict(status='HASHES_NATIVE_METADATA_AND_INDEPENDENT_SUMS_PASS',new_solves=6,reused_references=3,
        verified_receipt_files=hashed,checks=replays,contract_sha256=sha(a.root/'contract.json'),
        analysis_script_sha256=sha(Path(__file__)),new_fdtd_seconds=sum(v['new_official_fdtd_seconds'] for v in rows),
        new_build_solver_seconds=sum(v['new_build_solver_seconds'] for v in rows)))
    fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for col,s in enumerate([150.,162.,174.]):
        ref=refs[f'full_x{s:.1f}_y']
        for dl in [.0325,.0375]:
            row=next(v for v in rows if v['grid_m']==dl and v['station_m']==s)
            sf=results[row['case']];ratio=sf['H']/ref['H']
            axes[0,col].plot(sf['f']/1e6,20*np.log10(np.maximum(np.abs(ratio),1e-15)),label=f'{dl*100:g} cm')
            axes[1,col].plot(sf['f']/1e6,np.angle(ratio,deg=True),label=f'{dl*100:g} cm')
        axes[0,col].set(title=f'{s:g} m｜相对2.5cm复谱',ylabel='幅度比（dB）')
        axes[1,col].set(xlabel='频率（MHz）',ylabel='相位差（°，主值）')
    for ax in axes.flat:ax.axhline(0,color='black',lw=.7);ax.legend();ax.grid(alpha=.2)
    fig.suptitle('未经拟合/归一化的官方direct复响应｜同一20–170MHz/501点')
    fig.savefig(a.out/'complex_spectrum_comparison.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    for dl in [.0325,.0375]:
        rr=sorted([v for v in rows if v['grid_m']==dl],key=lambda v:v['station_m']);stations=[v['station_m'] for v in rr]
        axes[0].plot(stations,[v['complex_band_relative_l2']*100 for v in rr],'-o',label=f'{dl*100:g} cm')
        axes[1].plot(stations,[max(w['complex_relative_l2'] for w in v['windows'])*100 for v in rr],'-o',label=f'{dl*100:g} cm')
        axes[2].plot(stations,[v['build_solver_speedup'] for v in rr],'-o',label=f'{dl*100:g} cm')
    axes[0].axhline(.5,ls='--',color='black');axes[1].axhline(1,ls='--',color='black');axes[2].axhline(1,ls='--',color='black')
    for ax,title,ylabel in zip(axes,['全工作带复谱差','最差预定义窗口差','单次构建＋求解计时比'],['相对L2（%）','相对L2（%）','基准时间／候选时间']):
        ax.set(title=title,xlabel='站位（m）',ylabel=ylabel);ax.legend();ax.grid(alpha=.2)
    fig.suptitle('原冻结项目容差保持；旧基准本身未认证绝对收敛')
    fig.savefig(a.out/'errors_and_speed.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.3),constrained_layout=True)
    cmap=ListedColormap(['#e5eef3','#d5ad69','#87786c','#c0c7ba']);norm=BoundaryNorm(np.arange(-.5,4.5),4)
    geometries=[(.025,a.reference_inputs/next(c for c in prior['cases'] if c['id']=='full_x150.0_y')['geometry'])]
    geometries += [(c['grid_m'],a.inputs/c['geometry']) for c in contract['cases'] if c['station_m']==150.]
    for ax,(dl,path) in zip(axes,geometries):
        with h5py.File(path) as h:
            ox=float(h.attrs['global_profile_x_origin_m']);data=h['data'][:,:,0]
        ix0,ix1=np.floor((np.array([148.,152.])-ox)/dl).astype(int)
        iy0,iy1=np.floor((np.array([433.4,434.3])-410)/dl).astype(int)
        im=ax.imshow(data[ix0:ix1,iy0:iy1].T,origin='lower',aspect='auto',cmap=cmap,norm=norm,interpolation='nearest',
                     extent=[ox+ix0*dl,ox+ix1*dl,410+iy0*dl,410+iy1*dl])
        ax.set(title=f'{dl*100:g} cm｜原材料体素',xlabel='剖面里程（m）',ylabel='高程（m）')
    cb=fig.colorbar(im,ax=axes,ticks=[0,1,2,3],shrink=.8);cb.ax.set_yticklabels(['空气','粉质粘土','泥岩','砂岩'])
    fig.suptitle('150m接触面接近区｜同一连续界面重离散；不是波场或平滑后系数')
    fig.savefig(a.out/'thin_segment_geometry.png',dpi=180);plt.close(fig)
    print(summary)


if __name__=='__main__':main()
