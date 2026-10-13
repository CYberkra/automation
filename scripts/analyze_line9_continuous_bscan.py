"""Audit102 matched runs; compare8 frozen backgrounds without clean labels."""
from __future__ import annotations
import argparse
from pathlib import Path
import time
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from study_line9_v31_smoothing_domain import load,save,sha
from analyze_line9_v31_smoothing_domain import read_sf,independent_check,compare,coefficients,relative
from research_operator_contract import apply_configuration,local_mean_control,ConfigUnavailable


LABELS={'identity':'不处理','mean':'均值道扣除','sliding5':'滑动背景5道','sliding11':'滑动背景11道',
        'sliding21':'滑动背景21道','svd1':'SVD去前1阶','svd2':'SVD去前2阶','svd3':'SVD去前3阶'}


def derivative_metrics(candidate,reference,dx):
    stats={}
    for order in [1,2]:
        dref=np.diff(reference,n=order,axis=1)/dx**order
        dcand=np.diff(candidate,n=order,axis=1)/dx**order
        stats[f'order{order}_relative_l2']=relative(dcand,dref)
        stats[f'order{order}_reference_norm']=float(np.linalg.norm(dref))
        stats[f'order{order}_difference_norm']=float(np.linalg.norm(dcand-dref))
    return stats


def metric_checks():
    x=np.arange(51)*.6
    linear=np.tile(3*x+2,(5,1));quadratic=np.tile(x*x,(5,1))
    e1=float(np.abs(np.diff(linear,axis=1)/.6-3).max())
    e2=float(np.abs(np.diff(quadratic,n=2,axis=1)/.6**2-2).max())
    assert e1<1e-10 and e2<1e-9
    assert derivative_metrics(quadratic,quadratic,.6)['order2_relative_l2']==0
    ident=apply_configuration(quadratic,'B0_G1_BG')['output'];assert np.array_equal(ident,quadratic)
    return dict(status='PASS',linear_first_derivative_max_error=e1,quadratic_second_derivative_max_error=e2,
        identical_inputs_difference_zero=True,identity_exact=True,scope='Numerical metric checks, not physical acceptance')


def background(x,method):
    started=time.perf_counter()
    if 'width' in method:
        y=local_mean_control(x,method['width'])
        diag=dict(kind='sliding_mean',width=method['width'],strength=1,edges='truncated_count_normalized')
    else:
        result=apply_configuration(x,method['config_id']);y=result['output']
        assert np.array_equal(result['gain_curve'],np.ones(x.shape[0]))
        diag=result['steps'][0]['diagnostics']
    return y,dict(**diag,wall_seconds=time.perf_counter()-started)


def plot_panel(ax,x,stations,t,lim,title):
    im=ax.imshow(x,origin='upper',aspect='auto',cmap='gray',vmin=-lim,vmax=lim,interpolation='nearest',
        extent=[stations[0]-.3,stations[-1]+.3,t[-1]+(t[1]-t[0])/2,t[0]-(t[1]-t[0])/2])
    ax.set(title=title,xlabel='测线里程（m）',ylabel='时间（ns）')
    return im


def event_lines(ax,stations,fine_cases,tlo,thi):
    for selector,label,color in [('cover','覆盖层底垂向参考','#f0b43c'),('lower_top','下部泥岩顶参考','#90c2ff'),
                                 ('lower_bottom','下部泥岩底参考','#e3a1e5')]:
        values=[]
        for case in fine_cases:
            events=case['vertical_reference_events']
            ev=(events[1] if selector=='cover' and len(events)>1 else
                next((e for e in reversed(events) if e['label']==('3_to_2_vertical_reference' if selector=='lower_top' else '2_to_3_vertical_reference')),None))
            values.append(np.nan if ev is None else ev['time_ns'])
        ys=np.array(values);ys[(ys<tlo)|(ys>thi)]=np.nan
        ax.plot(stations,ys,ls='--',lw=.8,color=color,label=label)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    c=load(a.root/'contract.json');complete=load(a.root/'stage_A_complete.json');start=load(a.root/'execution_start.json')
    assert complete['status']=='COMPLETE' and complete['cases']==102==c['max_new_solves']
    assert sha(a.root/'contract.json')==start['contract_sha256'] and sha(a.root/'study.py')==start['script_sha256']
    assert start['gprMax']=='4.0.1' and '4090' in start['gpu']
    for r in load(a.root/'input_manifest.json')['files']:
        f=a.root/r['path'];assert sha(f)==r['sha256'] and f.stat().st_size==r['bytes']
    sf={};checks={};receipts={};hashed=0;tail=[]
    for case in c['cases']:
        cid=case['id'];dest=a.root/'cases'/cid;rec=load(dest/'receipt.json')
        assert rec==next(r for r in complete['receipts'] if r['case']==cid) and rec['status']=='PASS'
        for r in rec['files']:
            f=a.root/r['path'];assert sha(f)==r['sha256'] and f.stat().st_size==r['bytes'];hashed+=1
        assert sha(a.root/case['geometry'])==case['geometry_sha256']
        with h5py.File(dest/'impulse.h5') as h:
            y=h['rxs/rx1/Ez'];x=h['srcs/src1/excitation'];source=x['samples'][:]
            assert y.dtype==np.dtype('float32') and y.shape==(case['expected_samples'],)
            assert np.isfinite(y[:]).all() and np.count_nonzero(source)==1 and np.max(source)==1
            assert np.isclose(y.attrs['SampleInterval'],case['dt_s'],rtol=1e-12,atol=0)
            assert x.attrs['TimeSampleOffset']==float(h.attrs['dt'])/2 and y.attrs['TimeSampleOffset']==0
            for group,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:
                assert np.allclose(h[group].attrs['Position'],case[key],atol=1e-8,rtol=0)
        log=(dest/'solver.log').read_text(encoding='utf-8')
        assert 'CUDA backend using Device 0' in log and 'Interface averaging is enabled' in log
        sf[cid]=read_sf(dest/'sfcw.h5',frequency_spec=(20e6,110e6,301))
        checks[cid]=independent_check(dest/'impulse.h5',sf[cid])
        co=coefficients(a.root,cid);assert co['mixed_ez_edge_count']>0
        checks[cid]['mixed_ez_edge_count']=co['mixed_ez_edge_count']
        receipts[cid]=rec;tail.append(rec['tail_db'])
    stations=np.array(c['anchors_m']);assert len(stations)==51 and np.allclose(np.diff(stations),.6)
    groups={label:[next(v for v in c['cases'] if v['id']==f'{label}_x{s:.1f}') for s in stations] for label in ['fine','coarse']}
    rows=[]
    for j,s in enumerate(stations):
        fine,coarse=groups['fine'][j],groups['coarse'][j]
        for key in ['tx_global_x_elevation_m','rx_global_x_elevation_m','origin_global_x_elevation_m','domain_m']:
            assert fine[key]==coarse[key]
        r,q=sf[fine['id']],sf[coarse['id']]
        assert np.array_equal(r['t'],q['t']) and np.array_equal(r['f'],q['f'])
        stats=compare(q,r,fine,c['domain_screen'])
        rows.append(dict(station_m=s,**stats))
    t=sf[groups['fine'][0]['id']]['t'];freq=sf[groups['fine'][0]['id']]['f']
    real={g:np.column_stack([sf[v['id']]['real'] for v in cases]) for g,cases in groups.items()}
    env={g:np.column_stack([sf[v['id']]['env'] for v in cases]) for g,cases in groups.items()}
    H={g:np.column_stack([sf[v['id']]['H'] for v in cases]) for g,cases in groups.items()}
    pmask=(t>=0)&(t<1100);pt=t[pmask];input_arrays={g:real[g][pmask] for g in groups}
    outputs={};diagnostics={};arrays=dict(stations_m=stations,frequency_hz=freq,time_ns_full=t,time_ns_processing=pt)
    for g in groups:
        arrays[f'{g}_official_real_full']=real[g];arrays[f'{g}_official_envelope_full']=env[g];arrays[f'{g}_official_H']=H[g]
    process_rows=[];raw_windows=[]
    for lo,hi in c['evaluation_windows_ns']:
        m=(pt>=lo)&(pt<hi);q,r=input_arrays['coarse'][m],input_arrays['fine'][m]
        fullm=(t>=lo)&(t<hi)
        raw_windows.append(dict(lo_ns=lo,hi_ns=hi,real_relative_l2=relative(q,r),
            complex_envelope_relative_l2=relative(env['coarse'][fullm],env['fine'][fullm]),
            magnitude_relative_l2=relative(np.abs(env['coarse'][fullm]),np.abs(env['fine'][fullm])),
            **derivative_metrics(q,r,.6)))
    for method in c['background_methods']:
        name=method['name'];outputs[name]={};diagnostics[name]={}
        for g in groups:
            try:
                y,diag=background(input_arrays[g],method);assert y.shape==input_arrays[g].shape and np.isfinite(y).all()
                outputs[name][g]=y
                diagnostics[name][g]={**diag,'operator_status':diag.get('status'),'status':'AVAILABLE'}
                arrays[f'{g}_{name}']=y;arrays[f'{g}_{name}_removed']=input_arrays[g]-y
            except ConfigUnavailable as e:
                diagnostics[name][g]=dict(status='UNAVAILABLE',reason=str(e))
        if len(outputs[name])!=2:
            continue
        for lo,hi in c['evaluation_windows_ns']:
            m=(pt>=lo)&(pt<hi);q,r=outputs[name]['coarse'][m],outputs[name]['fine'][m]
            rawr=input_arrays['fine'][m];rawdelta=input_arrays['coarse'][m]-rawr;delta=q-r
            process_rows.append(dict(method=name,lo_ns=lo,hi_ns=hi,
                processed_reference_relative_l2=relative(q,r),
                difference_norm_on_raw_fine_scale=float(np.linalg.norm(delta)/np.linalg.norm(rawr)),
                cross_grid_difference_amplification=float(np.linalg.norm(delta)/np.linalg.norm(rawdelta)) if np.linalg.norm(rawdelta)>0 else None,
                fine_output_norm_over_raw=float(np.linalg.norm(r)/np.linalg.norm(rawr)),
                fine_removed_norm_over_raw=float(np.linalg.norm(rawr-r)/np.linalg.norm(rawr)),
                **derivative_metrics(q,r,.6)))
    np.savez_compressed(a.out/'bscan_arrays.npz',**arrays)
    # Timings are kept separately so metric/figure replays remain deterministic.
    timings={name:{g:diag.pop('wall_seconds',None) for g,diag in pair.items()} for name,pair in diagnostics.items()}
    save(a.out/'processing_timings.json',timings)
    save(a.out/'processing_metrics.json',dict(raw_windows=raw_windows,methods=process_rows,diagnostics=diagnostics,
        input_domain='time_real,0<=t<1100ns',amplitude_scale='official Ez/source response; no gain or amplitude fit',
        reference_status='finite2.5cm reference, not clean truth; no SNR/detection labels',
        operator_script_sha256=sha(Path(__file__).with_name('research_operator_contract.py'))))
    summary=dict(traces_per_grid=51,new_solves=102,all_original_screen_pass=all(r['passes_frozen_domain_screen'] for r in rows),
        passed_trace_count=sum(r['passes_frozen_domain_screen'] for r in rows),
        max_spectrum_relative_l2=max(r['complex_band_relative_l2'] for r in rows),
        max_window_relative_l2=max(w['complex_relative_l2'] for r in rows for w in r['windows']),
        tail_db_min=float(min(tail)),tail_db_max=float(max(tail)),tail_exceed_minus60_count=sum(v>-60 for v in tail),
        fdtd_seconds=sum(r['official_fdtd_seconds'] for r in receipts.values()),
        build_solver_seconds=sum(r['official_build_and_solver_seconds'] for r in receipts.values()),
        batch_wall_seconds=complete['elapsed_since_first_stage_s'])
    save(a.out/'grid_metrics.json',dict(summary=summary,traces=rows,thresholds=c['domain_screen']))
    save(a.out/'audit.json',dict(status='PASS',verified_receipt_files=hashed,checks=checks,
        contract_sha256=sha(a.root/'contract.json'),analysis_script_sha256=sha(__file__),
        source_dtype='float32',analysis_dtype='float64/complex128; does not upgrade FDTD precision',
        strict_pose_and_domain_equal_all51=True,metric_checks=metric_checks()))
    target=(t>=100)&(t<400);lim=float(np.quantile(np.abs(real['fine'][target]),.995))
    for lo,hi,name in [(0,1100,'raw_bscan_0_1100.png'),(20,450,'raw_bscan_20_450.png')]:
        m=(t>=lo)&(t<hi);fig,axes=plt.subplots(3,1,figsize=(13,11),sharex=True,constrained_layout=True)
        for ax,x,title in zip(axes,[real['fine'],real['coarse'],real['coarse']-real['fine']],
            ['2.5cm基准｜原始总场','精确1/30m｜原始总场','粗－细｜差量，不自动叫噪声']):
            im=plot_panel(ax,x[m],stations,t[m],lim,title);event_lines(ax,stations,groups['fine'],lo,hi)
        axes[0].legend(loc='lower right',fontsize=8)
        fig.colorbar(im,ax=axes,label='官方带符号响应；三行共用色标',shrink=.65)
        fig.suptitle(f'150–180m／51道／20–110MHz｜严格同位置／50m域／平滑y\n共同±{lim:.4g}色标；浅部强波可饱和显示，数值未截断、未加增益')
        fig.savefig(a.out/name,dpi=170);plt.close(fig)
    m=(pt>=100)&(pt<400)
    for row_scaled,filename in [(False,'background_methods_common_scale.png'),(True,'background_methods_row_scale.png')]:
        fig,axes=plt.subplots(8,3,figsize=(17,24),constrained_layout=True)
        for i,method in enumerate(c['background_methods']):
            name=method['name'];pair=outputs[name]
            if len(pair)!=2:
                for ax in axes[i]:ax.text(.1,.5,'该配置不可用；见diagnostics',transform=ax.transAxes)
                continue
            r,q=pair['fine'][m],pair['coarse'][m];rowlim=max(float(np.quantile(np.abs(r),.995)),np.finfo(float).tiny) if row_scaled else lim
            for ax,x,label in zip(axes[i],[r,q,q-r],['2.5cm','1/30m','粗－细']):
                plot_panel(ax,x,stations,pt[m],rowlim,f'{LABELS[name]}｜{label}｜±{rowlim:.3g}')
        fig.suptitle('100–400ns｜'+('每行独立色标，只能在同一行横比；无逐道归一化' if row_scaled else '所有面板同一色标，保留真实幅度差；无增益')+
            '\n两套各自独立执行相同参数；不以亮度或连续性判定成功')
        fig.savefig(a.out/filename,dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(12,8),constrained_layout=True)
    for lo,hi in c['evaluation_windows_ns']:
        mask=(t>=lo)&(t<hi);errors=[relative(env['coarse'][mask,j],env['fine'][mask,j])*100 for j in range(51)]
        axes[0].plot(stations,errors,'-o',ms=3,label=f'{lo}–{hi}ns复包络')
    axes[0].set(ylabel='逐道粗细相对L2（%）',xlabel='站位（m）');axes[0].legend();axes[0].grid(alpha=.2)
    for method in c['background_methods']:
        r=next((r for r in process_rows if r['method']==method['name'] and r['lo_ns']==100),None)
        if r:axes[1].bar(LABELS[method['name']],r['cross_grid_difference_amplification'])
    axes[1].axhline(1,color='black',ls='--');axes[1].set(ylabel='处理后粗细差范数／处理前粗细差范数',title='100–400ns｜粗细差放大系数，不是SNR或地层保真评分')
    fig.savefig(a.out/'trace_errors_and_processing_sensitivity.png',dpi=170);plt.close(fig)
    # Removed components are shown separately; no assumption that they are noise.
    fig,axes=plt.subplots(4,2,figsize=(14,16),constrained_layout=True)
    for ax,method in zip(axes.flat,c['background_methods']):
        name=method['name']
        if 'fine' not in outputs[name]:continue
        removed=input_arrays['fine']-outputs[name]['fine']
        plot_panel(ax,removed[m],stations,pt[m],lim,LABELS[name]+'｜细网格去除分量，统一色标')
    fig.suptitle('2.5cm去除分量｜100–400ns｜真实地层也可能被扣除，未经分离不标作噪声')
    fig.savefig(a.out/'removed_components.png',dpi=140);plt.close(fig)
    print(summary)


if __name__=='__main__':
    main()
