"""Paired full-complex attribution diagnostics; never generate physical labels."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from hs_capsule_identity import sha256
from hs4t2d_path_diagnostic import C0, ray_timings, cover_index
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime, OFFICIAL_PROCESSING_SHA256
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response, reconstruct_time_response


def verify_manifest(folder):
    for row in json.loads((folder/'manifest.json').read_text('utf-8')):
        if sha256(folder/row['file'])!=row['sha256']:
            raise ValueError('input capsule identity changed')


def load_response(folder):
    rows=json.loads((folder/'audit.json').read_text('utf-8'))
    traces=[]
    source=receiver=None
    for row in rows:
        p=folder/row['file']
        if sha256(p)!=row['sha256']:
            raise ValueError('raw output identity differs')
        s=load_source(p)
        r=load_receiver(p,receiver_path='/rxs/rx1',component='Ey')
        if source is None:
            source,receiver=s,r
        if s.dt!=source.dt or s.time_offset!=source.time_offset or not np.array_equal(s.samples,source.samples) or r.dt!=receiver.dt:
            raise ValueError('source/time axis changed across group')
        traces.append(r.samples)
    receiver=replace(receiver,samples=np.stack(traces,axis=1))
    taper=(round(200e-9/receiver.dt)-.25)/len(receiver.samples)
    response=direct_frequency_response(source,receiver,FREQ,tail_taper_fraction=taper)
    return response,rows


def weighted_comparison(reference,candidate):
    weights=np.hanning(len(FREQ))[:,None]
    a,b=weights*reference,weights*candidate
    norm=np.linalg.norm(a)
    if norm==0:
        return {'reason':'zero reference'}
    cross=np.vdot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)) if np.linalg.norm(b)>0 else 0j
    return {'hann_weighted_complex_relative_L2':float(np.linalg.norm(b-a)/norm),
            'hann_weighted_amplitude_L2_ratio':float(np.linalg.norm(b)/norm),
            'complex_coherence_real':float(np.real(cross)),'complex_coherence_imag':float(np.imag(cross)),
            'relative_L2_by_trace':(np.linalg.norm(b-a,axis=0)/np.linalg.norm(a,axis=0)).tolist()}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new analysis directory required')
    verify_official_runtime()
    c=json.loads(args.contract.read_text('utf-8'))
    study=args.contract.parent
    events=[json.loads(s) for s in (study/'execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('status')!='COMPLETED' or events[-1].get('traces')!=67 or events[0]['contract_sha256']!=sha256(args.contract):
        raise ValueError('complete frozen 67-trace study required')
    root=Path(__file__).resolve().parents[1]
    old=root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_results'
    old_raw=root/c['original_121_station_capsule']
    verify_manifest(old)
    verify_manifest(old_raw)
    with np.load(old/'comparison_arrays.npz') as data:
        take=np.arange(0,121,10)
        midpoint=data['midpoint_m'][take].copy()
        old_rough=data['relief2p0_frequency_complex'][:,take].copy()
        old_flat=data['flat_frequency_complex'][:,take].copy()
    responses={}
    identities={}
    for g in c['groups']:
        r,rows=load_response(study/g['id'])
        actual=np.asarray([(row['source_m'][0]+row['receiver_m'][0])/2-g['x_translation_m'] for row in rows])
        expected=midpoint if g['traces']>1 else midpoint[6:7]
        if not np.allclose(actual,expected,atol=1e-12,rtol=0):
            raise ValueError('paired station coordinates differ')
        responses[g['id']]=r
        identities[g['id']]={'audit_sha256':sha256(study/g['id']/'audit.json'),'source_quantity':r.source.quantity,
                             'source_units':r.source.units,'source_spatial_scale':r.source.spatial_scale,
                             'dt_s':r.source.dt,'traces':len(rows)}
    base=responses['halfspace12']
    spectra={'narrow_halfspace_difference':old_rough-base.response,
             'flat_halfspace_difference':old_flat-base.response,
             'old_flat_difference':old_rough-old_flat,
             'wide_halfspace_difference':responses['wide36_rough'].response-responses['wide36_halfspace'].response,
             'low_halfspace_difference':responses['low2_rough'].response-responses['low2_halfspace'].response,
             'zfine_halfspace_difference':responses['zfine_rough'].response-responses['zfine_halfspace'].response}
    air_shift=2*(15-2)/C0
    spectra['low_air_delay_aligned']=spectra['low_halfspace_difference']*np.exp(-2j*np.pi*FREQ[:,None]*air_shift)
    products={name:reconstruct_time_response(replace(base,response=value),zero_pad_factor=8,window='hann') for name,value in spectra.items()}
    time_ns=products['narrow_halfspace_difference'].time*1e9
    window=(time_ns>=160)&(time_ns<=220)
    report_metrics={}
    for name,p in products.items():
        if name=='low_halfspace_difference':
            selected=(time_ns>=160-air_shift*1e9)&(time_ns<=220-air_shift*1e9)
        else:
            selected=window
        envelope=2*np.abs(p.complex_envelope[selected])
        peaks=time_ns[selected][np.argmax(envelope,axis=0)]
        report_metrics[name]={'window_ns':[float(time_ns[selected][0]),float(time_ns[selected][-1])],
                              'signed_window_energy':float(np.sum(p.real_bandpass[selected]**2)),
                              'signed_window_peak_abs':float(np.max(np.abs(p.real_bandpass[selected]))),
                              'envelope_peak_time_ns':peaks.tolist(),'envelope_peak_span_ns':float(np.ptp(peaks)),
                              'envelope_max_on_window_boundary_traces':int(np.count_nonzero((np.argmax(envelope,axis=0)==0)|(np.argmax(envelope,axis=0)==len(envelope)-1)))}
    comparisons={'lateral_12_to_36':weighted_comparison(spectra['narrow_halfspace_difference'],spectra['wide_halfspace_difference']),
                 'z_grid_5cm_to_2p5cm_center':weighted_comparison(spectra['narrow_halfspace_difference'][:,6:7],spectra['zfine_halfspace_difference']),
                 'reference_flat_to_halfspace':weighted_comparison(spectra['narrow_halfspace_difference'],spectra['old_flat_difference'])}
    # Processing-only controls: same frequency samples and reconstruction, no solver.
    profile=np.genfromtxt(root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',delimiter=',',names=True)
    z=profile['grid_z_m']
    rays={str(h):ray_timings(z,midpoint,h) for h in (15,2)}
    toy_time=np.asarray(rays['15']['vertical_group_time_ns'])*1e-9
    toy_spectrum=np.exp(-2j*np.pi*FREQ[:,None]*toy_time)
    toy=reconstruct_time_response(replace(base,response=toy_spectrum),zero_pad_factor=8,window='hann')
    toy_peaks=time_ns[window][np.argmax(np.abs(toy.complex_envelope[window]),axis=0)]
    processing_control={'type':'array-only known delayed single echo; not FDTD truth',
                        'expected_peak_span_ns':float(np.ptp(toy_time)*1e9),
                        'reconstructed_peak_span_ns':float(np.ptp(toy_peaks)),
                        'maximum_peak_time_quantization_error_ns':float(np.max(np.abs(toy_peaks-toy_time*1e9)))}
    # Isolate late-record taper sensitivity without another simulation.
    raw_rough_path=old_raw/'relief2p0/profile61.h5'
    raw_half_path=study/'halfspace12'/json.loads((study/'halfspace12/audit.json').read_text('utf-8'))[6]['file']
    untapered=[]
    truncated=[]
    for path in (raw_rough_path,raw_half_path):
        s=load_source(path)
        r=load_receiver(path,receiver_path='/rxs/rx1',component='Ey')
        untapered.append(direct_frequency_response(s,r,FREQ,tail_taper_fraction=0).response[:,None])
        count=int(np.floor(300e-9/r.dt))+1
        truncated.append(direct_frequency_response(s,replace(r,samples=r.samples[:count]),FREQ,tail_taper_fraction=0).response[:,None])
    spectra['center_no_tail_taper']=untapered[0]-untapered[1]
    spectra['center_record_300ns_no_taper']=truncated[0]-truncated[1]
    comparisons['tail_taper_200ns_to_none_center']=weighted_comparison(spectra['narrow_halfspace_difference'][:,6:7],spectra['center_no_tail_taper'])
    comparisons['record_600ns_to_300ns_center']=weighted_comparison(spectra['center_no_tail_taper'],spectra['center_record_300ns_no_taper'])
    args.out.mkdir(parents=True)
    saved=(time_ns>=0)&(time_ns<=250)
    arrays={'time_ns':time_ns[saved],'midpoint_m':midpoint,'frequency_Hz':FREQ,
            'toy_vertical_single_echo_signed':toy.real_bandpass[saved],
            'toy_vertical_single_echo_complex_envelope':toy.complex_envelope[saved]}
    for name,p in products.items():
        arrays[name+'_signed']=p.real_bandpass[saved]
        arrays[name+'_complex_envelope']=p.complex_envelope[saved]
    for name,value in spectra.items():
        arrays[name+'_frequency_complex']=value
    np.savez_compressed(args.out/'diagnostic_arrays.npz',**arrays)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    columns=[('old_flat_difference','12 m域 / 15 m航高\n减平界面（旧参考）'),
             ('narrow_halfspace_difference','12 m域 / 15 m航高\n减全覆盖层（新参考）'),
             ('wide_halfspace_difference','36 m域 / 15 m航高\n减全覆盖层'),
             ('low_air_delay_aligned','12 m域 / 2 m航高\n减全覆盖层；补回空气时延')]
    fig,axes=plt.subplots(2,4,figsize=(18,8),constrained_layout=True)
    scales={}
    for j,(name,title) in enumerate(columns):
        p=products[name]
        scale=float(np.quantile(np.abs(p.real_bandpass[window]),.999)) or 1
        scales[name]=scale
        for i,(value,label,cmap) in enumerate([(p.real_bandpass,'带符号响应','gray'),(2*np.abs(p.complex_envelope),'复包络模','viridis')]):
            ax=axes[i,j]
            limit=1 if i==0 else 2
            im=ax.pcolormesh(midpoint,time_ns[window],value[window]/scale,shading='nearest',cmap=cmap,
                             vmin=-1 if i==0 else 0,vmax=limit,rasterized=True)
            ax.set_ylim(220,160)
            ax.set_xlabel('原模型收发中点 X（m）')
            ax.set_ylabel('重建时间（ns）')
            ax.set_title(title+'\n'+label)
            fig.colorbar(im,ax=ax,label='每个工况一个全矩阵标量归一化')
    fig.suptitle('V4.0.0 / CUDA FP64 / 0.8 m界面 / 13个匹配站位；官方20–170 MHz/Hann\n'
                 '形态诊断用单标量归一化，无逐道增益/SVD；幅度比较须读原数组与指标')
    fig.savefig(args.out/'controlled_bscan_shapes.png',dpi=140)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,5),constrained_layout=True)
    axes[0].stairs(12-z,np.r_[profile['x0_m'],profile['x1_m'][-1]],baseline=None)
    axes[0].set_xlim(midpoint[0],midpoint[-1]); axes[0].set_ylim(3.4,2.4)
    axes[0].set_title('输入深度剖面'); axes[0].set_ylabel('距地表深度（m）')
    for h,label in [('15','15 m航高'),('2','2 m航高')]:
        shift=air_shift*1e9 if h=='2' else 0
        ray=np.asarray(rays[h]['minimum_ray_group_time_ns'])+shift
        vertical=np.asarray(rays[h]['vertical_group_time_ns'])+shift
        axes[1].plot(midpoint,ray,label=label+' 最早路径')
        axes[1].plot(midpoint,vertical,linestyle='--',label=label+' 垂直估算')
    axes[1].invert_yaxis(); axes[1].set_title('95 MHz模型路径诊断\n最早路径不是最强包络峰真值')
    axes[1].set_ylabel('时间（ns，2 m补回空气时延）'); axes[1].legend(fontsize=8)
    im=axes[2].pcolormesh(midpoint,time_ns[window],2*np.abs(toy.complex_envelope[window]),shading='nearest',cmap='viridis')
    axes[2].set_ylim(220,160); axes[2].set_title('已知垂直到时的单回波数组\n同一官方频带/重建，可保留曲线')
    axes[2].set_ylabel('时间（ns）'); fig.colorbar(im,ax=axes[2],label='构造响应，任意幅度')
    for ax in axes:ax.set_xlabel('收发中点 X（m）')
    fig.savefig(args.out/'ray_and_processing_control.png',dpi=145)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    center=products['narrow_halfspace_difference'].real_bandpass[:,6]
    ax.plot(time_ns[window],center[window],label='Z网格5 cm')
    ax.plot(time_ns[window],products['zfine_halfspace_difference'].real_bandpass[window,0],label='Z网格2.5 cm')
    ax.set_xlabel('重建时间（ns）'); ax.set_ylabel('带符号场/源响应'); ax.legend()
    ax.set_title('中心道：同条件全覆盖层差分；无幅度归一化，只有竖向网格细化')
    fig.savefig(args.out/'vertical_grid_center.png',dpi=145); plt.close(fig)
    report={'status':'COMPLETED','contract_sha256':sha256(args.contract),'new_traces':67,
            'source_audit':identities,'official_processing_sha256':OFFICIAL_PROCESSING_SHA256,
            'metrics':report_metrics,'comparisons':comparisons,'ray_diagnostic':rays,
            'processing_only_control':processing_control,'low_height_fixed_air_delay_s':air_shift,
            'shape_plot_normalization_scalars':scales,
            'cells_per_cover_phase_wavelength_5cm':{str(f):float(C0/f/cover_index(f).real/.05) for f in (20e6,95e6,170e6)},
            'interpretation':'Controlled development diagnostics. No complete grid/domain convergence, PML-only attribution, clean labels or physical acceptance.',
            'code_identities':{Path(p).name:sha256(p) for p in [__file__,root/'scripts/hs4t2d_path_diagnostic.py']}}
    (args.out/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    # Reject concurrent changes to all raw inputs after processing.
    verify_manifest(old); verify_manifest(old_raw)
    for g in c['groups']:
        for row in json.loads((study/g['id']/'audit.json').read_text('utf-8')):
            if sha256(study/g['id']/row['file'])!=row['sha256']:
                raise ValueError('raw changed during analysis')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)}
        for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'metrics':report_metrics,'comparisons':comparisons},indent=2))


if __name__=='__main__':
    main()
