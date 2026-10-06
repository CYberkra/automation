"""Audit native Line9 traces, retain501 complex tones, render labelled SFCW panels."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from build_pdf_profile_geometry import digest, save_json
import run_line9_2d_first as runner

PROCESSING_SHA='adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b'
FREQUENCY=np.arange(501,dtype=float)*.3e6+20e6


def relative(a,b):
    denominator=np.linalg.norm(b)
    if denominator==0: raise ValueError('Zero spectral reference')
    return float(np.linalg.norm(a-b)/denominator)


def process_trace(raw):
    with h5py.File(raw) as h:
        if str(h['srcs/src1/excitation'].attrs['Polarisation'])!='z':
            raise ValueError('Expected native z polarization')
    source=sf.load_source(raw);receiver=sf.load_receiver(raw,'/rxs/rx1','Ez')
    if source.spatial_scale!=.025 or source.source_type.lower().find('hertzian')<0:
        raise ValueError('Expected native z Hertzian current element0.025m')
    if not np.isfinite(source.samples).all() or not np.isfinite(receiver.samples).all():
        raise ValueError('Nonfinite native histories')
    product=sf.direct_frequency_response(source,receiver,FREQUENCY,source_floor_db=-100,tail_taper_fraction=0)
    if not product.source_valid.all(): raise ValueError('At least one SFCW tone has insufficient source support')
    product=replace(product,response=product.response/source.spatial_scale)
    take=np.array([0,83,167,250,333,417,500])
    xf=source.dt*(np.exp(-2j*np.pi*FREQUENCY[take,None]*source.times)@source.samples)
    yf=receiver.dt*(np.exp(-2j*np.pi*FREQUENCY[take,None]*receiver.times)@receiver.samples)
    error=relative(product.response[take],yf/xf/source.spatial_scale)
    if error>1e-9: raise ValueError('Independent exact physical-frequency DFT mismatch')
    fraction=(round(20e-9/receiver.dt)-.25)/len(receiver.samples)
    tapered=sf.direct_frequency_response(source,receiver,FREQUENCY,tail_taper_fraction=fraction)
    tail=receiver.samples[-round(20e-9/receiver.dt):]
    peak=float(np.max(abs(receiver.samples)))
    if peak==0: raise ValueError('All-zero received field')
    return product,receiver,dict(independent_DFT_relative_L2=error,
        source_length_m=source.spatial_scale,source_time_offset_s=source.time_offset,
        receiver_time_offset_s=receiver.time_offset,
        sampled_peak_current_moment_Am=float(np.max(abs(source.samples))*source.spatial_scale),
        source_min_over_max_spectrum=float(np.min(abs(product.source_spectrum))/np.max(abs(product.source_spectrum))),
        last20ns_peak_over_record_peak=float(np.max(abs(tail))/peak),
        last20ns_RMS_over_record_peak=float(np.sqrt(np.mean(tail**2))/peak),
        last_sample_over_record_peak=float(abs(receiver.samples[-1])/peak),
        last20ns_taper_relative_spectral_change=relative(tapered.response/source.spatial_scale,product.response))


def main(study,out):
    if out.exists(): raise ValueError('Fresh analysis output required')
    if digest(Path(sf.__file__))!=PROCESSING_SHA: raise ValueError('Reviewed V4 SFCW implementation required')
    path=study/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    if c['stage']=='snapshot':
        from line9_large_domain_snapshot import audit as snapshot_audit
        current=snapshot_audit(path,True)
    else:
        current=runner.audit(path,True)
    stored=json.loads((study/'completed_verification.json').read_text('utf-8'))
    if current!=stored: raise ValueError('Completed native verification changed')
    rows=[]
    for g,r in zip(c['groups'],current['groups']):
        rows.append(dict(id=g['id'],s=g['acquisition_s_m'],x=g['profile_x_m'],
                         raw=r['raw_path'],sha=r['raw_sha256'],reused=False))
    if c['stage']=='preview':
        task=json.loads(runner.TASK.read_text('utf-8'))
        for name,item in c['reuse'].items():
            if digest(Path(item['raw_path']))!=item['raw_sha256']: raise ValueError('Reused output changed')
            index=int(name.rsplit('s',1)[1]);s=index*.5
            rows.append(dict(id=name,s=s,x=220-s,raw=item['raw_path'],sha=item['raw_sha256'],reused=True))
        if set(r['id'] for r in rows)!=set(task['acquisition']['preview_ids']): raise ValueError('Incomplete99-station preview')
    rows.sort(key=lambda r:r['s'])
    spectra=[];raws=[];records=[];reference_time=None
    for row in rows:
        raw=Path(row['raw'])
        if digest(raw)!=row['sha']: raise ValueError('Raw identity changed')
        with h5py.File(raw) as h:
            if h['rxs/rx1/Ez'].dtype!=np.float64: raise ValueError('Original dtype is not float64')
        product,receiver,record=process_trace(raw)
        if reference_time is None: reference_time=receiver.times
        np.testing.assert_array_equal(receiver.times,reference_time)
        records.append(dict(**row,**record));spectra.append(product.response);raws.append(receiver.samples)
    matrix=np.stack(spectra,axis=1);raw_matrix=np.stack(raws,axis=1)
    positions=np.asarray([r['s'] for r in rows]);arrays=dict(frequency_Hz=FREQUENCY,response_complex=matrix,
        acquisition_s_m=positions,profile_X_m=220-positions,native_time_s=reference_time,native_Ez_V_m=raw_matrix)
    profiles={};inverse_errors={}
    for window in ['rectangular','hann']:
        result=sf.reconstruct_time_response(replace(product,response=matrix),window=window,
            zero_pad_factor=8,normalise_window=True,time_shift=0)
        take=np.arange(0,len(result.time),53)
        manual=(np.exp(2j*np.pi*result.time[take,None]*FREQUENCY[None,:])@(result.weights[:,None]*matrix))/501
        error=relative(result.complex_bandpass[take],manual)
        if error>1e-9: raise ValueError('Independent inverse sum mismatch')
        inverse_errors[window]=error;profiles[window]=result.real_bandpass
        arrays[window+'_complex_bandpass']=result.complex_bandpass
        arrays[window+'_complex_envelope']=result.complex_envelope
        arrays[window+'_signed']=result.real_bandpass;arrays[window+'_weights']=result.weights
    arrays['sfcw_time_s']=result.time
    out.mkdir(parents=True);np.savez_compressed(out/'sfcw_arrays.npz',**arrays)
    # All panels share physical axes and each product family's color scale.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(3,1,figsize=(15,11),layout='constrained')
    panels=[(reference_time*1e9,raw_matrix,'原始Ricker总场 Ez；未去背景、未增益'),
            (result.time*1e9,profiles['rectangular'],'SFCW矩形频窗：源归一化总场；20–170MHz / 501点'),
            (result.time*1e9,profiles['hann'],'SFCW Hann频窗：源归一化总场；20–170MHz / 501点')]
    mask=(result.time<=1200e-9)
    sfcw_scale=max(np.max(abs(p[mask])) for p in profiles.values())
    for i,(t,v,label) in enumerate(panels):
        visible=t<=1200;scale=float(np.max(abs(v[visible]))) if i==0 else float(sfcw_scale)
        tv=t[visible];vv=v[visible]
        time_edges=np.r_[tv[0]-(tv[1]-tv[0])/2,(tv[:-1]+tv[1:])/2,tv[-1]+(tv[-1]-tv[-2])/2]
        if c['stage']!='preview':
            for j,s in enumerate(positions):
                axes[i].pcolormesh([s-.25,s+.25],time_edges,vv[:,j:j+1],cmap='gray',vmin=-scale,vmax=scale,rasterized=True)
            label+=f'；{len(positions)}站，白色间隔未计算'
        else:
            edges=np.r_[0,(positions[:-1]+positions[1:])/2,195]
            axes[i].pcolormesh(edges,time_edges,vv,cmap='gray',vmin=-scale,vmax=scale,rasterized=True)
            label+='；99站粗采样，未插值补道'
        axes[i].axhline(2*np.sqrt(15**2+.65**2)/299792458*1e9,color='#cf2525',ls='--',lw=.9,
            label='平地镜面地表参考约100ns；地形/色散事件需另解释')
        axes[i].set(xlim=(0,195),ylim=(1200,0),ylabel='时间 / ns',title=label)
        axes[i].legend(loc='lower right',fontsize=8)
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-scale,scale),cmap='gray'),ax=axes[i],
            label='Ez / (V/m)' if i==0 else '场响应 / 电流矩；共享色标')
    axes[-1].set_xlabel('采集距离 s / m（左端原X220，右端原X25）')
    fig.suptitle('九号线四材料：15m离地；2D沿线1.3m代理；网格2.5cm\n灰度为带符号数据，共享SFCW色标；材料是研究假设，未认证现场效果')
    fig.savefig(out/'line9_sfcw_panels.png',dpi=140);plt.close(fig)
    report=dict(status='PASS_PROCESSING_IDENTITIES_NOT_PHYSICAL_ACCEPTANCE',stage=c['stage'],
        calls_solver=False,contract_sha256=digest(path),completed_verification_sha256=digest(study/'completed_verification.json'),
        task_sha256=c['task_sha256'],script_sha256=digest(Path(__file__)),official_processing_sha256=PROCESSING_SHA,
        trace_count=len(rows),frequency_count=501,source_moment_normalisation=True,
        windows=['rectangular','hann'],gain=False,background_removal=False,extra_time_shift_s=0,
        tail_taper_fraction_primary=0,independent_inverse_relative_L2=inverse_errors,records=records,
        NPZ_sha256=digest(out/'sfcw_arrays.npz'),PNG_sha256=digest(out/'line9_sfcw_panels.png'),
        limits='Native and transform consistency only. Tail diagnostics do not prove unrecorded late paths absent.2m spacing is coarse morphology, no spatial-Nyquist/mesh/PML acceptance.')
    save_json(out/'analysis_report.json',report)
    if c['stage']=='pilots':
        save_json(out/'pilot_review.template.json',dict(status='REQUIRES_REVIEW',
            verification_sha256=digest(study/'completed_verification.json'),
            analysis_report_path=str((out/'analysis_report.json').resolve()),analysis_report_sha256=digest(out/'analysis_report.json'),
            native_and_processing='',record_tail='',boundary_and_mesh_limits='',resource_and_ETA='',
            instruction='Review the real panels, traces, tail sensitivity, elapsed/peak resources. Explain limitations; only then write pilot capsule/pilot_review.json APPROVE_PREVIEW. No claim of physical acceptance.'))
    print(f'Native/SFCW identity PASS: {len(rows)} traces,501complex tones,2windows; review physics separately')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.study.resolve(),a.out.resolve())
