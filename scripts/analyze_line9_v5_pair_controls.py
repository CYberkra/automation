"""Exact SFCW and FP32/top-boundary contrasts for completed v5 basal pairs."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from audit_line9_postprocessing import FREQ,inverse,weights,rel
from hs_capsule_identity import sha256 as sha
from analyze_line9_basal_pair import plotting,save
from diagnose_line9_layer_kinematics import primary_paths
from diagnose_line9_layer_stack import stack_response,admittance_response
from review_line9_result_packages import indices
from trace_line9_time_origin import correlation


def main(a):
    assert not a.out.exists() and not a.numerical.exists(),'Fresh outputs required'
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=c['study_manifest']
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json') and sha(a.prepared/'manifest.json')==c['file_identities'][c['package']+'\\manifest.json']
    groups=[*c.get('reused_groups',[]),*c['groups']]
    names=[g['id'] for g in groups]+['old_FP32'];assert names[:4]==['base_H1','base_H0','top20_H1','top20_H0']
    paths=[a.source/(name+'.h5') for name in names[:4]]+[a.prepared/'reference_fp32.h5']
    digests=[g['native_sha256'] for g in v['groups']]+[m['source_native_fp32_sha256']]
    responses=[];errors=[]
    for j,(p,digest) in enumerate(zip(paths,digests)):
        assert sha(p)==digest
        source=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez');assert source.spatial_scale==.025
        with h5py.File(p) as h:assert h['rxs/rx1/Ez'].dtype==(np.float32 if j==4 else np.float64)
        q=sf.direct_frequency_response(source,rx,FREQ,tail_taper_fraction=0);assert q.source_valid.all();response=q.response.reshape(501)/.025
        k=np.r_[0,np.arange(17,501,71),500]
        manual=rx.dt*(np.exp(-2j*np.pi*FREQ[k,None]*rx.times)@rx.samples)/(source.dt*(np.exp(-2j*np.pi*FREQ[k,None]*source.times)@source.samples))/.025
        err=rel(response[k],manual);assert err<1e-9;errors.append(err);responses.append(response)
    response=np.column_stack(responses)
    db=a.prepared/'base_H1/geometries/line9_research_materials_v1_smoothed.json';assert sha(db)==m['material_sha256']
    materials=read(db)['materials'];geometry=m['geometry_diagnostic'];paths=primary_paths(geometry,materials);primary=paths[-1]
    n=np.array([indices(materials,f) for f in FREQ]);thickness=np.r_[geometry['midpoint_agl_m'],np.diff([x['depth_m'] for x in geometry['boundaries']]),0.]
    full1,single1=stack_response(n,thickness,FREQ);full0,single0=stack_response(n[:,:3],thickness[:3],FREQ)
    scalar_check=max(rel(full1,admittance_response(n,thickness,FREQ)),rel(full0,admittance_response(n[:,:3],thickness[:3],FREQ)));assert scalar_check<1e-10
    correction=geometry['surface']['time95_ns']*1e-9-2*thickness[0]/299792458.
    phase=np.exp(-2j*np.pi*FREQ*correction)
    models=np.column_stack([primary,(full1-full0)*phase,(full0-single0)*phase])
    columns=names+['base_delta','top20_delta','FP64_minus_FP32','top20_minus_base_H1','top20_minus_base_H0','top20_minus_base_delta']
    response=np.column_stack([response,response[:,0]-response[:,1],response[:,2]-response[:,3],response[:,0]-response[:,4],response[:,2]-response[:,0],response[:,3]-response[:,1],(response[:,2]-response[:,3])-(response[:,0]-response[:,1])])
    planar=None
    if a.planar:
        with h5py.File(a.planar) as h:
            np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);planar=h['response'][:]
        assert planar.shape==(501,6)
    profiles={};model_profiles={};metrics={};inverse_errors={}
    for window in ['hann','blackman']:
        w=weights(window,501);z,t=inverse(response,FREQ,w);modes,_=inverse(models,FREQ,w);profiles[window]=z;model_profiles[window]=modes
        take=np.arange(0,len(t),103);manual=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        err=rel(z[take],manual);assert err<1e-9;inverse_errors[window]=err
        records={}
        for gate_name,bounds in [('basal',m['basal_gate_ns']),('deep',m['deep_gate_ns']),('early',m['early_gate_ns'])]:
            mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);b=z[mask];dn=np.linalg.norm(b[:,5]);h0n=np.linalg.norm(b[:,1]);h1n=np.linalg.norm(b[:,0])
            record=dict(gate_ns=bounds,base_H0_over_delta=float(h0n/dn),base_delta_over_H1=float(dn/h1n),
                top20_H0_over_delta=float(np.linalg.norm(b[:,3])/np.linalg.norm(b[:,6])),base_H0_delta_phase_deg=float(np.angle(np.vdot(b[:,1],b[:,5]),deg=True)),
                FP64_minus_FP32_over_base_delta=float(np.linalg.norm(b[:,7])/dn),FP64_minus_FP32_over_base_H1=float(np.linalg.norm(b[:,7])/h1n),
                top20_H1_change_over_base_H1=float(np.linalg.norm(b[:,8])/h1n),top20_H0_change_over_base_H0=float(np.linalg.norm(b[:,9])/h0n),
                top20_H1_change_over_base_delta=float(np.linalg.norm(b[:,8])/dn),top20_H0_change_over_base_delta=float(np.linalg.norm(b[:,9])/dn),
                top20_delta_change_over_base_delta=float(np.linalg.norm(b[:,10])/dn),delta_top_base_correlation=correlation(b[:,5],b[:,6]))
            if gate_name=='basal':
                record.update(base_H1_primary_correlation=correlation(b[:,0],modes[mask,0]),base_delta_primary_correlation=correlation(b[:,5],modes[mask,0]),top20_delta_primary_correlation=correlation(b[:,6],modes[mask,0]),
                    base_delta_full_stack_contrast_correlation=correlation(b[:,5],modes[mask,1]),base_H0_multiple_proxy_correlation=correlation(b[:,1],modes[mask,2]),
                    base_delta_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(b[:,5]))]]*1e9),top20_delta_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(b[:,6]))]]*1e9),
                    base_H0_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(b[:,1]))]]*1e9),multiple_proxy_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(modes[mask,2]))]]*1e9))
                if planar is not None:
                    pz,_=inverse(planar,FREQ,w);pd=pz[mask,5];p0=pz[mask,:5].sum(axis=1)
                    record.update(planar_bottom_unfitted_relative_L2=float(np.linalg.norm(b[:,5]-pd)/dn),actual_bottom_over_planar_L2=float(dn/np.linalg.norm(pd)),
                        planar_bottom_inner_phase_deg=float(np.angle(np.vdot(pd,b[:,5]),deg=True)),planar_bottom_correlation=correlation(b[:,5],pd),
                        planar_H0_unfitted_relative_L2=float(np.linalg.norm(b[:,1]-p0)/h0n),planar_H0_correlation=correlation(b[:,1],p0),planar_cover_twice_H0_correlation=correlation(b[:,1],pz[mask,3]))
            records[gate_name]=record
        metrics[window]=records
    a.out.mkdir(parents=True);plt=plotting()
    for window,z in profiles.items():
        fig,axes=plt.subplots(3,2,figsize=(14,13),layout='constrained')
        choices=[(axes[0,0],[0,1,2,3],0,600,'四个原生总场/剩余场：全时段共同绝对色标',None),
            (axes[0,1],[0,1,2,3],300,450,'深窗：四配置共同绝对色标',None),
            (axes[1,0],[0,1,2,3],*m['basal_gate_ns'],'底砂固定窗：四配置共同绝对色标',None),
            (axes[1,1],[5,6,7],*m['basal_gate_ns'],'底砂差场与精度差：以base差场为共同标尺，超出裁色',5)]
        for ax,ix,lo,hi,title,reference in choices:
            mask=(t*1e9>=lo)&(t*1e9<=hi);val=z[mask][:,ix].real;lim=float(abs(val).max()) if reference is None else float(abs(z[mask,reference].real).max())
            im=ax.imshow(val,extent=[-.5,len(ix)-.5,hi,lo],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
            labels={'base_H1':'原域H1','base_H0':'原域H0','top20_H1':'加20m空气H1','top20_H0':'加20m空气H0','base_delta':'原域H1−H0','top20_delta':'加空气H1−H0','FP64_minus_FP32':'原域H1：64−32位'}
            ax.set_xticks(range(len(ix)),[labels.get(columns[j],columns[j]) for j in ix],fontsize=9);ax.set(ylabel='SFCW时间 / ns',title=title)
            fig.colorbar(im,ax=ax,label='双极 / (V/m)/(A·m)')
        for ax,ix,title in [(axes[2,0],[0,1,2,3,5,6],'原场与底砂差场：共同绝对包络'),(axes[2,1],[5,6,7,10],'差场、精度差与边界差：共同绝对包络')]:
            show=(t*1e9>=330)&(t*1e9<=410)
            for j in ix:ax.plot(t[show]*1e9,abs(z[show,j]),lw=1,label=columns[j])
            ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.08);ax.axvline(m['basal_template_peak_ns'],color='black',ls=':',lw=.8)
            ax.set(xlabel='时间 / ns',ylabel='包络 / (V/m)/(A·m)',title=title);ax.legend(fontsize=8)
        fig.suptitle(f'v5第98道 / 剖面198.6m / AGL8m / {window}；同站配置列非连续测线\n精确501点20–170MHz/0.3MHz，原生FP64；旧FP32仅诊断；无AGC/去背景/拟合')
        fig.savefig(a.out/f'v5_pair_{window}_gray.png',dpi=140);plt.close(fig)
    # Show phase/shape comparison without fitting phase, delay or applying it to data.
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,window in zip(axes,['hann','blackman']):
        z=profiles[window];mp=model_profiles[window];show=(t*1e9>=340)&(t*1e9<=405)
        for value,label in [(z[:,5],'原域底砂差场'),(mp[:,0],'局部一次反射近似'),(mp[:,1],'完整1D层栈底砂对比'),(mp[:,2],'H0中1D多次候选')]:
            ax.plot(t[show]*1e9,abs(value[show])/max(abs(value[show])),label=label,lw=1)
        ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.08);ax.set(title=window+'：每曲线峰归一，仅形状诊断',xlabel='时间 / ns',ylabel='显示归一包络，非保幅');ax.legend(fontsize=9)
    fig.savefig(a.out/'v5_template_shapes.png',dpi=140);plt.close(fig)
    if planar is not None:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for k,window in enumerate(['hann','blackman']):
            z=profiles[window];pz,_=inverse(planar,FREQ,weights(window,501));show=(t*1e9>=330)&(t*1e9<=410)
            for value,label in [(z[:,5],'实际非平模型底砂差场'),(pz[:,5],'局部平层二维解析底砂对比')]:axes[k,0].plot(t[show]*1e9,abs(value[show]),label=label,lw=1)
            for value,label in [(z[:,1],'实际非平模型H0'),(pz[:,:5].sum(axis=1),'局部平层二维解析H0'),(pz[:,3],'解析覆盖层第二次往返')]:axes[k,1].plot(t[show]*1e9,abs(value[show]),label=label,lw=1)
            for j in range(2):
                axes[k,j].set(xlabel='SFCW时间 / ns',ylabel='包络 / (V/m)/(A·m)',title=window+'：同单位原幅度，无幅相/时间拟合');axes[k,j].axvspan(*m['basal_gate_ns'],color='red',alpha=.08);axes[k,j].legend(fontsize=9)
        fig.suptitle('解析平层二维与真实非平体素模型不同；高形状相关不等于幅相已验证')
        fig.savefig(a.out/'v5_planar_observed.png',dpi=140);plt.close(fig)
    if a.cache:
        review=read(a.review)
        with np.load(a.cache) as cache:
            np.testing.assert_array_equal(cache['frequency_Hz'],FREQ);np.testing.assert_array_equal(cache['ids'],[r['id'] for r in review['records']]);chain=cache['chainage_m'];old_response=cache['response']
        assert review['geometry_sha256']==m['original_geometry_sha256'] and old_response.shape==(501,194)
        z,tt=inverse(old_response,FREQ,weights('hann',501));show=(tt*1e9>=150)&(tt*1e9<=450);lim=float(np.percentile(abs(z[show].real),99.5));centers=[]
        for row in review['records']:
            prim=primary_paths(row['geometry'],materials)[-1];zz,_=inverse(prim,FREQ,weights('hann',501));centers.append(float(tt[np.argmax(abs(zz))]*1e9))
        fig,ax=plt.subplots(figsize=(13,6),layout='constrained');im=ax.imshow(z[show].real,extent=[chain[0],chain[-1],450,150],cmap='gray',vmin=-lim,vmax=lim,aspect='auto')
        ax.plot(chain,centers,color='cyan',lw=1,label='局部底砂一次反射近似，非实际检出');ax.axvline(m['profile_midpoint_m'],color='red',lw=1,label='此次FP64配对站位')
        ax.set(xlabel='剖面里程 / m（采集方向大→小）',ylabel='SFCW时间 / ns',title='既有v5/194道FP32总场：Hann/无AGC；仅红线站位新算四配置\n其他站位未做本批配对或边界控制，不插值H0');ax.legend();fig.colorbar(im,ax=ax,label='全图共同99.5百分位截色，显示裁切')
        fig.savefig(a.out/'existing_v5_bscan_station.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('time_s',data=t);h.create_dataset('response',data=response);h.create_dataset('template_response',data=models)
        if planar is not None:h.create_dataset('planar_response',data=planar)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']=','.join(columns);h.attrs['units']='(V/m)/(A*m)'
    result=dict(script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),
        columns=columns,native_sha256=digests,tones=501,frequency_MHz=[20,170],step_MHz=.3,AGC=False,tail_taper=False,independent_DFT_relative_L2=errors,independent_inverse_relative_L2=inverse_errors,
        independent_stack_admittance_relative_L2=scalar_check,metrics=metrics,basal_template_peak_ns=m['basal_template_peak_ns'],limits=m['limits'],
        extra_limits='FP64-FP32 is a total-field precision diagnostic, not paired FP32 difference error. Normal-incidence1D stack lacks2D spreading/lateral geometry. No inferred physical error floor/unique mechanism from correlation alone.',
        planar_numerical_sha256=sha(a.planar) if a.planar else None,existing_cache_sha256=sha(a.cache) if a.cache else None,existing_review_sha256=sha(a.review) if a.review else None)
    save(a.out/'analysis.json',result);print(json.dumps(dict(metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','prepared','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    for key in ['cache','review','planar']:p.add_argument('--'+key,type=Path)
    main(p.parse_args())
