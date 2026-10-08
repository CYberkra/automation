"""Exact SFCW comparison of frozen single-station boundary-distance controls; no solver."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from audit_line9_postprocessing import FREQ,inverse,weights,rel,sha,save
from trace_line9_time_origin import correlation
from analyze_line9_basal_pair import plotting


def main(source,reference,out,numerical,parent_manifest=None,materials=None):
    if out.exists() or numerical.exists():raise ValueError('Fresh products required')
    c=json.loads((source/'execution_contract.json').read_text('utf-8'))
    v=json.loads((source/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256']==sha(source/'execution_contract.json')
    groups=list(v['groups']);names=[g['id'] for g in groups]
    assert names==[g['id'] for g in c['groups']]
    reference_is_prior=bool(c['study_manifest'].get('reference_is_prior_H0'))
    if reference_is_prior:
        groups.insert(0,dict(id='base',native_sha256=sha(reference)))
        names.insert(0,'base')
    responses=[];native=[];errors=[]
    for g in groups:
        p=reference if reference_is_prior and g['id']=='base' else source/(g['id']+'.h5');assert sha(p)==g['native_sha256']
        src=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez')
        assert src.spatial_scale==.025
        product=sf.direct_frequency_response(src,rx,FREQ,tail_taper_fraction=0)
        assert product.source_valid.all()
        h=np.asarray(product.response).reshape(501)/src.spatial_scale
        k=np.array([0,83,167,250,333,417,500]);f=FREQ[k,None]
        expected=rx.dt*(np.exp(-2j*np.pi*f*rx.times)@rx.samples)/(src.dt*(np.exp(-2j*np.pi*f*src.times)@src.samples))/.025
        e=rel(h[k],expected);assert e<1e-9
        responses.append(h);native.append(rx.samples);errors.append(e)
    with h5py.File(reference) as ref:old=ref['rxs/rx1/Ez'][:]
    assert sha(reference)==c['study_manifest']['parent_H0_native_sha256']
    observer_error=rel(native[0],old) if not reference_is_prior else None
    if observer_error is not None and observer_error>1e-10:raise ValueError('No-observer repeat is not numerically matching reference')
    response=np.column_stack(responses);profiles={};checks={};metrics={}
    gates=c['study_manifest']['predefined_gates_ns']
    for window in ['hann','blackman']:
        w=weights(window,501);z,t=inverse(response,FREQ,w);profiles[window]=z
        take=np.arange(0,len(t),101);expected=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        e=rel(z[take],expected);assert e<1e-9;checks[window]=e
        metrics[window]={}
        for gate_name,(lo,hi) in gates.items():
            mask=(t*1e9>=lo)&(t*1e9<=hi);base=z[mask,0]
            item={}
            for j,name in enumerate(names):
                value=z[mask,j]
                item[name]=dict(L2_over_base=float(np.linalg.norm(value)/np.linalg.norm(base)),
                    change_L2_over_base=float(np.linalg.norm(value-base)/np.linalg.norm(base)),
                    complex_correlation_with_base=correlation(value,base),
                    envelope_peak_ns=float(t[np.argmax(np.where(mask,abs(z[:,j]),0))]*1e9))
            if not reference_is_prior:item['factor_interaction_L2_over_base']=float(np.linalg.norm((z[:,0]-z[:,1]-z[:,2]+z[:,3])[mask])/np.linalg.norm(base))
            metrics[window][gate_name]=item
    out.mkdir(parents=True);plt=plotting()
    labels=c['study_manifest'].get('case_labels',['原域：110×40m','只抬顶边界：110×80m','只移远侧边：190×40m','顶边与侧边均移远：190×80m'])
    fig,axs=plt.subplots(3,len(names),figsize=(3.8*len(names),10),layout='constrained');z=profiles['hann']
    for row,(lo,hi) in enumerate([(0,800),(200,410),(320,380)]):
        mask=(t*1e9>=lo)&(t*1e9<=hi);limit=float(np.max(abs(z[mask].real)))
        for j,ax in enumerate(axs[row]):
            im=ax.imshow(z[mask,j].real[:,None],extent=[-.5,.5,hi,lo],aspect='auto',cmap='gray',vmin=-limit,vmax=limit)
            ax.axhline(344.31137724550894,color='red',lw=.7);ax.set(title=labels[j],xticks=[0],xticklabels=['同一站位'],ylabel='时间 / ns')
        fig.colorbar(im,ax=axs[row],label='本行全部图共同色标 / (V/m)/(A·m)')
    title='地层/地表/全空气对照；基准复用原H0' if reference_is_prior else '边界距离对照：原域地质与源收相对位置不变'
    fig.suptitle('第8道 / X179.25m / AGL10.275m / V4 FP64：'+title+'\n各列是不同配置的同一站位，非不同空间道；红线为原底砂预测')
    fig.savefig(out/'boundary_controls_grayscale.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(13,10),layout='constrained')
    for ax,window,(lo,hi) in zip(axs,['hann','hann','blackman'],[(200,550),(320,380),(320,380)]):
        z=profiles[window];mask=(t*1e9>=lo)&(t*1e9<=hi)
        for j,label in enumerate(labels):ax.plot(t*1e9,abs(z[:,j]),label=label,lw=1.2)
        ax.axvspan(*gates['bottom'],color='red',alpha=.1)
        ax.set(xlim=(lo,hi),ylim=(0,float(abs(z[mask]).max())*1.08),xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=window+'窗；同一绝对幅度，不移时/不归一化');ax.legend(fontsize=8)
    fig.savefig(out/'boundary_controls_envelopes.png',dpi=145);plt.close(fig)
    with h5py.File(numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=response);h.create_dataset('time_s',data=t)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']=','.join(names);h.attrs['units']='(V/m)/(A*m)'
    proxy=None
    if parent_manifest and materials:
        from diagnose_line9_layer_stack import stack_response,admittance_response
        from review_line9_result_packages import indices,C0
        m=json.loads(parent_manifest.read_text('utf-8'));mat=json.loads(materials.read_text('utf-8'))['materials']
        assert sha(parent_manifest)==c['study_manifest']['parent_manifest_sha256'] and sha(materials)==m['material_sha256']
        g=m['geometry_diagnostic'];b=g['boundaries'][:-1]
        media=[p['above'] for p in b]+[b[-1]['below']]
        thick=np.r_[g['midpoint_agl_m'],np.diff([p['depth_m'] for p in b]),0.]
        n=np.array([indices(mat,f) for f in FREQ])[:,media]
        full,single=stack_response(n,thick,FREQ);other=admittance_response(n,thick,FREQ)
        phase=np.exp(-2j*np.pi*FREQ*(g['surface']['time95_ns']*1e-9-2*thick[0]/C0))
        z=profiles['hann'][:,0];v,_=inverse(full*phase,FREQ,weights('hann',501));s,_=inverse(single*phase,FREQ,weights('hann',501))
        shallow=abs(t*1e9-229.54091816367264)<12;deep=abs(t*1e9-344.31137724550894)<12
        coefficient=np.vdot(v[shallow],z[shallow])/np.vdot(v[shallow],v[shallow]);prediction=v*coefficient
        proxy=dict(independent_admittance_relative_L2=rel(full,other),shallow_correlation=correlation(v[shallow],z[shallow]),
            deep_correlation=correlation(v[deep],z[deep]),deep_predicted_over_actual_L2=float(np.linalg.norm(prediction[deep])/np.linalg.norm(z[deep])),
            deep_multiple_over_actual_L2=float(np.linalg.norm(((v-s)*coefficient)[deep])/np.linalg.norm(z[deep])),
            media=media,thickness_m=thick.tolist(),coefficient=[float(coefficient.real),float(coefficient.imag)],
            limitation='Local1D plane-wave proxy includes all layer multiples; shallow fit coefficient only. Omits2D spreading/angular spectrum/topography/boundaries; failure is not exclusion of all geological mechanisms.')
        fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
        ax.plot(t*1e9,abs(z),label='原域H0：实际二维FP64')
        ax.plot(t*1e9,abs(prediction),label='局部1D完整层栈：仅浅部拟合整体复比例')
        ax.plot(t*1e9,abs((v-s)*coefficient),label='该1D近似的层间多次分量')
        mask=(t*1e9>=320)&(t*1e9<=380)
        ax.set(xlim=(320,380),ylim=(0,float(abs(z[mask]).max())*1.1),xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title='局部1D理论核查：幅度来自浅部标定；非新正演、非生产去多次');ax.legend()
        fig.savefig(out/'local_layer_stack_proxy.png',dpi=145);plt.close(fig)
    report=dict(script_sha256=sha(__file__),calls_solver=False,contract_sha256=sha(source/'execution_contract.json'),verification_sha256=sha(source/'completed_verification.json'),
        reference_is_prior_H0=reference_is_prior,observer_repeat_native_relative_L2=observer_error,observer_repeat_bitwise_equal=bool(np.array_equal(native[0],old)) if not reference_is_prior else None,
        frequency_MHz=[20,170],step_MHz=.3,tones=501,tail_taper=False,AGC=False,source_spatial_scale_m=.025,
        native_dtype='float64',independent_DFT_relative_L2=errors,independent_inverse_relative_L2=checks,metrics=metrics,
        native_sha256=[g['native_sha256'] for g in groups],numerical_sha256=sha(numerical),
        limits=c['study_manifest']['limits'],gates_ns=gates,hardware=c['hardware_at_freeze'],local_H0_layer_stack_proxy=proxy)
    save(out/'analysis.json',report);print(json.dumps(metrics,ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','reference','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--parent-manifest',type=Path);p.add_argument('--materials',type=Path)
    a=p.parse_args();main(a.source,a.reference,a.out,a.numerical,a.parent_manifest,a.materials)
