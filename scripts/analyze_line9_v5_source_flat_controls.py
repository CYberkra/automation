"""Exact SFCW, source convolution, and unfitted planar comparisons for v5."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import fftconvolve
from gprMax.toolboxes.SFCW import processing as sf
from audit_line9_postprocessing import FREQ,inverse,weights,rel
from hs_capsule_identity import sha256 as sha
from analyze_line9_basal_pair import plotting,save
from trace_line9_time_origin import correlation


def comparison(actual,model):
    return dict(unfitted_relative_L2=float(np.linalg.norm(actual-model)/np.linalg.norm(actual)),actual_over_model_L2=float(np.linalg.norm(actual)/np.linalg.norm(model)),inner_phase_deg=float(np.angle(np.vdot(model,actual),deg=True)),complex_correlation=correlation(actual,model))


def main(a):
    assert not a.out.exists() and not a.numerical.exists(),'Fresh output required'
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=c['study_manifest']
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.prepared/'manifest.json')==c['file_identities'][c['package']+'\\manifest.json']
    pv=read(a.prior/'completed_verification.json')
    names=['nonflat_impulse_H1','nonflat_impulse_H0']+[g['id'] for g in c['groups']]
    assert names[2:]==['nonflat_ricker_H0','flat_ricker_H0','flat_ricker_H1']
    paths=[a.prior/'base_H1.h5',a.prior/'base_H0.h5']+[a.source/(name+'.h5') for name in names[2:]]
    hashes=[g['native_sha256'] for g in pv['groups'][:2]]+[g['native_sha256'] for g in v['groups']]
    responses=[];native=[];sources=[];DFT=[]
    for p,digest in zip(paths,hashes):
        assert sha(p)==digest
        with h5py.File(p) as h:assert h['rxs/rx1/Ez'].dtype==np.float64 and h.attrs['dt']==m['dt_s']
        s=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez');q=sf.direct_frequency_response(s,rx,FREQ,tail_taper_fraction=0)
        assert q.source_valid.all() and s.spatial_scale==.025
        response=q.response.reshape(501)/.025;k=np.unique(np.r_[np.arange(0,501,67),500])
        manual=rx.dt*(np.exp(-2j*np.pi*FREQ[k,None]*rx.times)@rx.samples)/(s.dt*(np.exp(-2j*np.pi*FREQ[k,None]*s.times)@s.samples))/.025
        error=rel(response[k],manual);assert error<1e-9;DFT.append(error)
        responses.append(response);native.append(rx.samples);sources.append(s.samples)
    np.testing.assert_array_equal(sources[2],sources[3]);np.testing.assert_array_equal(sources[3],sources[4])
    pred=fftconvolve(native[1],sources[2],mode='full')[:len(native[2])]/40.
    native_error=rel(pred,native[2]);assert native_error<1e-9
    response=np.column_stack(responses)
    response=np.column_stack([response,response[:,0]-response[:,1],response[:,4]-response[:,3],response[:,2]-response[:,1],response[:,3]-response[:,2]])
    names+=['nonflat_impulse_delta','flat_ricker_delta','source_change_H0','geometry_change_H0']
    with h5py.File(a.planar) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);planar=h['response'][:]
    assert planar.shape==(501,6)
    profiles={};planars={};metrics={};inverse_errors={}
    for window in ['hann','blackman']:
        w=weights(window,501);z,t=inverse(response,FREQ,w);pz,_=inverse(planar,FREQ,w);profiles[window]=z;planars[window]=pz
        take=np.arange(3,len(t),131);manual=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        err=rel(z[take],manual);assert err<1e-9;inverse_errors[window]=err;rows={}
        for gate,bounds in [('basal',m['basal_gate_ns']),('deep',m['deep_gate_ns']),('early',m['early_gate_ns']),('late',[1000,1200])]:
            mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);q=z[mask];p=pz[mask];dn=np.linalg.norm(q[:,5]);pn=np.linalg.norm(q[:,6])
            row=dict(gate_ns=bounds,source_change_over_nonflat_H0=float(np.linalg.norm(q[:,7])/np.linalg.norm(q[:,1])),source_change_over_nonflat_delta=float(np.linalg.norm(q[:,7])/dn),geometry_change_over_nonflat_H0=float(np.linalg.norm(q[:,8])/np.linalg.norm(q[:,2])),flat_H0_over_nonflat_H0=float(np.linalg.norm(q[:,3])/np.linalg.norm(q[:,2])),flat_H0_over_flat_delta=float(np.linalg.norm(q[:,3])/pn),flat_delta_over_nonflat_delta=float(pn/dn),flat_vs_nonflat_delta=comparison(q[:,6],q[:,5]),flat_delta_vs_planar=comparison(q[:,6],p[:,5]),flat_H0_vs_planar=comparison(q[:,3],p[:,:5].sum(axis=1)),nonflat_H0_vs_planar=comparison(q[:,2],p[:,:5].sum(axis=1)))
            if gate=='basal':
                row.update(flat_H0_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(q[:,3]))]]*1e9),flat_delta_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(q[:,6]))]]*1e9),flat_H0_vs_cover_twice=comparison(q[:,3],p[:,3]),flat_H0_minus_other_planar_components_vs_cover_twice=comparison(q[:,3]-p[:,[0,1,2,4]].sum(axis=1),p[:,3]))
            rows[gate]=row
        metrics[window]=rows
    k=np.array([0,83,250,417,500]);ratio=response[k,6]/planar[k,5]
    frequency_comparison=dict(frequency_MHz=(FREQ[k]/1e6).tolist(),flat_bottom_over_planar_complex_ratio_real=ratio.real.tolist(),flat_bottom_over_planar_complex_ratio_imag=ratio.imag.tolist(),amplitude_ratio=abs(ratio).tolist(),phase_deg=np.angle(ratio,deg=True).tolist(),note='Pointwise diagnostic, no fitted amplitude/phase/delay applied. Tiny spectral model values may be sensitive.')
    a.out.mkdir(parents=True);plt=plotting()
    labels=['原非平冲激H1','原非平冲激H0','非平Ricker H0','平层Ricker H0','平层Ricker H1','非平底砂差场','平层底砂差场','H0换源变化','H0拉平变化']
    for window,z in profiles.items():
        fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
        for ax,cols,bounds,title in [(axes[0,0],[0,1,2,3,4],[150,450],'总场/剩余场：五配置共用绝对灰度'),(axes[0,1],[0,1,2,3,4],m['basal_gate_ns'],'预冻结底砂窗：五配置共用绝对灰度'),(axes[1,0],[5,6,7],m['basal_gate_ns'],'底砂差场与换源差：以非平底砂差场为共同标尺，超出裁色')]:
            mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);data=z[mask][:,cols].real
            limit=float(abs(z[mask,5].real).max()) if cols==[5,6,7] else float(abs(data).max())
            im=ax.imshow(data,extent=[-.5,len(cols)-.5,bounds[1],bounds[0]],aspect='auto',cmap='gray',vmin=-limit,vmax=limit)
            ax.set_xticks(range(len(cols)),[labels[j] for j in cols],rotation=12,fontsize=8);ax.set(title=title,ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='双极 / (V/m)/(A·m)')
        ax=axes[1,1];mask=(t*1e9>=340)&(t*1e9<=405)
        for j in [1,2,3,5,6,7]:ax.plot(t[mask]*1e9,abs(z[mask,j]),lw=1,label=labels[j])
        ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.08);ax.set(xlabel='时间 / ns',ylabel='包络 / (V/m)/(A·m)',title='共同绝对幅度，无AGC/移时/幅相拟合');ax.legend(fontsize=8)
        fig.suptitle(f'v5 / 剖面198.6m / AGL8m / {window} / 精确501点20–170MHz\n同站诊断配置列，非空间B-scan；新算三项原生FP64，其余复用')
        fig.savefig(a.out/f'v5_source_flat_{window}_gray.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for row,window in enumerate(['hann','blackman']):
        z=profiles[window];p=planars[window];mask=(t*1e9>=340)&(t*1e9<=405)
        for value,label in [(z[:,6],'平层FDTD底砂差场'),(p[:,5],'平层二维解析底砂对比'),(z[:,5],'非平冲激底砂差场')]:axes[row,0].plot(t[mask]*1e9,abs(value[mask]),lw=1,label=label)
        for value,label in [(z[:,3],'平层FDTD H0'),(p[:,:5].sum(axis=1),'平层二维解析完整H0'),(p[:,3],'解析覆盖层第二次往返'),(z[:,2],'非平Ricker H0')]:axes[row,1].plot(t[mask]*1e9,abs(value[mask]),lw=1,label=label)
        for ax in axes[row]:ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.08);ax.set(title=window+'：原幅度、原相位，不拟合',xlabel='SFCW时间 / ns',ylabel='包络 / (V/m)/(A·m)');ax.legend(fontsize=8)
    fig.suptitle('连续平层二维解析与离散有限域FDTD直接核对；材料均为研究假设')
    fig.savefig(a.out/'v5_flat_unfitted_planar.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4),layout='constrained');tt=np.arange(len(pred))*m['dt_s']*1e9;mask=(tt>=340)&(tt<=420)
    ax.plot(tt[mask],native[2][mask],label='非平Ricker独立正演',lw=1);ax.plot(tt[mask],pred[mask],'--',label='旧冲激原生序列与Ricker源离散卷积/40',lw=1)
    ax.set(xlabel='原生时间 / ns（含源延迟）',ylabel='原生Ez / V/m',title='线性卷积核对，非SFCW或新空间测线');ax.legend();fig.savefig(a.out/'v5_native_convolution.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        for key,value in [('frequency_Hz',FREQ),('time_s',t),('response',response),('planar_response',planar)]:h.create_dataset(key,data=value)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']=','.join(names);h.attrs['units']='(V/m)/(A*m)'
    save(a.out/'analysis.json',dict(script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),prepared_manifest_sha256=sha(a.prepared/'manifest.json'),prior_verification_sha256=sha(a.prior/'completed_verification.json'),native_sha256=hashes,columns=names,numerical_sha256=sha(a.numerical),planar_numerical_sha256=sha(a.planar),tones=501,frequency_MHz=[20,170],step_MHz=.3,AGC=False,tail_taper=False,native_convolution_relative_L2=native_error,independent_DFT_relative_L2=DFT,independent_inverse_relative_L2=inverse_errors,metrics=metrics,frequency_comparison=frequency_comparison,limits=m['limits'],extra_limits='Early/late bottom-difference denominators can be tiny; ratios are not physical energy fractions or certified error floors. Nonflat basal pair retains impulse source; source-only control is H0.'))
    print(json.dumps(dict(native_convolution_relative_L2=native_error,metrics=metrics,frequency_comparison=frequency_comparison)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','prepared','prior','planar','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
