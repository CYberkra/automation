"""Offline exact SFCW, causal convolution and multi-station material contrasts."""
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
from diagnose_line9_layer_kinematics import primary_paths
from trace_line9_time_origin import correlation


def main(a):
    if a.out.exists() or a.numerical.exists():raise ValueError('Fresh results required')
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=c['study_manifest']
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    previous=read(a.prior/'completed_verification.json');paths=[a.prior/'H1.h5',a.prior/'H0.h5']+[a.source/(g['id']+'.h5') for g in v['groups']]
    names=['c0008_H1','c0008_H0']+[g['id'] for g in v['groups']]
    hashes=[g['raw_sha256'] for g in previous['groups']]+[g['native_sha256'] for g in v['groups']]
    assert [g['id'] for g in c['groups']]==[g['id'] for g in v['groups']]
    data={};sources={};responses=[];dft=[];support=[]
    for name,path,digest in zip(names,paths,hashes):
        assert sha(path)==digest
        src=sf.load_source(path);rx=sf.load_receiver(path,'/rxs/rx1','Ez');assert src.spatial_scale==.025
        with h5py.File(path) as h:assert h['rxs/rx1/Ez'].dtype==np.float64
        q=sf.direct_frequency_response(src,rx,FREQ,tail_taper_fraction=0);assert q.source_valid.all()
        response=np.asarray(q.response).reshape(501)/src.spatial_scale;k=np.unique(np.r_[np.arange(0,501,61),500])
        ss=src.dt*(np.exp(-2j*np.pi*FREQ[k,None]*src.times)@src.samples)
        manual=rx.dt*(np.exp(-2j*np.pi*FREQ[k,None]*rx.times)@rx.samples)/ss/.025
        error=rel(response[k],manual);assert error<1e-9
        responses.append(response);dft.append(error);data[name]=rx.samples;sources[name]=src.samples
        support.append(float(abs(q.source_spectrum).min()/abs(q.source_spectrum).max()))
    response=np.column_stack(responses);ix={name:j for j,name in enumerate(names)}
    # These comparisons are native sequences, without phase fitting or inverse processing.
    prefix={}
    for short,long in [('c0008_H0','impulse_1600'),('ricker_800','ricker_1600')]:
        n=len(data[short]);prefix[short+'_vs_'+long]=dict(bitwise_equal=bool(np.array_equal(data[short],data[long][:n])),relative_L2=rel(data[long][:n],data[short]))
    convolution={};predictions={}
    for name in ['ricker_800','ricker_1600']:
        n=len(data[name]);pred=fftconvolve(data['impulse_1600'],sources[name],mode='full')[:n]/40.
        predictions[name]=pred;dt=v['groups'][0]['dt_s'];mask=(np.arange(n)*dt>=250e-9)&(np.arange(n)*dt<=450e-9)
        convolution[name]=dict(full_native_relative_L2=rel(pred,data[name]),late_250_450ns_relative_L2=rel(pred[mask],data[name][mask]),normalization='Discrete convolution of native impulse40A sequence with native Ricker source, divided by40; no delay or coefficient fitting.')
    profiles={};metrics={};checks={}
    materials_path=a.parent/'H1/geometries/line9_research_materials_v1_smoothed.json'
    parent=read(a.parent/'manifest.json');assert sha(a.parent/'manifest.json')==m['parent_manifest_sha256'] and sha(materials_path)==parent['material_sha256']
    materials=read(materials_path)['materials']
    templates={s:primary_paths(g['geometry'],materials)[g['geometry']['boundaries'].index(g['geometry']['basal_sand'])] for s,g in m['stations'].items()}
    for window in ['hann','blackman']:
        w=weights(window,501);z,t=inverse(response,FREQ,w);profiles[window]=z
        take=np.arange(0,len(t),87);manual=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        checks[window]=rel(z[take],manual);assert checks[window]<1e-9
        record=dict(source_window={},stations={})
        for gate_name,(lo,hi) in [('basal',m['stations']['c0008']['basal_gate_ns']),('underground',[200,410]),('late',[600,800])]:
            mask=(t*1e9>=lo)&(t*1e9<=hi);base=z[mask,1];bn=np.linalg.norm(base);rows={}
            for name in ['ricker_800','impulse_1600','ricker_1600']:
                value=z[mask,ix[name]]
                rows[name]=dict(change_L2_over_impulse800=float(np.linalg.norm(value-base)/bn),L2_over_impulse800=float(np.linalg.norm(value)/bn),complex_correlation=correlation(value,base),envelope_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(value))]]*1e9))
            rows['factor_interaction_over_impulse800']=float(np.linalg.norm((z[:,ix['ricker_1600']]-z[:,ix['impulse_1600']]-z[:,ix['ricker_800']]+z[:,1])[mask])/bn)
            record['source_window'][gate_name]=dict(gate_ns=[lo,hi],models=rows)
        for station,g in m['stations'].items():
            lo,hi=g['basal_gate_ns'];mask=(t*1e9>=lo)&(t*1e9<=hi);j,k=ix[station+'_H1'],ix[station+'_H0'];delta=z[:,j]-z[:,k]
            template,_=inverse(templates[station],FREQ,w)
            row=dict(profile_midpoint_m=g['profile_midpoint_m'],midpoint_AGL_m=g['geometry']['midpoint_agl_m'],gate_ns=[lo,hi],
                total_template_correlation=correlation(z[mask,j],template[mask]),difference_template_correlation=correlation(delta[mask],template[mask]),
                H0_over_difference_L2=float(np.linalg.norm(z[mask,k])/np.linalg.norm(delta[mask])),difference_over_total_L2=float(np.linalg.norm(delta[mask])/np.linalg.norm(z[mask,j])),
                H0_difference_phase_deg=float(np.angle(np.vdot(z[mask,k],delta[mask]),deg=True)),difference_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(delta[mask]))]]*1e9))
            if station!='c0008':
                far=z[:,ix[station+'_far']];row.update(far_remaining_L2_over_H0=float(np.linalg.norm(far[mask])/np.linalg.norm(z[mask,k])),far_change_L2_over_H0=float(np.linalg.norm(far[mask]-z[mask,k])/np.linalg.norm(z[mask,k])))
            record['stations'][station]=row
        metrics[window]=record
    a.out.mkdir(parents=True);plt=plotting()
    fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    labels=['冲激800ns：复用原H0','Ricker100MHz/800ns','冲激1600ns','Ricker100MHz/1600ns'];cols=[1,2,3,4]
    for row,window in enumerate(['hann','blackman']):
        z=profiles[window]
        for col,(lo,hi) in enumerate([(200,410),(320,380)]):
            ax=axs[row,col];mask=(t*1e9>=lo)&(t*1e9<=hi)
            for j,label in zip(cols,labels):ax.plot(t*1e9,abs(z[:,j]),label=label,lw=1)
            ax.axvspan(*m['stations']['c0008']['basal_gate_ns'],alpha=.1,color='red')
            ax.set(xlim=(lo,hi),ylim=(0,abs(z[mask][:,cols]).max()*1.1),xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=window+'窗：源归一化后，同一绝对幅度；无移时或幅相拟合');ax.legend(fontsize=8)
    fig.savefig(a.out/'source_window_envelopes.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(2,4,figsize=(15,10),layout='constrained');z=profiles['hann']
    for row,(lo,hi) in enumerate([(0,800),(320,380)]):
        mask=(t*1e9>=lo)&(t*1e9<=hi);lim=abs(z[mask][:,cols].real).max()
        for ax,j,label in zip(axs[row],cols,labels):
            im=ax.imshow(z[mask,j].real[:,None],extent=[-.5,.5,hi,lo],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
            ax.axhline(m['stations']['c0008']['basal_template_peak_ns'],color='red',lw=1);ax.set(title=label,xticks=[0],xticklabels=['同一站位'],ylabel='时间 / ns')
        fig.colorbar(im,ax=axs[row],label='本行四图同一绝对灰度 / (V/m)/(A·m)')
    fig.suptitle('第8道H0 / AGL10.275m：源与时窗2×2控制，红线为原底砂预测；非空间B-scan')
    fig.savefig(a.out/'source_window_grayscale.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(13,7),layout='constrained')
    for ax,name in zip(axs,['ricker_800','ricker_1600']):
        n=len(data[name]);tt=np.arange(n)*dt*1e9;ax.plot(tt,data[name],label='独立Ricker正演 / 原生Ez');ax.plot(tt,predictions[name],'--',label='冲激正演与实际源的离散卷积预测')
        mask=(tt>=250)&(tt<=450);lim=max(abs(data[name][mask]).max(),abs(predictions[name][mask]).max())
        ax.set(xlim=(250,450),ylim=(-1.1*lim,1.1*lim),xlabel='仿真时间 / ns（包含Ricker源延迟）',ylabel='原生Ez / V/m',title=name+'：验证离散线性响应；不是SFCW图');ax.legend()
    fig.savefig(a.out/'native_convolution_check.png',dpi=145);plt.close(fig)
    stations=['c0008','c0058','c0108'];fig,axs=plt.subplots(3,4,figsize=(15,13),layout='constrained');mask=(t*1e9>=200)&(t*1e9<=420)
    lim=max(abs(profiles['hann'][mask,ix[s+'_H1']].real).max() for s in stations)
    for row,station in enumerate(stations):
        z=profiles['hann'];j,k=ix[station+'_H1'],ix[station+'_H0'];signals=[z[:,j],z[:,k],z[:,j]-z[:,k],None if station=='c0008' else z[:,ix[station+'_far']]]
        for ax,value,label in zip(axs[row],signals,['原H1','去底砂H0','底砂复数差场H1−H0','H0中只去x45–60m上覆夹层']):
            if value is None:ax.set_axis_off();ax.text(.1,.5,'原站的局部消融见上一批；\n此次不重复求解',transform=ax.transAxes);continue
            im=ax.imshow(value[mask].real[:,None],extent=[-.5,.5,420,200],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
            ax.axhline(m['stations'][station]['basal_template_peak_ns'],color='red',lw=1)
            ax.set(title=f"{station} / 剖面{m['stations'][station]['profile_midpoint_m']:.2f}m\n{label}",xticks=[0],xticklabels=['该独立站位'],ylabel='时间 / ns')
    fig.colorbar(im,ax=axs,label='全部图共同绝对灰度 / (V/m)/(A·m)',shrink=.7)
    fig.suptitle('三站有界配对：站距10m；原站复用、其余两站新算\n红线为各站底砂预测，图间其他站位未新增配对，不能当稠密新B-scan')
    fig.savefig(a.out/'three_station_grayscale.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    for row,station in enumerate(stations):
        for col,window in enumerate(['hann','blackman']):
            ax=axs[row,col];z=profiles[window];j,k=ix[station+'_H1'],ix[station+'_H0'];center=m['stations'][station]['basal_template_peak_ns'];mask=abs(t*1e9-center)<=45
            signals=[(z[:,j],'原H1'),(z[:,k],'去底砂H0'),(z[:,j]-z[:,k],'底砂差场')]
            if station!='c0008':signals.append((z[:,ix[station+'_far']],'只去ROI上覆夹层'))
            for value,label in signals:ax.plot(t*1e9,abs(value),label=label,lw=1)
            ax.axvspan(*m['stations'][station]['basal_gate_ns'],alpha=.1,color='red');ax.set(xlim=(center-45,center+45),ylim=(0,max(abs(value[mask]).max() for value,label in signals)*1.1),xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=f"{station} / 剖面{m['stations'][station]['profile_midpoint_m']:.2f}m / {window}；各面板独立绝对纵轴");ax.legend(fontsize=8)
    fig.savefig(a.out/'three_station_envelopes.png',dpi=145);plt.close(fig)
    if a.cache and a.review:
        review=read(a.review)
        with np.load(a.cache) as cache:
            r=cache['response'];chain=cache['chainage_m'];np.testing.assert_array_equal(cache['frequency_Hz'],FREQ)
            np.testing.assert_array_equal(cache['ids'],[row['id'] for row in review['records']])
            np.testing.assert_array_equal(chain,[row['chainage_m'] for row in review['records']])
        assert r.shape==(501,142) and len(review['records'])==142
        zz,tt=inverse(r,FREQ,weights('hann',501));mask=(tt*1e9>=150)&(tt*1e9<=450);lim=float(np.percentile(abs(zz[mask].real),99.5))
        fig,ax=plt.subplots(figsize=(13,6),layout='constrained');im=ax.imshow(zz[mask].real,extent=[chain[0],chain[-1],450,150],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
        centers=[]
        for row in review['records']:
            g=row['geometry'];p=primary_paths(g,materials)[g['boundaries'].index(g['basal_sand'])];v0,_=inverse(p,FREQ,weights('hann',501));centers.append(float(tt[np.argmax(abs(v0))]*1e9))
        ax.plot(chain,centers,color='cyan',lw=1,label='局部底砂一次反射近似，非实际检出')
        for station,g in m['stations'].items():ax.axvline(g['profile_midpoint_m'],color='red',lw=.8)
        ax.set(xlabel='剖面里程 / m（包内chainage；采集方向大→小）',ylabel='时间 / ns',title='原142道FP32总场：Hann/无AGC，无新插值；红线标此次三站FP64配对位置\n新配对仅三站，其余站位未计算H0/ROI消融');ax.legend();fig.colorbar(im,ax=ax,label='全图共同99.5百分位截色，仅显示裁切')
        fig.savefig(a.out/'existing_bscan_three_stations.png',dpi=145);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('time_s',data=t);h.create_dataset('response',data=response)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']=','.join(names);h.attrs['units']='(V/m)/(A*m)'
    report=dict(script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),
        columns=names,native_sha256=hashes,frequency_MHz=[20,170],step_MHz=.3,tones=501,AGC=False,tail_taper=False,independent_DFT_relative_L2=dft,independent_inverse_relative_L2=checks,source_band_min_over_max=support,native_prefix_checks=prefix,native_convolution_checks=convolution,metrics=metrics,limits=m['limits'],
        existing_Bscan_cache_sha256=sha(a.cache) if a.cache else None,existing_Bscan_review_sha256=sha(a.review) if a.review else None)
    save(a.out/'analysis.json',report);print(json.dumps(dict(prefix=prefix,convolution=convolution,metrics=metrics),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','prior','parent','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    for key in ['cache','review']:p.add_argument('--'+key,type=Path)
    main(p.parse_args())
