"""Exact complex local-ablation and same-domain height-pair diagnostics; no solver."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from audit_line9_postprocessing import FREQ,inverse,weights,rel
from analyze_line9_basal_pair import plotting,save
from hs_capsule_identity import sha256 as sha
from trace_line9_time_origin import correlation
from review_line9_result_packages import geometry_at,indices
from diagnose_line9_layer_kinematics import primary_paths


def main(a):
    if a.out.exists() or a.numerical.exists():raise ValueError('Fresh outputs required')
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=read(a.parent/'manifest.json')
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.parent/'manifest.json')==c['study_manifest']['parent_manifest_sha256']
    old=read(a.prior/'completed_verification.json')
    paths=[a.prior/'H1.h5',a.prior/'H0.h5']+[a.source/(g['id']+'.h5') for g in v['groups']]
    hashes=[g['raw_sha256'] for g in old['groups']]+[g['native_sha256'] for g in v['groups']]
    names=['high_H1','high_H0']+[g['id'] for g in v['groups']]
    assert names==['high_H1','high_H0','far_interbed_removed','near_interbed_removed','low_H1','low_H0']
    responses=[];dft=[]
    for path,expected_hash in zip(paths,hashes):
        assert sha(path)==expected_hash
        with h5py.File(path) as h:assert h['rxs/rx1/Ez'].dtype==np.float64
        src=sf.load_source(path);rx=sf.load_receiver(path,'/rxs/rx1','Ez')
        assert src.spatial_scale==.025
        p=sf.direct_frequency_response(src,rx,FREQ,tail_taper_fraction=0);assert p.source_valid.all()
        response=np.asarray(p.response).reshape(501)/src.spatial_scale;k=np.arange(0,501,71)
        manual=(np.exp(-2j*np.pi*FREQ[k,None]*rx.times)@rx.samples)/(np.exp(-2j*np.pi*FREQ[k,None]*src.times)@src.samples)/.025
        error=rel(response[k],manual);assert error<1e-9
        responses.append(response);dft.append(error)
    response=np.column_stack(responses);db=a.parent/'H1/geometries/line9_research_materials_v1_smoothed.json'
    assert sha(db)==m['material_sha256'];materials=read(db)['materials']
    assert sha(a.parent/'H1/geometries/full2d_compact.h5')==m['groups'][0]['geometry_sha256']
    with h5py.File(a.parent/'H1/geometries/full2d_compact.h5') as h:geo=h['data'][:,:,0]
    geometries={height:geometry_at(geo,np.array([.025]*3),np.array([78.6,y,0.]),np.array([79.9,y,0.]),indices(materials,95e6)) for height,y in [('high',35.),('low',25.725)]}
    templates={height:primary_paths(g,materials)[g['boundaries'].index(g['basal_sand'])] for height,g in geometries.items()}
    profiles={};metrics={};inverse_errors={}
    for window in ['hann','blackman']:
        w=weights(window,501);z,t=inverse(response,FREQ,w);profiles[window]=z
        take=np.arange(0,len(t),91);manual=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        inverse_errors[window]=rel(z[take],manual);assert inverse_errors[window]<1e-9
        lo,hi=c['study_manifest']['predefined_gates_ns']['high_basal'];gate=(t*1e9>=lo)&(t*1e9<=hi);base=z[gate,1]
        metrics[window]={'localized':{}}
        for j,name in [(2,'far_interbed_removed'),(3,'near_interbed_removed')]:
            metrics[window]['localized'][name]=dict(remaining_L2_over_high_H0=float(np.linalg.norm(z[gate,j])/np.linalg.norm(base)),change_L2_over_high_H0=float(np.linalg.norm(z[gate,j]-base)/np.linalg.norm(base)),correlation_with_high_H0=correlation(z[gate,j],base))
        for height,j in [('high',0),('low',4)]:
            template,_=inverse(templates[height],FREQ,w);center=float(t[np.argmax(abs(template))]*1e9);mask=abs(t*1e9-center)<=12;delta=z[:,j]-z[:,j+1]
            metrics[window][height]=dict(template_center_ns=center,gate_ns=[center-12,center+12],difference_peak_ns=float(t[np.argmax(np.where(mask,abs(delta),0))]*1e9),
                total_template_correlation=correlation(z[mask,j],template[mask]),difference_template_correlation=correlation(delta[mask],template[mask]),
                H0_over_difference_L2=float(np.linalg.norm(z[mask,j+1])/np.linalg.norm(delta[mask])),difference_over_total_L2=float(np.linalg.norm(delta[mask])/np.linalg.norm(z[mask,j])),
                H0_difference_inner_product_phase_deg=float(np.angle(np.vdot(z[mask,j+1],delta[mask]),deg=True)))
        metrics[window]['low_H0_same_air_shift_proxy_gate']=dict(gate_ns=c['study_manifest']['predefined_gates_ns']['low_basal_air_shift_proxy'])
    a.out.mkdir(parents=True);plt=plotting();labels=['高航高原H0','只去远处x45–60m夹层','只去近处x65–90m夹层']
    fig,axs=plt.subplots(2,1,figsize=(12,8),layout='constrained')
    for ax,window in zip(axs,['hann','blackman']):
        z=profiles[window]
        for j,label in zip([1,2,3],labels):ax.plot(t*1e9,abs(z[:,j]),label=label)
        ax.axvspan(*c['study_manifest']['predefined_gates_ns']['high_basal'],color='red',alpha=.1)
        mask=(t*1e9>=300)&(t*1e9<=380)
        ax.set(xlim=(300,380),ylim=(0,abs(z[mask,1:4]).max()*1.1),xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=window+'：同一高航高站位；绝对幅度，无移时/归一化');ax.legend()
    fig.savefig(a.out/'localized_interbed_envelopes.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(12,9),layout='constrained');z=profiles['hann']
    for row,(lo,hi) in enumerate([(0,800),(300,380)]):
        mask=(t*1e9>=lo)&(t*1e9<=hi);lim=abs(z[mask,1:4].real).max()
        for ax,j,label in zip(axs[row],[1,2,3],labels):
            im=ax.imshow(z[mask,j].real[:,None],extent=[-.5,.5,hi,lo],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
            ax.axhline(344.31137724550894,color='red',lw=1);ax.set(title=label,xticks=[0],xticklabels=['同一站位'],ylabel='时间 / ns')
        fig.colorbar(im,ax=axs[row],label='本行三图共同绝对色标 / (V/m)/(A·m)')
    fig.suptitle('第8道高航高局部夹层替换：各列为配置，非不同空间测点；红线为原底砂预测')
    fig.savefig(a.out/'localized_interbed_grayscale.png',dpi=145);plt.close(fig)
    z=profiles['hann'];signals=[z[:,0],z[:,1],z[:,0]-z[:,1],z[:,4],z[:,5],z[:,4]-z[:,5]]
    labels=['高航高10.275m：原H1','高航高10.275m：去底砂H0','高航高：复数差场H1−H0','低航高1m：原H1','低航高1m：去底砂H0','低航高：复数差场H1−H0']
    fig,axs=plt.subplots(2,6,figsize=(22,10),layout='constrained')
    for row,(lo,hi) in enumerate([(0,800),(230,380)]):
        mask=(t*1e9>=lo)&(t*1e9<=hi);lim=max(abs(s[mask].real).max() for s in signals)
        for j,ax in enumerate(axs[row]):
            im=ax.imshow(signals[j][mask].real[:,None],extent=[-.5,.5,hi,lo],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
            center=metrics['hann']['high' if j<3 else 'low']['template_center_ns']
            ax.axhline(center,color='red',lw=1);ax.set(title=labels[j],xticks=[0],xticklabels=['同一水平站位'],ylabel='时间 / ns')
        fig.colorbar(im,ax=axs[row],label='本行六图共同绝对色标 / (V/m)/(A·m)')
    fig.suptitle('相同110×40m域/原地质/2.5cm/V4 FP64：仅航高改变；不是完整B-scan\n高航高复用已完成两模型，低航高新算；红线为各自独立局部层序预测')
    fig.savefig(a.out/'same_domain_height_grayscale.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for row,window in enumerate(['hann','blackman']):
        z=profiles[window]
        for col,(height,j) in enumerate([('high',0),('low',4)]):
            ax=axs[row,col];center=metrics[window][height]['template_center_ns'];mask=abs(t*1e9-center)<=35
            for s,label in [(z[:,j],'原H1'),(z[:,j+1],'去底砂H0'),(z[:,j]-z[:,j+1],'底砂复数差场')]:ax.plot(t*1e9,abs(s),label=label)
            ax.axvspan(*metrics[window][height]['gate_ns'],alpha=.12,color='red');ax.set(xlim=(center-35,center+35),ylim=(0,1.1*max(abs(z[mask,j]).max(),abs(z[mask,j+1]).max(),abs((z[:,j]-z[:,j+1])[mask]).max())),
                xlabel='时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=f'{window} / '+('AGL10.275m' if height=='high' else 'AGL1m')+'；独立纵轴刻度，未归一化');ax.legend()
    fig.savefig(a.out/'same_domain_height_envelopes.png',dpi=145);plt.close(fig)
    from matplotlib.patches import Rectangle
    fig,ax=plt.subplots(figsize=(13,5),layout='constrained')
    from matplotlib.colors import ListedColormap,BoundaryNorm
    cmap=ListedColormap(['#eef7fc','#e5c36a','#bc7a9b','#b5b5b5']);norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4)
    im=ax.imshow(geo[::4,::4].T,origin='lower',extent=[0,110,0,40],aspect='equal',interpolation='nearest',cmap=cmap,norm=norm)
    cb=fig.colorbar(im,ax=ax,ticks=[0,1,2,3],shrink=.8);cb.ax.set_yticklabels(['空气','粉质黏土覆盖层','泥岩','砂岩'])
    for x,width,color,label in [(45,15,'red','远处替换ROI：只换其中上覆砂岩'),(65,25,'blue','近处替换ROI：只换其中上覆砂岩')]:ax.add_patch(Rectangle((x,0),width,40,fill=False,edgecolor=color,label=label))
    for y,label in [(35,'高航高收发'),(25.725,'1m收发')]:ax.plot([78.6,79.9],[y,y],'k.',label=label)
    ax.set(xlabel='局部x / m（剖面X=局部x+100m）',ylabel='局部y / m',title='原H1地质：空气/覆盖层/泥岩/砂岩ID0/1/2/3；替换对照使用已去底砂H0');ax.legend(fontsize=8)
    fig.savefig(a.out/'localized_geometry.png',dpi=145);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('time_s',data=t);h.create_dataset('response',data=response)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']=','.join(names);h.attrs['units']='(V/m)/(A*m)'
    report=dict(script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
        columns=names,native_sha256=hashes,numerical_sha256=sha(a.numerical),independent_DFT_relative_L2=dft,independent_inverse_relative_L2=inverse_errors,
        frequency_MHz=[20,170],step_MHz=.3,tones=501,tail_taper=False,AGC=False,geometries=geometries,metrics=metrics,limits=c['study_manifest']['limits'])
    save(a.out/'analysis.json',report);print(json.dumps(metrics,ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','prior','parent','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
