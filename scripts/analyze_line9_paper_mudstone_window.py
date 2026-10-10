"""Official direct SFCW: old material1200, new material crop1200/full2400."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import gprMax
from gprMax.toolboxes.SFCW import processing as sf
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import audit_source, save

FREQ=20e6+np.arange(501)*300000.


def relative(x,y):
    den=np.linalg.norm(y)
    return float(np.linalg.norm(x-y)/den) if den else None


def transfer(path, crop=None):
    s=sf.load_source(path); r=sf.load_receiver(path,'/rxs/rx1','Ez')
    if crop:
        s=replace(s,samples=s.samples[:crop]); r=replace(r,samples=r.samples[:crop])
    assert s.spatial_scale==.025 and s.time_offset==.5*s.dt and r.time_offset==0
    q=sf.direct_frequency_response(s,r,FREQ,tail_taper_fraction=0)
    assert q.source_valid.all()
    independent=np.empty(501,complex)
    for k in range(0,501,16):
        f=FREQ[k:k+16,None]
        independent[k:k+16]=(r.dt*(np.exp(-2j*np.pi*f*r.times)@r.samples))/(
            s.dt*(np.exp(-2j*np.pi*f*s.times)@s.samples))/.025
    q=replace(q,response=q.response/.025)
    error=relative(q.response,independent); assert error<1e-9
    return q,error


def main(a):
    assert not a.out.exists() and not a.numerical.exists() and gprMax.__version__=='4.0.1'
    c=json.loads((a.new/'execution_contract.json').read_text('utf-8'))
    v=json.loads((a.new/'completed_verification.json').read_text('utf-8'))
    oc=json.loads((a.old/'execution_contract.json').read_text('utf-8'))
    ov=json.loads((a.old/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256']==sha(a.new/'execution_contract.json')
    assert ov['completed'] and ov['contract_sha256']==sha(a.old/'execution_contract.json')
    assert sha(sf.__file__)==c['source_identities']['toolboxes/SFCW/processing.py']
    oldg={g['id']:g for g in oc['groups']}; oldr={r['id']:r for r in ov['groups']}
    qs=[[],[],[]]; raw=[]; rows=[]
    for g,record in zip(c['groups'],v['groups']):
        name=g['id']; p=a.new/(name+'.h5'); op=a.old/(name+'.h5'); og=oldg[name]
        assert sha(p)==record['native_sha256'] and sha(op)==oldr[name]['native_sha256']
        assert g['parent_geometry_sha256']==og['geometry_sha256']
        assert g['parent_material_sha256']==og['material_sha256']
        for key in ['native_shape','tx_m','rx_m','dt_s']:
            assert g[key]==og[key],key
        with h5py.File(p) as h, h5py.File(op) as oh:
            for key in ['gprMax','nx_ny_nz','dx_dy_dz','dt']:
                np.testing.assert_array_equal(h.attrs[key],oh.attrs[key])
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==40703
            x=h['rxs/rx1/Ez'][:]; src=h['srcs/src1/excitation/samples'][:]
            assert x.dtype==src.dtype==np.float64 and x.shape==src.shape==(40703,)
            assert np.isfinite(x).all() and np.isfinite(src).all() and np.any(x)
            audit_source(h['srcs/src1/excitation'],g,float(h.attrs['dt']))
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
            raw.append((x,oh['rxs/rx1/Ez'][:],float(h.attrs['dt'])))
        errors=[]; tails=[]
        for j,(path,crop) in enumerate([(op,None),(p,g['comparison_crop_samples']),(p,None)]):
            q,e=transfer(path,crop); qs[j].append(q); errors.append(e); tails.append(q.receiver_tail_relative_db)
        # Additivity of DTFT: appended late samples must explain full-minus-crop.
        r=sf.load_receiver(p,'/rxs/rx1','Ez'); s=sf.load_source(p)
        n=g['comparison_crop_samples']; late=np.empty(501,complex)
        for k in range(0,501,16):
            f=FREQ[k:k+16,None]
            late[k:k+16]=r.dt*(np.exp(-2j*np.pi*f*r.times[n:])@r.samples[n:])
        late/=qs[2][-1].source_spectrum*.025
        additivity=relative(qs[2][-1].response-qs[1][-1].response,late)
        assert additivity<1e-8
        rows.append(dict(id=name,native_sha256=sha(p),old_native_sha256=sha(op),samples=40703,
                         dtype='float64',independent_DFT_relative_L2=errors,
                         receiver_final_5pct_peak_relative_db=tails,
                         appended_tail_DTFT_additivity_relative_L2=additivity,
                         full_vs_crop_complex_spectrum_change=relative(qs[2][-1].response,qs[1][-1].response)))
    a.out.mkdir(parents=True); a.numerical.parent.mkdir(parents=True,exist_ok=True)
    profiles={}; metrics={}
    labels=['旧泥岩1200ns','论文泥岩前1200ns','论文泥岩完整2400ns']
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ)
        h.attrs['configurations']=json.dumps(labels,ensure_ascii=False)
        for w in ['hann','blackman']:
            pp=[]
            for j in range(3):
                z=np.column_stack([q.response for q in qs[j]])
                z=np.column_stack([z,z[:,1]-z[:,0]])
                profile=sf.reconstruct_time_response(replace(qs[j][0],response=z),window=w,zero_pad_factor=8,time_shift=0)
                k=np.arange(7,len(profile.time),101)
                independent=np.exp(2j*np.pi*profile.time[k,None]*FREQ)@(z*profile.weights[:,None])/profile.weights.sum()
                assert relative(profile.complex_bandpass[k],independent)<1e-9
                h.create_dataset(f'config{j}/{w}/complex_bandpass',data=profile.complex_bandpass)
                h.create_dataset(f'config{j}/{w}/real_bandpass',data=profile.real_bandpass)
                if w=='hann': h.create_dataset(f'config{j}/response',data=z)
                pp.append(profile)
            profiles[w]=pp; t=pp[0].time*1e9
            if w=='hann':h.create_dataset('time_s',data=pp[0].time)
            gates={}
            for gate,bounds in c['study_manifest']['predeclared_analysis_windows_ns'].items():
                mask=(t>=bounds[0])&(t<=bounds[1]); details={}
                for j,name in enumerate(['H0','H1','H1_minus_H0']):
                    z=[p.complex_bandpass[mask,j] for p in pp]
                    details[name]=dict(old_peak_ns=float(t[mask][np.argmax(abs(z[0]))]),
                        new_crop_peak_ns=float(t[mask][np.argmax(abs(z[1]))]),
                        new_full_peak_ns=float(t[mask][np.argmax(abs(z[2]))]),
                        material_change_at_equal1200ns=relative(z[1],z[0]),
                        full_vs_crop_complex_change=relative(z[2],z[1]),
                        new_crop_over_old_peak=float(abs(z[1]).max()/abs(z[0]).max()),
                        new_full_over_crop_peak=float(abs(z[2]).max()/abs(z[1]).max()))
                z=pp[2].complex_bandpass[mask]
                details['new_full_H0_over_delta_L2']=float(np.linalg.norm(z[:,0])/np.linalg.norm(z[:,2]))
                details['new_full_delta_over_total_L2']=float(np.linalg.norm(z[:,2])/np.linalg.norm(z[:,1]))
                gates[gate]=dict(bounds_ns=bounds,configurations=details)
            metrics[w]=gates
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    for w,pp in profiles.items():
        t=pp[0].time*1e9
        fig,axs=plt.subplots(2,2,figsize=(15,10),layout='constrained')
        for ax,bounds,title,indices in [(axs[0,0],[0,600],'总场：早段与深部',[0,1]),
                (axs[0,1],[250,550],'总场：深部放大窗',[0,1]),
                (axs[1,0],[250,550],'连通底砂配对复差場 H1−H0',[2])]:
            mask=(t>=bounds[0])&(t<=bounds[1])
            data=np.column_stack([p.complex_bandpass[mask][:,indices] for p in pp]).real
            lim=float(abs(data).max())
            im=ax.imshow(data,cmap='gray',aspect='auto',vmin=-lim,vmax=lim,
                         extent=[-.5,data.shape[1]-.5,t[mask][-1],t[mask][0]])
            ticks=[label+' '+role for label in labels for role in (['H0','H1'] if len(indices)==2 else ['差场'])]
            ax.set_xticks(range(len(ticks)),ticks,rotation=18,fontsize=7)
            ax.set(title=title+'；本面板共用绝对灰度',ylabel='SFCW 时间 / ns')
            ax.axhline(c['study_manifest']['vertical_95MHz_group_delay_guide_ns'],color='tab:green',ls=':',lw=1)
            fig.colorbar(im,ax=ax,label='Re(复带通) / (V/m)/(A·m)')
        mask=(t>=250)&(t<=550)
        for j,(label,p) in enumerate(zip(labels,pp)):
            axs[1,1].plot(t[mask],abs(p.complex_bandpass[mask,2]),label=label+' 差场',lw=1)
        axs[1,1].plot(t[mask],abs(pp[2].complex_bandpass[mask,0]),'--',label='新2400ns H0总场',lw=1)
        axs[1,1].set(title='深部差场与总场：共同幅度，无增益或拟合',xlabel='SFCW 时间 / ns',ylabel='包络 / (V/m)/(A·m)')
        axs[1,1].legend(fontsize=8)
        fig.suptitle(f'190m同站位两项对照，非空间测线 / 约8m离地 / 4.0.1原生FP64\n论文泥岩 εr18、σ0.006；覆盖层/砂岩不变；官方impulse/direct，20–170MHz、501点；{w}、补零8；绿虚线仅为垂直传播参考')
        fig.savefig(a.out/f'paper_mudstone_{w}_gray.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(13,8),layout='constrained')
    for ax,(x,old,dt),g in zip(axs,raw,c['groups']):
        t=np.arange(len(x))*dt*1e9
        ax.plot(t,20*np.log10(np.maximum(abs(x)/abs(x).max(),1e-15)),label='新泥岩2400ns原生场（相对本记录峰值）',lw=.5)
        ax.axvline(1200,color='tab:red',ls='--',label='1200ns比较截取位置')
        ax.axhline(-60,color='black',ls=':',label='官方末5% −60dB建议')
        ax.set(title=g['id']+' 原生响应衰减诊断',xlabel='FDTD 时间 / ns',ylabel='相对场幅 / dB',ylim=[-150,5])
        ax.legend(fontsize=8)
    fig.suptitle('只用于记录尾部诊断；未对求解结果加尾窗，未据此认证时窗收敛')
    fig.savefig(a.out/'paper_mudstone_native_tail.png',dpi=140);plt.close(fig)
    result=dict(status='COMPLETED_TWO_CASE_MATERIAL_AND_RECORD_LENGTH_COMPARISON_NOT_FIELD_VALIDATION',
                script_sha256=sha(__file__),official_processing_sha256=sha(sf.__file__),
                contract_sha256=sha(a.new/'execution_contract.json'),verification_sha256=sha(a.new/'completed_verification.json'),
                numerical_sha256=sha(a.numerical),configuration_order=labels,rows=rows,metrics=metrics,
                method='official_direct',tones=501,band_Hz=[20e6,170e6],step_Hz=300000.,
                source_reference='actual saved impulse samples; physical half-time-step retained; /0.025m current moment, not port S21',
                tail_taper=False,time_shift=False,AGC=False,fit=False,
                vertical_delay_guide_ns=c['study_manifest']['vertical_95MHz_group_delay_guide_ns'],
                limits=c['study_manifest']['limits']+' Crop1200 uses same2400 run, not a separate1200 solve. Full2400 compared only to1200; further convergence not certified.')
    save(a.out/'analysis.json',result)
    print(json.dumps(dict(rows=rows,metrics=metrics),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['new','old','out','numerical']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
