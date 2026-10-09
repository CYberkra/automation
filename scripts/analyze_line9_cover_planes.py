"""Native invariant, exact-tone SFCW and conditional lossy cover-plane closure."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from scipy.signal.windows import tukey
from hs_capsule_identity import sha256 as sha
from analyze_line9_native_probes import collocate
from analyze_line9_v401_version_controls import FREQ,inverse,response
from line9_cover_angular_transport import transport,cover_permittivity


def read(p):return json.loads(Path(p).read_text('utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def metrics(v):
    target,pred,pu,pd,au,ad,raw=[v[:,:,i] for i in range(7)];norm=np.linalg.norm
    return {'total_relative_error':float(norm(pred-target)/norm(target)),
            'up_relative_error':float(norm(pu-au)/norm(au)),
            'down_relative_error':float(norm(pd-ad)/norm(ad)),
            'up_complex_cosine':float(abs(np.vdot(pu,au))/(norm(pu)*norm(au))),
            'down_complex_cosine':float(abs(np.vdot(pd,ad))/(norm(pd)*norm(ad))),
            'destination_down_up_norm_ratio':float(norm(ad)/norm(au)),
            'sector_raw_relative_change':float(norm(target-raw)/norm(raw))}


def project_pair(spectra,pair,cfg):
    n=cfg['padding'];loc=slice((n-41)//2,(n-41)//2+41)
    q=np.zeros((501,6,n),complex);q[:,:,loc]=spectra*tukey(41,cfg['alpha'])
    i,j=pair;dy=(j-i)*1.5
    lower=transport(q[:,2*i],q[:,2*i+1],FREQ,.5,dy,cfg['sector'],cfg['rolloff'])
    upper=transport(q[:,2*j],q[:,2*j+1],FREQ,.5,0,cfg['sector'],cfg['rolloff'])
    return np.stack([upper['total'][:,loc],(lower['predicted_up']+lower['predicted_down'])[:,loc],
                     lower['predicted_up'][:,loc],lower['predicted_down'][:,loc],upper['up'][:,loc],upper['down'][:,loc],q[:,2*j,loc]],axis=-1)


def main(a):
    from audit_line9_cover_plane_inputs import check
    assert not a.out.exists() and not a.arrays.exists();check(a.package)
    m=read(a.package/'manifest.json');cert=read(a.source/'completed_verification.json');c=read(a.source/'execution_contract.json')
    assert cert['completed'] and cert['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.source/'profile.h5')==cert['groups'][0]['native_sha256']
    primary={'padding':512,'alpha':.25,'sector':.9,'rolloff':.25}
    configs=[dict(padding=n,alpha=al,sector=s,rolloff=.25) for n in [256,512,1024] for al in [0,.25,.5] for s in [.8,.9,.95]]
    variants=[('full',None),('native_40_230',[40,230]),('native_180_400',[180,400]),('native_220_380',[220,380])]
    contract={'status':'POST_OBSERVER_CONDITIONAL_COVER_ANALYSIS_FROZEN','manifest_sha256':sha(a.package/'manifest.json'),
              'raw_sha256':sha(a.source/'profile.h5'),'reference_sha256':sha(a.reference),'execution_contract_sha256':cert['contract_sha256'],
              'source_files_sha256':{p:sha(Path(__file__).with_name(p)) for p in ['analyze_line9_cover_planes.py','line9_cover_angular_transport.py','line9_air_smooth_transport.py','analyze_line9_native_probes.py','analyze_line9_v401_version_controls.py']},
              'primary':primary,'sensitivity':configs,'plane_y_m':m['cover_plane_y_m'],'x_m':list(np.arange(154,174.01,.5)),
              'pairs':[[0,1],[1,2],[0,2]],'core_x_m':[160,168],'source_reference_m':.025,'native_variants':variants,
              'native_window_alpha':.25,'gates_sfcw_ns':{'early':[50,220],'late':[220,330],'wide':[50,400]},
              'epsilon_formula':'11+0.5/(1+i*2*pi*f*6.4567e-9)-i*0.001/(2*pi*f*epsilon0)',
              'field_projection':'TM Ez/Hx, exp(+i omega t), outgoing complex ky; smooth air-radiating q sector',
              'limits':'Finite20m aperture, .5m sampling aliases non-air-radiating cover modes; exterior is nonflat. No unique ray/bounce count or energy fraction; gated native sensitivity is not formal SFCW; no fitted delay/amplitude/phase, no site or finite3D certification.'}
    a.out.mkdir(parents=True);save(a.out/'contract.json',contract)
    channels=0;ee=[];hh=[];yy=[];probes=m['probes'][287:]
    with h5py.File(a.source/'profile.h5') as h,h5py.File(a.reference) as old:
        assert h.attrs['nrx']==len(h['rxs'])==1231;dt=float(h.attrs['dt']);assert dt==m['dt_s']
        source=h['srcs/src1/excitation/samples'][:];assert source.tobytes()==old['srcs/src1/excitation/samples'][:].tobytes()
        for name,r in old['rxs'].items():
            for attr in ['Name','Position','GridPosition']:np.testing.assert_array_equal(h['rxs/'+name].attrs[attr],r.attrs[attr])
            for field in r:assert h['rxs/'+name+'/'+field][:].tobytes()==r[field][:].tobytes();channels+=1
        for p in probes:
            rr=[]
            for q in p['anchors']:
                r=h[f"rxs/rx{q['receiver_index']}"];assert r.attrs['Name']==q['name'] and set(r)==set(q['outputs'])
                np.testing.assert_array_equal(r.attrs['GridPosition'],q['coord']);rr.append({key:r[key][:] for key in q['outputs']})
                for key in q['outputs']:
                    v=rr[-1][key];assert v.dtype==np.float64 and v.shape==(20352,) and np.isfinite(v).all()
                    assert r[key].attrs['TimeSampleOffset']==(0 if key=='Ez' else -.5*dt)
            e,hx,hy=collocate(rr[0]['Ez'],rr[0]['Hx'],rr[1]['Hx'],rr[0]['Hy'],rr[2]['Hy']);ee.append(e);hh.append(hx);yy.append(hy)
        rx=h['rxs/rx1/Ez'][:]
    assert channels==1436
    e=np.array(ee);hx=np.array(hh);hy=np.array(yy);t=np.arange(e.shape[1])*dt;st_src=np.arange(len(source))*dt
    fields=np.stack([v[i*41:(i+1)*41] for i in range(3) for v in [e,hx]])
    full={};x=np.arange(154,174.01,.5);core=(x>=160)&(x<=168)
    for name,bounds in variants:
        win=np.ones(len(t))
        if bounds:
            win[:]=0;ix=np.flatnonzero((t*1e9>=bounds[0])&(t*1e9<=bounds[1]));win[ix]=tukey(len(ix),.25)
        z=np.empty((501,6,41),complex)
        for start in range(0,501,16):
            f=FREQ[start:start+16,None];ss=dt*(np.exp(-2j*np.pi*f*st_src)@source)
            z[start:start+16]=(np.exp(-2j*np.pi*f*t)@(fields*win).reshape(246,-1).T).reshape(len(f),6,41)*dt/ss[:,None,None]/.025
        full[name]=z
    zrx,err=response(a.source/'profile.h5',.025);zold,old_err=response(a.reference,.025);assert np.array_equal(zrx,zold)
    rows=[];profiles={};inverse_errors={}
    for name,spectra in full.items():
        for pair in contract['pairs']:
            for cfg in (configs if name=='full' else [primary]):
                z=project_pair(spectra,pair,cfg)
                for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
                    w=w/w.mean();v,st=inverse(z.reshape(501,-1),w);v=v.reshape(4008,41,7)
                    row={'variant':name,'pair':pair,**cfg,'window':window,'gates':{}}
                    for gate,bounds in contract['gates_sfcw_ns'].items():
                        mask=(st*1e9>=bounds[0])&(st*1e9<=bounds[1]);row['gates'][gate]=metrics(v[mask][:,core])
                    rows.append(row)
                    if cfg==primary:profiles[name+'_'+str(pair[0])+str(pair[1])+'_'+window]=v
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();v,st=inverse(zrx[:,None],w);pick=np.arange(7,4008,101)
        direct=np.exp(2j*np.pi*st[pick,None]*FREQ)@(zrx*w)/501
        inverse_errors[window]=float(np.linalg.norm(v[pick,0]-direct)/np.linalg.norm(direct));assert inverse_errors[window]<1e-9
        profiles['main_'+window]=v[:,0]
    a.arrays.parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.arrays,time_s=t,Ez=e,Hx=hx,Hy=hy,source=source,rx=rx,sfcw_time_s=st,frequency_Hz=FREQ,x_m=x,
             **profiles,**{'spectra_'+name:z for name,z in full.items()})
    result={'status':'IMMUTABLE_OBSERVER_EXACT_SFCW_CONDITIONAL_COVER_ANALYSIS_COMPLETE','contract_sha256':sha(a.out/'contract.json'),
            'arrays_sha256':sha(a.arrays),'previous_receiver_channels_bitwise_equal':1436,'source_bitwise_equal':True,
            'main_SFCW_bitwise_equal':True,'main_DFT_errors':[err,old_err],'inverse_errors':inverse_errors,
            'metrics':rows,'cover_epsilon_selected':[[float(f),float(z.real),float(z.imag)] for f,z in zip(FREQ[[0,250,500]],cover_permittivity(FREQ[[0,250,500]]))],
            'new_solver_runs':1,'limits':contract['limits']}
    save(a.out/'analysis.json',result);plot(a.out,t,x,e,hx,st,profiles,rows)
    print(json.dumps({'status':result['status'],'primary_full_late':[r for r in rows if r['variant']=='full' and r['window']=='hann' and all(r[k]==v for k,v in primary.items())]},ensure_ascii=False))


def plot(out,t,x,e,hx,st,profiles,rows):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    show=(t*1e9>=40)&(t*1e9<=400);fig,axs=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    lim=float(abs(e[:,show]).max());fluxlim=float(abs((e*hx)[:,show]).max())
    for i,y in enumerate([27.025,28.525,30.025]):
        for j,data in enumerate([e,e*hx]):
            im=axs[i,j].imshow(data[i*41:(i+1)*41,show].T,cmap='gray' if j==0 else 'RdBu_r',vmin=-[lim,fluxlim][j],vmax=[lim,fluxlim][j],aspect='auto',extent=[153.75,174.25,400,40])
            axs[i,j].set(title=f'覆盖层 y={y}m：'+['原始总场 Ez','总场 Sy（红上/蓝下）'][j],xlabel='模型 x / m',ylabel='原生时间 / ns');fig.colorbar(im,ax=axs[i,j],label=['Ez / (V/m)','Sy / (W/m²)'][j])
    fig.suptitle('190m高损耗H0：三条覆盖层水平观测线，各列统一物理色标\n固定一次激发的41列时空图；不是移动天线B-scan，无增益、去背景或拟合')
    fig.savefig(out/'native_three_cover_planes.png',dpi=140);plt.close(fig)
    show=(st*1e9>=220)&(st*1e9<=330);z=profiles['full_02_hann'];lim=float(abs(z[show,:,:6].real).max());fig,axs=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    titles=['上平面同扇区实际总场','下平面上+下行传播预测','下平面上行传播预测','下平面下行传播预测','上平面上行投影','上平面下行投影']
    for j,ax in enumerate(axs.flat):
        im=ax.imshow(z[show,:,j].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[153.75,174.25,330,220]);ax.set(title=titles[j],xlabel='模型 x / m',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='(V/m)/(A·m)；共同尺度')
    fig.suptitle('覆盖层27.025→30.025m：复杂介电谱分波/传播，晚窗精确501点 Hann\n20m有限孔径/.5m空间采样；仅平滑空气可辐射扇区，其他覆盖层模式可能混叠；无幅相/时延拟合')
    fig.savefig(out/'conditional_cover_split_gray.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained');mask=(st*1e9>=50)&(st*1e9<=400);ix=int(np.flatnonzero(x==164)[0])
    for pair,ax in zip(['01','12','02'],axs.flat):
        v=profiles['full_'+pair+'_hann']
        for j,label in [(4,'目标平面上行投影'),(2,'来源上行传播预测'),(5,'目标平面下行投影'),(3,'来源下行传播预测')]:ax.plot(st[mask]*1e9,v[mask,ix,j].real,lw=.9,label=label)
        ax.set(title='平面 '+pair+'；模型x164m，原幅值',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');ax.legend(fontsize=8)
    p=profiles['main_hann'];data=np.column_stack([p[mask].real,p[mask].real,np.zeros(mask.sum())]);lm=float(abs(data).max())
    im=axs[1,1].imshow(data,cmap='gray',vmin=-lm,vmax=lm,aspect='auto',extent=[-.5,2.5,400,50]);axs[1,1].set_xticks(range(3),['原H0','新增探针H0','新−原']);axs[1,1].set(title='单站配置列：主Rx精确SFCW不变',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axs[1,1],label='(V/m)/(A·m)')
    fig.suptitle('覆盖层三对平面的同位分波比较；单站配置图不是空间测线B-scan\n原高损耗粉质粘土假设配方 / 非平地H0，无底部连通砂岩；Hann与Blackman数字全部另存')
    fig.savefig(out/'cover_pair_waveforms_and_invariant.png',dpi=140);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['package','source','reference','arrays','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
