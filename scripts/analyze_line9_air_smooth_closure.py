"""Post-failure smooth angular transport and native late-window sensitivity."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.signal.windows import tukey
from line9_air_smooth_transport import smooth_transport,smooth_field
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import FREQ,inverse
from gprMax.toolboxes.SFCW import processing as sf


def read(p):return json.loads(p.read_text('utf-8'))
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(a):
    assert not a.out.exists() and not a.arrays.exists()
    old=read(a.hard/'analysis.json');c1=read(a.hard/'contract.json');parent=read(a.parent/'analysis.json');m=read(a.package/'manifest.json')
    assert sha(a.spectra)==old['spectra_sha256'] and sha(a.native_arrays)==parent['arrays_sha256'] and sha(a.raw)==parent['native_sha256']
    assert sha(a.package/'manifest.json')==parent['manifest_sha256']
    with np.load(a.spectra) as h:full=h['fields'];x=h['x_m'];f=h['frequency_Hz']
    assert np.array_equal(f,FREQ)
    variants=[('full',None)]+[(f'native_{lo}_{hi}',[lo,hi]) for lo,hi in [(190,450),(220,420),(240,410)]]
    configurations=[{'padding':n,'alpha':alpha,'sector':sector,'rolloff':.25} for n in [256,512,1024] for alpha in [0.,.25,.5] for sector in [.8,.9,.95]]
    primary={'padding':512,'alpha':.25,'sector':.9,'rolloff':.25}
    configurations += [dict(primary,rolloff=r) for r in [.1,.4]]
    contract={'status':'FROZEN_POST_FAILURE_SMOOTH_APERTURE_DIAGNOSTIC','script_sha256':sha(__file__),
              'helper_sha256':{n:sha(Path(__file__).with_name(n)) for n in ['line9_air_smooth_transport.py','line9_air_angular_transport.py','analyze_line9_v401_version_controls.py']},
              'hard_analysis_sha256':sha(a.hard/'analysis.json'),'native_parent_analysis_sha256':sha(a.parent/'analysis.json'),
              'manifest_sha256':sha(a.package/'manifest.json'),'raw_sha256':sha(a.raw),'input_spectra_sha256':sha(a.spectra),
              'native_arrays_sha256':sha(a.native_arrays),'configurations':configurations,'primary':primary,
              'variant_native_gates':variants,'native_gate_tukey_alpha':.25,'frequency_Hz':f.tolist(),'x_m':x.tolist(),
              'plane_y_m':c1['plane_y_m'],'spacing_m':.5,'height_difference_m':4.,'source_reference_m':.025,
              'gates_sfcw_ns':{'wide':[250,400],'late':[290,355]},'core_x_m':[160,168],
              'method_change_reason':'Hard angular cutoff creates excessive late tails; smooth cosine transition and explicitly windowed native sensitivity are post-observation followups, not preregistered physics gates',
              'limits':'20m aperture; smooth propagating sector, finite aperture/grazing/evanescent/exterior excluded. Native gate is diagnostic only, not changed formal SFCW or a causal event isolator; no unique cover ray/bounce or energy fraction.'}
    a.out.mkdir(parents=True);save(a.out/'contract.json',contract)
    with np.load(a.native_arrays) as h:
        t=h['time_s'];e=h['Ez'];hx=h['Hx']
    lookup={p['id']:i for i,p in enumerate(m['probes'])}
    band_indices={b:[lookup[p['id']] for p in sorted([q for q in m['probes'] if q['band']==b],key=lambda p:p['x_m'])] for b in ['air_34','air_38']}
    fields=np.array([v[band_indices[b]] for b in ['air_34','air_38'] for v in [e,hx]])
    source=sf.load_source(a.raw);all_spectra={'full':full}
    for name,bounds in variants[1:]:
        ix=np.flatnonzero((t*1e9>=bounds[0])&(t*1e9<=bounds[1]));window=np.zeros(len(t));window[ix]=tukey(len(ix),.25)
        native=fields*window;z=np.empty_like(full)
        for start in range(0,501,16):
            ff=f[start:start+16,None];ss=source.dt*(np.exp(-2j*np.pi*ff*source.times)@source.samples)
            z[start:start+16]=(np.exp(-2j*np.pi*ff*t)@native.reshape(164,-1).T).reshape(len(ff),4,41)*parent['dt_s']/ss[:,None,None]/.025
        all_spectra[name]=z
    metrics=[];profiles={};core=(x>=160)&(x<=168)
    for variant,spectra in all_spectra.items():
        configs=configurations if variant=='full' else [primary]
        for cfg in configs:
            n=cfg['padding'];left=(n-41)//2;loc=slice(left,left+41);q=np.zeros((501,4,n),complex);q[:,:,loc]=spectra*tukey(41,cfg['alpha'])
            lower=smooth_transport(q[:,0],q[:,1],f,.5,4,cfg['sector'],cfg['rolloff'])
            upper=smooth_transport(q[:,2],q[:,3],f,.5,0,cfg['sector'],cfg['rolloff'])
            actual=smooth_field(q[:,2],f,.5,cfg['sector'],cfg['rolloff'])
            values=np.stack([actual[:,loc],(lower['predicted_up']+lower['predicted_down'])[:,loc],
                lower['predicted_up'][:,loc],lower['predicted_down'][:,loc],upper['up'][:,loc],upper['down'][:,loc],q[:,2,loc]],axis=-1)
            for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
                w/=w.mean();z,st=inverse(values.reshape(501,-1),w);z=z.reshape(4008,41,7);row={'variant':variant,**cfg,'window':window,'gates':{}}
                for gate,bounds in contract['gates_sfcw_ns'].items():
                    mask=(st*1e9>=bounds[0])&(st*1e9<=bounds[1]);v=z[mask][:,core];norm=np.linalg.norm
                    target,pred,pu,pd,au,ad,raw=[v[:,:,j] for j in range(7)]
                    cosine=lambda a,b:float(abs(np.vdot(a,b))/(norm(a)*norm(b)))
                    row['gates'][gate]={'total_transport_relative_error':float(norm(pred-target)/norm(target)),
                        'up_only_vs_total_relative_error':float(norm(pu-target)/norm(target)),
                        'up_to_up_transport_relative_error':float(norm(pu-au)/norm(au)),
                        'up_to_up_complex_cosine':cosine(pu,au),'up_to_total_complex_cosine':cosine(pu,target),
                        'upper_down_over_up_norm':float(norm(ad)/norm(au)),
                        'predicted_down_over_up_norm':float(norm(pd)/norm(pu)),
                        'sector_over_raw_relative_change':float(norm(target-raw)/norm(raw))}
                metrics.append(row)
                if cfg==primary:profiles[variant+'_'+window]=z
    a.arrays.parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.arrays,time_s=st,x_m=x,**profiles,**{'spectra_'+name:v for name,v in all_spectra.items()})
    result={'status':'SMOOTH_FINITE_APERTURE_AIR_TRANSPORT_DIAGNOSTIC_COMPLETE','script_sha256':sha(__file__),
            'contract_sha256':sha(a.out/'contract.json'),'arrays_sha256':sha(a.arrays),'sensitivity_metrics':metrics,
            'new_solver_runs':0,'formal_SFCW_changed':False,'limitations':contract['limits']}
    save(a.out/'analysis.json',result);plot(a.out,st,x,profiles,metrics)
    print(json.dumps({'status':result['status'],'primary_metrics':[r for r in metrics if all(r[k]==v for k,v in primary.items())]},ensure_ascii=False))


def plot(out,t,x,p,m):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    ns=t*1e9;show=(ns>=250)&(ns<=400);z=p['full_hann'];lim=float(abs(z[show,:,:6].real).max())
    titles=['上平面实际同扇区总场','下平面上+下行传播预测','下平面上行传播预测','下平面下行传播预测','上平面上行投影','上平面下行投影']
    fig,axs=plt.subplots(3,2,figsize=(14,13),layout='constrained')
    for j,ax in enumerate(axs.flat):
        im=ax.imshow(z[show,:,j].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[153.75,174.25,400,250]);ax.set(title=titles[j],xlabel='模型 x / m',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='(V/m)/(A·m)；共同尺度')
    fig.suptitle('空气平面传播：原生全时段输入、平滑角谱投影；固定一次激发，非移动天线B-scan\n20m孔径/Tukey0.25，角度上限0.9、余弦过渡宽0.25，Hann精确501点；无幅相/时延拟合')
    fig.savefig(out/'smooth_air_plane_gray.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained');ix=int(np.flatnonzero(x==170)[0])
    for j,label in [(6,'上平面实际原场'),(0,'上平面平滑扇区总场'),(2,'下平面上行传播预测'),(4,'上平面上行投影')]:axs[0,0].plot(ns[show],z[show,ix,j].real,label=label,lw=1)
    axs[0,0].set(xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)',title='模型x170m；原幅值，无拟合');axs[0,0].legend(fontsize=8)
    for alpha,marker in [(0,'o'),(.25,'s'),(.5,'^')]:
        for window,style in [('hann','-'),('blackman','--')]:
            rows=[r for r in m if r['variant']=='full' and r['padding']==512 and r['alpha']==alpha and r['rolloff']==.25 and r['window']==window]
            axs[0,1].plot([r['sector'] for r in rows],[r['gates']['wide']['up_to_up_transport_relative_error'] for r in rows],style+marker,label=window+f'/Tukey{alpha}')
    axs[0,1].set(xlabel='角度扇区上限 |qx|/k0',ylabel='上行→上行相对L2误差',title='有限孔径/角度敏感性；核心x160–168m');axs[0,1].legend(fontsize=8)
    for variant,label in [('full','原生全时段'),('native_190_450','原生190–450ns窗'),('native_220_420','原生220–420ns窗'),('native_240_410','原生240–410ns窗')]:
        v=p[variant+'_hann'];axs[1,0].plot(ns[show],v[show,ix,2].real,label=label,lw=1)
    axs[1,0].set(xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)',title='下平面上行预测；时间窗仅作诊断');axs[1,0].legend(fontsize=8)
    rows=[r for r in m if r['padding']==512 and r['alpha']==.25 and r['sector']==.9 and r['rolloff']==.25 and r['window']=='hann']
    xx=np.arange(len(rows));axs[1,1].bar(xx-.15,[r['gates']['late']['up_to_up_transport_relative_error'] for r in rows],.3,label='上行传播误差')
    axs[1,1].bar(xx+.15,[r['gates']['late']['upper_down_over_up_norm'] for r in rows],.3,label='上平面下/上行范数比')
    axs[1,1].set_xticks(xx,['全时段','190–450','220–420','240–410'],rotation=15);axs[1,1].set(ylabel='无量纲比值',title='固定290–355ns晚窗；不是能量占比');axs[1,1].legend(fontsize=8)
    fig.suptitle('平滑方法为硬截断失败后的诊断修订；全配置保留，不认证唯一射线/反射阶次')
    fig.savefig(out/'smooth_air_plane_sensitivity.png',dpi=140);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['hard','spectra','native-arrays','raw','parent','package','out','arrays']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
