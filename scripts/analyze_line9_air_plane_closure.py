"""Reuse native air planes to test finite-aperture upward/downward transport."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from scipy.signal.windows import tukey
from hs_capsule_identity import sha256 as sha
from line9_air_angular_transport import project,sector_field
from analyze_line9_v401_version_controls import FREQ,inverse,response
from gprMax.toolboxes.SFCW import processing as sf


def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(a):
    assert not a.out.exists() and not a.spectra.exists() and not a.profiles.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(a.package/'manifest.json');parent=read(a.parent/'analysis.json');audit=read(a.parent/'independent_audit.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256']==sha(a.parent/'analysis.json')
    assert sha(a.arrays)==parent['arrays_sha256'] and sha(a.raw)==parent['native_sha256']
    probes=m['probes'];ids={p['id']:i for i,p in enumerate(probes)}
    bands={b:sorted([p for p in probes if p['band']==b],key=lambda p:p['x_m']) for b in ['air_34','air_38']}
    x=np.array([p['x_m'] for p in bands['air_34']]);assert np.array_equal(x,[p['x_m'] for p in bands['air_38']])
    assert np.array_equal(np.diff(x),np.full(40,.5))
    heights=[bands[b][0]['y_m'] for b in ['air_34','air_38']]
    assert all(p['y_m']==height for b,height in zip(['air_34','air_38'],heights) for p in bands[b])
    height=heights[1]-heights[0]; assert abs(height-4)<1e-12
    with np.load(a.arrays) as h:
        t=h['time_s']; e=h['Ez']; hx=h['Hx']; sy=h['Sy'];sx=h['Sx'];rx=h['rx']
    source=sf.load_source(a.raw)
    fields=np.array([v[[ids[p['id']] for p in bands[b]]] for b in ['air_34','air_38'] for v in [e,hx]])
    contract={'status':'FROZEN_CPU_FINITE_APERTURE_AIR_PLANE_CLOSURE','script_sha256':sha(__file__),
              'helper_sha256':{n:sha(Path(__file__).with_name(n)) for n in ['line9_air_angular_transport.py','analyze_line9_v401_version_controls.py']},
              'parent_analysis_sha256':sha(a.parent/'analysis.json'),'native_sha256':sha(a.raw),'arrays_sha256':sha(a.arrays),
              'x_m':x.tolist(),'plane_y_m':heights,'spacing_m':.5,'height_difference_m':height,
              'frequency_Hz':FREQ.tolist(),'sector_fractions':[.8,.9,.95],'aperture_tukey_alpha':[0.,.25,.5],
              'spatial_padding_samples':[256,512,1024],'primary':{'sector':.9,'tukey_alpha':.25,'padding':512},
              'windows':['hann','blackman'],'comparison_gate_sfcw_ns':[250,400],'interior_x_m':[160,168],
              'native_exploratory_up_windows_ns':{'cover_lower':[220,290],'cover_middle':[250,315],'cover_upper':[275,340],'air_surface':[280,350]},
              'native_flux_integral_relative_ns':[-4,4],'source_reference_m':.025,
              'limits':'20m sampled aperture; |qx|<sector*k0 only, grazing/evanescent/exterior contributions excluded; not full-space closure, unique exit/bounce or conserved-energy certificate; native local extrema post-observation exploratory'}
    a.out.mkdir(parents=True);save(a.out/'contract.json',contract)
    spectra=np.empty((501,4,41),complex)
    for start in range(0,501,16):
        f=FREQ[start:start+16,None]
        kernel=np.exp(-2j*np.pi*f*t)
        ss=source.dt*(np.exp(-2j*np.pi*f*source.times)@source.samples)
        spectra[start:start+16]=(kernel@fields.reshape(164,-1).T).reshape(len(f),4,41)*parent['dt_s']/ss[:,None,None]/.025
    primary_native_dft,err=response(a.raw,.025)
    dft_rx=np.array([np.exp(-2j*np.pi*freq*t)@rx for freq in FREQ])*parent['dt_s']
    ss=np.array([np.exp(-2j*np.pi*freq*source.times)@source.samples for freq in FREQ])*source.dt
    trim_error=float(np.linalg.norm(dft_rx/ss/.025-primary_native_dft)/np.linalg.norm(primary_native_dft))
    a.spectra.parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.spectra,frequency_Hz=FREQ,fields=spectra,x_m=x,plane_y_m=heights)
    metrics=[];primary={};core=(x>=160)&(x<=168)
    for padding in contract['spatial_padding_samples']:
        left=(padding-41)//2;loc=slice(left,left+41)
        for alpha in contract['aperture_tukey_alpha']:
            wspace=tukey(41,alpha);q=np.zeros((501,4,padding),complex);q[:,:,loc]=spectra*wspace
            for sector in contract['sector_fractions']:
                result=project(q[:,0],q[:,1],FREQ,.5,height,sector)
                actual=sector_field(q[:,2],FREQ,.5,sector)
                rawtarget=q[:,2,loc]
                values=np.stack([actual[:,loc],(result['predicted_up']+result['predicted_down'])[:,loc],
                    result['predicted_up'][:,loc],result['predicted_down'][:,loc],result['up'][:,loc],result['down'][:,loc],rawtarget],axis=-1)
                for window,wtime in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
                    wtime/=wtime.mean();z,st=inverse(values.reshape(501,-1),wtime);z=z.reshape(len(st),41,7)
                    gate=(st*1e9>=250)&(st*1e9<=400);v=z[gate][:,core];norm=np.linalg.norm
                    target=v[:,:,0];pred=v[:,:,1];up=v[:,:,2];down=v[:,:,3];raw=v[:,:,6]
                    metric={'padding':padding,'tukey_alpha':alpha,'sector':sector,'window':window,
                            'predicted_total_relative_error':float(norm(pred-target)/norm(target)),
                            'up_only_relative_error':float(norm(up-target)/norm(target)),
                            'prediction_target_complex_cosine':float(abs(np.vdot(target,pred))/(norm(target)*norm(pred))),
                            'up_target_complex_cosine':float(abs(np.vdot(target,up))/(norm(target)*norm(up))),
                            'predicted_down_over_up_norm':float(norm(down)/norm(up)),
                            'sector_target_over_raw_relative_change':float(norm(target-raw)/norm(raw))}
                    metrics.append(metric)
                    if padding==512 and alpha==.25 and sector==.9:
                        primary[window]=z
    np.savez(a.profiles,time_s=st,x_m=x,**primary)
    # Preserve all 41 native late-event rows including negative integrated Sy.
    lookup={(p['band'],p['x_m']):i for i,p in enumerate(probes)};events=[];ns=t*1e9
    for xx in x:
        row={'x_m':float(xx),'bands':{}}
        for band,bounds in contract['native_exploratory_up_windows_ns'].items():
            k=lookup[band,xx];ix=np.flatnonzero((ns>=bounds[0])&(ns<=bounds[1]));j=ix[np.argmax(sy[k,ix])]
            near=(ns>=ns[j]-4)&(ns<=ns[j]+4);sxsum=float(np.sum(sx[k,near]));sysum=float(np.sum(sy[k,near]))
            row['bands'][band]={'index':int(j),'time_ns':float(ns[j]),'Sy_peak_W_m2':float(sy[k,j]),
                               'signed8ns_Sy_J_m2':sysum*parent['dt_s'],'signed8ns_Sx_J_m2':sxsum*parent['dt_s'],
                               'direction_deg':float(np.degrees(np.arctan2(sysum,sxsum))),'peak_at_window_edge':bool(j in [ix[0],ix[-1]])}
        events.append(row)
    summary={'status':'FINITE_APERTURE_AIR_PLANE_DIAGNOSTIC_COMPLETE_NOT_UNIQUE_PATH','script_sha256':sha(__file__),
             'contract_sha256':sha(a.out/'contract.json'),'spectra_sha256':sha(a.spectra),'profiles_sha256':sha(a.profiles),
             'last_receiver_sample_trim_DFT_relative_L2':trim_error,'reference_original_DFT_relative_L2':err,
             'sensitivity_metrics':metrics,'native_late_events':events,'new_solver_runs':0,'limits':contract['limits']}
    save(a.out/'analysis.json',summary)
    plot(a.out,st,x,primary,metrics,events)
    print(json.dumps({'status':summary['status'],'primary_metrics':[v for v in metrics if v['padding']==512 and v['tukey_alpha']==.25 and v['sector']==.9],
                      'trim_error':trim_error},ensure_ascii=False))


def plot(out,t,x,profiles,metrics,events):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    ns=t*1e9;show=(ns>=250)&(ns<=400);z=profiles['hann'];lim=float(abs(z[show,:,:6].real).max())
    titles=['上平面实际场：同一角度扇区','下平面传播4m：上+下行预测','仅上行分量传播预测','仅下行分量传播预测','下平面上行投影','下平面下行投影']
    fig,axs=plt.subplots(3,2,figsize=(14,13),layout='constrained')
    for j,ax in enumerate(axs.flat):
        im=ax.imshow(z[show,:,j].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[153.75,174.25,400,250])
        ax.set(title=titles[j],xlabel='模型 x / m',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='(V/m)/(A·m)；六面板共用尺度')
    fig.suptitle('两条空气平面的有限孔径传播核查：固定一次激发，非移动天线B-scan\n20m孔径/Tukey0.25，保留|qx|≤0.9k0传播扇区，Hann精确501点；无AGC/幅相或时延拟合')
    fig.savefig(out/'air_plane_closure_gray.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    ix=int(np.flatnonzero(x==170)[0])
    for j,label in [(6,'上平面孔径内实际原场'),(0,'上平面同扇区场'),(1,'上下行传播预测'),(2,'上行传播预测')]:axs[0,0].plot(ns[show],z[show,ix,j].real,label=label,lw=1)
    axs[0,0].set(xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)',title='模型x170m原幅值；四曲线无拟合');axs[0,0].legend(fontsize=8)
    for window,style in [('hann','-'),('blackman','--')]:
        for alpha,marker in [(0,'o'),(.25,'s'),(.5,'^')]:
            rows=[r for r in metrics if r['padding']==512 and r['window']==window and r['tukey_alpha']==alpha]
            axs[0,1].plot([r['sector'] for r in rows],[r['up_only_relative_error'] for r in rows],style+marker,label=window+f' / Tukey{alpha}')
    axs[0,1].set(xlabel='保留的 |qx|/k0 最大值',ylabel='上行预测与同扇区实场相对L2',title='核心x160–168m；角度/孔径窗敏感性');axs[0,1].legend(fontsize=8)
    for band,label in [('cover_lower','覆盖层下部'),('cover_middle','覆盖层中部'),('cover_upper','覆盖层上部'),('air_surface','地表上方空气')]:
        axs[1,0].plot(x,[r['bands'][band]['time_ns'] for r in events],'.-',label=label)
        axs[1,1].plot(x,[r['bands'][band]['direction_deg'] for r in events],'.-',label=label)
    axs[1,0].set(xlabel='模型 x / m',ylabel='原生探索性Sy极值时间 / ns',title='时间先后不能独立认证传播链');axs[1,0].legend(fontsize=8)
    axs[1,1].set(xlabel='模型 x / m',ylabel='±4ns总场积分方向 / 度',title='反向/叠加结果保留；不强行串射线');axs[1,1].legend(fontsize=8)
    fig.suptitle('孔径/角度限制及层内方向歧义：全部配置另存，未挑最优替代主配置')
    fig.savefig(out/'air_plane_closure_sensitivity.png',dpi=140);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['package','parent','raw','arrays','out','spectra','profiles']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
