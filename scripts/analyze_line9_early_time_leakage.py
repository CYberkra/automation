"""Linear early/late partition of recorded fields; no solver or physical clean target."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from review_line9_result_packages import FREQ,sha,save,rel


def inverse(response,band):
    use=(FREQ>=band[0])&(FREQ<=band[1]);freq=FREQ[use]
    weights=sf.spectral_window('gaussian',len(freq))
    padded=np.zeros((len(freq)*8,response.shape[1]),complex)
    padded[:len(freq)]=response[use]*weights[:,None]
    return np.fft.ifft(padded,axis=0)*(len(freq)*8/len(freq)),np.arange(len(padded))/(len(padded)*300000)


def main(root,review,out):
    if out.exists(): raise ValueError('Fresh diagnostic output required')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    findings=[]
    for p in sorted(root.glob('line9_pkg*')):
        if not p.is_dir(): continue
        audit=json.loads((review/f'{p.name}_audit.json').read_text('utf-8'))
        rows=audit['records'];ids=[r['id'] for r in rows]
        raw=[]
        for r in rows:
            f=p/'cases'/r['id']/'profile.h5'
            if sha(f)!=r['native_sha256']: raise ValueError('Raw audit identity changed')
            with h5py.File(f) as h:
                raw.append(h['rxs/rx1/Ez'][:]);dt=float(h.attrs['dt'])
                source_offset=float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
        raw=np.column_stack(raw);native_t=np.arange(len(raw))*dt
        source=dt*40*np.exp(-2j*np.pi*FREQ*source_offset)
        total=sf.engineering_dft(raw,dt,FREQ)/source[:,None]/.025
        early={}
        for cutoff in [100,120]:
            start=(cutoff-20)*1e-9;stop=cutoff*1e-9
            gate=np.ones(len(raw));gate[native_t>=stop]=0
            m=(native_t>start)&(native_t<stop)
            gate[m]=.5*(1+np.cos(np.pi*(native_t[m]-start)/(stop-start)))
            early[cutoff]=sf.engineering_dft(raw*gate[:,None],dt,FREQ)/source[:,None]/.025
        entries={}
        for band in [(20e6,170e6),(20e6,40e6)]:
            full,time=inverse(total,band);wm=(time>=320e-9)&(time<=497e-9)
            key=f'{int(band[0]/1e6)}_{int(band[1]/1e6)}MHz';entries[key]={}
            for cutoff,spectrum in early.items():
                part,_=inverse(spectrum,band);late,_=inverse(total-spectrum,band)
                a=full[wm].ravel();b=part[wm].ravel()
                entries[key][str(cutoff)]=dict(early_over_total_complex_L2=float(np.linalg.norm(b)/np.linalg.norm(a)),
                    late_over_total_complex_L2=float(np.linalg.norm(late[wm])/np.linalg.norm(a)),
                    complex_correlation_magnitude=float(abs(np.vdot(a,b))/(np.linalg.norm(a)*np.linalg.norm(b))),
                    closure_relative_L2=rel(part+late,full))
                gates={}
                for name in ['cover_base','basal_sand']:
                    centers=np.asarray([r['geometry'][name]['time95_ns'] for r in rows])
                    selected=abs(time[:,None]*1e9-centers[None,:])<=10
                    a=full[selected];b=part[selected]
                    gates[name]=dict(early_over_total_complex_L2=float(np.linalg.norm(b)/np.linalg.norm(a)),
                        late_over_total_complex_L2=float(np.linalg.norm(late[selected])/np.linalg.norm(a)),
                        complex_correlation_magnitude=float(abs(np.vdot(a,b))/(np.linalg.norm(a)*np.linalg.norm(b))))
                entries[key][str(cutoff)]['actual_geometry_pm10ns_gates']=gates
            if band[1]==40e6:
                part,_=inverse(early[120],band);late,_=inverse(total-early[120],band)
                pos=np.array([r['chainage_m'] for r in rows]);order=np.argsort(pos)
                edges=np.r_[pos[order][0]-.1,(pos[order][:-1]+pos[order][1:])/2,pos[order][-1]+.1]
                keep=time<=600e-9;tt=time[keep]*1e9
                te=np.r_[tt[0]-(tt[1]-tt[0])/2,(tt[:-1]+tt[1:])/2,tt[-1]+(tt[-1]-tt[-2])/2]
                reference=float(np.max(abs(full)))
                fig,axes=plt.subplots(3,1,figsize=(12,11),layout='constrained')
                for ax,(v,title) in zip(axes,[(full,'完整记录重建'),(part,'仅早时分量：0–100ns保留，100–120ns平滑降至零'),(late,'剩余晚时分量：与上行复数相加等于完整结果')]):
                    db=20*np.log10(np.maximum(abs(v[keep][:,order])/reference,1e-10))
                    pic=ax.pcolormesh(edges,te,db,cmap='gray_r',vmin=-90,vmax=-30,shading='flat',rasterized=True)
                    basal=[r['geometry']['basal_sand']['time95_ns'] for r in rows]
                    ax.plot(pos[order],np.asarray(basal)[order],ls='--',color='#1ec766',lw=1,label='实际体素底部砂岩顶/95MHz近似')
                    ax.set(ylim=(600,0),ylabel='时间 / ns',title=title);ax.invert_xaxis()
                    ax.legend(fontsize=8,loc='lower right');fig.colorbar(pic,ax=ax,label='幅度dB；三行同一绝对参考')
                axes[-1].set_xlabel('测线里程 / m；只绘已完成区')
                fig.suptitle(f"{p.name.split('_')[1]}：20–39.8MHz/Gaussian；未AGC/叠道/尾渐消\n早时分量不含此模型深层回波；平滑时间分割本身也改变频谱，因此只认证带限泄漏机制")
                fig.savefig(out/f'{p.name}_early_partition.png',dpi=140);plt.close(fig)
        findings.append(dict(package=p.name,scope='early/late exact linear partition; ratios are not independent physical contribution percentages',windows_320_497ns=entries))
    psf={}
    for band in [(20e6,170e6),(20e6,40e6)]:
        use=(FREQ>=band[0])&(FREQ<=band[1]);freq=FREQ[use];weights=sf.spectral_window('gaussian',len(freq))
        t=np.linspace(0,200e-9,2001)
        value=abs(np.exp(2j*np.pi*t[:,None]*freq)@weights)/len(freq)
        crossing=int(np.flatnonzero(value<=.5)[0])
        psf[f'{int(band[0]/1e6)}_{int(band[1]/1e6)}MHz']=dict(realized_frequency_limits_Hz=[float(freq[0]),float(freq[-1])],frequency_count=len(freq),
            ideal_single_echo_full_width_at_half_amplitude_ns=float(2*t[crossing]*1e9))
    save(out/'early_partition_report.json',dict(calls_solver=False,findings=findings,gaussian_single_echo_psf=psf,
        script_sha256=sha(__file__),limits='No matched geological ablation. Early/late are time-windowed components, not pure direct/surface or pure target. Gates100/120ns include temporal partition effects. L2 ratios can exceed1 by coherent cancellation.'))
    print(json.dumps(dict(findings=findings,psf=psf),ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--review',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.root.resolve(),a.review.resolve(),a.out.resolve())
