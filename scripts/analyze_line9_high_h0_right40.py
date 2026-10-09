"""Exact501-tone original/right40 H0 comparison; boundary residual is not a target."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    events=[json.loads(s) for s in (a.source/'execution.jsonl').read_text('utf-8').splitlines()]
    assert c['max_runs']==1 and c['no_retry'] and v['completed']
    assert events[-1]['status']=='COMPLETED' and events[-1]['traces']==1
    assert events[-1]['verification_sha256']==sha(a.source/'completed_verification.json')
    assert v['contract_sha256']==sha(a.source/'execution_contract.json')
    m=c['study_manifest'];g=c['groups'][0]
    new=a.source/'profile.h5';old=a.package/'baseline/native_H0.h5'
    assert sha(new)==v['groups'][0]['native_sha256'] and sha(old)==m['baseline_sha256']['native_H0.h5']
    assert sha(a.package/'manifest.json')==next(h for p,h in c['file_identities'].items() if Path(p).name=='manifest.json')
    arrays=[];direct_errors=[];zs=[];metadata=[]
    for path,shape in [(old,[8400,1700,1]),(new,[10000,1700,1])]:
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1'
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],shape)
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert h.attrs['Iterations']==20352 and h.attrs['dt']==m['dt_s']
            for node,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                np.testing.assert_allclose(h[node].attrs['Position'],pos,atol=1e-12,rtol=0)
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
            assert x.dtype==s.dtype==np.float64 and x.shape==s.shape==(20352,)
            assert np.isfinite(x).all() and np.isfinite(s).all()
            if arrays:assert s.tobytes()==source.tobytes()
            else:source=s
            info=dict(source_offset=float(h['srcs/src1/excitation'].attrs['TimeSampleOffset']),
                receiver_offset=float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset']))
            metadata.append(info);arrays.append(x)
        z,error=response(path,.025);zs.append(z);direct_errors.append(error)
    assert metadata[0]==metadata[1]
    z=np.column_stack([*zs,zs[1]-zs[0]])
    dt=m['dt_s'];native_time=np.arange(20352)*dt+metadata[0]['receiver_offset']
    native_metrics={}
    for name,bounds in dict(early=[0,120],return_packet=[250,450],late=[450,1100]).items():
        k=(native_time*1e9>=bounds[0])&(native_time*1e9<=bounds[1]);oldx,newx=[x[k] for x in arrays]
        native_metrics[name]=dict(bitwise_equal=oldx.tobytes()==newx.tobytes(),
            relative_L2=float(np.linalg.norm(newx-oldx)/np.linalg.norm(oldx)),
            difference_maxabs_V_m=float(abs(newx-oldx).max()),original_maxabs_V_m=float(abs(oldx).max()))
    metrics={};profiles={}
    for window in ['hann','blackman']:
        w=spectral_window(window,501);x,t=inverse(z,w);profiles[window]=x;rows={}
        for gate,bounds in m['gates_ns'].items():
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);p,q=x[k,0],x[k,1];corr=np.vdot(p,q)/(np.linalg.norm(p)*np.linalg.norm(q))
            rows[gate]=dict(relative_complex_L2=float(np.linalg.norm(q-p)/np.linalg.norm(p)),
                peak_ns=[float(t[k][abs(j).argmax()]*1e9) for j in [p,q]],
                complex_correlation_abs=float(abs(corr)),correlation_phase_deg=float(np.angle(corr,deg=True)),
                maxabs_change=float(abs(q-p).max()),original_maxabs=float(abs(p).max()))
        j=np.arange(7,4008,97);direct=np.exp(2j*np.pi*t[j,None]*FREQ)@(z*w[:,None])/501
        err=float(np.linalg.norm(x[j]-direct)/np.linalg.norm(direct));assert err<1e-9
        metrics[window]=dict(gates=rows,independent_inverse_relative_L2=err)
    a.out.mkdir(parents=True)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z)
        h.attrs['columns']='original_H0,right40_H0,right40_minus_original_H0_boundary_residual'
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for ax,bounds in [(axes[0,0],[250,450]),(axes[0,1],[0,120])]:
        k=(native_time*1e9>=bounds[0])&(native_time*1e9<=bounds[1])
        for data,label in zip(arrays,['原域H0','右侧+40m H0']):ax.plot(native_time[k]*1e9,data[k],lw=.8,label=label)
        ax.set(xlabel='原生时间 / ns（含源时间参考）',ylabel='Ez / V/m',title=f'原生接收 / {bounds[0]}–{bounds[1]}ns');ax.legend(fontsize=8)
    k=(native_time*1e9>=250)&(native_time*1e9<=450)
    axes[1,0].plot(native_time[k]*1e9,(arrays[1]-arrays[0])[k],color='#604080',lw=.8)
    axes[1,0].set(xlabel='原生时间 / ns',ylabel='变化Ez / V/m',title='右侧扩域−原域：单独变化尺度，不是底砂差场')
    for q,label in zip(zs,['原域H0','右侧+40m H0']):axes[1,1].semilogy(FREQ/1e6,abs(q),label=label,lw=1)
    axes[1,1].set(xlabel='精确采样频率 / MHz',ylabel='复响应模 / (V/m)/(A·m)',title='完整501频点复响应模（复相位仍保留）');axes[1,1].legend(fontsize=8)
    fig.suptitle('190m高损耗H0右侧边界对照 / 约8m航高 / 4.0.1原生FP64\n只向右延续40m边缘列；原内部、源、收发、材料、网格、时窗不变；无新H1／快照。')
    fig.savefig(a.out/'native_and_frequency.png',dpi=130);plt.close(fig)
    for window,x in profiles.items():
        for bounds,label in [([300,450],'return'),([450,1100],'late')]:
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);ticks=t[k]*1e9
            fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
            lim=float(abs(x[k,:2].real).max());delta_lim=float(abs(x[k,2].real).max())
            im=axes[0,0].imshow(x[k,:2].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[-.5,1.5,ticks[-1],ticks[0]])
            axes[0,0].set_xticks([0,1],['原域H0','右侧+40m H0']);axes[0,0].set(title='同站两配置：共同物理灰度，非空间B-scan',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[0,0],label='(V/m)/(A·m)')
            im=axes[0,1].imshow(x[k,2:3].real,cmap='gray',vmin=-delta_lim,vmax=delta_lim,aspect='auto',interpolation='nearest',extent=[-.5,.5,ticks[-1],ticks[0]])
            axes[0,1].set_xticks([0],['右侧扩域−原域']);axes[0,1].set(title='边界变化：单独尺度；不是底砂／纯噪声',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[0,1],label='变化 / (V/m)/(A·m)')
            for j,name in enumerate(['原域H0','右侧+40m H0']):
                axes[1,0].plot(ticks,abs(x[k,j]),label=name,lw=1)
                axes[1,1].plot(ticks,x[k,j].real,label=name,lw=1)
            for ax,title in [(axes[1,0],'复包络：保持绝对尺度'),(axes[1,1],'带符号实部：不作幅相／到时拟合')]:
                ax.set(title=title,xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');ax.legend(fontsize=8);ax.grid(alpha=.15)
            fig.suptitle(f'190m高损耗H0 / 右侧地质+40m / 原生FP64 / 精确20–170MHz、0.3MHz、501点 / {window}\n无AGC或插值；只有一个站位的两个域配置，差列不是H1−H0。')
            fig.savefig(a.out/f'right40_{window}_{label}_gray.png',dpi=130);plt.close(fig)
    with h5py.File(a.package/m['groups'][0]['geometry']) as h:data=h['data'][::8,::8,0].T
    fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
    ax.imshow(data,origin='lower',extent=[0,250,0,42.5],aspect='auto',cmap=ListedColormap(['#d9eef6','#d7b594','#b6bbc4']),vmin=0,vmax=2,interpolation='nearest')
    ax.axvspan(248,250,color='red',alpha=.25,label='新右侧PML')
    ax.axvspan(208,210,facecolor='none',edgecolor='#168eaa',hatch='///',label='原右侧PML位置（本次已移走）')
    ax.axvline(210,color='#168eaa',ls='--',lw=.8)
    ax.axvspan(0,2,color='red',alpha=.15)
    for lo,hi in [(0,2),(40.5,42.5)]:ax.axhspan(lo,hi,color='red',alpha=.15)
    ax.scatter([g['tx_m'][0],g['rx_m'][0]],[g['tx_m'][1],g['rx_m'][1]],c=['orange','red'],s=25,label='Tx橙／Rx红（位置不变）')
    ax.set(xlabel='模型local x / m（剖面里程=x+20）',ylabel='模型y / m',title='唯一新配置：190m高损耗H0，250×42.5m\n右侧40m按原边缘列延续；浅蓝空气，黄色粉质黏土，灰色泥岩；H0没有底砂')
    ax.legend(fontsize=8,loc='lower left');fig.savefig(a.out/'right40_geometry.png',dpi=130);plt.close(fig)
    result=dict(status='PASS_ONE_RIGHT40_H0_NATIVE_EXACT501_AND_INVERSE',script_sha256=sha(__file__),
        contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
        numerical_sha256=sha(a.numerical),native_sha256=[sha(old),sha(new)],source_bitwise_equal=True,
        time_sample_offsets=metadata,direct_DFT_relative_L2=direct_errors,native_metrics=native_metrics,
        frequency_relative_complex_L2=float(np.linalg.norm(z[:,2])/np.linalg.norm(z[:,0])),metrics=metrics,
        solver_runs=1,physical_path_certified=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],native_metrics=native_metrics,metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
