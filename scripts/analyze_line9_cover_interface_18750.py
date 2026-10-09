"""Exact501-tone interface controls with common physical scales, no fitting."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=c['study_manifest']
    assert v['completed'] and c['max_runs']==2 and c['no_retry'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.package/'manifest.json')==next(h for p,h in c['file_identities'].items() if Path(p).name=='manifest.json')
    groups=c['groups'];names=['original_H0']+[g['id'] for g in groups]
    paths=[a.package/'baseline/native_H0.h5']+[a.source/g['id']/'profile.h5' for g in groups]
    hashes=[m['baseline_sha256']['native_H0.h5']]+[r['native_sha256'] for r in v['groups']]
    zs=[];raw=[];errors=[];positions=groups[0]
    for i,(path,hsh) in enumerate(zip(paths,hashes)):
        assert sha(path)==hsh
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==20352 and h.attrs['dt']==m['dt_s']
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            for node,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:np.testing.assert_allclose(h[node].attrs['Position'],positions[key],atol=1e-12,rtol=0)
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
            assert x.dtype==s.dtype==np.float64 and x.shape==s.shape==(20352,) and np.isfinite(x).all() and np.isfinite(s).all()
            offsets=[float(h['srcs/src1/excitation'].attrs['TimeSampleOffset']),float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset'])]
            if i:assert s.tobytes()==source.tobytes() and offsets==oldoffsets
            else:source=s;oldoffsets=offsets
            raw.append(x)
        z,err=response(path,.025);zs.append(z);errors.append(err)
    z=np.column_stack(zs);profiles={};metrics={};native={};times=np.arange(20352)*m['dt_s']+oldoffsets[1]
    for j,name in enumerate(names[1:],1):
        native[name]={}
        for gate,(lo,hi) in dict(early=[0,120],return_packet=[250,450],late=[450,1100]).items():
            k=(times*1e9>=lo)&(times*1e9<=hi);old,new=raw[0][k],raw[j][k]
            native[name][gate]=dict(new_over_original_L2=float(np.linalg.norm(new)/np.linalg.norm(old)),relative_change_L2=float(np.linalg.norm(new-old)/np.linalg.norm(old)),bitwise_equal=new.tobytes()==old.tobytes())
    for window in ['hann','blackman']:
        w=spectral_window(window,501);x,t=inverse(z,w);profiles[window]=x;rows={}
        for j,name in enumerate(names[1:],1):
            rows[name]={}
            for gate,(lo,hi) in m['gates_ns'].items():
                k=(t*1e9>=lo)&(t*1e9<=hi);old,new=x[k,0],x[k,j];corr=np.vdot(old,new)/(np.linalg.norm(old)*np.linalg.norm(new))
                jo=np.flatnonzero(k)[abs(old).argmax()];jn=np.flatnonzero(k)[abs(new).argmax()]
                rows[name][gate]=dict(new_over_original_complex_L2=float(np.linalg.norm(new)/np.linalg.norm(old)),relative_complex_change_L2=float(np.linalg.norm(new-old)/np.linalg.norm(old)),complex_correlation_abs=float(abs(corr)),phase_deg=float(np.angle(corr,deg=True)),peak_ns=[float(t[jo]*1e9),float(t[jn]*1e9)],new_over_old_at_original_peak=float(abs(x[jo,j])/abs(x[jo,0])),original_peak_time_coordinate_index=int(jo))
            k=np.arange(0,4008,67);direct=np.exp(2j*np.pi*t[k,None]*FREQ)@(z*w[:,None])/501
            assert np.linalg.norm(x[k]-direct)/np.linalg.norm(direct)<1e-9
        metrics[window]=rows
    # Unfitted local vertical group-delay hypotheses, not true nonlocal event times.
    material=a.package/'baseline/materials.json';assert sha(material)==m['baseline_sha256']['materials.json']
    cover=read(material)['materials']['material_001_cover'];assert cover['model']=='debye' and len(cover['poles'])==1
    base=cover['base'];pole=cover['poles'][0]
    epsinf=base['relative_permittivity'];sigma=base['electric_conductivity_s_per_m'];delta=pole['relative_permittivity_difference'];tau=pole['relaxation_time_s']
    assert base['relative_permeability']==1 and base['magnetic_conductivity_s_per_m']==0
    eps0=8.8541878128e-12;c0=299792458.
    def beta(f):
        omega=2*np.pi*f;epsilon=epsinf+delta/(1+1j*omega*tau)-1j*sigma/(eps0*omega)
        return (omega/c0*np.sqrt(epsilon)).real
    derivatives=[]
    for step in [100.,1000.,10000.]:derivatives.append((beta(FREQ+step)-beta(FREQ-step))/(4*np.pi*step))
    np.testing.assert_allclose(derivatives[0],derivatives[1],rtol=1e-8,atol=0)
    nominal=(2*.5*np.array(derivatives[1])*1e9)
    prediction=dict(basis='Local vertical extra0.5m cover round-trip group delay from declared Debye+DC; not nonlocal geometry prediction or fitted alignment',one_bottom_roundtrip_ns_range=[float(nominal.min()),float(nominal.max())],two_bottom_roundtrips_ns_range=[float(2*nominal.min()),float(2*nominal.max())],one_bottom_roundtrip_at95MHz_ns=float(nominal[250]),two_bottom_roundtrips_at95MHz_ns=float(2*nominal[250]),derivative_step_Hz=[100,1000,10000])
    a.out.mkdir(parents=True)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z);h.attrs['columns']=','.join(names)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    labels=['原域H0','覆盖层底对比移除H0','覆盖层底下移0.5m H0']
    for window,x in profiles.items():
        for bounds,tag in [([280,450],'return'),([0,120],'early'),([450,1100],'late')]:
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);tt=t[k]*1e9
            fig,axes=plt.subplots(2,1,figsize=(13,8),layout='constrained')
            lim=float(abs(x[k].real).max());im=axes[0].imshow(x[k].real,cmap='gray',vmin=-lim,vmax=lim,extent=[-.5,2.5,tt[-1],tt[0]],aspect='auto',interpolation='nearest')
            axes[0].set_xticks(range(3),labels);axes[0].set(ylabel='SFCW时间 / ns',title='同一站三个配置：共同灰度；非空间B-scan')
            fig.colorbar(im,ax=axes[0],label='(V/m)/(A·m)')
            for j,label in enumerate(labels):axes[1].plot(tt,abs(x[k,j]),lw=1,label=label)
            axes[1].set(xlabel='SFCW时间 / ns',ylabel='复包络 / (V/m)/(A·m)',title='绝对复包络：无AGC、逐道归一、移峰或幅相拟合');axes[1].legend();axes[1].grid(alpha=.15)
            fig.suptitle(f'187.5m高损耗H0／约8m AGL／4.0.1原生FP64／{window}\n完整20–170MHz、0.3MHz、501点；底部砂岩不在H0中；界面控制不是匹配底砂差场')
            fig.savefig(a.out/f'interface_{window}_{tag}.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(13,7),layout='constrained')
    for ax,bounds in zip(axes,[[0,120],[250,450]]):
        k=(times*1e9>=bounds[0])&(times*1e9<=bounds[1])
        for x,label in zip(raw,labels):ax.plot(times[k]*1e9,x[k],lw=.8,label=label)
        ax.set(xlabel='原生Ricker时间 / ns（含源时延）',ylabel='Ez / V/m');ax.legend();ax.grid(alpha=.15)
    fig.suptitle('原生宽频接收对照；与正式SFCW结果分域比较');fig.savefig(a.out/'native_controls.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,5),sharey=True,layout='constrained')
    with h5py.File(a.package/'baseline/geometry.h5') as h:original=h['data'][::8,::8,0].T
    geom=[original]
    for g in m['groups']:
        with h5py.File(a.package/g['geometry']) as h:geom.append(h['data'][::8,::8,0].T)
    for ax,x,label in zip(axes,geom,labels):
        ax.imshow(x,origin='lower',aspect='auto',extent=[0,210,0,42.5],interpolation='nearest',vmin=0,vmax=2,cmap=ListedColormap(['#d9eef6','#d7b594','#b6bbc4']))
        ax.scatter([positions['tx_m'][0],positions['rx_m'][0]],[positions['tx_m'][1],positions['rx_m'][1]],c=['orange','red'],s=15)
        ax.set(title=label,xlabel='模型x / m');ax.set_xlim(145,180);ax.set_ylim(15,42.5)
    axes[0].set_ylabel('模型y / m');fig.suptitle('真实体素局部图：浅蓝空气、黄色覆盖层、灰色泥岩；地表、源、接收不变')
    fig.savefig(a.out/'interface_geometry.png',dpi=140);plt.close(fig)
    result=dict(status='PASS_TWO_INTERFACE_CONTROLS_NATIVE_EXACT501',script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),native_sha256=hashes,source_bitwise_equal=True,direct_DFT_relative_L2=errors,metrics=metrics,native_metrics=native,nominal_delay_hypotheses=prediction,solver_runs=2,physical_path_certified=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],metrics=metrics,nominal_delay_hypotheses=prediction)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','package','out','numerical']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
