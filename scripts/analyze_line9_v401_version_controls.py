"""Audit exact-input 4.0.0/4.0.1 traces and compare phase-preserving SFCW."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha

FREQ = 20e6 + np.arange(501)*300000.


def relative(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def response(path, reference_scale=None):
    s=sf.load_source(path); r=sf.load_receiver(path,'/rxs/rx1','Ez')
    q=sf.direct_frequency_response(s,r,FREQ,tail_taper_fraction=0)
    assert q.source_valid.all() and s.spatial_scale>0
    scale=s.spatial_scale if reference_scale is None else reference_scale
    z=q.response.reshape(501)/scale
    direct=np.empty(501,complex)
    for k in range(0,501,16):
        f=FREQ[k:k+16,None]
        direct[k:k+16]=(r.dt*(np.exp(-2j*np.pi*f*r.times)@r.samples))/(
            s.dt*(np.exp(-2j*np.pi*f*s.times)@s.samples))/scale
    error=relative(z,direct); assert error<1e-9
    return z,error


def inverse(z,w):
    padded=np.zeros((4008,)+z.shape[1:],complex)
    padded[:501]=z*w[:,None]
    t=np.arange(4008)/(4008*300000.)
    return np.fft.ifft(padded,axis=0)*8*np.exp(2j*np.pi*20e6*t[:,None]),t


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    oldc=json.loads((a.old/'execution_contract.json').read_text('utf-8'))
    oldv=json.loads((a.old/'completed_verification.json').read_text('utf-8'))
    c=json.loads((a.new/'execution_contract.json').read_text('utf-8'))
    v=json.loads((a.new/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and oldv['completed']
    assert v['contract_sha256']==sha(a.new/'execution_contract.json')
    assert oldv['contract_sha256']==sha(a.old/'execution_contract.json')
    names=[g['id'] for g in c['groups']]
    assert names==[g['id'] for g in oldc['groups']]==['nonflat_ricker_H0','flat_ricker_H0','flat_ricker_H1']
    identities=[];zs=[[],[]];errors=[[],[]];native=[]
    for j,name in enumerate(names):
        for key in ['input_sha256','geometry_sha256','material_sha256']:
            assert c['groups'][j][key]==oldc['groups'][j][key], (name,key)
        pp=[a.old/(name+'.h5'),a.new/(name+'.h5')]
        xx=[]
        for version,p,record in zip(['4.0.0','4.0.1'],pp,[oldv['groups'][j],v['groups'][j]]):
            assert sha(p)==record['native_sha256']
            with h5py.File(p) as h:
                assert str(h.attrs['gprMax'])==version
                x=h['rxs/rx1/Ez'][:];source=h['srcs/src1/excitation/samples'][:]
                assert x.dtype==source.dtype==np.float64 and x.shape==(20352,) and np.isfinite(x).all()
                xx.append(x)
                if version=='4.0.0':
                    prior_source=source;meta={k:h.attrs[k] for k in ['dt','Iterations','dx_dy_dz','nx_ny_nz']}
                else:
                    np.testing.assert_array_equal(source,prior_source)
                    for k,value in meta.items():np.testing.assert_array_equal(h.attrs[k],value)
            z,error=response(p);i=0 if version=='4.0.0' else 1;zs[i].append(z);errors[i].append(error)
        native.append(dict(id=name,array_equal=bool(np.array_equal(*xx)),relative_L2=relative(xx[1],xx[0]),
                           native_sha256=[sha(p) for p in pp]))
        identities.append(dict(id=name,input_sha256=c['groups'][j]['input_sha256'],
                               geometry_sha256=c['groups'][j]['geometry_sha256'],material_sha256=c['groups'][j]['material_sha256']))
    z=[np.column_stack(x) for x in zs]
    z=[np.column_stack([x,x[:,2]-x[:,1]]) for x in z]
    metrics={};profiles={}
    for window in ['hann','blackman']:
        w=sf.spectral_window(window,501);old,t=inverse(z[0],w);new,_=inverse(z[1],w);profiles[window]=(old,new)
        k=np.arange(7,4008,101)
        direct=np.exp(2j*np.pi*t[k,None]*FREQ)@(z[1]*w[:,None])/501
        inv_error=relative(new[k],direct);assert inv_error<1e-9
        rows={}
        for gate,bounds in [('basal',oldc['study_manifest']['basal_gate_ns']),('deep',[300,450]),('early',[0,120])]:
            mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1])
            rows[gate]={name:relative(new[mask,j],old[mask,j]) for j,name in enumerate(names+['flat_delta'])}
        metrics[window]=dict(relative_version_changes=rows,independent_inverse_relative_L2=inv_error)
    ratio=z[1][:,3]/z[0][:,3];selected=[0,83,250,417,500]
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,(old,new) in profiles.items():
        fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,bounds,title in [(axs[0,0],[150,450],'三个配置的原始总场/剩余场'),(axs[0,1],oldc['study_manifest']['basal_gate_ns'],'预冻结底砂时间窗')]:
            mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);data=np.column_stack([old[mask,:3],new[mask,:3]]).real;lim=float(abs(data).max())
            im=ax.imshow(data,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[-.5,5.5,bounds[1],bounds[0]])
            ax.set_xticks(range(6),['旧非平H0','旧平层H0','旧平层H1','新非平H0','新平层H0','新平层H1'],rotation=15,fontsize=8)
            ax.set(title=title+'；六列共用绝对灰度',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='(V/m)/(A·m)')
        mask=(t*1e9>=340)&(t*1e9<=405)
        for value,label,style in [(old[:,3],'4.0.0 平层底砂差场','-'),(new[:,3],'4.0.1 平层底砂差场','--')]:axs[1,0].plot(t[mask]*1e9,value[mask].real,style,label=label,lw=1)
        axs[1,0].axvline(371.756487,color='gray',ls=':',label='既有一次反射模型到时')
        axs[1,0].set(xlabel='SFCW时间 / ns',ylabel='双极 / (V/m)/(A·m)',title='原幅度与相位直接对照，无拟合');axs[1,0].legend(fontsize=8)
        diff=new-old
        for j,label in [(0,'非平H0'),(1,'平层H0'),(2,'平层H1'),(3,'平层底砂差场')]:
            axs[1,1].plot(t[mask]*1e9,abs(diff[mask,j]),label=label,lw=1)
        axs[1,1].set(xlabel='SFCW时间 / ns',ylabel='复数版本差的包络 / (V/m)/(A·m)',title='版本差单独标尺，不作为新地下反射');axs[1,1].legend(fontsize=8)
        fig.suptitle(f'v5同站198.6m / AGL8m / {window} / 20–170MHz精确501频点\n4.0.0与4.0.1，输入/体素/材料逐字节一致；配置列不是连续空间测线')
        fig.savefig(a.out/f'v401_version_{window}_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('time_s',data=t)
        h.create_dataset('v400_response',data=z[0]);h.create_dataset('v401_response',data=z[1])
        for w,(old,new) in profiles.items():h.create_dataset('v400_'+w,data=old);h.create_dataset('v401_'+w,data=new)
    result=dict(status='PASS_EXACT_INPUT_VERSION_COMPARISON_NOT_FIELD_VALIDATION',script_sha256=sha(__file__),
                old_contract_sha256=sha(a.old/'execution_contract.json'),new_contract_sha256=sha(a.new/'execution_contract.json'),
                old_verification_sha256=sha(a.old/'completed_verification.json'),new_verification_sha256=sha(a.new/'completed_verification.json'),
                scientific_identities=identities,native=native,independent_DFT_relative_L2=errors,metrics=metrics,
                flat_delta_version_ratio=dict(frequency_MHz=(FREQ[selected]/1e6).tolist(),amplitude=abs(ratio[selected]).tolist(),phase_deg=np.angle(ratio[selected],deg=True).tolist()),
                numerical_sha256=sha(a.numerical),AGC=False,tail_taper=False,phase_amplitude_time_fit=False,
                limits='Three selected station configurations only; no whole-line, finite-3D antenna, field/material or mesh-convergence validation.')
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(native=native,metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['old','new','out','numerical']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
