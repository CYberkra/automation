"""Phase-preserving SFCW sensitivity of low-loss v5 to side/bottom extension."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window

IDS=['base','sides40','bottom20']
NAMES=['原域','左右各扩40m','底部扩20m']


def metrics(x,t,gates):
    result={}
    for gate,(lo,hi) in gates.items():
        k=(t*1e9>=lo)&(t*1e9<=hi);b=x[k,:3];rows={}
        assert np.linalg.norm(b[:,2])>0
        for j,name in enumerate(IDS):
            q=x[k,3*j:3*j+3];norm=np.linalg.norm
            rows[name]=dict(H0_over_delta=float(norm(q[:,0])/norm(q[:,2])),
                H1_vs_delta_correlation=float(abs(np.vdot(q[:,1],q[:,2]))/(norm(q[:,1])*norm(q[:,2]))),
                delta_norm_over_original=float(norm(q[:,2])/norm(b[:,2])),
                delta_peak_ns=float(t[k][np.argmax(abs(q[:,2]))]*1e9),
                relative_complex_change={r:float(norm(q[:,i]-b[:,i])/norm(b[:,i])) for i,r in enumerate(['H0','H1','delta'])},
                changes_over_original_delta={r:float(norm(q[:,i]-b[:,i])/norm(b[:,2])) for i,r in enumerate(['H0','H1','delta'])})
        result[gate]=rows
    return result


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=read(a.package/'manifest.json')
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    paths={r['id']:(a.package/r['path'],r['native_sha256']) for r in m['reused']}
    for g,r in zip(c['groups'],v['groups']):
        assert g['id']==r['id'];paths[g['id']]=(a.source/(g['id']+'.h5'),r['native_sha256'])
    vals=[];native=[]
    for name in IDS:
        pair=[]
        for role in ['H0','H1']:
            label=name+'_'+role;p,digest=paths[label];assert sha(p)==digest
            with h5py.File(p) as h:assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64
            z,err=response(p,.025);pair.append(z);native.append(dict(id=label,native_sha256=digest,independent_DFT_relative_L2=err))
        vals.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(vals);profiles={};allmetrics={}
    gates={key:m[key+'_gate_ns'] for key in ['basal','deep','early']};gates['late']=[450.,1100.]
    for window in ['hann','blackman']:
        x,t=inverse(z,spectral_window(window,501));profiles[window]=x;allmetrics[window]=metrics(x,t,gates)
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,x in profiles.items():
        for bounds,suffix in [([300,450],'target'),([450,1100],'late')]:
            lo,hi=bounds;k=(t*1e9>=lo)&(t*1e9<=hi)
            fig,axes=plt.subplots(2,2,figsize=(16,9),layout='constrained')
            lim=float(abs(x[k].real).max())
            for ax,offset,title in [(axes[0,0],1,'H1原始总场'),(axes[0,1],2,'H1−H0底砂替换差场')]:
                im=ax.imshow(x[k,offset::3].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[-.5,2.5,hi,lo])
                ax.set_xticks(range(3),NAMES);ax.set(title=title+'；两面板共用绝对灰度',ylabel='SFCW时间 / ns')
                if suffix=='target':
                    for v0 in gates['basal']:ax.axhline(v0,color='green',ls=':',lw=1)
                fig.colorbar(im,ax=ax,label='共同参考 (V/m)/(A·m)')
            changes=np.column_stack([x[k,3*j+i]-x[k,i] for j in [1,2] for i in [0,1,2]])
            limit=float(abs(changes.real).max());assert limit>0
            im=axes[1,0].imshow(changes.real,cmap='gray',vmin=-limit,vmax=limit,aspect='auto',interpolation='nearest',extent=[-.5,5.5,hi,lo])
            axes[1,0].set_xticks(range(6),['侧扩H0变化','侧扩H1变化','侧扩差场变化','底扩H0变化','底扩H1变化','底扩差场变化'],rotation=15,fontsize=8)
            axes[1,0].set(title='扩域−原域：变化分量另设共同标尺',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[1,0],label='变化 / (V/m)/(A·m)')
            for j,label in enumerate(NAMES):axes[1,1].plot(t[k]*1e9,abs(x[k,3*j+2]),lw=1,label=label+' 底砂差场')
            for j,label in [(1,'侧扩−原域'),(2,'底扩−原域')]:axes[1,1].plot(t[k]*1e9,abs(x[k,3*j+2]-x[k,2]),ls='--',lw=1,label=label+' 差场变化')
            axes[1,1].set(title='绝对包络，无AGC/到时或幅相拟合',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');axes[1,1].legend(fontsize=8)
            fig.suptitle('原非平v5 / X198.6m / AGL8m / 4.0.1 FP64 / '+window+'\n低跨度0.3、DC0.3mS/m；同站配置列非空间B-scan；原区域完全保留；正式材料未换')
            fig.savefig(a.out/f'low_loss_boundaries_{window}_{suffix}_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z)
        h.attrs['column_ids']=json.dumps([name+'_'+r for name in IDS for r in ['H0','H1','delta']])
    result=dict(status='COMPLETED_LOW_LOSS_SIDE_BOTTOM_SENSITIVITY_NOT_GLOBAL_CONVERGENCE',script_sha256=sha(__file__),
        contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
        numerical_sha256=sha(a.numerical),native=native,metrics=allmetrics,gate_ns=gates,new_solves=4,reused_native=2,
        AGC=False,fit=False,production_materials_changed=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(allmetrics))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','package','out','numerical']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
