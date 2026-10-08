"""Actual nonflat matched-pair DC-loss sensitivity, exact phase-preserving SFCW."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response, inverse, FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    c=json.loads((a.source/'execution_contract.json').read_text('utf-8'))
    ver=json.loads((a.source/'completed_verification.json').read_text('utf-8'))
    assert ver['completed'] and ver['contract_sha256']==sha(a.source/'execution_contract.json')
    m=c['study_manifest'];reused=a.package/'reused_nonflat_H0.h5'
    assert sha(reused)==m['reused_H0']['native_sha256']
    values=[];identities=[]
    paths=[('baseline_H0',reused,m['reused_H0']['native_sha256'])]
    for g,r in zip(c['groups'],ver['groups']):
        assert g['id']==r['id']
        paths.append((g['id'],a.source/(g['id']+'.h5'),r['native_sha256']))
    for name,p,digest in paths:
        assert sha(p)==digest
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1'
            assert h['rxs/rx1/Ez'].dtype==np.float64
        z,err=response(p,reference_scale=.025);values.append(z)
        identities.append(dict(id=name,native_sha256=digest,independent_DFT_relative_L2=err))
    z=np.column_stack(values)
    z=np.column_stack([z,z[:,1]-z[:,0],z[:,3]-z[:,2]])
    labels=['baseline_H0','baseline_H1','low_H0','low_H1','baseline_delta','low_delta']
    metrics={};profiles={}
    for window in ['hann','blackman']:
        x,t=inverse(z,spectral_window(window,501));profiles[window]=x
        metrics[window]={}
        for gate in ['basal','deep']:
            lo,hi=m[gate+'_gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi)
            norm=lambda j:np.linalg.norm(x[keep,j])
            corr=lambda j,k:float(abs(np.vdot(x[keep,j],x[keep,k]))/(norm(j)*norm(k)))
            metrics[window][gate]=dict(baseline_H0_over_delta=float(norm(0)/norm(4)),
                                      changed_H0_over_delta=float(norm(2)/norm(5)),
                                      changed_delta_over_baseline_delta=float(norm(5)/norm(4)),
                                      changed_H0_over_baseline_H0=float(norm(2)/norm(0)),
                                      changed_H1_vs_delta_correlation=corr(3,5),
                                      baseline_H1_vs_delta_correlation=corr(1,4))
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,x in profiles.items():
        keep=(t*1e9>=340)&(t*1e9<=405)
        fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
        for ax,indices,title,names in [(axes[0],[0,1,2,3],'原始H0/H1总场，共同灰度',['σ.003 H0','σ.003 H1','σ.0003 H0','σ.0003 H1']),
                                      (axes[1],[4,5],'配对底砂差场，共同灰度',['σ.003 差场','σ.0003 差场'])]:
            data=x[keep][:,indices].real;limit=float(abs(data).max())
            im=ax.imshow(data,cmap='gray',vmin=-limit,vmax=limit,aspect='auto',extent=[-.5,len(indices)-.5,405,340])
            ax.set_xticks(range(len(indices)),names,rotation=18,fontsize=8)
            ax.set(title=title,ylabel='SFCW时间 / ns')
            for v in m['basal_gate_ns']:ax.axhline(v,color='green',ls=':',lw=1)
            fig.colorbar(im,ax=ax,label='共同参考 (V/m)/(A·m)')
        for j,label,style in [(1,'σ.003 总场','-'),(4,'σ.003 底砂差场','--'),(3,'σ.0003 总场','-'),(5,'σ.0003 底砂差场','--')]:
            axes[2].plot(t[keep]*1e9,abs(x[keep,j]),style,lw=1,label=label)
        for v in m['basal_gate_ns']:axes[2].axvline(v,color='green',ls=':',lw=1)
        axes[2].set(title='共同绝对幅度，无AGC或拟合',xlabel='SFCW时间 / ns',ylabel='复轮廓包络');axes[2].legend(fontsize=8)
        fig.suptitle('原非平v5模型 / X198.6m站位，AGL8m / 4.0.1 FP64 / '+window+'\n仅泥岩DC电导率改变；同站配置列非空间测线；绿色虚线为预声明窗；正式材料未换')
        fig.savefig(a.out/f'v401_nonflat_conductivity_{window}_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z)
        h.attrs['column_ids']=json.dumps(labels)
    result=dict(status='COMPLETED_NONFLAT_SINGLE_STATION_DC_SENSITIVITY_NOT_SITE_CALIBRATION',
                script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),
                verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),
                native=identities,metrics=metrics,gate_ns={k:m[k+'_gate_ns'] for k in ['basal','deep']},
                new_solves=3,reused_H0=1,AGC=False,fit=False,production_materials_changed=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(metrics))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
