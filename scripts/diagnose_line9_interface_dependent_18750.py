"""Exploratory shared-background differences; not pure echoes or clean labels."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import inverse,FREQ
from audit_line9_v401_delivery import native_response


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    pub=read(a.public/'analysis.json');audit=read(a.public/'independent_audit.json')
    assert audit['analysis_sha256']==sha(a.public/'analysis.json') and pub['numerical_sha256']==sha(a.numerical)
    with h5py.File(a.numerical) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);z=h['response'][:]
    paths=[a.package/'baseline/native_H0.h5',a.source/'no_cover_contrast_H0/profile.h5',a.source/'cover_bottom_down050_H0/profile.h5']
    independently=np.column_stack([native_response(p,.025,h) for p,h in zip(paths,pub['native_sha256'])])
    # Same homogeneous-cover reference for both; model differences include interactions.
    q=np.column_stack([z[:,0]-z[:,1],z[:,2]-z[:,1]])
    directq=np.column_stack([independently[:,0]-independently[:,1],independently[:,2]-independently[:,1]])
    error=float(np.linalg.norm(q-directq)/np.linalg.norm(directq));assert error<1e-8
    metrics={};profiles={}
    gates=dict(main_fixed=[300,450],first_replication=[120,250])
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();x,t=inverse(q,w);profiles[window]=x;rows={}
        for gate,(lo,hi) in gates.items():
            k=(t*1e9>=lo)&(t*1e9<=hi);direct=np.exp(2j*np.pi*t[k,None]*FREQ)@(directq*w[:,None])/501
            assert np.linalg.norm(x[k]-direct)/np.linalg.norm(direct)<1e-7
            peaks=[float(t[k][abs(x[k,j]).argmax()]*1e9) for j in range(2)]
            dp=[float(t[k][abs(direct[:,j]).argmax()]*1e9) for j in range(2)]
            np.testing.assert_allclose(peaks,dp,atol=1e-9,rtol=0)
            rows[gate]=dict(peak_ns=peaks,peak_change_ns=peaks[1]-peaks[0],complex_L2_ratio=float(np.linalg.norm(x[k,1])/np.linalg.norm(x[k,0])),direct_inverse_relative_L2=float(np.linalg.norm(x[k]-direct)/np.linalg.norm(direct)))
        metrics[window]=rows
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,x in profiles.items():
        for lo,hi,tag in [(120,450,'full'),(300,450,'late')]:
            k=(t*1e9>=lo)&(t*1e9<=hi);tt=t[k]*1e9;fig,axes=plt.subplots(2,1,figsize=(13,8),layout='constrained')
            lim=float(abs(x[k].real).max());im=axes[0].imshow(x[k].real,cmap='gray',vmin=-lim,vmax=lim,extent=[-.5,1.5,tt[-1],tt[0]],aspect='auto',interpolation='nearest')
            labels=['原域−去对比参考','底下移0.5m−同一去对比参考'];axes[0].set_xticks([0,1],labels);axes[0].set_ylabel('SFCW时间 / ns');fig.colorbar(im,ax=axes[0],label='(V/m)/(A·m)；共同尺度')
            for j,label in enumerate(labels):axes[1].plot(tt,abs(x[k,j]),label=label,lw=1)
            axes[1].legend();axes[1].set(xlabel='SFCW时间 / ns',ylabel='差分复包络 / (V/m)/(A·m)');axes[1].grid(alpha=.15)
            fig.suptitle(f'187.5m高损耗H0／{window}／完整501点\n界面对比相关差分包含传播与交互变化，不是纯多次波；无移峰／缩放拟合；120–250ns为预先固定复核窗')
            fig.savefig(a.out/f'interface_dependent_{window}_{tag}.png',dpi=140);plt.close(fig)
    result=dict(status='PASS_EXPLORATORY_INTERFACE_DIFFERENCES_INDEPENDENT_DFT_AND_INVERSE',script_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),independent_audit_sha256=sha(a.public/'independent_audit.json'),numerical_sha256=sha(a.numerical),independent_difference_DFT_relative_L2=error,metrics=metrics,
        residual_definition=['original_H0_minus_no_cover_contrast_H0','cover_bottom_down050_H0_minus_same_no_cover_contrast_H0'],gates_ns=gates,first_gate_status='Predeclared replication of190m exploratory120to250ns window; not an acceptance threshold',solver_runs=0,
        limits='Shared background differences contain material propagation and scattering interactions. Peak changes do not prove unique bounce count; no alignment, scaling or independent field/3D validation.')
    (a.out/'diagnostic.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','public','package','numerical','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
