"""Exact501-tone boundary comparisons for the worst-clearance220m station."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window

NAMES=['base','top20','right40']
LABELS=['原210×42.5m域','顶部增加20m空气','右侧延续40m地质']

def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read=lambda p:json.loads(p.read_text('utf-8'));m=read(a.package/'manifest.json');c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    records={r['id']:r for r in v['groups']};groups={g['id']:g for g in m['groups']};z=[];native=[];prefix={};raw={}
    for name in NAMES:
        pair=[]
        for role in ['H0','H1']:
            sid=name+'_'+role;p=a.source/(sid+'.h5');assert sha(p)==records[sid]['native_sha256']
            with h5py.File(p) as h:
                assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64 and h['rxs/rx1/Ez'].shape==(20352,)
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],groups[sid]['native_shape'])
                for node,r in [('srcs/src1','tx'),('rxs/rx1','rx')]:np.testing.assert_allclose(h[node].attrs['Position'],groups[sid][r+'_m'],rtol=0,atol=1e-11)
                samples=h['srcs/src1/excitation/samples'][:]
                if not native:source=samples
                else:np.testing.assert_array_equal(source,samples)
                raw[sid]=h['rxs/rx1/Ez'][:]
            q,error=response(p,.025);pair.append(q);native.append(dict(id=sid,native_sha256=sha(p),direct_DFT_relative_L2=error))
        z.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(z);profiles={};metrics={}
    for name in NAMES[1:]:
        for role in ['H0','H1']:
            k=int(120e-9/m['dt_s']);a0=raw['base_'+role][:k];b0=raw[name+'_'+role][:k]
            prefix[name+'_'+role]=dict(first120ns_bitwise_equal=bool(np.array_equal(a0,b0)),relative_L2=float(np.linalg.norm(b0-a0)/np.linalg.norm(a0)))
    for window in ['hann','blackman']:
        x,t=inverse(z,spectral_window(window,501));profiles[window]=x;metrics[window]={}
        for gate,(lo,hi) in m['gates_ns'].items():
            k=(t*1e9>=lo)&(t*1e9<=hi);b=x[k,:3];norm=np.linalg.norm;rows={}
            for j,name in enumerate(NAMES):
                q=x[k,3*j:3*j+3]
                rows[name]=dict(H0_over_delta=float(norm(q[:,0])/norm(q[:,2])),H1_vs_delta_complex_correlation=float(abs(np.vdot(q[:,1],q[:,2]))/(norm(q[:,1])*norm(q[:,2]))),
                    delta_peak_ns=float(t[k][np.argmax(abs(q[:,2]))]*1e9),H1_peak_ns=float(t[k][np.argmax(abs(q[:,1]))]*1e9),
                    relative_complex_change={r:float(norm(q[:,i]-b[:,i])/norm(b[:,i])) for i,r in enumerate(['H0','H1','delta'])},
                    change_over_original_delta={r:float(norm(q[:,i]-b[:,i])/norm(b[:,2])) for i,r in enumerate(['H0','H1','delta'])})
            metrics[window][gate]=rows
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,x in profiles.items():
        for bounds,label in [([180,450],'target'),([450,1100],'late')]:
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);ticks=t[k]*1e9;lim=float(abs(x[k].real).max())
            fig,axes=plt.subplots(2,2,figsize=(16,9),layout='constrained')
            for ax,i,title in [(axes[0,0],1,'H1原始总场'),(axes[0,1],2,'H1−H0底砂替换差场')]:
                im=ax.imshow(x[k,i::3].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[-.5,2.5,ticks[-1],ticks[0]])
                ax.set_xticks(range(3),LABELS,fontsize=9);ax.set(title=title+'：共同绝对灰度',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='(V/m)/(A·m)')
                if label=='target':
                    for edge in m['gates_ns']['basal']:ax.axhline(edge,c='green',ls=':',lw=.8)
            changes=np.column_stack([x[k,3*j+i]-x[k,i] for j in [1,2] for i in [0,1,2]]);change_lim=float(abs(changes.real).max())
            im=axes[1,0].imshow(changes.real,cmap='gray',vmin=-change_lim,vmax=change_lim,aspect='auto',interpolation='nearest',extent=[-.5,5.5,ticks[-1],ticks[0]])
            axes[1,0].set_xticks(range(6),['顶扩H0','顶扩H1','顶扩差场','右扩H0','右扩H1','右扩差场'],fontsize=9)
            axes[1,0].set(title='扩域−原域：变化分量单独共同标尺',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[1,0],label='变化 / (V/m)/(A·m)')
            for j,name in enumerate(LABELS):axes[1,1].plot(ticks,abs(x[k,3*j+2]),lw=1,label=name+' 底砂差场')
            axes[1,1].set(title='绝对包络，不做延时或幅相拟合',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');axes[1,1].legend(fontsize=8)
            fig.suptitle('220m批前边界核查 / 约8m航高 / 4.0.1原生FP64 / '+window+'\n低损耗诊断配方；同站配置列非空间测线；原域顶部PML入口距接收器0.325m；批量尚未启动')
            fig.savefig(a.out/f'line_prebatch_{window}_{label}_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z)
    result=dict(status='COMPLETED_SIX_PREBATCH_NOT_BULK_NOT_GLOBAL_CONVERGENCE',script_sha256=sha(__file__),contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),native=native,metrics=metrics,native_early_prefix=prefix,new_solves=6,bulk_started=False,production_materials_changed=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],metrics=metrics)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
