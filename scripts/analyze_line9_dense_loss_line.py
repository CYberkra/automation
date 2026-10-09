"""Exact-band paired-loss line views; explicitly mask uncomputed positions."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window

VARIANTS=['high','low'];NAMES={'high':'原高损耗泥岩','low':'低损耗诊断泥岩'}
def read(p):return json.loads(p.read_text('utf-8'))

def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    m=read(a.package/'manifest.json');c=read(a.source/'execution_contract.json');v=read(a.source/'snapshot_manifest.json')
    assert c['study_manifest']==m and v['contract_sha256']==sha(a.source/'execution_contract.json')
    expected={r['id']:r for r in v['records']};expected.update({r['id']:r for r in m['reused']})
    items={r['id']:r for r in m['groups']+m['reused']};spectra={};native=[];source=None
    assert set(p.stem for p in a.source.glob('*.h5'))==set(expected)
    for sid,record in expected.items():
        path=a.source/(sid+'.h5');g=items[sid];assert sha(path)==record['native_sha256']
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64 and h['rxs/rx1/Ez'].shape==(20352,)
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
            for node,role in [('srcs/src1','tx'),('rxs/rx1','rx')]:np.testing.assert_allclose(h[node].attrs['Position'],g[role+'_m'],rtol=0,atol=1e-11)
            s=h['srcs/src1/excitation/samples'][:]
            if source is None:source=s
            else:np.testing.assert_array_equal(source,s)
        spectra[sid],err=response(path,.025);native.append(dict(id=sid,native_sha256=sha(path),direct_DFT_relative_L2=err))
    ids=list(spectra);z=np.column_stack([spectra[k] for k in ids]);profiles={};peaks={};anchor_metrics={};availability={}
    for variant in VARIANTS:
        availability[variant]=[f'{variant}_x{round(x*100):05d}_H1' in spectra for x in m['chainage_m']]
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    cmap=plt.get_cmap('gray').copy();cmap.set_bad('#d6e5ee')
    for window in ['hann','blackman']:
        x,t=inverse(z,spectral_window(window,501));profiles[window]=x;lookup={sid:x[:,j] for j,sid in enumerate(ids)}
        peaks[window]={};anchor_metrics[window]={};lines={};deltas={}
        for variant in VARIANTS:
            line=np.full((len(t),21),np.nan+1j*np.nan);delta=np.full_like(line,np.nan+1j*np.nan);pks=[];metrics=[]
            for j,s in enumerate(m['stations']):
                prefix=f'{variant}_x{round(s["chainage_m"]*100):05d}_';sid=prefix+'H1';h0=prefix+'H0';lo,hi=s['templates'][variant]['basal_gate_ns'];k=(t*1e9>=lo)&(t*1e9<=hi)
                if sid in lookup:
                    line[:,j]=lookup[sid];q=lookup[sid][k];pks.append(dict(chainage_m=s['chainage_m'],H1_peak_ns=float(t[k][np.argmax(abs(q))]*1e9),local_template_peak_ns=s['templates'][variant]['peak_ns']))
                if sid in lookup and h0 in lookup:
                    delta[:,j]=lookup[sid]-lookup[h0];q=lookup[sid][k];b=lookup[h0][k];d=q-b;norm=np.linalg.norm
                    metrics.append(dict(chainage_m=s['chainage_m'],H0_over_delta=float(norm(b)/norm(d)),H1_vs_delta_complex_correlation=float(abs(np.vdot(q,d))/(norm(q)*norm(d))),delta_peak_ns=float(t[k][np.argmax(abs(d))]*1e9),H1_peak_ns=float(t[k][np.argmax(abs(q))]*1e9),delta_norm=float(norm(d)),H0_norm=float(norm(b)),H1_norm=float(norm(q))))
            lines[variant]=line;deltas[variant]=delta;peaks[window][variant]=pks;anchor_metrics[window][variant]=metrics
        for bounds,tag in [([180,450],'main'),([450,1100],'late')]:
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);ticks=t[k]*1e9
            common=max(float(np.nanmax(abs(lines[v][k].real))) if np.isfinite(lines[v][k]).any() else 0. for v in VARIANTS)
            fig,axes=plt.subplots(3,2,figsize=(15,13),layout='constrained')
            for j,variant in enumerate(VARIANTS):
                for row,data,title,lim in [(0,lines[variant],'H1原始总场：两种损耗共同物理灰度',common),(1,lines[variant],'H1原始总场：该配方单独放大，不能跨图比较强弱',None),(2,deltas[variant],'H1−H0底砂差场：仅已计算配对锚点，其余留白',None)]:
                    block=data[k].real
                    if lim is None:lim=float(np.nanmax(abs(block))) if np.isfinite(block).any() else common
                    if lim==0:lim=1.
                    im=axes[row,j].imshow(np.ma.masked_invalid(block),cmap=cmap,vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[190.125,184.875,ticks[-1],ticks[0]])
                    axes[row,j].set(title=NAMES[variant]+' / '+title,xlabel='剖面里程 / m（190→185）',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[row,j],label='(V/m)/(A·m)')
                    if tag=='main':axes[row,j].plot(m['chainage_m'],[s['templates'][variant]['peak_ns'] for s in m['stations']],c='#008a47',lw=1,ls='--',label='局部底砂一次模板（非真值）');axes[row,j].legend(fontsize=8)
            fig.suptitle(f'5m连续高低损耗对照 / {window} / 4.0.1原生FP64 / AGL约8m / 0.25m间距\n新增完成{v["completed_new"]}/46，复用2份；蓝灰色为未计算；无插值、AGC、去背景或幅相拟合')
            fig.savefig(a.out/f'dense_loss_{window}_{tag}_gray.png',dpi=130);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z);h.create_dataset('ids',data=np.array(ids,dtype=h5py.string_dtype('utf-8')))
    result=dict(status='COMPLETE' if v['status']=='COMPLETE' else 'PARTIAL_NOT_FULL_LINE',script_sha256=sha(__file__),manifest_sha256=sha(a.package/'manifest.json'),contract_sha256=sha(a.source/'execution_contract.json'),snapshot_sha256=sha(a.source/'snapshot_manifest.json'),numerical_sha256=sha(a.numerical),native=native,availability=availability,H1_peaks=peaks,anchor_metrics=anchor_metrics,completed_new=v['completed_new'],expected_new=46,reused=2,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],completed_new=v['completed_new'],available_H1={k:sum(r) for k,r in availability.items()},anchor_metrics=anchor_metrics)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
