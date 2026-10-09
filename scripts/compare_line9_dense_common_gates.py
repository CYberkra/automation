"""Shared predeclared target windows: separate actual H0 change from target growth."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha

def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    assert not a.out.exists()
    m=read(a.package/'manifest.json');p=read(a.public/'analysis.json');audit=read(a.public/'independent_audit.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256']==sha(a.public/'analysis.json') and p['numerical_sha256']==sha(a.numerical)
    with h5py.File(a.numerical) as h:
        f=h['frequency_Hz'][:];z=h['response'][:];ids=[x.decode() if isinstance(x,bytes) else x for x in h['ids'][:]]
    np.testing.assert_array_equal(f,20e6+np.arange(501)*300000.)
    t=np.arange(4008)/(4008*300000.);gates={s['chainage_m']:[min(s['templates'][v]['basal_gate_ns'][0] for v in ['low','high']),max(s['templates'][v]['basal_gate_ns'][1] for v in ['low','high'])] for s in m['stations']}
    view=[min(lo for lo,hi in gates.values()),max(hi for lo,hi in gates.values())];keep=(t*1e9>=view[0])&(t*1e9<=view[1]);e=np.exp(2j*np.pi*t[keep,None]*f);a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False;cmap=plt.get_cmap('gray').copy();cmap.set_bad('#d6e5ee')
    allmetrics={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();x=e@(z*w[:,None])/501;lookup={sid:x[:,j] for j,sid in enumerate(ids)};metrics=[];lines={};bg={}
        for v in ['high','low']:
            line=np.full((keep.sum(),21),np.nan+1j*np.nan);background=np.full_like(line,np.nan+1j*np.nan)
            for j,s in enumerate(m['stations']):
                prefix=f'{v}_x{round(s["chainage_m"]*100):05d}_'
                if prefix+'H1' in lookup:line[:,j]=lookup[prefix+'H1']
                if prefix+'H0' in lookup:background[:,j]=lookup[prefix+'H0']
            lines[v]=line;bg[v]=background
        for s in m['stations']:
            xpos=s['chainage_m'];q={};lo,hi=gates[xpos];k=(t[keep]*1e9>=lo)&(t[keep]*1e9<=hi)
            for v in ['high','low']:
                prefix=f'{v}_x{round(xpos*100):05d}_'
                if prefix+'H1' not in lookup or prefix+'H0' not in lookup:break
                h1=lookup[prefix+'H1'][k];h0=lookup[prefix+'H0'][k];d=h1-h0;norm=np.linalg.norm
                q[v]=dict(H0_norm=float(norm(h0)),delta_norm=float(norm(d)),H1_norm=float(norm(h1)),H0_over_delta=float(norm(h0)/norm(d)),H1_vs_delta_correlation=float(abs(np.vdot(h1,d))/(norm(h1)*norm(d))))
            if len(q)==2:metrics.append(dict(chainage_m=xpos,common_gate_ns=[lo,hi],variants=q,low_over_high_delta_norm=q['low']['delta_norm']/q['high']['delta_norm'],low_over_high_H0_norm=q['low']['H0_norm']/q['high']['H0_norm']))
        allmetrics[window]=metrics
        common=max(float(np.nanmax(abs(lines[v].real))) if np.isfinite(lines[v]).any() else 0 for v in ['high','low'])
        bgcommon=max(float(np.nanmax(abs(bg[v].real))) if np.isfinite(bg[v]).any() else 0 for v in ['high','low'])
        fig,axes=plt.subplots(3,2,figsize=(14,10),layout='constrained')
        for j,v in enumerate(['high','low']):
            for row,data,lim,label in [(0,lines[v],common,'H1：高低损耗共同物理灰度'),(1,lines[v],None,'H1：该配方单独放大，色标不同'),(2,bg[v],bgcommon,'H0：仅已计算锚点，同一物理灰度，非噪声真值')]:
                if lim is None:lim=float(np.nanmax(abs(data.real))) if np.isfinite(data).any() else common
                if not lim:lim=1.
                im=axes[row,j].imshow(np.ma.masked_invalid(data.real),cmap=cmap,vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[190.125,184.875,t[keep][-1]*1e9,t[keep][0]*1e9])
                axes[row,j].set(title=('高损耗' if v=='high' else '低损耗诊断')+' / '+label,xlabel='剖面里程 / m',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axes[row,j],label='(V/m)/(A·m)')
                axes[row,j].plot(m['chainage_m'],[s['templates'][v]['peak_ns'] for s in m['stations']],c='#008a47',ls='--',lw=1,label='局部底砂一次模板（非真值）');axes[row,j].legend(fontsize=8)
        fig.suptitle(f'底砂窗细看 / {window} / 新增完成{p["completed_new"]}/46 / 原生FP64 / 约8m航高\n窗口由求解前两种材料模板窗并集确定；蓝灰色未计算，禁止插值、AGC、幅相或延时拟合')
        fig.savefig(a.out/f'dense_loss_{window}_basal_zoom.png',dpi=130);plt.close(fig)
    result=dict(status=p['status'],script_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),independent_audit_sha256=sha(a.public/'independent_audit.json'),numerical_sha256=sha(a.numerical),view_gate_ns=view,metrics=allmetrics,limits='Same predeclared union gate within each high/low pair. Absolute H0 change distinguished from delta growth; no claim H0 is noise,phase-matched clean,continuous H0 or unique mechanism.')
    (a.out/'common_gate_comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],metrics=allmetrics)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','public','numerical','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
