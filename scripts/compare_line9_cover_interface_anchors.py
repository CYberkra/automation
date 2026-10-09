"""Two measured simulation anchors on a sparse grid; no interpolated traces."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import inverse,FREQ


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    data=[];references=[]
    for pub,num,station in [(a.pub18750,a.num18750,187.5),(a.pub190,a.num190,190.)]:
        r=read(pub/'analysis.json');audit=read(pub/'independent_audit.json');d=read(pub/'interface_dependent_diagnostic/diagnostic.json')
        assert audit['analysis_sha256']==d['analysis_sha256']==sha(pub/'analysis.json') and d['independent_audit_sha256']==sha(pub/'independent_audit.json')
        assert r['numerical_sha256']==d['numerical_sha256']==sha(num)
        with h5py.File(num) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);z=h['response'][:]
        assert z.shape==(501,3);data.append(z)
        references.append(dict(chainage_m=station,analysis_sha256=sha(pub/'analysis.json'),audit_sha256=sha(pub/'independent_audit.json'),difference_diagnostic_sha256=sha(pub/'interface_dependent_diagnostic/diagnostic.json'),numerical_sha256=sha(num),contract_sha256=r['contract_sha256']))
    data=np.stack(data,axis=1);stations=np.arange(187.5,190.0001,.25);rows={}
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    cmap=plt.get_cmap('gray').copy();cmap.set_bad('#d8e1e8')
    labels=['原高损耗H0','覆盖层底介电对比移除','覆盖层底下移0.5m']
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();x,t=inverse(data.reshape(501,6),w);x=x.reshape(4008,2,3)
        k=(t*1e9>=300)&(t*1e9<=450);tt=t[k]*1e9
        mat=np.full((k.sum(),len(stations),3),np.nan);mat[:,[0,-1],:]=x[k].real
        difference=x[:,:, [0,2]]-x[:,:,1,None];dmat=np.full((k.sum(),len(stations),2),np.nan);dmat[:,[0,-1],:]=difference[k].real
        lim=float(np.max(abs(x[k].real)));dlim=float(np.max(abs(difference[k].real)))
        fig,axes=plt.subplots(2,3,figsize=(17,9),layout='constrained')
        for j in range(3):
            im=axes[0,j].imshow(mat[:,:,j],aspect='auto',extent=[stations[0]-.125,stations[-1]+.125,tt[-1],tt[0]],interpolation='nearest',cmap=cmap,vmin=-lim,vmax=lim)
            axes[0,j].set(title=labels[j],xlabel='测线里程 / m',ylabel='SFCW时间 / ns');axes[0,j].set_xticks([187.5,190])
            for s in [187.5,190]:axes[0,j].axvline(s,color='tab:blue',ls='--',lw=.5)
        fig.colorbar(im,ax=list(axes[0]),label='原场共用尺度 / (V/m)/(A·m)',shrink=.8)
        for j,label in enumerate(['原H0−本站去对比参考','下移H0−本站同一参考']):
            im=axes[1,j].imshow(dmat[:,:,j],aspect='auto',extent=[stations[0]-.125,stations[-1]+.125,tt[-1],tt[0]],interpolation='nearest',cmap=cmap,vmin=-dlim,vmax=dlim)
            axes[1,j].set(title=label,xlabel='测线里程 / m',ylabel='SFCW时间 / ns');axes[1,j].set_xticks([187.5,190])
        fig.colorbar(im,ax=list(axes[1,:2]),label='界面相关差分共同尺度 / (V/m)/(A·m)',shrink=.8)
        axes[1,2].axis('off');axes[1,2].text(.02,.92,'只有两个站位，蓝灰色为未计算\n不插值、不按站位归一\n\nH0无底部连通砂岩\n差分含传播与交互变化\n不能直接称纯多次波/clean\n\n虚线为实际求解站位锚点\n两个站位不认证整线或三维',va='top',fontsize=13)
        fig.suptitle(f'187.5m与190m覆盖层控制复核／4.0.1原生FP64／精确501点／{window}\n每站原场三个配置共享同一灰度；差分单独声明尺度；原地表/采集保持')
        fig.savefig(a.out/f'two_anchor_{window}_sparse.png',dpi=140);plt.close(fig)
        rows[window]=dict(total_common_limit=lim,difference_common_limit=dlim,computed_chainage_m=[187.5,190],missing_chainage_m=stations[1:-1].tolist())
    result=dict(status='PASS_TWO_ANCHOR_SPARSE_COMMON_SCALE_NO_INTERPOLATION',script_sha256=sha(__file__),sources=references,display=rows,solver_runs=0,limits='Two simulated anchor positions, not complete line/3D/field validation. Common-background differences include interactions.')
    (a.out/'comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],computed_stations=2,missing_columns=9)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['pub18750','num18750','pub190','num190','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
