"""Coordinate-labelled raw key details; shared scale excludes the source instant."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from hs_capsule_identity import sha256 as sha


def render(root,package,out):
    c=json.loads((root/'execution_contract.json').read_text('utf-8'));receipts=[]
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    with h5py.File(package/'geometries/H1.h5') as h:geo=h['data'][:,:,0]
    for x in [170.,190.]:
        pair={r:next(g for g in c['groups'] if g['role']==r and g['kind']=='local' and g['profile_x_m']==x) for r in ['H1','H0']}
        h1={p.name:p for p in (root/'keys'/pair['H1']['id']).glob('*.h5')};rows=[];vmax=0.
        for name,p in sorted(h1.items()):
            q=root/'keys'/pair['H0']['id']/name
            with h5py.File(p) as a,h5py.File(q) as b:
                assert a.attrs['time']==b.attrs['time'];e,f=a['Ez'][:,:,0],b['Ez'][:,:,0]
                rows.append((float(a.attrs['time'])*1e9,e,f,p,q));vmax=max(vmax,float(abs(e).max()),float(abs(f).max()),float(abs(e-f).max()))
        g=pair['H1'];roi=g['roi'];left=roi['bounds'][0];xx=left+(np.arange(800)+.5)*.025;columns=np.floor(xx/.025).astype(int)
        idx=np.arange(1700)[None,:]
        lines=[np.max(np.where(geo[columns]!=0,idx+1,0),axis=1)*.025,np.min(np.where(geo[columns]==1,idx,1700),axis=1)*.025,np.max(np.where(geo[columns]==3,idx+1,0),axis=1)*.025]
        for target in [240,360]:
            tt,e,f,p,q=min(rows,key=lambda r:abs(r[0]-target));fig,axes=plt.subplots(1,3,figsize=(18,10),constrained_layout=True)
            for ax,z,title in zip(axes,[e,f,e-f],['H1 原始总场','H0 去底砂反差','H1−H0 场差（含交互）']):
                im=ax.imshow(z.T,origin='lower',extent=[left+20,left+40,12,40],aspect='equal',interpolation='none',cmap='RdBu_r',norm=SymLogNorm(vmax*1e-6,vmin=-vmax,vmax=vmax))
                for line,color in zip(lines,['green','darkorange','black']):ax.plot(xx+20,line,color=color,linewidth=1)
                for pos,tag in [(g['tx_m'],'Tx'),(g['rx_m'],'Rx')]:ax.plot(pos[0]+20,pos[1],'o',color='purple',markersize=3);ax.annotate(tag,(pos[0]+20,pos[1]),xytext=(3,3),textcoords='offset points',fontsize=9)
                ax.set(title=title,xlabel='剖面x / m',ylabel='模型内部y / m')
            fig.suptitle(f'{x:g}m · {tt:.3f}ns 原始宽带冲激场（非SFCW）\n2.5cm原生数组；两情形和场差共用此站8个关键时刻的绝对色标；绿地表/橙覆底/黑底砂参考',fontsize=13)
            fig.colorbar(im,ax=axes,shrink=.8,label='Ez / V/m；固定对数阈值，相同色标，未改变数值')
            dest=out/f'x{int(x)}_local_detail_{target}ns.png';fig.savefig(dest,dpi=180);plt.close(fig)
            receipts.append(dict(file=dest.name,sha256=sha(dest),native_input_sha256=[sha(p),sha(q)],time_ns=tt,
                raw_shape=[800,1120],shared_vmax_V_m=vmax,scale_reference='same eight retrieved late key instants, not source instant',display_only=True))
    (out/'key_detail_render_receipts.json').write_text(json.dumps(dict(script_sha256=sha(__file__),files=receipts),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['root','package','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();render(a.root,a.package,a.out)
