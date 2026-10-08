"""Phase-preserving four-site sparse spatial check; no filled missing traces."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window
from diagnose_line9_layer_kinematics import primary_paths
from review_line9_result_packages import geometry_at,indices

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/'configs/research/line9_v401_spatial_diagnostic_v1.json'


def correlation(a,b):
    return float(abs(np.vdot(a,b))/(np.linalg.norm(a)*np.linalg.norm(b)))


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(a.package/'manifest.json');c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    group={g['id']:g for g in m['groups']};rec={g['id']:g for g in v['groups']}
    anchor=m['reused'][1]['parent_group'];material=a.package/group['full2d_r0141_H1']['material']
    db=read(material)['materials'];n95=indices(db,95e6)
    with h5py.File(a.package/group['full2d_r0141_H1']['geometry']) as h:data=h['data'][:,:,0];dl=h.attrs['dx_dy_dz']
    positions=[198.6,190.,180.,146.4];z=[];templates=[];rows=[];native=[]
    for position in positions:
        if position==198.6:
            g=anchor;sid='anchor';paths={r['id'].rsplit('_',1)[1]:(a.package/r['path'],r['native_sha256']) for r in m['reused']}
            basal_gate=m['basal_gate_ns'];geometry=geometry_at(data,dl,np.array(g['tx_m']),np.array(g['rx_m']),n95)
            center=(basal_gate[0]+basal_gate[1])/2;wide_gate=[center-60,center+60]
        else:
            s=next(s for s in m['station_design']['stations'] if s['chainage_m']==position);sid=s['id'];g=group[sid+'_H1'];geometry=s['geometry']
            paths={role:(a.source/(sid+'_'+role+'.h5'),rec[sid+'_'+role]['native_sha256']) for role in ['H0','H1']}
            basal_gate=g['basal_gate_ns'];wide_gate=g['wide_gate_ns']
        pair=[]
        for role in ['H0','H1']:
            p,digest=paths[role];assert sha(p)==digest
            with h5py.File(p) as h:assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64
            q,err=response(p,.025);pair.append(q);native.append(dict(id=sid+'_'+role,native_sha256=digest,independent_DFT_relative_L2=err))
        z.extend([*pair,pair[1]-pair[0]])
        template=primary_paths(geometry,db)[geometry['boundaries'].index(geometry['basal_sand'])];templates.append(template)
        rows.append(dict(id=sid,chainage_m=position,tx_m=g['tx_m'],rx_m=g['rx_m'],geometry=geometry,basal_gate_ns=basal_gate,wide_gate_ns=wide_gate))
    z=np.column_stack(z);templates=np.column_stack(templates);profiles={};modelprofiles={};metrics={}
    for window in ['hann','blackman']:
        w=spectral_window(window,501);x,t=inverse(z,w);q,_=inverse(templates,w);profiles[window]=x;modelprofiles[window]=q;metrics[window]=[]
        common=(t*1e9>=180)&(t*1e9<=450)
        for j,row in enumerate(rows):
            model_peak=float(t[np.argmax(abs(q[:,j]))]*1e9);gates={}
            for label,bounds in [('basal',row['basal_gate_ns']),('wide',row['wide_gate_ns'])]:
                k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);h0,h1,delta=x[k,3*j:3*j+3].T
                imax=int(np.argmax(abs(delta)));hmax=int(np.argmax(abs(h1)));ticks=t[k]*1e9
                gates[label]=dict(H0_over_delta_L2=float(np.linalg.norm(h0)/np.linalg.norm(delta)),H1_vs_delta_complex_correlation=correlation(h1,delta),
                    delta_vs_local_primary_complex_correlation=correlation(delta,q[k,j]),H1_vs_local_primary_complex_correlation=correlation(h1,q[k,j]),
                    delta_peak_ns=float(ticks[imax]),H1_peak_ns=float(ticks[hmax]),peak_on_gate_edge=bool(imax in [0,len(ticks)-1]),
                    actual_minus_template_peak_ns=float(ticks[imax]-model_peak),H1_minus_template_peak_ns=float(ticks[hmax]-model_peak),
                    delta_L2=float(np.linalg.norm(delta)),H0_L2=float(np.linalg.norm(h0)))
            metrics[window].append(dict(chainage_m=row['chainage_m'],template_peak_ns=model_peak,gates=gates,
                common_180_450_delta_peak_ns=float(t[common][np.argmax(abs(x[common,3*j+2]))]*1e9),
                common_180_450_H1_peak_ns=float(t[common][np.argmax(abs(x[common,3*j+1]))]*1e9),
                common_180_450_H0_peak_ns=float(t[common][np.argmax(abs(x[common,3*j]))]*1e9)))
        anchor_result=metrics[window][0]
        for row in metrics[window]:
            row['delta_anchor_relative_time_change_minus_model_time_change_ns']=float((row['gates']['basal']['delta_peak_ns']-anchor_result['gates']['basal']['delta_peak_ns'])-(row['template_peak_ns']-anchor_result['template_peak_ns']))
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.patches import Patch,Rectangle
    from matplotlib.colors import ListedColormap,BoundaryNorm
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(figsize=(16,6),layout='constrained')
    colors=['#e8f4fa','#e0b992','#9c8c87','#e6cf70'];cmap=ListedColormap(colors)
    # Display every fourth voxel; all numerical checks use the full original map.
    ax.imshow(data[::4,::4].T,origin='lower',extent=[20,230,417.5,460],aspect='auto',interpolation='nearest',cmap=cmap,norm=BoundaryNorm(np.arange(5)-.5,4))
    for rect in [(20,417.5,2,42.5),(228,417.5,2,42.5),(20,417.5,210,2),(20,458,210,2)]:
        ax.add_patch(Rectangle(rect[:2],rect[2],rect[3],facecolor='gray',alpha=.25,hatch='//',edgecolor='gray'))
    for row in rows:
        tx=np.array(row['tx_m']);rx=np.array(row['rx_m']);mid=(tx+rx)/2
        ax.scatter([tx[0]+20,rx[0]+20],[tx[1]+417.5,rx[1]+417.5],c=['#c0392b','#2366a8'],s=25,marker='v')
        ax.text(row['chainage_m'],mid[1]+416.9,f'{row["chainage_m"]:g}m',ha='center',va='top',fontsize=8,rotation=35)
    ax.annotate('测线方向：220→25m',xy=(30,452),xytext=(85,452),arrowprops=dict(arrowstyle='->'),ha='center')
    handles=[Patch(facecolor=c,label=s) for c,s in zip(colors,['空气','粉质粘土','泥岩','底连通砂岩'])]
    handles.append(Patch(facecolor='gray',alpha=.25,hatch='//',label='80格/2m PML'))
    ax.legend(handles=handles,loc='lower right',ncol=5,fontsize=8)
    ax.set(xlabel='剖面里程 / m',ylabel='模型参考高程 / m',title='原v5全域H1几何与四个实际站位（红Tx、蓝Rx）；AGL约8m\nH0把所有底连通砂岩换为泥岩；图仅每4格显示，计算与审核均用原2.5cm网格；纵横比例非1:1')
    fig.savefig(a.out/'spatial_pairs_geometry_stations.png',dpi=140);plt.close(fig)
    labels=[f'{p:g}m'+('（复用）' if j==0 else '') for j,p in enumerate(positions)]
    keep=(t*1e9>=180)&(t*1e9<=450);ticks=t[keep]*1e9
    grid=146.4+np.arange(262)*.2;columns=np.rint((np.array(positions)-146.4)/.2).astype(int)
    np.testing.assert_allclose(grid[columns],positions,atol=1e-10,rtol=0)
    for window,x in profiles.items():
        lim=float(abs(x[keep].real).max());fig,axes=plt.subplots(1,3,figsize=(18,7),layout='constrained')
        for ax,offset,title in zip(axes,[1,0,2],['H1原始总场','H0底连通砂岩换泥岩','H1−H0底砂替换差场']):
            pic=ax.imshow(x[keep,offset::3].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[-.5,3.5,ticks[-1],ticks[0]])
            ax.set_xticks(range(4),labels,rotation=15);ax.set(title=title,ylabel='SFCW时间 / ns')
            ax.scatter(range(4),[r['template_peak_ns'] for r in metrics[window]],marker='_',s=150,c='green',label='局部一次反射预测')
            ax.legend(fontsize=8);fig.colorbar(pic,ax=ax,label='共同绝对灰度 (V/m)/(A·m)')
        fig.suptitle('原非平v5 / AGL约8m / 4.0.1 FP64 / '+window+' / 四站配置列，间隔不按距离\n低跨度0.3、DC0.3mS/m机制对照；精确501点，无AGC/拟合；正式材料未换')
        fig.savefig(a.out/f'spatial_pairs_{window}_station_gray.png',dpi=140);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(18,7),layout='constrained');cmap=plt.get_cmap('gray').copy();cmap.set_bad('#cbd7e4')
        for ax,offset,title in zip(axes,[1,0,2],['H1原始总场','H0剩余地质场','底砂替换差场']):
            sparse=np.ma.masked_invalid(np.full((keep.sum(),len(grid)),np.nan));sparse[:,columns]=x[keep,offset::3].real
            pic=ax.imshow(sparse,cmap=cmap,vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[146.3,198.7,ticks[-1],ticks[0]])
            ax.scatter(positions,[r['template_peak_ns'] for r in metrics[window]],marker='_',s=80,c='green')
            ax.set_xlim(198.7,146.3);ax.set(title=title,xlabel='实际剖面里程 / m',ylabel='SFCW时间 / ns')
            ax.set_xticks(positions,[f'{p:g}' for p in positions],rotation=20)
            ax.legend(handles=[Patch(facecolor='#cbd7e4',label='该0.2m位置未计算')],fontsize=8);fig.colorbar(pic,ax=ax,label='共同 (V/m)/(A·m)')
        fig.suptitle('四站稀疏物理B-scan / '+window+' / 198.6m复用、其余三站新算\n仅显示146.4–198.6m范围，全部缺道保留空白；不插值、不连成已算测线')
        fig.savefig(a.out/f'spatial_pairs_{window}_sparse_bscan.png',dpi=140);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
        for j,ax in enumerate(axes.flat):
            for offset,title in [(1,'H1总场'),(0,'H0剩余场'),(2,'H1−H0差场')]:ax.plot(ticks,abs(x[keep,3*j+offset]),lw=1,label=title)
            ax.axvline(metrics[window][j]['template_peak_ns'],color='green',ls=':',label='局部一次反射预测')
            for edge in rows[j]['basal_gate_ns']:ax.axvline(edge,color='gray',ls='--',lw=.7)
            ax.set(title=labels[j]+'；同站绝对包络',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');ax.legend(fontsize=8)
        fig.suptitle('四站包络 / '+window+' / 无逐道归一化；子图纵轴刻度各自标明\n虚线为求解前固定窗，点线为模型近似；不能当盲测检出或完整地层恢复')
        fig.savefig(a.out/f'spatial_pairs_{window}_envelopes.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z);h.create_dataset('templates',data=templates)
        h.create_dataset('chainage_m',data=positions);h.create_dataset('sparse_chainage_m',data=grid)
        mask=np.zeros(len(grid),bool);mask[columns]=True;h.create_dataset('completed_mask',data=mask)
    result=dict(status='COMPLETED_SIX_NEW_TWO_REUSED_FOUR_SPARSE_STATIONS_NOT_FULL_LINE',script_sha256=sha(__file__),protocol_sha256=sha(PROTOCOL),
        contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),numerical_sha256=sha(a.numerical),
        native=native,stations=rows,metrics=metrics,new_solves=6,reused_native=2,AGC=False,fit=False,production_materials_changed=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(metrics))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
