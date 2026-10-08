"""Render audited native Ricker snapshots on the solver host; no solver mutations."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.patches import Rectangle
    from PIL import Image
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    read=lambda p:json.loads(p.read_text('utf-8'))
    cp=a.execution/'execution_contract.json';vp=a.execution/'completed_verification.json'
    c=read(cp);v=read(vp);m=c['study_manifest']
    assert v['completed'] and v['contract_sha256']==sha(cp) and len(v['groups'])==2
    assert not a.out.exists(),'Fresh render directory required'
    rows={g['id']:g for g in v['groups']};raw={};files={}
    for name,row in rows.items():
        p=Path(row['native_path']);assert sha(p)==row['native_sha256']
        with h5py.File(p) as h:
            raw[name]=h['rxs/rx1/Ez'][:];dt=float(h.attrs['dt'])
            assert raw[name].dtype==np.float64
        files[name]={r['iteration']:r for r in row['snapshots']}
        assert len(files[name])==m['snapshot_count']==799
    with h5py.File(a.reference) as h:reference=h['rxs/rx1/Ez'][:]
    assert np.array_equal(reference,raw['H0']),'Observers changed native receiver sequence'
    selected=[j for k,j in enumerate(m['snapshot_iterations']) if k%3==0 and j*dt<=500e-9]
    arrays=[];digests=[]
    for j in selected:
        frame=[]
        for name in ['H0','far_removed']:
            row=files[name][j];p=Path(row['file']);assert sha(p)==row['sha256']
            with h5py.File(p) as h:
                assert h.attrs['iteration']==j and abs(h.attrs['time']-j*dt)<1e-20
                x=h['Ez'][:,:,0];assert x.dtype==np.float64 and x.shape==(450,270)
            frame.append(x.T)
            digests.append(dict(group=name,iteration=j,sha256=row['sha256']))
        arrays.append([frame[0],frame[1],frame[0]-frame[1]])
    scales=[max(float(abs(frame[k]).max()) for frame in arrays for k in [0,1])]
    scales.append(max(float(abs(frame[2]).max()) for frame in arrays))
    assert min(scales)>0
    a.out.mkdir(parents=True)
    # Exact H0 material boundaries in observer coordinates. These are guides only.
    gp=next(g for g in c['groups'] if g['id']=='H0')
    with h5py.File(Path(c['package'])/gp['geometry']) as h:geo=h['data'][1600:3400:4,400:1480:4,0].T
    assert geo.shape==(270,450)
    xx=40+np.arange(450)*.1;yy=10+np.arange(270)*.1
    fig=plt.figure(figsize=(13,6),layout='constrained')
    grid=fig.add_gridspec(2,3,height_ratios=[3,1]);axes=[fig.add_subplot(grid[0,k]) for k in range(3)]
    names=['原H0：保留上覆夹层，已去底砂','只将x45–60m夹层换泥岩','材料响应差场：原H0−替换场']
    images=[]
    for k,ax in enumerate(axes):
        limit=scales[0 if k<2 else 1]
        im=ax.imshow(arrays[0][k],origin='lower',extent=[40,85,10,37],aspect='equal',cmap='gray',norm=SymLogNorm(limit*1e-5,vmin=-limit,vmax=limit,base=10))
        images.append(im)
        ax.contour(xx,yy,geo,levels=[.5,1.5,2.5],colors=['lime'],linewidths=.45,alpha=.7)
        ax.add_patch(Rectangle((45,10),15,27,fill=False,edgecolor='magenta',lw=.7))
        ax.scatter([78.6,79.9],[35,35],c=['orange','red'],s=14,zorder=4)
        ax.set(xlabel='模型local x / m',ylabel='模型y / m',title=names[k],xlim=(40,85),ylim=(10,37))
        fig.colorbar(im,ax=ax,shrink=.8,label='原生Ez / V/m；固定对称对数显示')
    graph=fig.add_subplot(grid[1,:]);tn=np.arange(len(raw['H0']))*dt*1e9
    for name,label in [('H0','原H0'),('far_removed','替换夹层')]:graph.plot(tn,raw[name],lw=.7,label=label)
    graph.plot(tn,raw['H0']-raw['far_removed'],lw=.7,label='原H0−替换场')
    graph.set(xlim=(180,450),xlabel='原生时间 / ns（含Ricker源延迟14.142ns）',ylabel='接收Ez / V/m')
    gate=(tn>=180)&(tn<=450)
    limit=max(float(abs(raw['H0'][gate]).max()),float(abs(raw['far_removed'][gate]).max()))
    graph.set_ylim(-1.05*limit,1.05*limit);graph.legend(loc='upper right',ncol=3)
    cursor=graph.axvline(0,color='red',lw=1)
    title=fig.suptitle('原生Ricker波场，非501频点SFCW；绿线为原H0层界，紫框为材料替换范围')
    frames=[];static=[]
    targets={min(selected,key=lambda j:abs(j*dt*1e9-t)) for t in [100,150,200,250,300,350,400,450]}
    for j,values in zip(selected,arrays):
        for im,value in zip(images,values):im.set_data(value)
        cursor.set_xdata([j*dt*1e9]*2)
        title.set_text(f'原生Ricker波场 t={j*dt*1e9:.2f}ns，非SFCW；绿线原层界，紫框材料替换ROI')
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()).convert('P',palette=Image.Palette.ADAPTIVE,colors=128))
        if j in targets:
            p=a.out/f'wavefield_{j:05d}.png';fig.savefig(p,dpi=120);static.append(dict(file=p.name,iteration=j,time_ns=j*dt*1e9,sha256=sha(p)))
    movie=a.out/'interbed_wavefield.gif'
    frames[0].save(movie,save_all=True,append_images=frames[1:],duration=90,loop=0,optimize=True)
    plt.close(fig)
    report=dict(status='AUDITED_NATIVE_RICKER_MOVIE_NOT_SFCW_OR_UNIQUE_PATH_CLASSIFICATION',script_sha256=sha(__file__),contract_sha256=sha(cp),verification_sha256=sha(vp),reference_native_sha256=sha(a.reference),
        H0_receiver_bitwise_equal_to_no_observers=True,total_snapshots_per_case=799,rendered_frames=len(frames),rendered_snapshot_hashes=digests,
        sampling_time_step_ns=51*dt*1e9,display_end_ns=selected[-1]*dt*1e9,original_and_control_shared_abs_limit_V_per_m=scales[0],difference_abs_limit_V_per_m=scales[1],symlog_linear_threshold_fraction=1e-5,
        gif_sha256=sha(movie),gif_bytes=movie.stat().st_size,static_frames=static,
        limits='ROI changes edges and all interactions. Spatial observer sampling .1m; full FDTD remains .025m. Difference has its own fixed color scale. No claim of unique scattering point/order or real antenna/field equivalence.')
    (a.out/'render_receipt.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['status','rendered_frames','gif_bytes','H0_receiver_bitwise_equal_to_no_observers']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['execution','reference','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
