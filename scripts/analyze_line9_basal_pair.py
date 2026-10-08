"""Phase-preserving paired single-station diagnostics and fixed-scale wavefield movies."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def save(p,v):
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def plotting():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    return plt


def traces(a):
    from gprMax.toolboxes.SFCW import processing as sf
    from diagnose_line9_layer_kinematics import primary_paths
    from audit_line9_postprocessing import inverse,weights,rel
    from trace_line9_time_origin import correlation
    m=json.loads(a.manifest.read_text('utf-8'))
    verification=json.loads(a.verification.read_text('utf-8'))
    assert verification['completed'] and len(verification['groups'])==2
    paths=[a.h1,a.h0,a.fp32]
    for p,g in zip(paths,verification['groups']):assert sha(p)==g['raw_sha256']
    assert sha(a.fp32)==m['source_native_sha256']
    with h5py.File(paths[0]) as h:
        dt=float(h.attrs['dt']);shape=h.attrs['nx_ny_nz'].copy()
        assert h['rxs/rx1/Ez'].dtype==np.float64
    responses=[];dft_errors=[]
    f=20e6+np.arange(501)*300000.
    for p in paths:
        with h5py.File(p) as h:
            assert float(h.attrs['dt'])==dt
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],shape)
        src=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez')
        assert src.spatial_scale==.025
        product=sf.direct_frequency_response(src,rx,f,tail_taper_fraction=0)
        assert product.source_valid.all()
        response=np.asarray(product.response).reshape(501)/src.spatial_scale
        k=np.array([0,83,167,250,333,417,500])
        expected=(rx.dt*(np.exp(-2j*np.pi*f[k,None]*rx.times)@rx.samples))/(src.dt*(np.exp(-2j*np.pi*f[k,None]*src.times)@src.samples))/src.spatial_scale
        e=rel(response[k],expected);assert e<1e-9;dft_errors.append(e);responses.append(response)
    responses=np.column_stack(responses)
    profiles={};inverse_errors={};direct_errors={}
    for window in ['hann','blackman']:
        w=weights(window,501)
        profile,t=inverse(responses,f,w)
        official=sf.reconstruct_time_response(replace(product,response=responses),window=window,zero_pad_factor=8)
        e=rel(profile,official.complex_bandpass);assert e<1e-9
        take=np.arange(0,len(t),97)
        expected=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*responses)/len(f)
        direct_errors[window]=rel(profile[take],expected);assert direct_errors[window]<1e-9
        profiles[window]=profile;inverse_errors[window]=e
    materials=json.loads(a.materials.read_text('utf-8'))['materials']
    assert sha(a.materials)==m['material_sha256']
    g=m['geometry_diagnostic'];primaries=primary_paths(g,materials)
    columns={name:g['boundaries'].index(g[name]) for name in ['cover_base','first_sand','basal_sand']}
    plt=plotting();a.out.mkdir(parents=True,exist_ok=False)
    report=dict(calls_solver=False,calls_training=False,script_sha256=sha(__file__),verification_sha256=sha(a.verification),
        native_sha256=[sha(p) for p in paths],dt_s=dt,frequency_MHz=[20,170],step_MHz=.3,tones=501,
        spatial_normalisation_m=.025,tail_taper=False,AGC=False,per_trace_normalisation=False,
        station_chainage_m=m['chainage_m'],midpoint_AGL_m=m['midpoint_agl_m'],
        independent_DFT_relative_L2=dft_errors,official_inverse_relative_L2=inverse_errors,independent_inverse_relative_L2=direct_errors,
        limitations='One informed selected station; lower-material contrast includes propagation, multiples and changed PML continuation. H1-H0 is not pure single reflection or field clean truth. FP64-FP32 includes build/rounding/snapshot configuration differences; no universal precision floor.',windows={})
    for window,z in profiles.items():
        delta=z[:,0]-z[:,1];precision=z[:,0]-z[:,2]
        record={}
        for name,index in columns.items():
            model,_=inverse(primaries[index],f,weights(window,501))
            center=float(t[np.argmax(abs(model))]*1e9)
            gate=abs(t*1e9-center)<=12
            def ratio(x,y):return float(np.linalg.norm(x[gate])/np.linalg.norm(y[gate]))
            coefficient=np.vdot(model[gate],delta[gate])/np.vdot(model[gate],model[gate])
            record[name]=dict(primary_peak_ns=center,gate_ns=[center-12,center+12],
                difference_over_H1_L2=ratio(delta,z[:,0]),precision_difference_over_H1_L2=ratio(precision,z[:,0]),
                H0_over_material_difference_L2=ratio(z[:,1],delta),
                H0_difference_inner_product_phase_deg=float(np.angle(np.vdot(z[gate,1],delta[gate]),deg=True)),
                precision_difference_over_material_difference_L2=ratio(precision,delta),
                H1_H0_complex_correlation=correlation(z[gate,0],z[gate,1]),
                H1_primary_complex_correlation=correlation(z[gate,0],model[gate]),
                difference_primary_complex_correlation=correlation(delta[gate],model[gate]),
                difference_primary_fit_relative_L2=rel(coefficient*model[gate],delta[gate]),
                difference_peak_ns=float(t[np.argmax(np.where(gate,abs(delta),0))]*1e9))
        report['windows'][window]=record
    z=profiles['hann'];signals=np.column_stack([z[:,:2],z[:,0]-z[:,1],z[:,0]-z[:,2]])
    labels=['H1：原地质 / FP64','H0：底部砂岩改泥岩 / FP64','H1−H0：先作复数差分','FP64−原FP32：数值配置差异']
    fig,axs=plt.subplots(3,4,figsize=(15,11),layout='constrained')
    for row,(lo,hi) in enumerate([(0,800),(200,410),(320,380)]):
        gate=(t*1e9>=lo)&(t*1e9<=hi)
        scale=float(np.max(abs(signals[gate,:2].real)))
        for j in range(4):
            ax=axs[row,j];im=ax.imshow(signals[gate,j].real[:,None],extent=[-.5,.5,hi,lo],aspect='auto',cmap='gray',vmin=-scale,vmax=scale)
            for name,item in report['windows']['hann'].items():
                if lo<item['primary_peak_ns']<hi:ax.axhline(item['primary_peak_ns'],color={'cover_base':'#e69100','first_sand':'#1565c0','basal_sand':'#b80058'}[name],lw=1)
            ax.set(title=labels[j],ylabel='双程时间 / ns',xticks=[0],xticklabels=['第8道'])
        fig.colorbar(im,ax=axs[row],label='带符号实部；本行四图共用色标 / (V/m)/(A·m)')
    fig.suptitle('X179.25m，AGL10.275m：单站对照（非完整B-scan）\n三行分别为全时窗/地下段/底砂段，各行四图共同色标；橙/蓝/红线：覆盖层底/首砂层/底砂的模型预测')
    fig.savefig(a.out/'single_station_grayscale.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(13,10),layout='constrained')
    for ax,window,(lo,hi) in zip(axs,['hann','hann','blackman'],[(200,410),(320,380),(320,380)]):
        z=profiles[window]
        values=[z[:,0],z[:,1],z[:,0]-z[:,1],z[:,0]-z[:,2]]
        for v,label in zip(values,labels):ax.plot(t*1e9,abs(v),label=label,lw=1)
        for name,item in report['windows'][window].items():ax.axvspan(*item['gate_ns'],alpha=.12,color={'cover_base':'orange','first_sand':'blue','basal_sand':'red'}[name])
        gate=(t*1e9>=lo)&(t*1e9<=hi)
        ax.set(xlim=(lo,hi),ylim=(0,1.1*max(float(abs(v[gate]).max()) for v in values)),xlabel='双程时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=window+'窗；同一绝对幅度；着色区为模型预测±12ns，未经移动搜索');ax.legend(fontsize=8)
    fig.savefig(a.out/'paired_complex_envelopes.png',dpi=150);plt.close(fig)
    if a.cache:
        with np.load(a.cache) as cache:
            response=cache['response'];chain=cache['chainage_m']
        zz,tt=inverse(response,f,weights('hann',501));gg=(tt*1e9>=150)&(tt*1e9<=420)
        vmax=float(np.percentile(abs(zz[gg].real),99.5))
        fig,ax=plt.subplots(figsize=(13,6),layout='constrained')
        im=ax.imshow(zz[gg].real,extent=[chain[0],chain[-1],420,150],cmap='gray',aspect='auto',vmin=-vmax,vmax=vmax)
        ax.axvline(m['chainage_m'],color='red',label='此次对照：X179.25m / 第8道');ax.legend()
        ax.set(title='原142道高航高B-scan：Hann窗，统一幅度，无去背景/AGC；红线为此次单站验证',xlabel='剖面横坐标 / m（采集方向从大到小）',ylabel='双程时间 / ns')
        fig.colorbar(im,ax=ax,label='带符号实部；全图共同99.5百分位截色，仅显示裁切')
        fig.savefig(a.out/'existing_bscan_station_marked.png',dpi=150);plt.close(fig)
    if a.numerical_out:
        if a.numerical_out.exists():raise ValueError('Fresh private numerical output required')
        with h5py.File(a.numerical_out,'x') as h:
            h.create_dataset('frequency_Hz',data=f);h.create_dataset('response',data=responses)
            h.create_dataset('time_s',data=t)
            for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
            h.attrs['columns']='H1 float64; H0 float64; original float32'
            h.attrs['units']='(V/m)/(A*m); 2D field proxy, not port S21'
            h.attrs['source_spatial_scale_m']=.025;h.attrs['script_sha256']=sha(__file__)
        report['private_numerical_sha256']=sha(a.numerical_out)
    save(a.out/'analysis.json',report)
    print(json.dumps(report['windows'],ensure_ascii=False))


def wavefield(a):
    from matplotlib.colors import SymLogNorm
    from matplotlib.patches import Rectangle
    from PIL import Image
    plt=plotting();a.out.mkdir(parents=True,exist_ok=False)
    verification=json.loads(a.verification.read_text('utf-8'))
    assert verification['completed']
    sets=[sorted(a.h1snap.glob('snap*.h5')),sorted(a.h0snap.glob('snap*.h5'))]
    assert len(sets[0])==len(sets[1])==256
    for files,row in zip(sets,verification['groups']):
        assert len(files)==row['snapshot_count']
        for file,item in zip(files,row['snapshots']):assert sha(file)==item['sha256']
    with h5py.File(a.geometry) as h:geo=h['data'][::8,::8,0].T
    frames=[];times=[];maxima=[];differences=[]
    for p,q in zip(*sets):
        with h5py.File(p) as h1,h5py.File(q) as h0:
            for key in ['iteration','time','origin','dx_dy_dz']:np.testing.assert_array_equal(h1.attrs[key],h0.attrs[key])
            x=h1['Ez'][::2,::2,0].T;y=h0['Ez'][::2,::2,0].T
            assert x.dtype==y.dtype==np.float64
            times.append(float(h1.attrs['time'])*1e9);frames.append((x,y,x-y))
            maxima.append(float(max(abs(x).max(),abs(y).max())));differences.append(float(abs(x-y).max()))
    limit=max(maxima);dlimit=max(differences)
    norms=[SymLogNorm(limit*1e-5,vmin=-limit,vmax=limit)]*2+[SymLogNorm(dlimit*1e-5,vmin=-dlimit,vmax=dlimit)]
    fig,axs=plt.subplots(1,3,figsize=(16,4.4),layout='constrained');images=[]
    titles=['H1：原始地质 / Ez','H0：底部砂岩改泥岩 / Ez','H1−H0：底部材料对照差场 / Ez']
    for j,ax in enumerate(axs):
        im=ax.imshow(frames[0][j],origin='lower',extent=[0,110,0,40],aspect='equal',cmap='gray',norm=norms[j]);images.append(im)
        ax.contour(.0125+np.arange(geo.shape[1])*.2,.0125+np.arange(geo.shape[0])*.2,geo,levels=[.5,1.5,2.5],colors=['#49a9d4','#db9e24','#b54862'],linewidths=.55)
        ax.add_patch(Rectangle((2,2),106,36,fill=False,edgecolor='purple',lw=.8,linestyle='--'))
        ax.plot([78.6,79.9],[35,35],'r.',ms=4)
        ax.set(title=titles[j],xlabel='局部x / m',ylabel='局部y / m')
        cb=fig.colorbar(im,ax=ax,orientation='horizontal',label='Ez / (V/m)；固定对称对数色标',shrink=.9)
        bound=limit if j<2 else dlimit
        ticks=[-bound,-bound*1e-5,0,bound*1e-5,bound]
        cb.set_ticks(ticks);cb.set_ticklabels([f'{value:.1e}' if value else '0' for value in ticks]);cb.ax.tick_params(labelsize=8)
    # Keep both total fields on one fixed range; difference uses its own fixed range.
    selected=np.unique(np.rint(np.linspace(0,255,96)).astype(int))
    movie=[]
    for k in selected:
        for j,im in enumerate(images):im.set_data(frames[k][j])
        fig.suptitle(f'第8道 / X179.25m / V4 FP64：t={times[k]:.2f}ns\n两总场全程同色标；差场单独固定色标；紫虚线内侧为PML接口；地质线按H1叠加')
        fig.canvas.draw();movie.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()).quantize(colors=128))
    movie[0].save(a.out/'wavefield_pair.gif',save_all=True,append_images=movie[1:],duration=160,loop=0,optimize=True)
    for wanted in [120,200,280,350,450,650,790]:
        k=int(np.argmin(abs(np.array(times)-wanted)))
        for j,im in enumerate(images):im.set_data(frames[k][j])
        fig.suptitle(f'第8道 / X179.25m / t={times[k]:.2f}ns；同一套固定色标；差场色标独立')
        fig.savefig(a.out/f'wavefield_{wanted}ns.png',dpi=130)
    plt.close(fig)
    save(a.out/'wavefield_analysis.json',dict(script_sha256=sha(__file__),verification_sha256=sha(a.verification),
        frames_audited=256,GIF_frames=len(selected),GIF_frame_indices=selected.tolist(),time_ns=times,
        sampled_space_step_m=.2,max_abs_total_Ez_by_time=maxima,max_abs_difference_Ez_by_time=differences,
        total_fixed_abs_limit=limit,difference_fixed_abs_limit=dlimit,symlog_linear_threshold_fraction=1e-5,
        limits='0.1m stored sparse snapshots displayed at0.2m; nonuniform time sampling and96 selected animation frames. Observer movie is diagnostic, not full-band phase validation or isolated PML control. Difference has an independent fixed color scale.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['traces','wavefield'])
    for key in ['h1','h0','fp32','manifest','verification','materials','out','cache','h1snap','h0snap','geometry','numerical-out']:p.add_argument('--'+key,type=Path)
    a=p.parse_args();globals()[a.action](a)
