"""Exact-band native probes and native-resolution snapshot delivery.

Snapshot movies are raw impulse fields, never a substitute for the exact
501-tone probe calculation. Total-field flux does not identify unique rays.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha

FREQ=20e6+np.arange(501)*300000.
FIELDS=['Ez','Hx','Hy']
BANDS=['air','cover_upper','cover_middle','cover_lower','mud','basal_below']
LABELS=['地表上方空气','粉质粘土上部','粉质粘土中部','粉质粘土下部','覆底下方泥岩','底砂参考面下方']
# Known matched-anchor events from the completed parent batch. Avoid choosing
# the 190 m first-reflection sidelobe at 222 ns as the later cover return.
UPPER_GATES={170.:(220,270),190.:(300,350)}


def save(p,value):
    p.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def rel(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def read(p):
    return json.loads(p.read_text('utf-8'))


def spectral_tail(raw,dt):
    tail=raw[np.arange(len(raw))*dt>=1800e-9]
    f=np.fft.rfftfreq(len(tail),dt);p=abs(np.fft.rfft(tail*np.hanning(len(tail))))**2
    return dict(diagnostic_window_ns=[1800,2400],diagnostic_window='Hann, not a formal SFCW tail taper',
        peak_Hz=float(f[p.argmax()]),nyquist_Hz=.5/dt,
        fraction_20_170_MHz=float(p[(f>=20e6)&(f<=170e6)].sum()/p.sum()),
        fraction_above_1GHz=float(p[f>1e9].sum()/p.sum())),f,p


def extract(root,g,private):
    raw=Path(g['input']).with_suffix('.h5');v=read(raw.parent/'observer_verification.json')
    assert v['main_source_receiver_bitwise_equal'] and v['native_sha256']==sha(raw)
    source=sf.load_source(raw);receiver=sf.load_receiver(raw,'/rxs/rx1','Ez')
    q=sf.direct_frequency_response(source,receiver,FREQ,tail_taper_fraction=0)
    q=replace(q,response=q.response/.025)
    summary=dict(id=g['id'],native_sha256=sha(raw),verified_snapshots=v['snapshot_count'],
        source_receiver_bitwise_equal=True,snapshot_probe_max_error=v['snapshot_probe_max_normalised_abs_error'],
        full_2400ns_dt_s=source.dt,official_receiver_tail_db=q.receiver_tail_relative_db)
    tail,ff,pp=spectral_tail(receiver.samples,source.dt);summary['late_diagnostic_spectrum']=tail
    values={};native={};band={};checks={}
    if g['probes']:
        with h5py.File(raw) as h:
            E=[];Hx=[];Hy=[]
            for p in g['probes']:
                rr=[h[f"rxs/rx{a['receiver_index']}"] for a in p['anchors'][:3]]
                E.append(rr[0]['Ez'][:]);Hx.append(.5*(rr[0]['Hx'][:]+rr[1]['Hx'][:]))
                Hy.append(.5*(rr[0]['Hy'][:]+rr[2]['Hy'][:]))
        values=dict(Ez=np.array(E).T,Hx=np.array(Hx).T,Hy=np.array(Hy).T)
        native['Ez']=values['Ez'][:-1];native['Hx']=.5*(values['Hx'][:-1]+values['Hx'][1:])
        native['Hy']=.5*(values['Hy'][:-1]+values['Hy'][1:])
        native['Sx']=-native['Ez']*native['Hy'];native['Sy']=native['Ez']*native['Hx']
        np.testing.assert_array_equal(values['Ez'][:,-1],receiver.samples)
        for component,samples in values.items():
            offset=0 if component=='Ez' else -.5*source.dt
            r=sf.SampledSignal(path='spatially_collocated_native_'+component,samples=samples,dt=source.dt,
                               time_offset=offset,quantity=component,units='V/m' if component=='Ez' else 'A/m')
            fq=sf.direct_frequency_response(source,r,FREQ,tail_taper_fraction=0)
            fq=replace(fq,response=fq.response/.025)
            picks=np.unique(np.r_[np.arange(0,501,61),500]);independent=[]
            for f in FREQ[picks]:
                independent.append((np.exp(-2j*np.pi*f*r.times)@samples)/
                                   (np.exp(-2j*np.pi*f*source.times)@source.samples)/.025)
            err=rel(fq.response[picks],np.array(independent));assert err<1e-9
            checks[component+'_independent_DFT_relative_L2']=err
            for window in ['hann','blackman']:
                p=sf.reconstruct_time_response(fq,window=window,zero_pad_factor=8,time_shift=0)
                band[window,component]=p.complex_bandpass
                k=np.arange(13,len(p.time),173)
                z=np.exp(2j*np.pi*p.time[k,None]*FREQ)@(fq.response*p.weights[:,None])/501
                er=rel(z,p.complex_bandpass[k]);assert er<1e-9
                checks[component+'_'+window+'_independent_inverse_relative_L2']=er
        for window in ['hann','blackman']:
            e,hx,hy=[band[window,k] for k in FIELDS]
            # This index is twice the usual analytic-field cycle-average
            # expression. Retain an explicit source-normalised direction
            # index; do not call it individual path energy.
            band[window,'Sx_index']=-np.real(e*np.conj(hy))
            band[window,'Sy_index']=np.real(e*np.conj(hx))
    profiles={w:sf.reconstruct_time_response(q,window=w,zero_pad_factor=8,time_shift=0) for w in ['hann','blackman']}
    path=private/(g['id']+'_arrays.h5');assert not path.exists()
    with h5py.File(path,'x') as h:
        h.attrs['source_native_sha256']=summary['native_sha256'];h.attrs['method']='official_direct'
        h.attrs['limits']='Snapshot=raw impulse; band transport=source-normalised total-field index, not separated pulse energy'
        h['native_time_s']=np.arange(len(receiver.samples))*source.dt;h['native_main_Ez']=receiver.samples
        h['frequency_Hz']=FREQ;h['main_response']=q.response;h['sfcw_time_s']=profiles['hann'].time
        h['tail_diagnostic_frequency_Hz']=ff;h['tail_diagnostic_power']=pp
        if g['probes']:
            h['probe_positions_m']=[[p['x_m'],p['y_m']] for p in g['probes']]
            h['probe_names']=np.array([p['id'] for p in g['probes']],dtype=h5py.string_dtype())
            for name,z in values.items():h.create_dataset('native_staggered/'+name,data=z,compression='gzip',compression_opts=2)
            for name,z in native.items():h.create_dataset('native_aligned/'+name,data=z,compression='gzip',compression_opts=2)
            for (window,name),z in band.items():h.create_dataset(window+'/probes/'+name,data=z)
        for w,p in profiles.items():h[w+'/main_complex_bandpass']=p.complex_bandpass
    summary.update(independent_checks=checks,private_arrays_sha256=sha(path))
    return summary,path


def plot_pair(x,groups,arrays,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    g={role:next(g for g in groups if g['profile_x_m']==x and g['kind']=='local' and g['role']==role) for role in ['H1','H0']}
    with h5py.File(arrays[g['H1']['id']]) as a,h5py.File(arrays[g['H0']['id']]) as b:
        t=a['sfcw_time_s'][:]*1e9;p=a['hann/main_complex_bandpass'][:];q=b['hann/main_complex_bandpass'][:];d=p-q
        fig,axes=plt.subplots(3,1,figsize=(12,11),constrained_layout=True)
        for ax in axes[:2]:
            ax.plot(t,abs(p),label='H1 总场');ax.plot(t,abs(q),label='H0 保留覆盖层、去底砂反差');ax.plot(t,abs(d),label='同站复差 H1−H0')
            ax.set(xlabel='官方 direct SFCW 时间 / ns',ylabel='复包络 / (V/m)/(A·m)');ax.legend()
        structure=(t>=100)&(t<=550)
        axes[0].set(xlim=[100,550],ylim=[0,1.15*max(float(abs(z[structure]).max()) for z in [p,q,d])],title=f'剖面{x:g}m：覆盖层首次返回、后续强返回和底砂响应；按本显示窗设纵轴')
        axes[1].set(xlim=[220,450],ylim=[0,1.15*max(abs(p[(t>=220)&(t<=450)]))],title=f'{x:g}m：后续覆盖层返回和底砂响应；未经AGC或拟合')
        for h,name in [(a,'H1'),(b,'H0')]:
            ff=h['tail_diagnostic_frequency_Hz'][:]/1e9;pp=h['tail_diagnostic_power'][:]
            axes[2].plot(ff,10*np.log10(np.maximum(pp/pp.max(),1e-30)),label=name)
        axes[2].axvspan(.02,.17,color='green',alpha=.2,label='20–170MHz 工作带')
        axes[2].set(xlabel='频率 / GHz',ylabel='诊断谱 / 相对本谱峰 dB',xlim=[0,8.5],ylim=[-160,5],
                    title='1800–2400ns原生残留：Hann诊断谱；未给正式SFCW添加尾窗');axes[2].legend()
        fig.savefig(out/f'x{int(x)}_main_and_tail.png',dpi=160);plt.close(fig)
        ids=[[i for i,p in enumerate(g['H1']['probes']) if p['band']==band] for band in BANDS]
        fig,axes=plt.subplots(6,3,figsize=(15,18),constrained_layout=True)
        h1=a['hann/probes/Sy_index'][:];h0=b['hann/probes/Sy_index'][:]
        for row,(band,label,index) in enumerate(zip(BANDS,LABELS,ids)):
            vals=[h1[:,index],h0[:,index],h1[:,index]-h0[:,index]]
            mask=(t>=0)&(t<=550);scale=max(float(abs(z[mask]).max()) for z in vals)
            norm=SymLogNorm(max(scale*1e-6,1e-30),vmin=-scale,vmax=scale)
            for col,(z,title) in enumerate(zip(vals,['H1净方向指标','H0净方向指标','指标差（含交叉项）'])):
                im=axes[row,col].imshow(z[mask],aspect='auto',origin='lower',extent=[-.5,4.5,t[mask][0],t[mask][-1]],cmap='RdBu_r',norm=norm)
                axes[row,col].set(title=label+' / '+title,xlabel='同次激发的观察列（站位偏移−12、−8、−4、0、+4m）',ylabel='SFCW 时间 / ns')
                axes[row,col].invert_yaxis();fig.colorbar(im,ax=axes[row,col],fraction=.035)
        fig.suptitle(f'{x:g}m：完整原生时序→官方501点→Hann；红上行、蓝下行的总场指标\n每行三图同尺度；固定对数显示，非移动测线B-scan，非单次反射能量',fontsize=13)
        fig.savefig(out/f'x{int(x)}_band_direction.png',dpi=150);plt.close(fig)
        event_rows=[]
        for role,h in [('H1',a),('H0',b)]:
            sy=h['hann/probes/Sy_index'][:];sx=h['hann/probes/Sx_index'][:]
            for i,probe in enumerate(g[role]['probes']):
                windows=[('early',[40,220]),('cover_return',[220,340]),('basal',[330,430])]
                events={}
                for name,(lo,hi) in windows:
                    ix=np.flatnonzero((t>=lo)&(t<=hi));j=ix[np.argmax(abs(sy[ix,i]))];near=(t>=t[j]-4)&(t<=t[j]+4)
                    vx=float(sx[near,i].sum());vy=float(sy[near,i].sum())
                    events[name]=dict(time_ns=float(t[j]),Sy_index=float(sy[j,i]),sum_Sx_index=vx,sum_Sy_index=vy,
                                      direction_deg=float(np.degrees(np.arctan2(vy,vx))),extremum_at_gate_edge=bool(j in [ix[0],ix[-1]]))
                event_rows.append(dict(role=role,id=probe['id'],band=probe['band'],position_m=[probe['x_m'],probe['y_m']],events=events))
        main_directions=[]
        sx=a['hann/probes/Sx_index'][:,-1];sy=a['hann/probes/Sy_index'][:,-1]
        for label,z,(lo,hi) in [('first_cover',q,(100,220)),('upper_return',q,UPPER_GATES[x]),('basal_delta',d,(330,430))]:
            ix=np.flatnonzero((t>=lo)&(t<=hi));j=ix[np.argmax(abs(z[ix]))];near=(t>=t[j]-4)&(t<=t[j]+4)
            main_directions.append(dict(event=label,gate_ns=[lo,hi],time_ns=float(t[j]),peak_abs=float(abs(z[j])),
                H1_net_band_direction_deg=float(np.degrees(np.arctan2(sy[near].sum(),sx[near].sum())))))
        candidates=[]
        for role,h in [('H1',a),('H0',b)]:
            sy=h['hann/probes/Sy_index'][:];sx=h['hann/probes/Sx_index'][:]
            for i,probe in enumerate(g[role]['probes']):
                if probe['band']!='air':continue
                lo,hi=(150,245) if x==170 else (245,350)
                ix=np.flatnonzero((t>=lo)&(t<=hi));j=ix[np.argmax(sy[ix,i])];near=(t>=t[j]-4)&(t<=t[j]+4)
                delta=np.array(g[role]['rx_m'][:2])-np.array([probe['x_m'],probe['y_m']])
                delay=float(np.linalg.norm(delta)/299792458*1e9)
                candidates.append(dict(role=role,id=probe['id'],profile_x_m=probe['x_m']+20,
                    up_Sy_extremum_ns=float(t[j]),geometric_air_delay_ns=delay,sum_arrival_ns=float(t[j]+delay),
                    measured_direction_deg=float(np.degrees(np.arctan2(sy[near,i].sum(),sx[near,i].sum()))),
                    direct_to_rx_angle_deg=float(np.degrees(np.arctan2(delta[1],delta[0]))),
                    extremum_at_gate_edge=bool(j in [ix[0],ix[-1]])))
        itinerary=[]
        if x==170:
            fig,axes=plt.subplots(4,1,figsize=(12,12),constrained_layout=True)
            for row,band in enumerate(['cover_upper','cover_middle','cover_lower']):
                i=next(i for i,p in enumerate(g['H1']['probes']) if p['band']==band and p['x_m']==150.)
                for role,h in [('H1',a),('H0',b)]:
                    s=h['hann/probes/Sy_index'][:,i]
                    axes[row].plot(t,s,label=role,linestyle='-' if role=='H1' else '--',linewidth=1.5)
                    for lo,hi,sign,name in [(20,80,-1,'首次下行'),(80,125,1,'首次上行'),(125,165,-1,'再次下行'),(175,215,1,'再次上行')]:
                        k=np.flatnonzero((t>=lo)&(t<=hi));j=k[np.argmax(sign*s[k])]
                        itinerary.append(dict(role=role,probe=g[role]['probes'][i]['id'],band=band,event=name,
                            gate_ns=[lo,hi],time_ns=float(t[j]),Sy_index=float(s[j]),extremum_at_gate_edge=bool(j in [k[0],k[-1]])))
                        if role=='H1':axes[row].annotate(f'{name}\n{t[j]:.2f}ns',(t[j],s[j]),xytext=(4,14 if sign<0 else -32),textcoords='offset points',fontsize=9)
                axes[row].set_yscale('symlog',linthresh=float(abs(a['hann/probes/Sy_index'][:,i]).max())*1e-5)
                axes[row].set(xlim=[20,220],ylabel='Sy 指标（固定对数），正上行 / 负下行');axes[row].set_title(LABELS[BANDS.index(band)]+'：同一位置的工作带总场方向指标',pad=18);axes[row].axhline(0,color='grey',linewidth=.5);axes[row].legend()
            axes[3].plot(t,abs(p),label='H1 总场');axes[3].plot(t,abs(q),'--',label='H0 去底砂反差');axes[3].plot(t,abs(d),label='复差包络')
            axes[3].set(xlim=[220,440],ylim=[0,1.2],xlabel='官方 direct SFCW 时间 / ns',ylabel='接收复包络 / (V/m)/(A·m)',title='约232ns强返回在H0中保留；约356ns出现底砂相关复差');axes[3].legend()
            fig.suptitle('170m覆盖层往返证据：原生完整时序→官方501频点→Hann\nH0/H1往返次序相同；指标为总场，未拆分独立光线或反射能量')
            fig.savefig(out/'x170_cover_roundtrip.png',dpi=170);plt.close(fig)
    save(out/f'x{int(x)}_probe_events.json',dict(status='EXPLORATORY_TOTAL_FIELD_BAND_EXTREMA',events=event_rows,
        main_peak_directions=main_directions,air_return_candidates=candidates,cover_itinerary=itinerary,
        limits='Absolute Sy maximum may select standing/interfering components. Gate edges kept; not unique path, direction uncertainty or bounce certification.'))


def movies(groups,out,package):
    from PIL import Image,ImageDraw,ImageFont
    import matplotlib
    lut=(matplotlib.colormaps['RdBu_r'](np.linspace(0,1,256))[:,:3]*255).astype('uint8')
    font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',22)
    with h5py.File(package/'geometries/H1.h5') as h:geo=h['data'][:,:,0]
    receipts=[]
    for x,kind in [(170.,'local'),(190.,'local'),(170.,'global')]:
        pair={role:next(g for g in groups if g['profile_x_m']==x and g['kind']==kind and g['role']==role) for role in ['H1','H0']}
        entries={role:read(Path(g['input']).parent/'observer_verification.json')['snapshots'] for role,g in pair.items()}
        nlist=pair['H1']['snapshot_iterations'];roi=pair['H1']['roi'];dt=pair['H1']['dt_s'];nx,ny,_=roi['shape']
        def load(role,k):
            row=entries[role][k];p=Path(pair[role]['input']).parent/'profile_snaps'/row['file']
            assert sha(p)==row['sha256']
            with h5py.File(p) as h:v=h['Ez'][:,:,0]
            receipts.append(dict(group=pair[role]['id'],iteration=row['iteration'],sha256=row['sha256']))
            return v
        # A common physical scale from both observed fields, reused throughout
        # all clips. Fixed log display does not alter the numeric arrays.
        scale=0.
        for k in range(len(nlist)):
            h1,h0=load('H1',k),load('H0',k)
            scale=max(scale,float(abs(h1).max()),float(abs(h0).max()),float(abs(h1-h0).max()))
        linear=scale*1e-6
        xcentres=roi['bounds'][0]+(np.arange(nx)+.5)*roi['spacing'][0]
        columns=np.floor(xcentres/.025).astype(int)
        surface=np.max(np.where(geo[columns]!=0,np.arange(1700)[None,:]+1,0),axis=1)*.025
        cover=np.min(np.where(geo[columns]==1,np.arange(1700)[None,:],1700),axis=1)*.025
        basal=np.max(np.where(geo[columns]==3,np.arange(1700)[None,:]+1,0),axis=1)*.025
        def pixel(xx,yy):return ((xx-roi['bounds'][0])/roi['spacing'][0]-.5,ny-.5-(yy-roi['bounds'][1])/roi['spacing'][1])
        def frame(k):
            h1,h0=load('H1',k),load('H0',k);values=[h1,h0,h1-h0]
            canvas=Image.new('RGB',(3*nx+40,ny+140),'white');draw=ImageDraw.Draw(canvas)
            draw.text((12,8),f'剖面{x:g}m / {kind} / 原始impulse时刻 {nlist[k]*dt*1e9:.3f}ns / 非SFCW场图',fill='black',font=font)
            for j,(z,title) in enumerate(zip(values,['H1 原始总场 Ez','H0 去底砂反差 Ez','H1−H0 场差（含交互）'])):
                q=np.sign(z)*np.log1p(abs(z)/linear)/np.log1p(scale/linear)
                idx=np.clip(np.rint((q+1)*127.5),0,255).astype('uint8')
                img=Image.fromarray(lut[idx.T[::-1]],'RGB');left=10+j*(nx+10);top=70
                canvas.paste(img,(left,top));draw.text((left,40),title,fill='black',font=font)
                for yy,color in [(surface,'green'),(cover,'darkorange'),(basal,'black')]:
                    line=[(left+i,top+pixel(xcentres[i],yy[i])[1]) for i in range(nx)]
                    draw.line(line,fill=color,width=2)
                for pos,tag in [(pair['H1']['tx_m'],'Tx'),(pair['H1']['rx_m'],'Rx')]:
                    px,py=pixel(pos[0],pos[1]);px+=left;py+=top
                    draw.ellipse((px-4,py-4,px+4,py+4),fill='yellow',outline='black');draw.text((px+5,py-12),tag,fill='black',font=font)
                if kind=='global':
                    for xx in [2.,208.]:
                        px=left+pixel(xx,0)[0];draw.line((px,top,px,top+ny-1),fill='purple',width=2)
                    for yy in [2.,40.5]:
                        py=top+pixel(0,yy)[1];draw.line((left,py,left+nx-1,py),fill='purple',width=2)
                    draw.text((left+15,top+15),'紫线：80格PML内沿',fill='purple',font=font)
            draw.text((12,ny+78),f'空间采样{roi["spacing"][0]:g}m；三图共用 ±{scale:.3g} V/m，固定对数阈值 {linear:.3g}；红/蓝为正/负',fill='black',font=font)
            draw.text((12,ny+108),f'x为剖面{roi["bounds"][0]+20:g}→{roi["bounds"][3]+20:g}m；y为模型内部{roi["bounds"][1]:g}→{roi["bounds"][4]:g}m；H0黑线仅为底砂参考位置',fill='black',font=font)
            return canvas
        clips=[('first',120,220,1),('return',220,340,1),('basal',340,440,1)] if kind=='local' else [('global_early',0,600,5),('global_late',600,2400,17)]
        for name,lo,hi,skip in clips:
            selected=[k for k,n in enumerate(nlist) if lo<=n*dt*1e9<=hi][::skip]
            images=[frame(k) for k in selected];assert images
            path=out/f'x{int(x)}_{kind}_{name}.gif'
            images[0].save(path,save_all=True,append_images=images[1:],duration=70,loop=0,optimize=False)
            receipts.append(dict(movie=path.name,frames=len(images),bytes=path.stat().st_size,sha256=sha(path),
                native_pixel_resolution=[nx,ny],shared_vmax_V_m=scale,shared_log_threshold_V_m=linear,
                temporal_display_decimation=skip,frame_step_ns=9*dt*1e9 if kind=='local' else skip*17*dt*1e9))
            preview=[]
            for im in images[::2]:
                small=im.resize((im.width//2,im.height//2),resample=Image.Resampling.BOX)
                view=Image.new('RGB',(small.width,small.height+35),'white');view.paste(small,(0,35))
                ImageDraw.Draw(view).text((8,4),'预览：图像缩小2倍，每2保存帧展示；原生图和动画另存',fill='black',font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',16))
                preview.append(view);small.close()
            pp=out/(path.stem+'_preview.gif')
            preview[0].save(pp,save_all=True,append_images=preview[1:],duration=140,loop=0,optimize=False)
            receipts.append(dict(movie=pp.name,frames=len(preview),bytes=pp.stat().st_size,sha256=sha(pp),
                source_movie=path.name,spatial_display_decimation=2,temporal_display_decimation=2*skip,display_only=True))
            for im in preview:im.close()
            for im in images:im.close()
            print(json.dumps(dict(rendered=path.name,bytes=path.stat().st_size,preview_bytes=pp.stat().st_size)),flush=True)
        wanted=[150,200,240,280,320,360,400,440] if kind=='local' else [100,200,400,600,1200,1800,2300]
        for target in wanted:
            k=int(np.argmin(abs(np.array(nlist)*dt*1e9-target)))
            frame(k).save(out/f'x{int(x)}_{kind}_{target}ns.png')
    save(out/'snapshot_render_receipts.json',dict(status='NATIVE_SPATIAL_RESOLUTION_FIXED_DISPLAY',script_sha256=sha(__file__),receipts=receipts,
        limits='Raw broadband snapshots have sampling limitations. Native-grid pixels are retained in local clips; arrows/unique-ray attribution are not invented. Exact-band claims use full-time probe data.'))


def main(a):
    assert not a.out.exists();a.out.mkdir(parents=True);a.private.mkdir(parents=True,exist_ok=True)
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert len(v['groups'])==6
    rows=[];arrays={}
    for g in c['groups']:
        row,path=extract(a.source,g,a.private);rows.append(row);arrays[g['id']]=path
        print(json.dumps(dict(extracted=g['id'])),flush=True)
    for x in [170.,190.]:plot_pair(x,c['groups'],arrays,a.out)
    save(a.out/'analysis.json',dict(status='SIX_NATIVE_FP64_OBSERVER_CASES_ANALYZED',solver_runs=6,
        script_sha256=sha(__file__),official_processing_sha256=sha(sf.__file__),contract_sha256=sha(a.source/'execution_contract.json'),
        completed_certificate_sha256=sha(a.source/'completed_verification.json'),groups=rows,
        limits=c['study_manifest']['limits'],frequency_band_Hz=[20e6,170e6],tones=501))
    if not a.skip_movies:movies(c['groups'],a.out,Path(c['package']))
    print(json.dumps(dict(status='DELIVERED',out=str(a.out))),flush=True)


def refresh_verified_plots(source,arrays,out):
    """Correct diagnostic plot/gate metadata; retain numerical and movie hashes."""
    c=read(source/'execution_contract.json');v=read(source/'completed_verification.json')
    assert v['completed'] and len(v['groups'])==6
    old=read(out/'analysis.json')
    assert old['contract_sha256']==sha(source/'execution_contract.json')
    mapping={g['id']:arrays/(g['id']+'_arrays.h5') for g in c['groups']}
    for row in old['groups']:assert sha(mapping[row['id']])==row['private_arrays_sha256']
    for x in [170.,190.]:plot_pair(x,c['groups'],mapping,out)
    old['plot_revision']=dict(status='R2_DIAGNOSTIC_GATE_AND_VISIBLE_AXIS_CORRECTION',
        script_sha256=sha(__file__),upper_gates_ns={str(x):list(g) for x,g in UPPER_GATES.items()},
        reason='190m broad 220–350ns gate selected the first-return sidelobe at222ns. Parent matched event is327.678ns. Main display axis now uses only displayed samples.',
        native_and_private_arrays_unchanged=True,formal_SFCW_unchanged=True,new_solver_runs=0,
        movie_render_script_sha256=read(out/'snapshot_render_receipts.json')['script_sha256'])
    save(out/'analysis.json',old)
    write_gallery(out)


def write_gallery(out):
    """Static scientific report; animations load only on deliberate clicks."""
    panels=[]
    for x in [170,190]:
        for kind,stamp,label in [('return',240,'覆盖层返回阶段'),('basal',360,'底砂返回阶段')]:
            gif=f'x{x}_local_{kind}_preview.gif';png=f'x{x}_local_detail_{stamp}ns.png'
            assert (out/gif).exists() and (out/png).exists()
            panels.append(f'<article><h2>{x}m · {label}</h2><p>原始宽带冲激场，三图同尺度：H1、H0、场差。原生图每面板800×1120；动画为明确标注的紧凑预览。</p><button onclick="play(this,\'{gif}\')">播放预览</button> <button onclick="play(this,\'{png}\')">显示静图</button> <a href="{png}" target="_blank">打开原生静图并放大</a><div class="picture"><img loading="lazy" src="{png}" alt="{x}m {label}三幅配对波场"></div></article>')
    content='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Line9 精细波场诊断</title>
<style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#182231;background:#f3f5f7;max-width:1350px;margin:auto;padding:25px}article,.lead{background:white;padding:20px;margin:20px 0;border-radius:10px}h1,h2{line-height:1.3}button,a{font:inherit;margin-right:12px}button{padding:7px 14px;cursor:pointer}.picture{overflow:auto;background:#fff;margin-top:15px}img{max-width:100%;height:auto}.warning{border-left:4px solid #ba7b19;padding-left:16px}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:7px 15px}</style>
<h1>覆盖层多次返回与底砂响应</h1><div class="lead"><p>六项FP64求解、8122帧。170m的约232ns强带在去底砂反差H0中保留；约356ns出现底砂相关复差。两轮覆盖层上下传播在H0/H1中保持。</p><p class="warning">动画是原始宽带impulse场，含密集高频数值纹理；它不是20–170MHz波场。工作频带结论以完整原生探针和官方501点direct处理为依据。预览缩小不改变数值记录，播放速度不是实际传播时间。</p><a href="x170_cover_roundtrip.png" target="_blank">放大主要工作带证据</a><img src="x170_cover_roundtrip.png" alt="覆盖层两轮往返与接收对照"><p><a href="x190_main_and_tail.png">190m工作带与尾段谱</a> <a href="x170_main_and_tail.png">170m工作带与尾段谱</a> <a href="mechanism_checks.json">无拟合平层估算与边界贡献诊断</a> <a href="delivery_verification.json">取回及数值验收</a></p></div>'''
    content+=''.join(panels)
    content+='''<article><h2>全域晚期 · 边界观察</h2><p>紫线为PML内沿。输出20cm；晚期预览还抽帧，未将亮纹直接称为边界反射。</p><button onclick="play(this,'x170_global_global_late_preview.gif')">播放全域晚期预览</button> <button onclick="play(this,'x170_global_1800ns.png')">显示1800ns静图</button><div class="picture"><img loading="lazy" src="x170_global_1800ns.png" alt="170m全域晚期配对场"></div></article><article><h2>证据边界</h2><p>支持覆盖层重复返回及侧偏传播；底砂对应回波已可见。末尾约8.465GHz纹理与工作带强条纹须分开。未认证真实天线/横向收发布局、现场复电性或3D/网格/PML/时窗全面收敛。</p><p>完整原生动画与约90.6GiB快照保留在ROG：E:\\djh\\line9_fine_wavefield_20261011_r1\\delivery_r1 和 execution\\prepared\\cases。完整路径及哈希见public_retrieval_manifest.json、snapshot_render_receipts.json。</p></article><script>function play(button,path){button.closest('article').querySelector('img').src=path}</script></html>'''
    (out/'wavefield_report.html').write_text(content,encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','out','private']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--skip-movies',action='store_true');p.add_argument('--refresh-plots',action='store_true')
    a=p.parse_args()
    if a.refresh_plots:refresh_verified_plots(a.source,a.private,a.out)
    else:main(a)
