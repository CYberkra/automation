"""Read-only audit and comparable reconstruction of the three supplied Line9 result packages."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

C0=299792458.
EPS0=8.8541878128e-12
FREQ=20e6+300000*np.arange(501)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def rel(a,b):
    denominator=np.linalg.norm(b)
    if denominator==0: raise ValueError('Zero reference')
    return float(np.linalg.norm(a-b)/denominator)


def indices(materials,f):
    values=[]
    for material in materials.values():
        base=material['base']
        er=complex(base['relative_permittivity'],-base['electric_conductivity_s_per_m']/(2*np.pi*f*EPS0))
        for pole in material.get('poles',[]):
            er+=pole['relative_permittivity_difference']/(1+2j*np.pi*f*pole['relaxation_time_s'])
        values.append(np.sqrt(er))
    return np.asarray(values)


def ray_time(thickness,index,separation):
    # Local horizontal layers, symmetric bistatic ray, phase indices at95MHz.
    keep=thickness>0;d=thickness[keep];n=index[keep].real
    if not len(d): raise ValueError('Empty path')
    low=0.;high=float(n.min())*(1-1e-12)
    for _ in range(60):
        p=(low+high)/2
        offset=float(np.sum(2*d*p/np.sqrt(n*n-p*p)))
        if offset<separation: low=p
        else: high=p
    p=(low+high)/2
    return float(np.sum(2*d*n*n/np.sqrt(n*n-p*p))/C0*1e9)


def geometry_at(data,spacing,tx,rx,index):
    midpoint=(tx+rx)/2;ix=round(midpoint[0]/spacing[0]);column=data[ix]
    top=int(np.flatnonzero(column!=0)[-1]+1)
    ground=top*spacing[1];agl=float(midpoint[1]-ground)
    changes=np.flatnonzero(np.diff(column)!=0)+1
    boundaries=[]
    for j in changes[::-1]:
        if j>top: continue
        above=int(column[j]);below=int(column[j-1])
        thickness=np.bincount(column[j:top],minlength=4)*spacing[1]
        thickness[0]+=agl
        boundaries.append(dict(y_m=float(j*spacing[1]),depth_m=float((top-j)*spacing[1]),above=above,below=below,
            time95_ns=ray_time(thickness,index,float(abs(rx[0]-tx[0]))),
            vertical_absorption95_db=float(-40/np.log(10)*np.sum(thickness[1:]*(-2*np.pi*95e6/C0*index[1:].imag)))))
    def pick(rows,which):
        if not rows: return None
        return rows[which]
    cover=pick([b for b in boundaries if b['above']==1 and b['below']!=1],0)
    sand=[b for b in boundaries if b['below']==3 and b['above']!=3]
    return dict(midpoint_agl_m=agl,surface_y_m=float(ground),boundaries=boundaries,
        surface=boundaries[0],cover_base=cover,first_sand=pick(sand,0),basal_sand=pick(sand,-1))


def package(path,out):
    manifest=json.loads((path/'package_manifest.json').read_text('utf-8'))
    cases={c['id']:c for c in manifest['cases']}
    if len(cases)!=manifest['station_count']: raise ValueError('Duplicate manifest cases')
    rawpaths=sorted((path/'cases').glob('*/profile.h5'))
    storedpaths={p.name.replace('_sfcw.h5',''):p for p in (path/'sfcw').glob('*_sfcw.h5')}
    events=[json.loads(line) for line in (path/'sweep_progress.jsonl').read_text('utf-8').splitlines()]
    completed=[r['station'] for r in events if r['ok'] and r['returncode']==0]
    ids=[p.parent.name for p in rawpaths]
    if len(set(completed))!=len(completed) or set(completed)!=set(ids) or set(ids)!=set(storedpaths):
        raise ValueError('Raw/SFCW/unique completion mismatch')
    if any(not e['ok'] or e['returncode'] for e in events): raise ValueError('Failed rows require separate audit')
    material_path=next((path/'geometries').glob('*.json'))
    materials=json.loads(material_path.read_text('utf-8'))['materials']
    n95=indices(materials,95e6)
    geometry=next((path/'geometries').glob('*.h5'))
    with h5py.File(geometry) as h:
        data=h['data'][:,:,0];spacing=np.array(h.attrs['dx_dy_dz']);shape=h['data'].shape
        keys=[k.decode() for k in h['material_keys'][:]]
    if keys!=list(materials): raise ValueError('Material voxel mapping mismatch')
    histories=[];spectra=[];source_spectra=[];stored_band=[];records=[];sources=[]
    source_offsets=set();receiver_offsets=set();iterations=set();dtypes=set();dt_values=set();tapers=set()
    for raw in rawpaths:
        sid=raw.parent.name;case=cases[sid]
        card=(raw.parent/'profile.in').read_text('utf-8')
        commands={s.split(':',1)[0][1:]:s.split(':',1)[1].strip().split() for s in card.splitlines() if s.startswith('#') and ':' in s}
        with h5py.File(raw) as h:
            rx=h['rxs/rx1/Ez'][:];src=h['srcs/src1/excitation/samples'][:]
            ea=h['srcs/src1/excitation'].attrs
            dt=float(h.attrs['dt']);ni=int(h.attrs['Iterations'])
            tx=np.array(h['srcs/src1'].attrs['Position']);receiver=np.array(h['rxs/rx1'].attrs['Position'])
            if str(h.attrs['gprMax'])!='4.0.0' or rx.shape!=(ni,) or src.shape!=(ni,) or not np.isfinite(rx).all():
                raise ValueError('Native schema/version/finiteness')
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],shape)
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],spacing)
            np.testing.assert_allclose(tx[:2],np.array(commands['hertzian_dipole'][1:3],float),atol=1e-10,rtol=0)
            np.testing.assert_allclose(receiver[:2],np.array(commands['rx'][:2],float),atol=1e-10,rtol=0)
            if commands['waveform']!=['impulse','40','1','pulse'] or np.count_nonzero(src)!=1 or src[0]!=40:
                raise ValueError('Unexpected impulse source')
            if str(ea['Polarisation'])!='z' or float(ea['SpatialScale'])!=.025: raise ValueError('Source identity')
            if not np.isclose(dt,.025/(C0*np.sqrt(2)),rtol=1e-12): raise ValueError('CFL step differs')
            if abs((ni-1)*dt-float(commands['time_window'][0]))>dt: raise ValueError('Native record length')
            for pos in [tx,receiver]:
                cell=np.rint(pos[:2]/spacing[:2]).astype(int)
                if data[tuple(cell)]!=0: raise ValueError('Antenna not in air')
            source_offset=float(ea['TimeSampleOffset']);receiver_offset=float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset'])
            source_offsets.add(source_offset);receiver_offsets.add(receiver_offset)
            dtypes.add(str(rx.dtype));iterations.add(ni);dt_values.add(dt)
        with h5py.File(storedpaths[sid]) as h:
            np.testing.assert_array_equal(h['frequency'][:],FREQ)
            if not h['source_valid'][:].all(): raise ValueError('Invalid source frequency')
            spectrum=h['response'][:];ss=h['source_spectrum'][:]
            np.testing.assert_allclose(h['I'][:],spectrum.real,rtol=1e-12,atol=1e-25)
            np.testing.assert_allclose(h['Q'][:],spectrum.imag,rtol=1e-12,atol=1e-25)
            fraction=float(h.attrs['TailTaperFraction']);tapers.add(fraction)
            np.testing.assert_allclose(h['time_response/weights'][:],sf.spectral_window('gaussian',501),rtol=1e-12)
            if h['time_response'].attrs['Window']!='gaussian': raise ValueError('Unexpected package window')
            stored_band.append(h['time_response/complex_bandpass'][:]);times=h['time_response/time'][:]
        geom=geometry_at(data,spacing,tx,receiver,n95)
        top=float(shape[1]*spacing[1]-2)
        mirror_top=np.sqrt((receiver[0]-tx[0])**2+(2*top-tx[1]-receiver[1])**2)/C0*1e9
        side_times=[]
        for xface in [2,shape[0]*spacing[0]-2]:
            side_times.append(float(np.hypot(2*xface-tx[0]-receiver[0],tx[1]-receiver[1])/C0*1e9))
        records.append(dict(id=sid,chainage_m=case['chainage_m'],native_sha256=sha(raw),sfcw_sha256=sha(storedpaths[sid]),
            input_sha256=sha(raw.parent/'profile.in'),geometry=geom,
            top_buffer_to_pml_m=float(top-max(tx[1],receiver[1])),top_pml_air_mirror_ns=float(mirror_top),
            closest_side_pml_air_mirror_ns=min(side_times),right_side_pml_air_mirror_ns=side_times[1],
            raw_peak=float(np.max(abs(rx))),raw_last5percent_peak_over_peak=float(np.max(abs(rx[-max(8,int(np.ceil(.05*ni))):]))/np.max(abs(rx)))))
        histories.append(rx);sources.append(src);spectra.append(spectrum);source_spectra.append(ss)
    if any(len(v)!=1 for v in [iterations,dt_values,source_offsets,receiver_offsets,tapers,dtypes]):
        raise ValueError('Inconsistent native sampling/processing')
    rx=np.column_stack(histories);dt=dt_values.pop();ni=iterations.pop();fraction=tapers.pop()
    source=np.asarray(sources[0],float)
    if any(not np.array_equal(source,s) for s in sources): raise ValueError('Source histories differ')
    ss=dt*40*np.exp(-2j*np.pi*FREQ*next(iter(source_offsets)))
    stored=np.column_stack(spectra)
    np.testing.assert_allclose(np.column_stack(source_spectra),ss[:,None]*np.ones((1,len(ids))),rtol=1e-9,atol=1e-25)
    package_samples=min(ni,int(np.floor(1200e-9/dt))+1)
    tapered=sf.apply_tail_taper(rx[:package_samples],fraction)
    recomputed=sf.engineering_dft(tapered,dt,FREQ,time_offset=next(iter(receiver_offsets)))/ss[:,None]
    raw_response=sf.engineering_dft(rx,dt,FREQ,time_offset=next(iter(receiver_offsets)))/ss[:,None]
    error=rel(recomputed,stored)
    if error>1e-9: raise ValueError('Packaged spectrum does not reproduce raw+declared taper')
    take=[0,83,167,250,333,417,500]
    direct=dt*(np.exp(-2j*np.pi*FREQ[take,None]*(np.arange(package_samples)*dt+next(iter(receiver_offsets))))@tapered)/ss[take,None]
    independent=rel(recomputed[take],direct)
    if independent>1e-9: raise ValueError('Independent direct transform mismatch')
    product=sf.direct_frequency_response(sf.load_source(rawpaths[0]),sf.load_receiver(rawpaths[0],'/rxs/rx1','Ez'),FREQ)
    profiles={}
    inverse_errors={}
    for window in ['gaussian','hann','rectangular']:
        tr=sf.reconstruct_time_response(replace(product,response=raw_response/.025),window=window,zero_pad_factor=8)
        profiles[window]=tr
        ticks=np.arange(0,4008,97)
        direct_inverse=(np.exp(2j*np.pi*tr.time[ticks,None]*FREQ)@(tr.weights[:,None]*raw_response/.025))/501
        inverse_errors[window]=rel(tr.complex_bandpass[ticks],direct_inverse)
        if inverse_errors[window]>1e-9: raise ValueError('Independent inverse mismatch')
    old=sf.reconstruct_time_response(replace(product,response=stored/.025),window='gaussian',zero_pad_factor=8)
    inverse_packaged=rel(old.complex_bandpass*.025,np.column_stack(stored_band))
    if inverse_packaged>1e-9: raise ValueError('Stored Gaussian reconstruction mismatch')
    for j,r in enumerate(records):
        r['timing_gates']={}
        for name in ['cover_base','first_sand','basal_sand']:
            boundary=r['geometry'][name]
            if boundary is None: continue
            center=boundary['time95_ns'];t=profiles['gaussian'].time*1e9
            gate=abs(t-center)<=10;guard=((abs(t-center)>=25)&(abs(t-center)<=45))
            value=abs(profiles['gaussian'].complex_envelope[:,j]);base=float(np.mean(value[guard]))
            r['timing_gates'][name]=dict(center95_ns=center,mean_magnitude=float(np.mean(value[gate])),
                peak_magnitude=float(np.max(value[gate])),local_guard_contrast_db=float(20*np.log10(np.mean(value[gate])/base)))
    times=profiles['gaussian'].time
    mask=(times>=320e-9)&(times<=497e-9)
    sensitivity={}
    for lo,hi in [(100,300),(300,600),(600,800)]:
        use=(times>=lo*1e-9)&(times<hi*1e-9)
        sensitivity[f'{lo}_{hi}ns']=rel(old.complex_bandpass[use],profiles['gaussian'].complex_bandpass[use])
    ranges={}
    for key in ['surface','cover_base','first_sand','basal_sand']:
        values=[r['geometry'][key] for r in records if r['geometry'][key] is not None]
        ranges[key]=dict(depth_m=[min(v['depth_m'] for v in values),max(v['depth_m'] for v in values)],
            time95_ns=[min(v['time95_ns'] for v in values),max(v['time95_ns'] for v in values)])
    stats=dict(package=path.name,raw_count=len(ids),planned_count=len(cases),recorded_completed_count=len(completed),
        coverage_chainage_m=[records[0]['chainage_m'],records[-1]['chainage_m']],native_dtype=next(iter(dtypes)),
        domain_m=[float(shape[0]*spacing[0]),float(shape[1]*spacing[1])],native_record_ns=(ni-1)*dt*1e9,
        cell_count=int(np.prod(shape)),dt_s=dt,iterations=ni,source='impulse40A on first half-step; z Hertzian0.025m',
        geometry_sha256=sha(geometry),materials_sha256=sha(material_path),geometry_ranges=ranges,
        midpoint_agl_m=[min(r['geometry']['midpoint_agl_m'] for r in records),max(r['geometry']['midpoint_agl_m'] for r in records)],
        top_buffer_to_pml_m=[min(r['top_buffer_to_pml_m'] for r in records),max(r['top_buffer_to_pml_m'] for r in records)],
        closest_side_air_path_ns=[min(r['closest_side_pml_air_mirror_ns'] for r in records),max(r['closest_side_pml_air_mirror_ns'] for r in records)],
        packaged_truncated_sample_count=package_samples,packaged_taper_ns=int(np.ceil(package_samples*fraction))*dt*1e9,package_reproduction_relative_L2=error,
        independent_7_frequency_all_trace_DFT_L2=independent,independent_inverse_L2=inverse_errors,
        packaged_inverse_L2=inverse_packaged,taper_sensitivity_complex_L2=sensitivity,
        logged_seconds_median=float(np.median([e['seconds'] for e in events])),logged_seconds_sum=float(sum(e['seconds'] for e in events)),
        records=records)
    save(out/f'{path.name}_audit.json',stats)
    print(json.dumps({k:v for k,v in stats.items() if k!='records'},ensure_ascii=False),flush=True)
    return dict(stats=stats,records=records,profiles=profiles,packaged=old,frequency=FREQ,response=raw_response/.025,
        chainage=np.array([r['chainage_m'] for r in records]),raw_matrix=rx,native_time=np.arange(ni)*dt,
        geometry=data,spacing=spacing,path=path)


def plots(results,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    colors={'surface':'#e34a33','cover_base':'#00a6d6','first_sand':'#c25be5','basal_sand':'#1ec766'}
    labels={'surface':'地表','cover_base':'覆盖层底','first_sand':'首个砂岩顶','basal_sand':'连续底部砂岩顶'}
    for number,result in enumerate(results,1):
        stats=result['stats'];pos=result['chainage'];idx=np.argsort(pos)
        edges=np.r_[pos[idx][0]-.1,(pos[idx][:-1]+pos[idx][1:])/2,pos[idx][-1]+.1]
        fig,axes=plt.subplots(3,1,figsize=(13,12),layout='constrained')
        for ax,window in zip(axes,['gaussian','hann','rectangular']):
            profile=result['profiles'][window];keep=profile.time<=800e-9
            mag=abs(profile.complex_envelope[keep][:,idx]);reference=float(np.max(abs(profile.complex_envelope)))
            db=20*np.log10(np.maximum(mag/reference,1e-9))
            tt=profile.time[keep]*1e9;te=np.r_[tt[0]-(tt[1]-tt[0])/2,(tt[:-1]+tt[1:])/2,tt[-1]+(tt[-1]-tt[-2])/2]
            pic=ax.pcolormesh(edges,te,db,cmap='gray_r',vmin=-90,vmax=-30,shading='flat',rasterized=True)
            for key,color in colors.items():
                t=np.array([r['geometry'][key]['time95_ns'] if r['geometry'][key] else np.nan for r in result['records']])
                ax.plot(pos[idx],t[idx],color=color,lw=1,ls='--',label=labels[key]+'：H5分层/95MHz局部射线近似')
            ax.set(ylim=(800,0),ylabel='时间 / ns',title=f'{window}：未叠道、未AGC、未去背景、未尾渐消；本面板全记录峰值为0dB')
            ax.invert_xaxis();ax.legend(fontsize=7,loc='lower right')
            fig.colorbar(pic,ax=ax,label='复包络幅度 / 本面板共同参考（dB）')
        axes[-1].set_xlabel('测线里程 / m（从大到小为飞行方向；仅绘已完成区）')
        fig.suptitle(f"包{number}：{stats['raw_count']}/{stats['planned_count']}站，原生{stats['native_dtype']}，域{stats['domain_m']}m\n曲线来自实际体素，非固定7m+7.3m；不同频窗独立共同参考，不能凭亮度比较绝对增益")
        fig.savefig(out/f'package{number}_corrected_bscan.png',dpi=140);plt.close(fig)
    # Direct comparison uses common v3 stations, absolute shared magnitude scale, no AGC.
    high,low=results[1:]
    n=len(low['chainage'])
    np.testing.assert_allclose(high['chainage'][:n],low['chainage'],rtol=0,atol=1e-10)
    np.testing.assert_array_equal(high['geometry'],low['geometry'][:high['geometry'].shape[0]])
    fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    for column,(flo,fhi) in enumerate([(20e6,170e6),(20e6,40e6)]):
        products=[]
        for r in [high,low]:
            take=(FREQ>=flo)&(FREQ<=fhi);freq=FREQ[take];spectrum=r['response'][take,:n]
            weights=sf.spectral_window('gaussian',len(freq));padded=np.zeros((len(freq)*8,n),complex)
            padded[:len(freq)]=spectrum*weights[:,None]
            env=np.fft.ifft(padded,axis=0)*(len(freq)*8/len(freq))
            time=np.arange(len(env))/(len(env)*300000)
            products.append((env,time))
        reference=max(float(np.max(abs(v))) for v,t in products)
        for row,(r,(env,time)) in enumerate(zip([high,low],products)):
            pos=r['chainage'][:n];sort=np.argsort(pos);edges=np.r_[pos[sort][0]-.1,(pos[sort][:-1]+pos[sort][1:])/2,pos[sort][-1]+.1]
            keep=time<=600e-9;db=20*np.log10(np.maximum(abs(env[keep][:,sort])/reference,1e-10))
            tt=time[keep]*1e9;te=np.r_[tt[0]-(tt[1]-tt[0])/2,(tt[:-1]+tt[1:])/2,tt[-1]+(tt[-1]-tt[-2])/2]
            ax=axes[row,column];pic=ax.pcolormesh(edges,te,db,cmap='gray_r',vmin=-90,vmax=-30,shading='flat',rasterized=True)
            for key,color in colors.items():
                ts=np.array([v['geometry'][key]['time95_ns'] if v['geometry'][key] else np.nan for v in r['records'][:n]])
                ax.plot(pos[sort],ts[sort],color=color,lw=1,ls='--',label=labels[key])
            ax.set(ylim=(600,0),title=f"{'定高35m，实际AGL约10m' if row==0 else '地形跟随AGL1m'}；{flo/1e6:g}–{freq[-1]/1e6:g}MHz；未叠道/AGC/尾渐消",ylabel='时间 / ns',xlabel='测线里程 / m')
            ax.invert_xaxis();ax.legend(fontsize=7,loc='lower right');fig.colorbar(pic,ax=ax,label='dB；同列两图同一绝对参考')
    fig.suptitle('同一v3地质的70个共同站位：降低高度对照\n两域右边界相差10m，因此不是严格单因素；彩线为实际H5/95MHz近似，不代表已经检出')
    fig.savefig(out/'height_comparison_70stations.png',dpi=140);plt.close(fig)


def main(root,out):
    if out.exists(): raise ValueError('Fresh review output required')
    out.mkdir(parents=True)
    results=[package(p,out) for p in sorted(root.glob('line9_pkg*')) if p.is_dir()]
    if len(results)!=3: raise ValueError('Exactly3 input packages required')
    plots(results,out)
    summaries=[]
    for r in results:
        stats={k:v for k,v in r['stats'].items() if k!='records'}
        stats['gate_contrast_proxy_median_db']={key:float(np.median([v['timing_gates'][key]['local_guard_contrast_db'] for v in r['records'] if key in v['timing_gates']])) for key in ['cover_base','first_sand','basal_sand']}
        summaries.append(stats)
    save(out/'summary.json',dict(status='AUDITED_RECONSTRUCTION_NOT_INTERFACE_DETECTION',calls_solver=False,packages=summaries,
        input_inventory_sha256=sha(root/'input_inventory.json'),script_sha256=sha(__file__),official_processing_sha256=sha(sf.__file__),
        limits='FP32 only. Model-informed95MHz local horizontal-layer phase-ray curves are approximate, not full dispersive waveform arrival truth. Local gate contrast is not SNR/detection or physical validation. Height controls also change right domain boundary.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.root.resolve(),a.out.resolve())
