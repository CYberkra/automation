"""Source-normalized native and dense wavefield diagnostics in the actual band."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import h5py
import numpy as np

from hs_capsule_identity import sha256
from hs4_height_wavefield_v0_2 import receivers, save, REFERENCE
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime
from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
    direct_frequency_response, engineering_dft, apply_tail_taper,
    reconstruct_time_response, spectral_window)


def relative(a,b):
    return float(np.linalg.norm(b-a)/np.linalg.norm(a)) if np.linalg.norm(a)>0 else None


def response(p):
    s=load_source(p); r=load_receiver(p,receiver_path='/rxs/rx1',component='Ey')
    fraction=(round(200e-9/r.dt)-.25)/len(r.samples)
    result=direct_frequency_response(s,r,FREQ,tail_taper_fraction=fraction)
    if not result.source_valid.all():
        raise ValueError('source contains invalid band bins')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists(): raise ValueError('new result directory required')
    verify_official_runtime()
    c=json.loads((args.study/'execution_contract.json').read_text('utf-8'))
    verification=json.loads((args.study/'independent_verification.json').read_text('utf-8'))
    if verification['status']!='PASS' or verification['contract_sha256']!=sha256(args.study/'execution_contract.json'):
        raise ValueError('independently verified wavefields required')
    args.out.mkdir(parents=True)
    fields_folder=args.out/'fields'; fields_folder.mkdir()
    dt=c['dt_s']; n=c['iterations']; iterations=np.asarray(c['snapshot_iterations'])
    taper=apply_tail_taper(np.ones(n),(round(200e-9/dt)-.25)/n)
    weights=spectral_window('hann',len(FREQ))
    results={}; native={}; field_files={}; profiles={}
    geometry=np.genfromtxt(REFERENCE.parents[0]/'2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',delimiter=',',names=True)
    x=c['snapshot_extent_m'][0]+.075+.15*np.arange(c['snapshot_shape_xyz'][0])
    z=c['snapshot_extent_m'][2]+.075+.15*np.arange(c['snapshot_shape_xyz'][2])
    interface=geometry['grid_z_m'][np.clip(np.floor((x-12)/.25).astype(int),0,47)]
    frame_times=np.arange(0,241,10)*1e-9
    synthesis=np.exp(2j*np.pi*frame_times[:,None]*FREQ[None,:])*weights[None,:]/len(FREQ)
    frame_data={}
    for g in c['groups']:
        directory=Path(g['input']).parent; name=g['id']
        raw=directory/'profile.h5'; r=response(raw); native[name]=r
        source_spectrum=engineering_dft(r.source.samples,r.source.dt,FREQ,time_offset=r.source.time_offset)
        with h5py.File(raw) as h:
            rx=receivers(h); probe_native=[]
            for probe in g['probes']:
                k=probe['id']
                probe_native.extend([rx[k+'c']['Ey'][:],
                    .5*(rx[k+'c']['Hx'][:]+rx[k+'z']['Hx'][:]),
                    .5*(rx[k+'c']['Hz'][:]+rx[k+'x']['Hz'][:])])
        probe_native=np.stack(probe_native,axis=1)*taper[:,None]
        full=engineering_dft(probe_native,dt,FREQ)
        decimation={}
        for hop in (5,10):
            sampled=engineering_dft(probe_native[::hop],dt*hop,FREQ)
            per_channel=np.linalg.norm((sampled-full)*weights[:,None],axis=0)/np.linalg.norm(full*weights[:,None],axis=0)
            decimation[str(hop)]={'hann_complex_relative_L2':relative(full*weights[:,None],sampled*weights[:,None]),
                'maximum_per_probe_component_relative_L2':float(per_channel.max()),'per_probe_component_relative_L2':per_channel.tolist()}
        chosen=np.arange(0,501,31)
        direct=dt*np.exp(-2j*np.pi*FREQ[chosen,None]*dt*np.arange(n)[None,:])@probe_native
        direct_error=relative(direct,full[chosen])
        results[name]={'raw_sha256':sha256(raw),'source_all501_valid':True,'decimation':decimation,
            'native_czt_vs_selected_direct_DFT_relative_L2':direct_error,'H_time_offset_s':-.5*dt}
        histories={f:np.empty((len(iterations),len(x)*len(z)),dtype=np.float64) for f in ('Ey','Hx','Hz')}
        for k,j in enumerate(iterations):
            p=directory/'profile_snaps'/f'snap{j:05d}.h5'
            with h5py.File(p) as h:
                for f in histories: histories[f][k]=h[f][:,0,:].reshape(-1)
        files={}; frames={}
        for f,history in histories.items():
            history*=taper[iterations,None]
            target=fields_folder/(name+'_'+f+'.npy')
            transfer=np.lib.format.open_memmap(target,mode='w+',dtype=np.complex128,shape=(501,len(x)*len(z)))
            for start in range(0,transfer.shape[1],512):
                stop=min(start+512,transfer.shape[1])
                transfer[:,start:stop]=engineering_dft(history[:,start:stop],dt*10,FREQ,
                    time_offset=0 if f=='Ey' else -.5*dt)/source_spectrum[:,None]
            transfer.flush()
            # Full-native interpolated probes independently close selected cells.
            errors=[]
            for k,probe in enumerate(g['probes']):
                channel=k*3+('Ey','Hx','Hz').index(f)
                offset_phase=1 if f=='Ey' else np.exp(2j*np.pi*FREQ*dt*.5)
                expected=full[:,channel]*offset_phase/source_spectrum
                cell=probe['snapshot_ix']*len(z)+probe['snapshot_iz']
                errors.append(relative(expected*weights,transfer[:,cell]*weights))
            results[name][f+'_band_probe_closure_max_relative_L2']=max(errors)
            frames[f]=(synthesis@transfer).reshape(len(frame_times),len(x),len(z))
            if f=='Ey':
                ix=np.arange(len(x)); iz=np.clip(np.rint((interface+.15-z[0])/.15).astype(int),0,len(z)-1)
                profiles[name]=np.sqrt(np.mean(abs(weights[:,None]*transfer[:,ix*len(z)+iz])**2,axis=0))
            files[f]={'file':target.relative_to(args.out).as_posix(),'sha256':sha256(target),'shape':list(transfer.shape),'dtype':str(transfer.dtype)}
            del transfer
        del histories
        field_files[name]=files; frame_data[name]=frames
        print('Analyzed',name,flush=True)
    reference={role:response(REFERENCE/('centre_'+role)/'profile.h5') for role in ('rough','halfspace')}
    source_equivalence={}
    for role in reference:
        a,b=reference[role].response,native['high_'+role].response
        source_equivalence[role]={'hann_complex_relative_L2':relative(weights*a,weights*b)}
    ref_delta=reference['rough'].response-reference['halfspace'].response
    high=native['high_rough'].response-native['high_halfspace'].response
    low=native['low_rough'].response-native['low_halfspace'].response
    source_equivalence['paired_difference']={'hann_complex_relative_L2':relative(weights*ref_delta,weights*high)}
    product_high=reconstruct_time_response(replace(native['high_rough'],response=high),window='hann',zero_pad_factor=8)
    product_low=reconstruct_time_response(replace(native['low_rough'],response=low),window='hann',zero_pad_factor=8)
    product_ref=reconstruct_time_response(replace(reference['rough'],response=ref_delta),window='hann',zero_pad_factor=8)
    mask=(product_high.time>=160e-9)&(product_high.time<=220e-9)
    source_equivalence['paired_difference']['signed160_220ns_relative_L2']=relative(product_ref.real_bandpass[mask],product_high.real_bandpass[mask])
    arrays={'x_m':x,'z_m':z,'interface_z_m':interface,'frame_times_ns':frame_times*1e9,
        'time_ns':product_high.time*1e9,'receiver_high_signed':product_high.real_bandpass,
        'receiver_low_signed':product_low.real_bandpass,'receiver_high_complex':product_high.complex_envelope,
        'receiver_low_complex':product_low.complex_envelope,'frequency_Hz':FREQ,'receiver_high_spectrum':high,'receiver_low_spectrum':low}
    for name,frame in frame_data.items():
        for f,v in frame.items(): arrays[name+'_'+f+'_complex_band_frames']=v
        arrays[name+'_above_interface_spectral_RMS']=profiles[name]
    np.savez_compressed(fields_folder/'band_frames.npz',**arrays)
    # Independent selected-time direct synthesis versus the official IFFT.
    checked_times=product_high.time[np.arange(0,240,13)]
    synthesis_check=(np.exp(2j*np.pi*checked_times[:,None]*FREQ[None,:])*weights[None,:]/len(FREQ))@high
    official=product_high.complex_bandpass[np.arange(0,240,13)]
    synth_error=relative(official,synthesis_check)
    summary={'status':'COMPLETED_DIAGNOSTIC_NOT_UNIQUE_INTERFACE_ATTRIBUTION','code_sha256':sha256(__file__),
        'contract_sha256':sha256(args.study/'execution_contract.json'),'independent_verification_sha256':sha256(args.study/'independent_verification.json'),
        'source_equivalence_to_same_grid_impulse':source_equivalence,'field_sampling_checks':results,'field_files':field_files,
        'direct_field_synthesis_vs_official_IFFT_relative_L2':synth_error,
        'coordinate_scope':'ROI X13.5-22.5 Z8-28.1m sampled at .15m centres; native FDTD .025m. No outer-PML observation. Invariant Y line source.',
        'processing':'501 exact20-170MHz tones; saved full-native source spectrum; sampled native200ns taper; H(-dt/2) included in DTFT; shared normalized Hann/8x; complex pair difference before magnitude.',
        'scope':'Native probe alias checks cover12 points, not every output cell. Source-normalized band fields have sidelobes; no pre-ringing interpreted as causal arrival. Field difference and direction diagnostic are not physical energy percentages.'}
    save(args.out/'summary.json',summary)
    plot(args.out,arrays)
    print('Source equivalence:',source_equivalence,flush=True)


def plot(out,a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    x,z=a['x_m'],a['z_m']
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    t=a['time_ns']; air_delay=2*13/299792458*1e9
    axes[0].plot(t,a['receiver_high_signed'],label='15m')
    axes[0].plot(t+air_delay,a['receiver_low_signed'],label='2m（仅诊断：增加固定空气双程时延）')
    axes[0].set_xlim(145,225); axes[0].set_ylabel('源归一化带符号响应'); axes[0].legend()
    for name,label in [('high_halfspace','15m匹配全覆盖层'),('low_halfspace','2m匹配全覆盖层')]:
        axes[1].plot(x-12,a[name+'_above_interface_spectral_RMS'],label=label)
    axes[1].set_xlabel('原模型X（m）'); axes[1].set_ylabel('界面上方入射场带内RMS代理'); axes[1].legend()
    axes[0].set_title('中心配对响应：20–170MHz / 相同2.5cm网格 / 不逐道归一化')
    axes[1].set_title('全覆盖层场在各界面附近采样：不是独立反射贡献或功率')
    fig.savefig(out/'height_receiver_and_incidence.png',dpi=150); plt.close(fig)
    chosen=[6,10,14,18]
    delta={h:a[h+'_rough_Ey_complex_band_frames']-a[h+'_halfspace_Ey_complex_band_frames'] for h in ('high','low')}
    limit=max(float(np.max(abs(2*v[chosen].real))) for v in delta.values())
    fig,axes=plt.subplots(2,4,figsize=(14,7),constrained_layout=True)
    for row,(height,name) in enumerate([(15,'high'),(2,'low')]):
        for col,k in enumerate(chosen):
            ax=axes[row,col]; image=ax.pcolormesh(x-12,z,2*delta[name][k].real.T,cmap='RdBu_r',vmin=-limit,vmax=limit,shading='nearest')
            ax.plot(x-12,a['interface_z_m'],'k-',lw=.8); ax.axhline(12,color='k',ls='--',lw=.7)
            ax.scatter([5.6,6.9],[12+height]*2,c=['red','black'],s=10)
            ax.set_title(f'{height}m / {a["frame_times_ns"][k]:.0f}ns'); ax.set_ylim(8,28.1); ax.set_aspect('equal')
            ax.set_xlabel('原模型X（m）'); ax.set_ylabel('Z（m）')
    fig.colorbar(image,ax=axes.ravel().tolist(),label='配对差分带符号场/源，共享线性色标')
    fig.suptitle('同频带空间响应：模型起伏不等于B-scan的垂直投影\n带限波列有前后振铃；此图不能单独签认事件来源')
    fig.savefig(out/'height_band_wavefield_atlas.png',dpi=160); plt.close(fig)


if __name__=='__main__': main()
