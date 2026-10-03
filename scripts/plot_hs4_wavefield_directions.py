"""Shared-scale field movies, flat-depth incident profiles and direct DFT checks."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs4_height_wavefield_v0_2 import receivers, save
from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import FREQ
from gprMax.toolboxes.SFCW.processing import apply_tail_taper, spectral_window


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--study',type=Path,required=True)
    p.add_argument('--results',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new visualization required')
    a.out.mkdir(parents=True)
    c=json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    s=json.loads((a.results/'summary.json').read_text('utf-8'))
    arrays=np.load(a.results/'fields/band_frames.npz'); x,z=arrays['x_m'],arrays['z_m']
    weights=spectral_window('hann',501); dt=c['dt_s']; n=c['iterations']
    taper=apply_tail_taper(np.ones(n),(round(200e-9/dt)-.25)/n)
    def oracle(values,hop):
        times=np.arange(0,n,hop)*dt
        result=np.empty((501,values.shape[1]),np.complex128)
        for start in range(0,501,32):
            stop=min(start+32,501)
            result[start:stop]=dt*hop*np.exp(-2j*np.pi*FREQ[start:stop,None]*times[None,:])@values[::hop]
        return result
    direct={}
    for g in c['groups']:
        with h5py.File(Path(g['input']).with_suffix('.h5')) as h:
            r=receivers(h); v=[]
            for probe in g['probes']:
                k=probe['id']; v.extend([r[k+'c']['Ey'][:],.5*(r[k+'c']['Hx'][:]+r[k+'z']['Hx'][:]),.5*(r[k+'c']['Hz'][:]+r[k+'x']['Hz'][:])])
        values=np.stack(v,axis=1)*taper[:,None]
        native=oracle(values,1)
        direct[g['id']]={}
        for hop in (5,10):
            sampled=oracle(values,hop)
            error=np.linalg.norm(weights[:,None]*(sampled-native),axis=0)/np.linalg.norm(weights[:,None]*native,axis=0)
            direct[g['id']][str(hop)]={'maximum_probe_component_relative_L2':float(error.max()),'all_probe_component_relative_L2':error.tolist()}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.animation import FuncAnimation, PillowWriter
    plt.rcParams['font.family']='Microsoft YaHei'
    incident={}; angles={}
    iz=int(np.argmin(abs(z-10.025)))
    for height in ('high','low'):
        name=height+'_halfspace'
        f={k:np.load(a.results/s['field_files'][name][k]['file'],mmap_mode='r')[:,np.arange(len(x))*len(z)+iz] for k in ('Ey','Hx','Hz')}
        incident[height]=np.sqrt(np.mean(abs(weights[:,None]*f['Ey'])**2,axis=0))
        sx=np.mean(weights[:,None]**2*np.real(f['Ey']*f['Hz'].conj()),axis=0)
        sz=-np.mean(weights[:,None]**2*np.real(f['Ey']*f['Hx'].conj()),axis=0)
        angles[height]=np.rad2deg(np.arctan2(sx,-sz))
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for height,label in [('high','15m'),('low','2m')]:
        axes[0].plot(x-12,incident[height],label=label)
        axes[1].plot(x-12,angles[height],label=label)
    axes[0].set_ylabel('源归一化Ey带内RMS代理'); axes[1].set_ylabel('向下场方向代理（度）')
    for ax in axes:
        ax.set_xlabel('原模型X（m）'); ax.legend(); ax.axvspan(5,6.5,color='gray',alpha=.1)
    fig.suptitle(f'匹配全覆盖层在固定Z={z[iz]:.3f}m上的入射场\n固定深度比较排除界面深度采样变化；方向代理不是能量归因')
    fig.savefig(a.out/'fixed_depth_incidence.png',dpi=150); plt.close(fig)
    delta={height:arrays[height+'_rough_Ey_complex_band_frames']-arrays[height+'_halfspace_Ey_complex_band_frames'] for height in ('high','low')}
    limit=max(float(np.max(abs(2*v.real))) for v in delta.values())
    norm=SymLogNorm(linthresh=limit*1e-4,vmin=-limit,vmax=limit)
    fig,axes=plt.subplots(1,2,figsize=(8,7),constrained_layout=True)
    images=[]
    for ax,(height,h) in zip(axes,[('high',15),('low',2)]):
        images.append(ax.pcolormesh(x-12,z,2*delta[height][0].real.T,cmap='RdBu_r',norm=norm,shading='nearest'))
        ax.plot(x-12,arrays['interface_z_m'],'k-',lw=1); ax.axhline(12,color='k',ls='--')
        ax.scatter([5.6,6.9],[12+h]*2,c=['red','black'],s=18); ax.set_aspect('equal')
        ax.set_xlabel('原模型X（m）'); ax.set_ylabel('Z（m）')
    fig.colorbar(images[0],ax=axes.tolist(),label='带符号差分场/源，共享对数色标')
    def update(k):
        for ax,img,(height,h) in zip(axes,images,[('high',15),('low',2)]):
            img.set_array((2*delta[height][k].real.T).ravel()); ax.set_title(f'{h}m / {arrays["frame_times_ns"][k]:.0f}ns')
        return images
    fig.suptitle('20–170MHz配对差分波场\n有前后振铃；不按帧归一化，不是因果波前证明')
    animation=FuncAnimation(fig,update,frames=len(arrays['frame_times_ns']),interval=150,blit=False)
    animation.save(a.out/'band_difference.gif',writer=PillowWriter(fps=6),dpi=100)
    update(18); fig.savefig(a.out/'band_difference_180ns.png',dpi=150); plt.close(fig)
    # Causal actual Ricker histories: same physical times and one global scale.
    times=np.arange(0,241,6); frames=[]
    native_j=np.asarray(c['snapshot_iterations'])
    for height in ('high','low'):
        matrix=[]
        for t in times:
            j=int(native_j[np.argmin(abs(native_j*dt*1e9-t))])
            with (h5py.File(a.study/(height+'_rough')/'profile_snaps'/f'snap{j:05d}.h5') as r,
                  h5py.File(a.study/(height+'_halfspace')/'profile_snaps'/f'snap{j:05d}.h5') as b):
                matrix.append(r['Ey'][:,0,:]-b['Ey'][:,0,:])
        frames.append(np.asarray(matrix))
    limit_raw=max(float(np.max(abs(v))) for v in frames)
    norm_raw=SymLogNorm(linthresh=limit_raw*1e-4,vmin=-limit_raw,vmax=limit_raw)
    fig,axes=plt.subplots(1,2,figsize=(8,7),constrained_layout=True); images=[]
    for ax,h in zip(axes,[15,2]):
        images.append(ax.pcolormesh(x-12,z,frames[len(images)][0].T,cmap='RdBu_r',norm=norm_raw,shading='nearest'))
        ax.plot(x-12,arrays['interface_z_m'],'k-',lw=1); ax.axhline(12,color='k',ls='--')
        ax.scatter([5.6,6.9],[12+h]*2,c=['red','black'],s=18); ax.set_aspect('equal'); ax.set_xlabel('原模型X（m）'); ax.set_ylabel('Z（m）')
    fig.colorbar(images[0],ax=axes.tolist(),label='实际Ricker差分Ey（V/m），共享对数色标')
    def update_raw(k):
        for row,(ax,img,h) in enumerate(zip(axes,images,[15,2])):
            img.set_array(frames[row][k].T.ravel()); actual=native_j[np.argmin(abs(native_j*dt*1e9-times[k]))]*dt*1e9
            ax.set_title(f'{h}m / {actual:.2f}ns')
        return images
    fig.suptitle('实际Ricker因果传播记录：起伏减匹配全覆盖层\n观察窗口不包含外侧PML；差分不是独立能量份额')
    FuncAnimation(fig,update_raw,frames=len(times),interval=120,blit=False).save(a.out/'raw_ricker_difference.gif',writer=PillowWriter(fps=8),dpi=100)
    plt.close(fig)
    save(a.out/'summary.json',{'status':'COMPLETED_DIAGNOSTIC','code_sha256':sha256(__file__),
        'input_summary_sha256':sha256(a.results/'summary.json'),'direct_DFT_sampling_checks':direct,
        'fixed_depth_z_m':float(z[iz]),'x_m':x.tolist(),'incident_RMS':{k:v.tolist() for k,v in incident.items()},
        'incident_direction_degrees':{k:v.tolist() for k,v in angles.items()},
        'band_shared_limit':limit,'raw_shared_limit_V_per_m':limit_raw,'log_linear_threshold_fraction':1e-4,
        'scope':'Whole-array scale shared by heights/frames. Direct DTFT validates12probe points, not entire-grid field accuracy. Allcover incident field on one flat plane shows spatial weight/angle differences, not identified reflected contributions.'})
    print('Direct DFT oracle and paired movies completed',flush=True)


if __name__=='__main__': main()
