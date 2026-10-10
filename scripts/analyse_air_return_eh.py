"""E/H return diagnostic on native 2D Yee fields, then official direct SFCW.

Inverse Laplace contour displacement evaluates the causal discrete propagator
away from real-axis grazing singularities. It is undone in the reconstructed
raw time history, not a taper of the official SFCW input. Padding/contour
convergence is reported; no production acceptance or fitted gain/delay.
"""
import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.fft import next_fast_len
from gprMax.toolboxes.SFCW import processing as sf


def fields(filename,count):
    with h5py.File(filename) as h:
        groups={str(g.attrs['Name']):g for g in h['rxs'].values()}
        keys=sorted((k for k in groups if k.startswith('return')),key=lambda k:int(k[6:]))
        return (np.column_stack([groups[k]['Ez'][:count] for k in keys]),
                np.column_stack([groups[k]['Hx'][:count] for k in keys]))


def propagate(e,h,dt,dl,distance,ix,sp,tp,period_decay):
    nt,nx=e.shape
    ns=next_fast_len(nx*sp);nf=next_fast_len(nt*tp)
    freq=np.fft.rfftfreq(nf,d=dt)
    sigma=period_decay/(nf*dt)
    damping=np.exp(-sigma*np.arange(nt)*dt)
    es=np.fft.rfft(e*damping[:,None],n=nf,axis=0)
    hs=np.fft.rfft(h*damping[:,None],n=nf,axis=0)
    omega=2*np.pi*freq-1j*sigma
    wt=2*np.sin(omega*dt/2)/dt
    q=2*np.pi*np.fft.fftfreq(ns,d=dl)
    qt=2*np.sin(q*dl/2)/dl
    xphase=np.exp(1j*q*ix*dl)
    out=np.zeros(len(freq),complex)
    mu=4*np.pi*1e-7
    for start in range(0,len(freq),64):
        stop=min(start+64,len(freq))
        root=np.sqrt((wt[start:stop,None]/299792458.)**2-qt[None,:]**2+0j)
        ky=2/dl*np.arcsin(dl*root/2)
        ky=ky.real-1j*np.abs(ky.imag)
        qy=2*np.sin(ky*dl/2)/dl
        fe=np.fft.fft(es[start:stop],n=ns,axis=1)
        # Hx sample n is at (n-1/2)dt and y_plane+dy/2.
        fh=np.fft.fft(hs[start:stop],n=ns,axis=1)*np.exp(1j*omega[start:stop,None]*dt/2)
        denominator=2*np.cos(ky*dl/2)
        up=(fe*np.exp(1j*ky*dl/2)+mu*wt[start:stop,None]/qy*fh)/denominator
        out[start:stop]=np.sum(up*np.exp(-1j*ky*distance)*xphase,axis=1)/ns
    raw=np.fft.irfft(out,n=nf)[:nt]/damping
    assert np.isfinite(raw).all()
    return raw


def analyse(public,private):
    start=time.perf_counter()
    c=json.loads((public/'contract.json').read_text(encoding='utf-8'))
    source=sf.load_source(private/'full_air/model.h5')
    dt=source.dt;dl=c['mesh_m'];n=4240
    xs=np.arange(round(2.05/dl),round((c['width_m']-2.05)/dl)+1)*dl
    ix=int(np.argmin(abs(xs-c['original_rx_logical_m'][0])))
    height=c['original_rx_logical_m'][1]-c['return_plane_y_m']
    values={name:fields(private/name/'model.h5',n) for name in ['full_air','full_cover','compact_air','compact_cover']}
    truth=sf.load_receiver(private/'full_cover/model.h5','rxs/rx1','Ez').samples[:n]-sf.load_receiver(private/'full_air/model.h5','rxs/rx1','Ez').samples[:n]
    f=20e6+.3e6*np.arange(501)
    sig=lambda v:sf.SampledSignal(path='EH_return',samples=v,dt=dt,time_offset=0)
    target=sf.direct_frequency_response(source,sig(truth),f).response
    rows=[];arrays=dict(frequency_hz=f,dt_s=dt,truth_native=truth,truth_direct=target)
    for sp,tp,decay in [(4,4,8),(8,8,12),(16,16,16)]:
        for label in ['full','compact']:
            e=values[label+'_cover'][0]-values[label+'_air'][0]
            h=values[label+'_cover'][1]-values[label+'_air'][1]
            raw=propagate(e,h,dt,dl,height,ix,sp,tp,decay)
            response=sf.direct_frequency_response(source,sig(raw),f).response
            early=np.arange(n)*dt<80e-9;low=f<=40e6
            rows.append(dict(plane=label,space_pad=sp,time_pad=tp,period_decay=decay,
                             fullband_relative_l2=float(np.linalg.norm(response-target)/np.linalg.norm(target)),
                             low20_40_relative_l2=float(np.linalg.norm((response-target)[low])/np.linalg.norm(target[low])),
                             time_early80ns_relative_l2=float(np.linalg.norm((raw-truth)[early])/np.linalg.norm(truth[early])),
                             time_full_relative_l2=float(np.linalg.norm(raw-truth)/np.linalg.norm(truth))))
            arrays[f'{label}_space{sp}_native']=raw
            arrays[f'{label}_space{sp}_direct']=response
    wall=time.perf_counter()-start
    (public/'eh_return_checks.json').write_text(json.dumps(dict(rows=rows,analysis_wall_s=wall,
        description='Exact-staggered E/H outgoing separation and causal discrete propagation; internal Laplace weight undone before official direct',
        incomplete=['Actual terrain/deep interface','Aperture convergence beyond current domain','Time-window convergence','Production speedup']),indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(public/'eh_return_checks.npz',**arrays)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',type=Path,required=True)
    p.add_argument('--private',type=Path,required=True)
    a=p.parse_args();analyse(a.public,a.private)
