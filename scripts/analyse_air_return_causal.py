"""Existing native histories: causal return followed by official direct SFCW.

Full-band numerical propagation is used only to reconstruct raw time histories.
Then crop to the original receiver's record before the unchanged 501-tone direct.
No fitted gain/delay, tail taper, hard angle cutoff, or new solver run.
"""
import argparse
import json
from pathlib import Path
import time

import h5py
import numpy as np
from scipy.fft import next_fast_len
from gprMax.toolboxes.SFCW import processing as sf


def plane(filename, label, count=4240):
    with h5py.File(filename) as f:
        gs={str(g.attrs['Name']):g for g in f['rxs'].values()}
        order=sorted((k for k in gs if k.startswith(label)),key=lambda k:int(k[len(label):]))
        return np.column_stack([gs[k]['Ez'][:count] for k in order])


def propagate(data,dt,dl,distance,ix,space_pad,time_pad):
    nt,nx=data.shape
    ns=next_fast_len(space_pad*nx)
    nfft=next_fast_len(time_pad*nt)
    freq=np.fft.rfftfreq(nfft,d=dt)
    q=2*np.pi*np.fft.fftfreq(ns,d=dl)
    qt=2*np.sin(q*dl/2)/dl
    wt=2*np.sin(np.pi*freq*dt)/dt
    spectrum=np.fft.rfft(data,n=nfft,axis=0)
    result=np.zeros(len(freq),dtype=np.complex128)
    xphase=np.exp(1j*q*ix*dl)
    # Blocked evaluation bounds memory even for padded time/space transforms.
    for start in range(0,len(freq),128):
        stop=min(start+128,len(freq))
        root=np.sqrt((wt[start:stop,None]/299792458.)**2-qt[None,:]**2+0j)
        ky=2/dl*np.arcsin(dl*root/2)
        ky=ky.real-1j*np.abs(ky.imag)
        kernel=np.exp(-1j*ky*distance)*xphase
        field=np.fft.fft(spectrum[start:stop],n=ns,axis=1)
        result[start:stop]=np.sum(field*kernel,axis=1)/ns
    return np.fft.irfft(result,n=nfft)[:nt]


def analyse(public,private):
    start=time.perf_counter()
    c=json.loads((public/'contract.json').read_text(encoding='utf-8'))
    source=sf.load_source(private/'full_air/model.h5')
    dl=c['mesh_m'];dt=source.dt;n=4240
    xs=np.arange(round(2.05/dl),round((c['width_m']-2.05)/dl)+1)*dl
    ix=int(np.argmin(abs(xs-c['original_rx_logical_m'][0])))
    distance=c['original_rx_logical_m'][1]-c['return_plane_y_m']
    full=plane(private/'full_cover/model.h5','return')-plane(private/'full_air/model.h5','return')
    compact=plane(private/'compact_cover/model.h5','return')-plane(private/'compact_air/model.h5','return')
    truth=sf.load_receiver(private/'full_cover/model.h5','rxs/rx1','Ez').samples[:n]-sf.load_receiver(private/'full_air/model.h5','rxs/rx1','Ez').samples[:n]
    f=20e6+.3e6*np.arange(501)
    sig=lambda v:sf.SampledSignal(path='causal_return',samples=v,dt=dt,time_offset=0)
    target=sf.direct_frequency_response(source,sig(truth),f).response
    arrays={'frequency_hz':f,'truth_scattered':truth,'dt_s':dt}
    rows=[]
    for sp,tp in [(4,2),(8,4),(16,8)]:
        for label,values in [('full',full),('compact',compact)]:
            raw=propagate(values,dt,dl,distance,ix,sp,tp)
            response=sf.direct_frequency_response(source,sig(raw),f).response
            low=f<=40e6;early=np.arange(n)*dt<80e-9
            rows.append(dict(plane=label,space_pad=sp,time_pad=tp,
                             fullband_relative_l2=float(np.linalg.norm(response-target)/np.linalg.norm(target)),
                             low20_40_relative_l2=float(np.linalg.norm((response-target)[low])/np.linalg.norm(target[low])),
                             time_early80ns_relative_l2=float(np.linalg.norm((raw-truth)[early])/np.linalg.norm(truth[early])),
                             time_full_relative_l2=float(np.linalg.norm(raw-truth)/np.linalg.norm(truth))))
            arrays[f'{label}_space{sp}_time{tp}_native']=raw
            arrays[f'{label}_space{sp}_time{tp}_direct']=response
    (public/'causal_return_checks.json').write_text(json.dumps(dict(rows=rows,analysis_wall_s=time.perf_counter()-start,scope='Existing raw histories only; upgoing, finite-aperture numerical transport remains diagnostic'),indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(public/'causal_return_checks.npz',**arrays)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',type=Path,required=True)
    p.add_argument('--private',type=Path,required=True)
    a=p.parse_args();analyse(a.public,a.private)
