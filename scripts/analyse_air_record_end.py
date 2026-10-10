"""Existing impulse data only: test whether the late profile peak follows crop.

All crops are diagnostic, not a change to the frozen 500ns official direct path.
Same source, tones, phase, Hann display; no solver, tail taper or amplitude fit.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf


def analyse(public,private):
    source=sf.load_source(private/'full_air/model.h5')
    h0=sf.load_receiver(private/'full_air/model.h5','rxs/rx1','Ez')
    h1=sf.load_receiver(private/'full_cover/model.h5','rxs/rx1','Ez')
    data=h1.samples-h0.samples[:len(h1.samples)]
    f=20e6+.3e6*np.arange(501)
    rows=[];arrays={};profiles=[]
    for duration in [200,300,400,500]:
        n=min(4240,round(duration*1e-9/source.dt))
        actual=n*source.dt*1e9
        signal=sf.SampledSignal(path=f'diagnostic_crop_{duration}',samples=data[:n],dt=source.dt,time_offset=0)
        response=sf.direct_frequency_response(source,signal,f)
        profile=sf.reconstruct_time_response(response,window='hann',zero_pad_factor=32)
        t=profile.time*1e9;amp=abs(profile.complex_envelope)
        near=(t>=actual-20)&(t<=actual+20)
        indices=np.flatnonzero(near);index=indices[np.argmax(amp[near])]
        primary=(t>30)&(t<80)
        level=20*np.log10(amp[index]/amp[primary].max())
        rows.append(dict(crop_requested_ns=duration,actual_record_length_ns=actual,
                         endpoint_peak_time_ns=float(t[index]),endpoint_peak_relative_primary_db=float(level),
                         original_tail_relative_db=float(response.receiver_tail_relative_db)))
        profiles.append((t,amp,duration));arrays[f'crop{duration}_response']=response.response
        arrays[f'crop{duration}_envelope']=profile.complex_envelope
    arrays['display_time_ns']=t;arrays['frequency_hz']=f
    np.savez_compressed(public/'record_end_checks.npz',**arrays)
    (public/'record_end_checks.json').write_text(json.dumps(dict(rows=rows,
        scope='Diagnostic crops of existing same-source same-geometry histories, not new formal windows',
        claim_boundary='Peak tracking supports finite-record contamination, not a unique explanation of all ringing or Line9 events'),indent=2)+'\n',encoding='utf-8')
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(1,2,figsize=(15,5),constrained_layout=True)
    for t,amp,duration in profiles:
        for axis in ax:axis.plot(t,amp,label=f'同一原始记录，截取{duration}ns')
    ax[0].set(xlim=(30,80),title='地表主回波位置保持',xlabel='时间（ns）',ylabel='源归一化包络幅度（共同尺度）')
    ax[1].set(xlim=(100,550),title='晚峰是否跟随记录截断位置',xlabel='时间（ns）',ylabel='源归一化包络幅度（共同尺度）')
    for axis in ax:axis.grid(alpha=.2);axis.legend(fontsize=9)
    fig.suptitle('8米离地平地单站位：仅裁取已有冲激输出，未启动新正演\n官方direct、同一501频点、Hann显示；原正式诊断仍为500ns',fontsize=13)
    fig.savefig(public/'record_end_comparison.png',dpi=150)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',type=Path,required=True)
    p.add_argument('--private',type=Path,required=True)
    a=p.parse_args();analyse(a.public,a.private)
