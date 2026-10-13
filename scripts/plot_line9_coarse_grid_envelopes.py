"""Separate phase-sensitive response differences from envelope changes; no filtering."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze_line9_v31_smoothing_domain import read_sf,relative
from study_line9_v31_smoothing_domain import load,save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['root','reference','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();contract=load(a.root/'contract.json')
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    rows=[];fig,axes=plt.subplots(2,3,figsize=(16,8),constrained_layout=True)
    for col,s in enumerate([150.,162.,174.]):
        ref=read_sf(a.reference/'cases'/f'full_x{s:.1f}_y'/'sfcw.h5')
        t=ref['t'];mask=(t>=100)&(t<400)
        axes[0,col].plot(t[mask],ref['real'][mask],color='black',lw=1.2,label='2.5cm 基准')
        axes[1,col].plot(t[mask],np.abs(ref['env'][mask]),color='black',lw=1.2,label='2.5cm 基准')
        for case in [v for v in contract['cases'] if v['station_m']==s]:
            sf=read_sf(a.root/'cases'/case['id']/'sfcw.h5')
            rows.append(dict(case=case['id'],complex_window_relative_l2=relative(sf['env'][mask],ref['env'][mask]),
                             envelope_only_window_relative_l2=relative(np.abs(sf['env'][mask]),np.abs(ref['env'][mask]))))
            label=f"{case['grid_m']*100:g}cm"
            axes[0,col].plot(t[mask],sf['real'][mask],lw=.85,label=label)
            axes[1,col].plot(t[mask],np.abs(sf['env'][mask]),lw=.85,label=label)
        axes[0,col].set(title=f'{s:g}m｜带符号总场',ylabel='官方实数带通重建幅值')
        axes[1,col].set(title=f'{s:g}m｜同一结果的复包络模',ylabel='官方复包络幅值',xlabel='时间（ns）')
    out=a.out/'envelope_diagnostic.json'
    if out.exists():assert load(out)==rows
    else:save(out,rows)
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('三站原幅值叠加｜无增益、归一化、滤波或时间对齐；包络相近不代表复相位一致')
    dest=a.out/'waveform_and_envelope.png'
    if dest.exists():raise FileExistsError(dest)
    fig.savefig(dest,dpi=170);plt.close(fig)


if __name__=='__main__':main()
