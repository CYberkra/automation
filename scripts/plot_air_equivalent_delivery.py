"""Chinese single-station delivery figure from completed diagnostic evidence."""
import argparse
import json
from pathlib import Path
from dataclasses import replace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf


def plot(public,private):
    report=json.loads((public/'results.json').read_text(encoding='utf-8'))
    d=np.load(public/'continuum_flat_checks.npz')
    f=d['frequency_hz']
    source=sf.load_source(private/'full_air/model.h5')
    rx=sf.load_receiver(private/'full_cover/model.h5','rxs/rx1','Ez')
    rx=replace(rx,samples=rx.samples[:4240])
    base=sf.direct_frequency_response(source,rx,f)
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(2,2,figsize=(16,10),constrained_layout=True)
    colours=['#2864a7','#d35d35','#3e8662']
    series=[('full_native_500ns_scattered','完整域原接收位置：地表差场'),
            ('compact_plane_EH_return','缩域等效源＋电磁场回传'),
            ('continuum_scattered','无限平地解析参照（含近场）')]
    for (key,label),colour in zip(series,colours):
        fr=replace(base,response=d[key],receiver_spectrum=d[key]*base.source_spectrum)
        profile=sf.reconstruct_time_response(fr,window='hann',zero_pad_factor=8)
        ax[0,0].plot(profile.time*1e9,abs(profile.complex_envelope),label=label,color=colour)
        ax[1,0].plot(f/1e6,abs(d[key]),label=label,color=colour)
    arrival=np.sqrt(16**2+1.3**2)/299792458.*1e9
    ax[0,0].axvline(arrival,color='gray',linestyle=':',label=f'几何地表双程约{arrival:.1f}ns')
    ax[0,0].set(xlim=(0,650),title='单站位地表差场：同一Hann重建包络，补零8',xlabel='从原始冲激起算的时间（ns）',ylabel='源归一化幅度（共同尺度）')
    ax[0,0].legend(fontsize=9);ax[0,0].grid(alpha=.2)
    ax[1,0].set(title='地表差场复谱幅度：有限记录尚未通过收敛',xlabel='频率（MHz）',ylabel='|Ez / 原冲激电流|（V/m/A）')
    ax[1,0].legend(fontsize=9);ax[1,0].grid(alpha=.2)
    reference=d['continuum_scattered']
    for (key,label),colour in zip(series[:2],colours[:2]):
        ax[0,1].plot(f/1e6,np.angle(d[key]*reference.conj(),deg=True),label=label,color=colour)
    ax[0,1].set(title='相对无限平地解析参照的相位差（无拟合）',xlabel='频率（MHz）',ylabel='相位差（度）')
    ax[0,1].legend(fontsize=9);ax[0,1].grid(alpha=.2)
    records={r['case']:r for r in report['runs']}
    names=['full_air','compact_air','full_cover','compact_cover']
    labels=['完整空气域','缩域空气重放','完整平地','缩域平地重放']
    solve=np.array([records[k]['native_solve_seconds'] for k in names])
    total=np.array([records[k]['wall_seconds'] for k in names])
    x=np.arange(4)
    ax[1,1].bar(x,solve,label='官方求解计时',color='#4f7396')
    ax[1,1].bar(x,total-solve,bottom=solve,label='准备／启动／写出等',color='#d8a05a')
    for i,t in enumerate(total):ax[1,1].text(i,t+1,f'{t:.1f}s',ha='center')
    ax[1,1].set(xticks=x,xticklabels=labels,title='串行实测耗时：缩域当前仍较慢',ylabel='秒',ylim=(0,total.max()*1.2))
    ax[1,1].legend(fontsize=9);ax[1,1].grid(axis='y',alpha=.2)
    fig.suptitle('二维平地省算验收：8米离地、收发距1.3米、官方impulse与direct SFCW\n单站位A-scan；不是移动测线B-scan，实际地形与深层界面未计算',fontsize=15)
    fig.savefig(public/'delivery_summary.png',dpi=150)
    (public/'figure_settings.json').write_text(json.dumps(dict(
        window='official Hann for display only',zero_pad_factor=8,time_shift_s=0,
        metrics='Unwindowed exact 501-tone complex spectra',scope='Single stationary Tx/Rx flat half-space, no measured line'),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',type=Path,required=True)
    p.add_argument('--private',type=Path,required=True)
    a=p.parse_args();plot(a.public,a.private)
