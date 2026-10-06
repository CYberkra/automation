"""Readable one-station display from already verified SFCW arrays; no signal changes."""
import argparse
from pathlib import Path
import numpy as np
from build_pdf_profile_geometry import digest, save_json


def render(arrays, out):
    if out.exists():
        raise ValueError('Fresh single-station display required')
    with np.load(arrays) as z:
        if z['native_Ez_V_m'].shape[1]!=1 or not np.array_equal(z['frequency_Hz'],20e6+np.arange(501)*.3e6):
            raise ValueError('Expected one verified20–170MHz/501point station')
        panels=[(z['native_time_s']*1e9,z['native_Ez_V_m'][:,0],'原始Ricker脉冲总场','Ez / (V/m)'),
                (z['sfcw_time_s']*1e9,z['rectangular_signed'][:,0],'SFCW矩形频窗总场','源归一化场响应（非端口S21）'),
                (z['sfcw_time_s']*1e9,z['hann_signed'][:,0],'SFCW Hann频窗总场','源归一化场响应（非端口S21）')]
        x=float(z['profile_X_m'][0]);s=float(z['acquisition_s_m'][0])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    out.mkdir(parents=True)
    common=max(float(np.max(abs(v[t<=1200]))) for t,v,_,_ in panels[1:])
    fig,axes=plt.subplots(1,3,figsize=(12,8),layout='constrained')
    for i,(t,v,label,unit) in enumerate(panels):
        mask=t<=1200;scale=float(np.max(abs(v[mask]))) if i==0 else common
        pic=axes[i].imshow(v[mask,None],extent=[s-.25,s+.25,t[mask][-1],t[mask][0]],
                           aspect='auto',cmap='gray',vmin=-scale,vmax=scale,interpolation='nearest')
        axes[i].set(title=label,ylabel='时间 / ns',xlabel=f'仅X{x:g}；s={s:g}m',xticks=[s])
        axes[i].axhline(100.16,color='red',ls='--',lw=.8)
        fig.colorbar(pic,ax=axes[i],label=unit)
    fig.suptitle('完整400×75m域；15m离地；2.5cm/FP64；只有一站，没有相邻道\n横向展开仅方便看单道，不能解释为连续B-scan；未去背景、未增益\n红线为平地约100ns气程参考；原始脉冲还包含源延迟')
    fig.savefig(out/'single_station_gray.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(15,10),layout='constrained')
    for row,(t,v,label,unit) in enumerate(panels):
        for col,(lo,hi) in enumerate([(0,150),(150,600),(600,1200)]):
            mask=(t>=lo)&(t<hi)
            axes[row,col].plot(t[mask],v[mask],lw=.8)
            axes[row,col].set(title=label+f'；{lo}–{hi}ns',xlabel='时间 / ns',ylabel=unit)
            axes[row,col].grid(alpha=.15)
    fig.suptitle('同一X220原数值A-scan；各窗独立纵轴便于看弱回波，无数值缩放或增益\nSFCW20–170MHz/501点；保留符号、相位与源归一化定义')
    fig.savefig(out/'single_station_windowed_ascans.png',dpi=140);plt.close(fig)
    save_json(out/'display_report.json',dict(status='READONLY_DISPLAY_FROM_VERIFIED_ARRAYS',calls_solver=False,
        arrays_sha256=digest(arrays),script_sha256=digest(Path(__file__)),signal_modified=False,
        scope='One station width expanded for visibility only; no interpolated traces. A-scan window axes are independent physical scales.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--arrays',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();render(a.arrays.resolve(),a.out.resolve())
