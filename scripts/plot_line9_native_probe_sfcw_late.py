"""Shared physical-scale late-window single-station observer acceptance panels."""
import argparse,json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    s=json.loads((a.results/'analysis.json').read_text('utf-8'))
    assert sha(a.arrays)==s['arrays_sha256'] and s['native_source_rx_bitwise_equal'] and s['SFCW_complex_response_equal']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    with np.load(a.arrays) as h:
        nt=h['time_s']*1e9;rx=h['rx'];t=h['sfcw_time_s']*1e9;z={w:h[w] for w in ['hann','blackman']}
    fig,axs=plt.subplots(2,2,figsize=(13,10),layout='constrained');peaks={}
    mask=(nt>=300)&(nt<=370)
    axs[0,0].plot(nt[mask],rx[mask],label='原主Rx',lw=1.8,color='black');axs[0,0].plot(nt[mask],rx[mask],'--',label='新增探针后主Rx',lw=1,color='tab:orange')
    axs[0,0].set(xlabel='原生Ricker时间 / ns',ylabel='Ez / (V/m)',title='原生晚波逐位一致；不拟合时延或振幅');axs[0,0].legend()
    mask=(t>=300)&(t<=450)
    for ax,w in zip([axs[0,1],axs[1,0]],['hann','blackman']):
        real=z[w][mask].real;data=np.column_stack([real,real,np.zeros_like(real)]);lim=float(abs(data).max())
        im=ax.imshow(data,cmap='gray',vmin=-lim,vmax=lim,extent=[-.5,2.5,450,300],aspect='auto')
        ax.set_xticks(range(3),['原H0','探针H0','新−原']);ax.set(ylabel='SFCW时间 / ns',title=w+'：晚窗独立物理色标，三配置共用')
        fig.colorbar(im,ax=ax,label='(V/m)/(A·m)')
        j=np.flatnonzero(mask)[np.argmax(abs(z[w][mask]))]
        peaks[w]={'late_envelope_peak_ns':float(t[j]),'late_envelope_peak_amplitude':float(abs(z[w][j])),'signed_grayscale_limit':lim}
        axs[1,1].plot(t[mask],z[w][mask].real,lw=1,label=w)
    axs[1,1].set(xlabel='SFCW时间 / ns',ylabel='双极 / (V/m)/(A·m)',title='两种固定频窗的晚波原幅值');axs[1,1].legend()
    fig.suptitle('190m / ~8m AGL / 高损耗H0：原模型与原接收器不变，仅新增观察器\n20–170MHz、0.3MHz、精确501频点；配置列不是移动天线B-scan；原生和SFCW时标分开')
    out=a.results/'native_observer_sfcw_late_acceptance.png';assert not out.exists();fig.savefig(out,dpi=145);plt.close(fig)
    result={'status':'PASS_OBSERVER_LATE_WINDOW_VISUAL_ACCEPTANCE','script_sha256':sha(__file__),
            'analysis_sha256':sha(a.results/'analysis.json'),'figure_sha256':sha(out),'window_ns':[300,450],
            'peaks':peaks,'gain':False,'per_trace_normalization':False,'time_or_amplitude_fit':False}
    p=a.results/'late_visual_acceptance.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True);p.add_argument('--arrays',type=Path,required=True)
    main(p.parse_args())
