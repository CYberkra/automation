"""Show the sole calculated station on the requested220-to25 m line; blanks stay blank."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    analysis=json.loads(a.analysis.read_text('utf-8'))
    assert sha(a.numerical)==analysis['numerical_sha256']
    m=json.loads((a.package/'manifest.json').read_text('utf-8')); g=m['groups'][1]
    with h5py.File(a.package/g['geometry']) as h:
        x=(g['tx_m'][0]+g['rx_m'][0])/2+float(h.attrs['profile_x_offset_m'])
    assert x==190.
    stations=np.arange(220.,24.9,-.5); ix=int(np.flatnonzero(stations==x)[0])
    with h5py.File(a.numerical) as h:
        t=h['time_s'][:]*1e9
        traces=[h[f'config{j}/hann/complex_bandpass'][:,1].real for j in range(3)]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    cmap=plt.get_cmap('gray').copy(); cmap.set_bad('white')
    fig,axs=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    labels=['旧泥岩1200ns H1总场','论文泥岩前1200ns H1总场','论文泥岩2400ns H1总场']
    for row,bounds in enumerate([[0,600],[250,550]]):
        mask=(t>=bounds[0])&(t<=bounds[1]); lim=max(float(abs(z[mask]).max()) for z in traces)
        for ax,z,label in zip(axs[row],traces,labels):
            data=np.full((int(mask.sum()),len(stations)),np.nan); data[:,ix]=z[mask]
            im=ax.imshow(data,cmap=cmap,vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',
                extent=[220.25,24.75,t[mask][-1],t[mask][0]])
            ax.set(title=label,xlabel='剖面横坐标 x / m（220 → 25）',ylabel='SFCW 时间 / ns')
            ax.text(.48,.5,'其余站位未计算',transform=ax.transAxes,color='gray',fontsize=11)
            ax.annotate('唯一已算站位 x=190m',xy=(190,t[mask][0]),xytext=(.02,1.015),
                textcoords='axes fraction',fontsize=8,arrowprops=dict(arrowstyle='->'))
            if row==1 and label.startswith('论文'):
                ax.plot(190,analysis['vertical_delay_guide_ns'],'o',mfc='none',mec='tab:green',ms=6)
        fig.colorbar(im,ax=list(axs[row]),label='本行共用绝对灰度：Re(复带通) / (V/m)/(A·m)')
    fig.suptitle('稀疏测线验收：仅 x=190m 已算，其余白色区域未计算、未插值\n约8m离地；gprMax4.0.1原生FP64；官方impulse/direct 20–170MHz、501点；Hann/补零8，无增益或拟合\n0.5m列宽仅用于显示已有单站，不代表新增391道；绿圈是该站垂直到时参考')
    fig.savefig(a.out,dpi=140);plt.close(fig)
    a.out.with_suffix('.json').write_text(json.dumps(dict(script_sha256=sha(__file__),
        numerical_sha256=sha(a.numerical),analysis_sha256=sha(a.analysis),figure_sha256=sha(a.out),
        profile_x_m=190,travel_distance_from220_m=30,computed_stations=1,
        blanks_not_interpolated=True,display_column_width_m=.5),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['numerical','analysis','package','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
