"""Supplemental raw late fields and sector-loss visual audit; no data gain."""
import argparse,json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    result=json.loads((a.out/'analysis.json').read_text('utf-8'));assert sha(a.arrays)==result['arrays_sha256']
    with np.load(a.arrays) as h:d={k:h[k] for k in ['time_s','Ez','Hx','sfcw_time_s','full_01_hann','full_12_hann','full_02_hann']}
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    t=d['time_s']*1e9;gate=(t>=220)&(t<=350);fields=[d['Ez'],d['Ez']*d['Hx']];limits=[.75,.005]
    fig,axs=plt.subplots(3,2,figsize=(14,12),layout='constrained');clip=[]
    for i,y in enumerate([27.025,28.525,30.025]):
        for j,v in enumerate(fields):
            z=v[i*41:(i+1)*41,gate].T;fraction=float(np.mean(abs(z)>limits[j]));clip.append({'y_m':y,'field':['Ez','Sy'][j],'clipped_fraction':fraction})
            im=axs[i,j].imshow(z,cmap=['gray','RdBu_r'][j],vmin=-limits[j],vmax=limits[j],aspect='auto',extent=[153.75,174.25,350,220])
            axs[i,j].set(title=f'y={y}m：'+['总场Ez','总场Sy（红上/蓝下）'][j]+f'；截色{fraction*100:.2f}%',xlabel='模型 x / m',ylabel='原生时间 / ns');fig.colorbar(im,ax=axs[i,j],label=['Ez / (V/m)','Sy / (W/m²)'][j])
    fig.suptitle('三条覆盖层观测线的弱晚场：各列固定共同尺度，原幅值，无增益\n一次固定激发的41列时空图，非移动天线B-scan；不能依据总场局部正峰串成唯一射线')
    paths=['native_cover_late_shared_scale.png','cover_sector_excludes_late_gray.png']
    for name in paths:assert not (a.out/name).exists()
    fig.savefig(a.out/paths[0],dpi=140);plt.close(fig)
    t=d['sfcw_time_s']*1e9;gate=(t>=220)&(t<=330);fig,axs=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    limit=max(float(abs(d['full_'+pair+'_hann'][gate][:,:,[0,6]].real).max()) for pair in ['01','12','02'])
    for i,pair in enumerate(['01','12','02']):
        z=d['full_'+pair+'_hann']
        for j,index in enumerate([6,0]):
            im=axs[i,j].imshow(z[gate,:,index].real,cmap='gray',vmin=-limit,vmax=limit,aspect='auto',extent=[153.75,174.25,330,220]);axs[i,j].set(title=f'平面对{pair}目标：'+['全横向样本总场','空气可辐射平滑扇区总场'][j],xlabel='模型 x / m',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=axs[i,j],label='(V/m)/(A·m)；全图共同尺度')
    fig.suptitle('覆盖层晚波：原采样场与受限角谱的差别；Hann精确501点\n同一空间Tukey0.25；未拟合幅度/相位/时延，扇区场不能代替全部覆盖层响应')
    fig.savefig(a.out/paths[1],dpi=140);plt.close(fig)
    proof={'status':'SUPPLEMENTAL_RAW_LATE_SHARED_SCALE_AND_SECTOR_AUDIT','script_sha256':sha(__file__),'analysis_sha256':sha(a.out/'analysis.json'),
           'raw_late_shared_limits':dict(Ez_V_m=.75,Sy_W_m2=.005),'clipping':clip,'sfcw_common_limit':limit,
           'figures':{name:sha(a.out/name) for name in paths},'new_solver_runs':0}
    assert not (a.out/'supplemental_figures.json').exists();(a.out/'supplemental_figures.json').write_text(json.dumps(proof,indent=2)+'\n',encoding='utf-8');print(proof['status'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--arrays',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
