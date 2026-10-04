"""Actual input geometry and paired B-scan use different vertical quantities."""
import argparse
import json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.out.exists(): raise ValueError('new output required')
    root=Path(__file__).resolve().parents[1]
    profile=root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv'
    patch=root/'artifacts/research_checks/2026-10-04_hs4_local_patch_results/patch_arrays.npz'
    geom=np.genfromtxt(profile,delimiter=',',names=True)
    with np.load(patch) as h:
        x=h['high_midpoint_m']; t=h['time_ns']; mask=(t>=155)&(t<=220)
        matrices=[h[k+'_signed'][mask].copy() for k in ('high_base','high_crest_change','high_slope_change')]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes=plt.subplots(4,1,figsize=(10,12),constrained_layout=True,sharex=True,gridspec_kw={'height_ratios':[1,2,2,2]})
    edges=np.r_[geom['x0_m'],geom['x1_m'][-1]]
    axes[0].stairs(geom['cover_depth_m'],edges,color='black',linewidth=2)
    axes[0].axvspan(5,6.5,color='tab:blue',alpha=.2,label='拱顶扰动范围')
    axes[0].axvspan(8,9.5,color='tab:orange',alpha=.2,label='坡面扰动范围')
    axes[0].set_ylim(3.6,2.3); axes[0].set_ylabel('真实模型深度（m）'); axes[0].legend(loc='upper left',ncol=2)
    axes[0].set_title('同一份 0.8 m 全剖面起伏模型，截取实际扫描范围')
    shared_change=max(float(abs(v).max()) for v in matrices[1:])
    for i,(ax,m,title) in enumerate(zip(axes[1:],matrices,('原界面减匹配全覆盖层：上部较平缓','只抬高拱顶：整个扫描范围的上部波列改变','等面积抬高坡面：后续波列改变较明显'))):
        limit=float(abs(m).max()) if i==0 else shared_change
        image=ax.pcolormesh(x,t[mask],m,cmap='RdBu_r',vmin=-limit,vmax=limit,shading='nearest')
        ax.invert_yaxis(); ax.set_ylabel('接收时间（ns）'); ax.set_title(title)
        fig.colorbar(image,ax=ax,label='带符号响应' if i==0 else '变化量（两图共享尺度）')
    axes[-1].set_xlabel('收发中点原模型 X（m）')
    axes[-1].set_xlim(3,9.5)
    fig.suptitle('模型几何与 B-scan 的对应：15 m 航高 / 20–170 MHz / 原生 2.5 cm 网格\nB-scan 每列含多个界面位置的相干响应；纵轴是时间，不能直接当深度剖面')
    a.out.mkdir(parents=True)
    fig.savefig(a.out/'model_and_causal_bscan.png',dpi=150); plt.close(fig)
    result={'status':'COMPLETED_EXPLANATORY_PLOT','code_sha256':sha256(__file__),'profile_sha256':sha256(profile),'patch_arrays_sha256':sha256(patch),
        'image_sha256':sha256(a.out/'model_and_causal_bscan.png'),'scope':'Unmigrated paired response and localized sensitivity. Geometry depth and reception time are distinct axes. Native main responses; no phase adjustment, gain or per-trace normalization. Perturbation differences are not isolated energy contributions.'}
    (a.out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__': main()
