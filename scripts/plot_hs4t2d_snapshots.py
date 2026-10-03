"""Plot audited broadband Ey snapshots with common signed logarithmic scales."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--profile',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new plot directory required')
    c=json.loads((args.run/'execution_contract.json').read_text('utf-8'))
    execution=json.loads((args.run/'execution.json').read_text('utf-8'))
    if execution['status']!='COMPLETED' or execution['contract_sha256']!=sha256(args.run/'execution_contract.json'):
        raise ValueError('completed frozen run required')
    source=Path(c['source_input']).read_bytes()
    actual=(args.run/'profile.in').read_bytes()
    appended=actual[len(source):].decode('ascii').splitlines()
    if not actual.startswith(source) or len(appended)!=9 or any(not line.startswith('#snapshot:') for line in appended):
        raise ValueError('input must only append frozen passive snapshots')
    profile_digest=sha256(args.profile)
    profile=np.genfromtxt(args.profile,delimiter=',',names=True)
    edges=np.r_[profile['x0_m'],profile['x1_m'][-1]]
    rows=sorted(execution['snapshots'],key=lambda r:r['actual_time_ns'])
    fields=[]
    for row in rows:
        p=args.run/row['file']
        if sha256(p)!=row['sha256']:
            raise ValueError('snapshot identity changed')
        with h5py.File(p,'r') as h:
            if not np.array_equal(h.attrs['dx_dy_dz'],[.05,.05,.05]) or not np.array_equal(h.attrs['origin'],[0,0,0]):
                raise ValueError('snapshot physical coordinates differ')
            if h['Ey'].shape!=(240,1,660):
                raise ValueError('unexpected X,Y,Z snapshot layout')
            fields.append(h['Ey'][:,0,:].T)
    limit=max(float(np.max(np.abs(v))) for v in fields)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.patches import Rectangle
    from PIL import Image
    plt.rcParams['font.family']='Microsoft YaHei'
    norm=SymLogNorm(linthresh=limit*1e-4,vmin=-limit,vmax=limit,base=10)
    args.out.mkdir(parents=True)
    def draw(ax,k,zoom=False):
        im=ax.imshow(fields[k],origin='lower',extent=(0,12,0,33),cmap='RdBu_r',norm=norm,
                     interpolation='nearest',aspect='equal')
        ax.axhline(12,color='black',linewidth=.8,linestyle='--')
        ax.stairs(profile['grid_z_m'],edges,baseline=None,color='black',linewidth=.8)
        ax.add_patch(Rectangle((1,1),10,31,fill=False,linestyle=':',edgecolor='gray'))
        ax.plot([5.6],[27],'r^',label='Tx')
        ax.plot([6.9],[27],'kv',label='Rx')
        ax.set_xlim(0,12)
        ax.set_ylim((6,13) if zoom else (0,33))
        ax.set_xlabel('X（m）')
        ax.set_ylabel('Z（m）')
        ax.set_title(f"Ey：{rows[k]['actual_time_ns']:.2f} ns")
        return im
    frame_paths=[]
    for k in range(len(rows)):
        fig,ax=plt.subplots(figsize=(6,9),constrained_layout=True)
        im=draw(ax,k)
        fig.colorbar(im,ax=ax,label='Ey（V/m，共享带符号对数色标）',shrink=.75)
        ax.legend(loc='upper left')
        fig.suptitle('V4.0.0 / CUDA FP64；0.8 m 起伏中心道\n'
                     '宽带 impulse 瞬时波场；非 20–170 MHz 合成波场\n'
                     '黑虚线地表；黑实线界面；灰点线 PML 内缘')
        p=args.out/f"wave_{rows[k]['requested_time_ns']:03.0f}.png"
        fig.savefig(p,dpi=115)
        plt.close(fig)
        frame_paths.append(p)
    frames=[Image.open(p).convert('RGB') for p in frame_paths]
    frames[0].save(args.out/'wavefield.gif',save_all=True,append_images=frames[1:],duration=800,loop=0)
    for frame in frames:
        frame.close()
    chosen=[2,3,4,5,6,7]
    fig,axes=plt.subplots(2,3,figsize=(13,6.5),constrained_layout=True)
    for ax,k in zip(axes.flat,chosen):
        im=draw(ax,k,zoom=True)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Ey（V/m，同动画共享带符号对数色标）',shrink=.8)
    fig.suptitle('地下局部：同一中心道瞬时宽带波场；黑实线为实际台阶界面\n'
                 '同色标、无归一化/AGC；快照不足以单独认定 SFCW B-scan 形态根因')
    fig.savefig(args.out/'underground_snapshots.png',dpi=145)
    plt.close(fig)
    report={'status':'PASS','snapshot_axes':'stored [X,Y,Z]; plotted Ey[:,0,:].T at cell centres',
            'shared_norm':'SymLogNorm','symmetric_limit_V_per_m':limit,'linear_threshold_V_per_m':limit*1e-4,
            'source_directives_unchanged':True,'snapshot_count':9,'profile_sha256':profile_digest,
            'execution_sha256':sha256(args.run/'execution.json'),'code_sha256':sha256(__file__),
            'physical_interpretation':'broadband impulse snapshots only; no geometry-to-SFCW acceptance'}
    (args.out/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)}
        for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
