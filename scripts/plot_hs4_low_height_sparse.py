"""Show actual2m/15m raw responses at3 stations; gaps are uncomputed.

CPU only. No spatial interpolation, background subtraction or time correction.
"""
import argparse
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256, read_manifest, verify_file
from hs4_cross_pc import ROOT, save, read
from analyze_hs4_cross_pc import response, reconstruct
from analyze_hs4_local_patch_controls import baseline_path
from sfcw_official_loader_v0_2 import verify_official_runtime


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True); args=parser.parse_args()
    if args.out.exists(): raise ValueError('new output required')
    verify_official_runtime()
    checks=ROOT/'artifacts/research_checks'
    study=checks/'2026-10-04_hs4_local_patch_controls'
    verification=read(study/'completed_verification.json')
    if verification['status']!='PASS' or verification['contract_sha256']!=sha256(study/'execution_contract.json'):
        raise ValueError('verified low-height controls required')
    low_records=next(g for g in verification['groups'] if g['group']=='low_base')['outputs']
    records=[]; products={}; spectra={}; poses=[]
    for height in (2,15):
        loaded=[]
        for j,index in enumerate((1,61,121),1):
            raw=study/'low_base'/f'profile{j}.h5' if height==2 else baseline_path('rough',index)
            if height==2:
                if sha256(raw)!=low_records[j-1]['sha256']: raise ValueError('raw low identity changed')
            else:
                base=raw.parents[1]; verify_file(base,read_manifest(base),raw.relative_to(base).as_posix())
            with h5py.File(raw) as h:
                dx=np.asarray(h.attrs['dx_dy_dz']); tx=h['srcs/src1'].attrs['GridPosition']*dx
                rx=h['rxs/rx1'].attrs['GridPosition']*dx
                if str(h.attrs['gprMax'])!='4.0.0' or h['rxs/rx1/Ey'].dtype!=np.float64 or not np.array_equal(dx,[.025,.05,.025]):
                    raise ValueError('native runtime/precision/grid mismatch')
                if not np.allclose([tx[2],rx[2]],[12+height]*2,rtol=0,atol=1e-12): raise ValueError('height differs')
                midpoint=(tx[0]+rx[0])/2-12
                if not np.isclose(midpoint,[3.25,6.25,9.25][j-1],rtol=0,atol=1e-12): raise ValueError('station differs')
            records.append({'height_m':height,'station_index':index,'file':raw.relative_to(ROOT).as_posix(),'sha256':sha256(raw),'tx_m':tx.tolist(),'rx_m':rx.tolist()})
            r=response(raw,'Ey'); loaded.append(r)
        spectra[str(height)]=np.stack([r.response for r in loaded],axis=1)
        products[str(height)]=reconstruct(loaded[0],spectra[str(height)])
    cached=checks/'2026-10-04_hs4_local_patch_results/patch_arrays.npz'
    with np.load(cached) as h:
        if not np.array_equal(h['low_base_signed'],products['2'].real_bandpass) or not np.array_equal(h['low_base_spectrum'],spectra['2']):
            raise ValueError('low raw reconstruction differs from archived arrays')
    profile=checks/'2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv'
    geom=np.genfromtxt(profile,delimiter=',',names=True)
    t=products['2'].time*1e9; x=np.array([3.25,6.25,9.25])
    # Absolute receiver-time windows; low height arrives sooner. No alignment.
    windows={'2':(55,145),'15':(140,230)}
    masks={k:(t>=lo)&(t<=hi) for k,(lo,hi) in windows.items()}
    maximum=max(float(np.max(abs(products[k].real_bandpass[masks[k]]))) for k in products)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    plt.rcParams['font.family']='Microsoft YaHei'
    norm=SymLogNorm(linthresh=maximum*1e-3,vmin=-maximum,vmax=maximum)
    fig,axes=plt.subplots(3,1,figsize=(11,10),constrained_layout=True,sharex=True,
                          gridspec_kw={'height_ratios':[1,2.4,2.4]})
    axes[0].stairs(geom['cover_depth_m'],np.r_[geom['x0_m'],geom['x1_m'][-1]],color='black',linewidth=2)
    axes[0].set_ylim(3.6,2.3); axes[0].set_ylabel('界面深度（m）')
    axes[0].set_title('同一份 0.8 m 起伏模型；地表水平，收发间距 1.3 m')
    for ax,k in zip(axes[1:],('2','15')):
        mask=masks[k]; times=t[mask]; step=t[1]-t[0]
        time_edges=np.r_[times-step/2,times[-1]+step/2]
        ax.set_facecolor('#dedede')
        for j,position in enumerate(x):
            im=ax.pcolormesh([position-.3,position+.3],time_edges,
                            products[k].real_bandpass[mask,j,None],cmap='RdBu_r',norm=norm,shading='flat')
        ax.set_ylim(windows[k][1],windows[k][0]); ax.set_ylabel('实际接收时间（ns）')
        ax.set_title(f'{k} 米航高：原始总响应的地下时间窗（仅左、中、右三道）')
    axes[-1].set_xticks(x,labels=['左 3.25','中 6.25','右 9.25']); axes[-1].set_xlim(2.8,9.7)
    axes[-1].set_xlabel('收发中点原模型 X（m）；灰色区域未计算，不做插值')
    fig.colorbar(im,ax=axes[1:].tolist(),label='源归一化带符号响应；两图共用 SymLog 色标')
    fig.suptitle('降低航高的现有结果：2 m 与 15 m 同站位稀疏 B-scan\n20–170 MHz / 原生 2.5 cm 网格 / 无 AGC、逐道归一化或背景扣除')
    args.out.mkdir(parents=True)
    fig.savefig(args.out/'low_height_sparse_comparison.png',dpi=150); plt.close(fig)
    np.savez_compressed(args.out/'raw_height_arrays.npz',time_ns=t,midpoint_m=x,
        low_signed=products['2'].real_bandpass,high_signed=products['15'].real_bandpass,
        low_spectrum=spectra['2'],high_spectrum=spectra['15'],
        low_complex=products['2'].complex_envelope,high_complex=products['15'].complex_envelope)
    save(args.out/'summary.json',{'status':'COMPLETED_EXISTING_3_STATION_DISPLAY_NOT_FULL_SCAN',
        'solver_executed':False,'script_sha256':sha256(__file__),'raw_files':records,
        'profile_sha256':sha256(profile),'low_cached_arrays_bitwise_equal':True,
        'windows_ns':windows,'color_norm':{'type':'SymLog','shared':True,'vmax':maximum,'linthresh':maximum*1e-3},
        'column_display_width_m':.6,'column_width_scope':'Drawing only, not antenna footprint or spatial sampling coverage.',
        'processing':'Same501 source-normalized complex tones20-170MHz;200ns tail;Hann/8x; raw total response without background removal; physical reception time, no delay alignment.',
        'image_sha256':sha256(args.out/'low_height_sparse_comparison.png'),
        'arrays_sha256':sha256(args.out/'raw_height_arrays.npz'),
        'limits':'Three stations cannot establish a continuous interface curve. Gray gaps are uncomputed. Raw response includes direct/surface/band sidelobes; this is not the matched whole-cover difference used for regional attribution. No migration, numerical convergence or real-antenna acceptance.'})
    print('PASS:6 verified native files;2m3 traces match archived signed arrays exactly')


if __name__=='__main__': main()
