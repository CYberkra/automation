"""Confirm source equivalence and spatial shape sensitivity on paired spectra."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np

from analyze_hs4t2d_cause_controls import load_response, verify_manifest, weighted_comparison
from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import verify_official_runtime, OFFICIAL_PROCESSING_SHA256
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def array_comparison(a,b):
    norm=np.linalg.norm(a)
    if norm==0:
        return {'reason':'zero reference'}
    return {'relative_L2':float(np.linalg.norm(b-a)/norm),
            'amplitude_L2_ratio':float(np.linalg.norm(b)/norm),
            'relative_L2_by_trace':(np.linalg.norm(b-a,axis=0)/np.linalg.norm(a,axis=0)).tolist()}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--first-results',type=Path,required=True)
    ap.add_argument('--first-controls',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new output required')
    verify_official_runtime()
    verify_manifest(args.first_results)
    c=json.loads(args.contract.read_text('utf-8'))
    events=[json.loads(s) for s in (args.contract.parent/'execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('status')!='COMPLETED' or events[-1].get('traces')!=32 or events[0]['contract_sha256']!=sha256(args.contract):
        raise ValueError('complete frozen second-stage run required')
    spectra={}
    sources={}
    responses={}
    for g in c['groups']:
        response,rows=load_response(args.contract.parent/g['id'])
        if not response.source_valid.all():
            raise ValueError('all 501 normalized source bins required')
        responses[g['id']]=response
        sources[g['id']]={'quantity':response.source.quantity,'units':response.source.units,
                          'spatial_scale':response.source.spatial_scale,'source_all_501_bins_valid':True,
                          'source_spectrum_min_relative_db':float(20*np.log10(np.min(np.abs(response.source_spectrum))/np.max(np.abs(response.source_spectrum))))}
    base=responses['ricker12_rough']
    for width in (12,24,36):
        spectra[f'ricker{width}_difference']=responses[f'ricker{width}_rough'].response-responses[f'ricker{width}_halfspace'].response
    spectra['zfine13_difference']=responses['zfine13_rough'].response-responses['zfine13_halfspace'].response
    with np.load(args.first_results/'diagnostic_arrays.npz') as data:
        t=data['time_ns'].copy()
        x=data['midpoint_m'].copy()
        old_signed={key:data[key+'_signed'].copy() for key in ['narrow_halfspace_difference','wide_halfspace_difference','old_flat_difference','flat_halfspace_difference','low_air_delay_aligned']}
        old_complex={key:data[key+'_complex_envelope'].copy() for key in old_signed}
        old_spectrum=data['narrow_halfspace_difference_frequency_complex'].copy()
    products={key:reconstruct_time_response(replace(base,response=value),window='hann',zero_pad_factor=8) for key,value in spectra.items()}
    if not np.array_equal(t,products['ricker12_difference'].time[:len(t)]*1e9):
        raise ValueError('reconstruction time axes differ')
    window=(t>=160)&(t<=220)
    new_signed={key:p.real_bandpass[:len(t)] for key,p in products.items()}
    new_complex={key:p.complex_envelope[:len(t)] for key,p in products.items()}
    comparisons={'source_impulse_to_ricker_center_spectral':weighted_comparison(old_spectrum[:,6:7],spectra['ricker12_difference']),
                 'source_impulse_to_ricker_center_underground':array_comparison(old_signed['narrow_halfspace_difference'][window,6:7],new_signed['ricker12_difference'][window]),
                 'ricker_domain12_to24_center_spectral':weighted_comparison(spectra['ricker12_difference'],spectra['ricker24_difference']),
                 'ricker_domain24_to36_center_spectral':weighted_comparison(spectra['ricker24_difference'],spectra['ricker36_difference']),
                 'ricker_domain12_to24_center_underground':array_comparison(new_signed['ricker12_difference'][window],new_signed['ricker24_difference'][window]),
                 'ricker_domain24_to36_center_underground':array_comparison(new_signed['ricker24_difference'][window],new_signed['ricker36_difference'][window]),
                 'impulse_domain12_to36_13stations_underground':array_comparison(old_signed['narrow_halfspace_difference'][window],old_signed['wide_halfspace_difference'][window]),
                 'z_refinement_13stations_spectral':weighted_comparison(old_spectrum,spectra['zfine13_difference']),
                 'z_refinement_13stations_underground_signed':array_comparison(old_signed['narrow_halfspace_difference'][window],new_signed['zfine13_difference'][window]),
                 'z_refinement_13stations_underground_envelope':array_comparison(2*np.abs(old_complex['narrow_halfspace_difference'][window]),2*np.abs(new_complex['zfine13_difference'][window]))}
    closure=old_signed['narrow_halfspace_difference']-old_signed['flat_halfspace_difference']-old_signed['old_flat_difference']
    flat_energy=np.sum(old_signed['flat_halfspace_difference'][window]**2)
    rough_energy=np.sum(old_signed['narrow_halfspace_difference'][window]**2)
    # Match the centre in the new series against the previously completed fine single trace.
    replay={}
    for role in ('rough','halfspace'):
        with h5py.File(args.first_controls/f'zfine_{role}/profile.h5','r') as a,h5py.File(args.contract.parent/f'zfine13_{role}/profile7.h5','r') as b:
            replay[role]=bool(np.array_equal(a['rxs/rx1/Ey'][:],b['rxs/rx1/Ey'][:]))
            if not replay[role]:
                raise ValueError('fine centre replay not bit-identical')
    peaks={}
    for name,values in [('coarse15',old_complex['narrow_halfspace_difference']),('fine15',new_complex['zfine13_difference']),('low2_aligned',old_complex['low_air_delay_aligned'])]:
        values=2*np.abs(values[window])
        positions=t[window][np.argmax(values,axis=0)]
        peaks[name]={'peak_time_ns':positions.tolist(),'span_ns':float(np.ptp(positions)),
                     'scope':'window maxima can switch lobes; not certified interface arrivals'}
    root=Path(__file__).resolve().parents[1]
    original=root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_results'
    verify_manifest(original)
    with np.load(original/'comparison_arrays.npz') as data:
        raw=data['relief2p0_raw_signed']
        raw_t=data['time_ns']
        raw_window=(raw_t>=160)&(raw_t<=220)
        peak=float(np.max(np.abs(raw[raw_window])))
    full_scale=json.loads((original/'results.json').read_text('utf-8'))['plot_color_limits']['原始带符号响应']
    display={'raw_0_250ns_shared_scale':full_scale,'underground_raw_window_peak':peak,
             'underground_peak_to_full_plot_scale_ratio':peak/full_scale,
             'scope':'display dynamic range only; not SNR or hardware sensitivity'}
    args.out.mkdir(parents=True)
    payload={'time_ns':t,'midpoint_m':x,'frequency_Hz':base.frequency}
    for key,p in products.items():
        payload[key+'_signed']=new_signed[key]
        payload[key+'_complex_envelope']=new_complex[key]
        payload[key+'_frequency_complex']=spectra[key]
    np.savez_compressed(args.out/'confirmation_arrays.npz',**payload)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    axes[0].plot(t[window],old_signed['narrow_halfspace_difference'][window,6],label='impulse / 12 m域')
    for width in (12,24,36):
        axes[0].plot(t[window],new_signed[f'ricker{width}_difference'][window,0],linestyle='--',label=f'Ricker95 / {width} m域')
    axes[0].set_title('中心道地下窗：都先除以各自保存的复源谱；无幅度归一化')
    axes[0].set_xlabel('重建时间（ns）'); axes[0].set_ylabel('带符号场/电流响应'); axes[0].legend()
    for width in (12,24,36):
        residual=new_signed[f'ricker{width}_difference'][window,0]-old_signed['narrow_halfspace_difference'][window,6]
        axes[1].plot(t[window],residual,label=f'Ricker {width} m − impulse 12 m')
    axes[1].set_xlabel('重建时间（ns）'); axes[1].set_ylabel('响应差'); axes[1].legend()
    axes[1].set_title('差值独立色标；不能据差值更亮认定原信号更强')
    fig.savefig(args.out/'source_and_domain_center.png',dpi=145); plt.close(fig)
    columns=[('coarse15','15 m航高 / Z网格5 cm',old_signed['narrow_halfspace_difference']),
             ('fine15','15 m航高 / Z网格2.5 cm',new_signed['zfine13_difference']),
             ('low2','2 m航高 / Z网格5 cm\n仅补回固定空气时延',old_signed['low_air_delay_aligned'])]
    profile=np.genfromtxt(root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',delimiter=',',names=True)
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True,gridspec_kw={'height_ratios':[1,2]})
    for j,(name,title,matrix) in enumerate(columns):
        ax=axes[0,j]
        ax.stairs(12-profile['grid_z_m'],np.r_[profile['x0_m'],profile['x1_m'][-1]],baseline=None)
        ax.set_xlim(x[0],x[-1]); ax.set_ylim(3.4,2.4); ax.set_ylabel('距地表深度（m）')
        ax.set_title('同一实际输入台阶剖面')
        scale=float(np.quantile(np.abs(matrix[window]),.999)) or 1
        ax=axes[1,j]
        im=ax.pcolormesh(x,t[window],matrix[window]/scale,shading='nearest',cmap='gray',vmin=-1,vmax=1,rasterized=True)
        ax.set_ylim(220,160); ax.set_title(title); ax.set_ylabel('重建时间（ns）'); ax.set_xlabel('原模型收发中点 X（m）')
        fig.colorbar(im,ax=ax,label='一个全矩阵标量归一化；非幅度比较')
    fig.suptitle('V4.0.0 / CUDA FP64；同条件全覆盖层参考；13匹配站位；20–170 MHz / Hann\n'
                 '未SVD/AGC/成像；降低航高是原因对照，不代表更改实际采集条件')
    fig.savefig(args.out/'geometry_height_grid_confirmation.png',dpi=145); plt.close(fig)
    # The old flat reference contains its own interface echo. Show all terms
    # at one physical scale, rather than conflating the contrast with geometry.
    terms=[('narrow_halfspace_difference','起伏模型 − 全覆盖层参考'),
           ('flat_halfspace_difference','平界面模型 − 全覆盖层参考'),
           ('old_flat_difference','起伏模型 − 平界面模型\n等于左项减去中项')]
    limit=float(np.quantile(np.abs(np.stack([old_signed[key][window] for key,_ in terms])),.999)) or 1
    fig,axes=plt.subplots(1,3,figsize=(13,5),constrained_layout=True)
    for ax,(key,title) in zip(axes,terms):
        im=ax.pcolormesh(x,t[window],old_signed[key][window],shading='nearest',cmap='gray',vmin=-limit,vmax=limit,rasterized=True)
        ax.set_ylim(220,160); ax.set_title(title); ax.set_xlabel('收发中点 X（m）'); ax.set_ylabel('重建时间（ns）')
        fig.colorbar(im,ax=ax,label='场/电流响应，三图共用物理色标')
    fig.suptitle('平界面不是无界面背景：相减会引入负的水平界面回波\n'
                 '先复频谱相减再重建；无归一化/SVD/AGC；不把差分当 clean 真值')
    fig.savefig(args.out/'reference_subtraction_mechanism.png',dpi=145); plt.close(fig)
    report={'status':'COMPLETED','new_traces':32,'contract_sha256':sha256(args.contract),
            'first_results_manifest_sha256':sha256(args.first_results/'manifest.json'),
            'official_processing_sha256':OFFICIAL_PROCESSING_SHA256,'source_audit':sources,
            'comparisons':comparisons,'fine_center_replay_bit_identical':replay,'window_peak_diagnostic':peaks,
            'reference_algebra':{'signed_complex_closure_max_abs':float(np.max(np.abs(closure))),
                                 'flat_to_rough_contrast_window_squared_norm_ratio':float(flat_energy/rough_energy),
                                 'reference_mechanism_shared_color_limit':limit,
                                 'scope':'exact contrast decomposition, not additive physical energy contributions'},
            'display_diagnostic':display,'code_sha256':sha256(__file__),
            'interpretation':'Evidence separates reference, propagation height, in-band source and vertical-grid effects. No 3D/field extrapolation or complete numerical convergence.'}
    (args.out/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    verify_manifest(args.first_results); verify_manifest(original)
    for g in c['groups']:
        for row in json.loads((args.contract.parent/g['id']/'audit.json').read_text('utf-8')):
            if sha256(args.contract.parent/g['id']/row['file'])!=row['sha256']:
                raise ValueError('raw output changed during analysis')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)} for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
