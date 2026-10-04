"""Paired local geometry sensitivities; do not interpret L2 as energy shares."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256
from hs4_local_patch_controls import save, ROOT, BASE
from analyze_hs4_height_wavefield import response
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def stack(paths):
    loaded=[response(p) for p in paths]
    return replace(loaded[0],response=np.stack([r.response for r in loaded],axis=1))


def baseline_path(role,index):
    if index==61: return BASE/('centre_'+role)/'profile.h5'
    remaining=ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_remaining'
    segment,first=('left',1) if index<61 else ('right',71)
    return remaining/(segment+'_'+role)/f'profile{(index-first)//10+1}.h5'


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--study',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.out.exists(): raise ValueError('new results required')
    c=json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    check=json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if check['status']!='PASS' or check['contract_sha256']!=sha256(a.study/'execution_contract.json'):
        raise ValueError('verified controls required')
    for g in check['groups']:
        for row in g['outputs']:
            if sha256(a.study/g['group']/row['file'])!=row['sha256']: raise ValueError('raw identity changed')
    new={g['id']:stack([a.study/g['id']/f'profile{k}.h5' for k in range(1,g['traces']+1)]) for g in c['groups']}
    high_indices=list(range(1,122,10)); low_indices=[1,61,121]
    high_rough=stack([baseline_path('rough',i) for i in high_indices])
    high_halfspace=stack([baseline_path('halfspace',i) for i in high_indices])
    high_base=high_rough.response-high_halfspace.response
    spectra={k:r.response for k,r in new.items()}
    spectra['high_base']=high_base
    spectra['high_crest_pair']=new['high_crest'].response-high_halfspace.response
    spectra['high_slope_pair']=new['high_slope'].response-high_halfspace.response
    for patch in ('crest','slope'):
        spectra['high_'+patch+'_change']=new['high_'+patch].response-high_rough.response
        spectra['low_'+patch+'_change']=new['low_'+patch].response-new['low_base'].response
    products={k:reconstruct_time_response(replace(high_rough,response=v),window='hann',zero_pad_factor=8) for k,v in spectra.items()}
    t=products['high_base'].time*1e9
    delay=2*13/299792458*1e9
    metrics={}; arrays={'time_ns':t,'high_midpoint_m':3.25+.5*np.arange(13),'low_midpoint_m':[3.25,6.25,9.25],'frequency_Hz':high_rough.frequency}
    for k,product in products.items():
        arrays[k+'_spectrum']=spectra[k]; arrays[k+'_signed']=product.real_bandpass; arrays[k+'_complex_envelope']=product.complex_envelope
    for height,shift in [('high',0),('low',delay)]:
        windows={}
        for label,(lo,hi) in {'early':(160,180),'later':(180,220),'full':(160,220)}.items():
            mask=(t+shift>=lo)&(t+shift<=hi)
            effects={patch:np.linalg.norm(products[height+'_'+patch+'_change'].complex_envelope[mask],axis=0) for patch in ('crest','slope')}
            record={patch+'_complex_change_L2':v.tolist() for patch,v in effects.items()}
            record['crest_to_slope_L2_ratio']=(effects['crest']/effects['slope']).tolist()
            record['physical_time_window_ns']=[lo-shift,hi-shift]
            if height=='high':
                base=np.linalg.norm(products['high_base'].complex_envelope[mask],axis=0)
                record['crest_change_relative_to_baseline_pair_L2']=(effects['crest']/base).tolist()
                record['slope_change_relative_to_baseline_pair_L2']=(effects['slope']/base).tolist()
            windows[label]=record
        metrics[height]=windows
    peaks={}
    # Maxima are explicitly window diagnostics; no claim that all are one event.
    mask=(t>=160)&(t<=185)
    for k in ('high_base','high_crest_pair','high_slope_pair'):
        peaks[k]=(t[mask][np.argmax(abs(products[k].complex_envelope[mask]),axis=0)]).tolist()
    with np.load(ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_results/joint_grid_arrays.npz') as old:
        delta=float(np.max(abs(old['fine_frequency_complex']-high_base)))
        if delta!=0: raise ValueError('baseline reconstruction differs from archived fine scan')
    low_ricker=response(ROOT/'artifacts/research_checks/2026-10-04_hs4_height_wavefield_continuation/low_rough/profile.h5')
    weights=np.hanning(501)
    low_source_eq=float(np.linalg.norm(weights*(new['low_base'].response[:,1]-low_ricker.response))/np.linalg.norm(weights*new['low_base'].response[:,1]))
    a.out.mkdir(parents=True)
    np.savez_compressed(a.out/'patch_arrays.npz',**arrays)
    summary={'status':'COMPLETED_LOCAL_CAUSAL_SENSITIVITY_NOT_ENERGY_DECOMPOSITION','code_sha256':sha256(__file__),
        'contract_sha256':sha256(a.study/'execution_contract.json'),'completed_verification_sha256':sha256(a.study/'completed_verification.json'),
        'new_traces':35,'high_station_indices':high_indices,'low_station_indices':low_indices,
        'same_cross_section_perturbation_m2':.125,'original_patch_ranges_m':{'crest':[5,6.5],'right_slope':[8,9.5]},
        'metrics':metrics,'high_window_envelope_maxima_ns':peaks,
        'cached_baseline_max_abs_spectrum_difference':delta,'low_Ricker_impulse_hann_complex_source_equivalence_relative_L2':low_source_eq,
        'scope':'Actual geometry perturbation at fixed grid/material/source identifies influence, not isolated reflection contribution. Ratios are waveform-change norms, not power percentages. Low windows use fixed air-time offset only; no high/low amplitude normalization or imaging. Whole-domain grid convergence and finite3D antenna remain separate.'}
    save(a.out/'summary.json',summary)
    plot(a.out,arrays,summary)
    print(json.dumps({'high_early_crest_to_slope_ratio':metrics['high']['early']['crest_to_slope_L2_ratio'],
        'low_early_crest_to_slope_ratio':metrics['low']['early']['crest_to_slope_L2_ratio'],'low_source_equivalence':low_source_eq},indent=2))


def plot(out,a,s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    t=a['time_ns']; mask=(t>=155)&(t<=220); x=a['high_midpoint_m']
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True,gridspec_kw={'height_ratios':[1,1]})
    base=[a[k+'_signed'][mask] for k in ('high_base','high_crest_pair','high_slope_pair')]
    limit=max(float(np.max(abs(v))) for v in base)
    for ax,m,label in zip(axes[0],base,['原界面配对','局部抬高拱顶','等面积抬高右坡面']):
        im=ax.pcolormesh(x,t[mask],m,cmap='RdBu_r',vmin=-limit,vmax=limit,shading='nearest')
        ax.invert_yaxis(); ax.set_title(label); ax.set_xlabel('收发中点原X（m）'); ax.set_ylabel('时间（ns）')
    change=[a[k+'_signed'][mask] for k in ('high_crest_change','high_slope_change')]
    limit_delta=max(float(np.max(abs(v))) for v in change)
    for ax,m,label in zip(axes[1,:2],change,['拱顶扰动导致的响应变化','坡面扰动导致的响应变化']):
        im2=ax.pcolormesh(x,t[mask],m,cmap='RdBu_r',vmin=-limit_delta,vmax=limit_delta,shading='nearest')
        ax.invert_yaxis(); ax.set_title(label); ax.set_xlabel('收发中点原X（m）'); ax.set_ylabel('时间（ns）')
    for patch in ('crest','slope'):
        axes[1,2].plot(x,s['metrics']['high']['early'][patch+'_complex_change_L2'],label=patch)
    axes[1,2].set_title('160–180ns响应变化范数'); axes[1,2].set_xlabel('原X（m）'); axes[1,2].legend()
    fig.colorbar(im,ax=axes[0].tolist(),label='配对带符号响应，共享尺度')
    fig.colorbar(im2,ax=axes[1,:2].tolist(),label='扰动响应变化，两图共享尺度')
    fig.suptitle('15m / 原始2.5cm格 / 等面积局部扰动 / 20–170MHz\n不改变边界或地表，不将变化量当能量占比')
    fig.savefig(out/'localized_patch_bscan.png',dpi=150); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,height,xx in zip(axes,['high','low'],[a['high_midpoint_m'],a['low_midpoint_m']]):
        for patch in ('crest','slope'):
            ax.plot(xx,s['metrics'][height]['early'][patch+'_complex_change_L2'],'o-',label=patch)
        ax.set_title('15m' if height=='high' else '2m（固定空气时延窗口对照）'); ax.set_xlabel('原X（m）'); ax.set_ylabel('早期窗复波形变化L2'); ax.legend()
    fig.savefig(out/'patch_spatial_sensitivity.png',dpi=150); plt.close(fig)


if __name__=='__main__': main()
