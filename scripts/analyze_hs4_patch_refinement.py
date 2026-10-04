"""Finer-grid and reverse-geometry robustness checks for regional attribution."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256
from analyze_hs4_height_wavefield import response, relative
from hs4_local_patch_controls import ROOT, BASE, save
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response, spectral_window


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--study',type=Path,required=True)
    p.add_argument('--patch-results',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new result required')
    check=json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if check['status']!='PASS' or check['contract_sha256']!=sha256(a.study/'execution_contract.json'): raise ValueError('verified refinement required')
    for g in check['groups']:
        for r in g['outputs']:
            if sha256(a.study/g['group']/r['file'])!=r['sha256']: raise ValueError('raw changed')
    third=ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_third'
    fb=response(third/'centre_rough/profile.h5'); fh=response(third/'centre_halfspace/profile.h5')
    fine={k:response(a.study/k/'profile.h5').response-fb.response for k in ('fine_crest','fine_slope')}
    coarse=response(BASE/'centre_rough/profile.h5'); ch=response(BASE/'centre_halfspace/profile.h5')
    reverse=response(a.study/'reverse_crest/profile.h5').response-ch.response
    with np.load(a.patch_results/'patch_arrays.npz') as old:
        coarse_change={k:old['high_'+k+'_change_spectrum'][:,6].copy() for k in ('crest','slope')}
        pair={k:old[k+'_spectrum'][:,6].copy() for k in ('high_base','high_crest_pair','high_slope_pair')}
    spectra={'coarse_'+k:v for k,v in coarse_change.items()}|fine|pair|{'reverse_crest_pair':reverse,'fine_base_pair':fb.response-fh.response}
    products={k:reconstruct_time_response(replace(coarse,response=v),window='hann',zero_pad_factor=8) for k,v in spectra.items()}
    t=products['high_base'].time*1e9; early=(t>=160)&(t<=180); full=(t>=160)&(t<=220)
    metrics={}
    for k in ('crest','slope'):
        c=products['coarse_'+k]; f=products['fine_'+k]
        metrics[k]={'early_complex_relative_L2':relative(c.complex_envelope[early],f.complex_envelope[early]),
            'full_signed_relative_L2':relative(c.real_bandpass[full],f.real_bandpass[full]),
            'early_coarse_change_norm':float(np.linalg.norm(c.complex_envelope[early])),
            'early_fine_change_norm':float(np.linalg.norm(f.complex_envelope[early]))}
    ratios={grid:float(np.linalg.norm(products[grid+'_crest'].complex_envelope[early])/np.linalg.norm(products[grid+'_slope'].complex_envelope[early])) for grid in ('coarse','fine')}
    # Diagnostic fit retains carrier/phase and positive real amplitude. It does
    # not identify an isolated event or fit a free complex phase to hide delay.
    frequencies=coarse.frequency; weights=spectral_window('hann',501)
    times=products['high_base'].time[early]
    lag_grid=np.arange(-8.,8.001,.01)
    target_times=times[None,:]-lag_grid[:,None]*1e-9
    template=(np.exp(2j*np.pi*target_times[:,:,None]*frequencies[None,None,:])@(weights*pair['high_base']))/501
    fits={}
    for k in ('high_crest_pair','reverse_crest_pair'):
        target=products[k].complex_bandpass[early]
        denominator=np.sum(abs(template)**2,axis=1)
        amplitude=np.maximum(0,np.sum((template.conj()*target).real,axis=1)/denominator)
        error=np.linalg.norm(amplitude[:,None]*template-target[None,:],axis=1)/np.linalg.norm(target)
        i=int(np.argmin(error))
        fits[k]={'lag_ns':float(lag_grid[i]),'positive_real_amplitude':float(amplitude[i]),'residual_relative_L2':float(error[i]),
            'scope':'fixed160-180ns composite-waveform fit over[-8,8]ns; not certified isolated-interface arrival'}
    a.out.mkdir(parents=True)
    arrays={'time_ns':t}
    for k,product in products.items(): arrays[k+'_signed']=product.real_bandpass; arrays[k+'_complex']=product.complex_envelope; arrays[k+'_spectrum']=spectra[k]
    np.savez_compressed(a.out/'refinement_arrays.npz',**arrays)
    summary={'status':'COMPLETED_REGIONAL_ATTRIBUTION_ROBUSTNESS_NOT_FULL_CONVERGENCE','code_sha256':sha256(__file__),
        'contract_sha256':sha256(a.study/'execution_contract.json'),'verification_sha256':sha256(a.study/'completed_verification.json'),
        'patch_input_summary_sha256':sha256(a.patch_results/'summary.json'),'metrics':metrics,
        'centre_early_crest_to_slope_change_L2_ratio':ratios,'reverse_geometry_template_checks':fits,
        'scope':'Equal-area crest/slope controls retain regional dominance at finer1/60m centre. Reverse crest changes waveform in opposite fitted-delay direction. Fine test covers centre only; no absolute FDTD error bound/full spatial convergence, energy partition or3D truth.'}
    save(a.out/'summary.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for k,label in [('high_base','原界面'),('high_crest_pair','局部抬高拱顶'),('reverse_crest_pair','局部降低拱顶')]:
        axes[0].plot(t[full],products[k].real_bandpass[full],label=label)
    for k,label in [('coarse_crest','2.5cm拱顶变化'),('fine_crest','1/60m拱顶变化'),('coarse_slope','2.5cm坡面变化'),('fine_slope','1/60m坡面变化')]:
        axes[1].plot(t[full],products[k].real_bandpass[full],label=label)
    for ax in axes: ax.set_xlabel('时间（ns）'); ax.set_ylabel('源归一化带符号响应'); ax.legend()
    fig.suptitle('15m中心道：反向扰动与区域敏感性的网格复核\n波形仍有有限网格差异，不声明绝对收敛')
    fig.savefig(a.out/'regional_robustness.png',dpi=150); plt.close(fig)
    print(json.dumps({'ratios':ratios,'reverse_checks':fits,'grid_metrics':metrics},indent=2))


if __name__=='__main__': main()
