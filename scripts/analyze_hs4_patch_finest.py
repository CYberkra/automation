"""Three native grids: verify regional attribution without an absolute error claim."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256
from hs4_local_patch_controls import ROOT, BASE, save
from analyze_hs4_height_wavefield import response, relative
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True)
    p.add_argument('--refinement',type=Path,required=True)
    p.add_argument('--patch-results',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new output directory required')
    checks=[]
    for folder in (a.study,a.refinement):
        check=json.loads((folder/'completed_verification.json').read_text('utf-8'))
        if check['status']!='PASS' or check['contract_sha256']!=sha256(folder/'execution_contract.json'):
            raise ValueError('verified capsule required')
        for group in check['groups']:
            for record in group['outputs']:
                if sha256(folder/group['group']/record['file'])!=record['sha256']:
                    raise ValueError('raw output identity changed')
        checks.append(sha256(folder/'completed_verification.json'))
    reference=response(BASE/'centre_rough/profile.h5')
    third=ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_third'
    midbase=response(third/'centre_rough/profile.h5')
    spectra={}
    with np.load(a.patch_results/'patch_arrays.npz') as old:
        spectra['2p5cm_base']=old['high_base_spectrum'][:,6].copy()
        for k in ('crest','slope'):
            spectra['2p5cm_'+k]=old['high_'+k+'_change_spectrum'][:,6].copy()
    spectra['1over60m_base']=midbase.response-response(third/'centre_halfspace/profile.h5').response
    for k in ('crest','slope'):
        spectra['1over60m_'+k]=response(a.refinement/('fine_'+k)/'profile.h5').response-midbase.response
    finestbase=response(a.study/'base/profile.h5')
    spectra['1p25cm_base']=finestbase.response-response(a.study/'halfspace/profile.h5').response
    for k in ('crest','slope'):
        spectra['1p25cm_'+k]=response(a.study/k/'profile.h5').response-finestbase.response
    products={k:reconstruct_time_response(replace(reference,response=v),window='hann',zero_pad_factor=8) for k,v in spectra.items()}
    t=next(iter(products.values())).time*1e9
    early=(t>=160)&(t<=180); full=(t>=160)&(t<=220)
    metrics={}; grids=('2p5cm','1over60m','1p25cm')
    for left,right in zip(grids[:-1],grids[1:]):
        metrics[left+'_to_'+right]={}
        for k in ('base','crest','slope'):
            c=products[left+'_'+k]; f=products[right+'_'+k]
            metrics[left+'_to_'+right][k]={
                'early_complex_relative_L2':relative(c.complex_envelope[early],f.complex_envelope[early]),
                'full_signed_relative_L2':relative(c.real_bandpass[full],f.real_bandpass[full]),
                'full_envelope_magnitude_relative_L2':relative(abs(c.complex_envelope[full]),abs(f.complex_envelope[full]))}
    ratios={g:float(np.linalg.norm(products[g+'_crest'].complex_envelope[early])/np.linalg.norm(products[g+'_slope'].complex_envelope[early])) for g in grids}
    a.out.mkdir(parents=True)
    arrays={'time_ns':t}
    for k,v in products.items():
        arrays[k+'_spectrum']=spectra[k]; arrays[k+'_signed']=v.real_bandpass; arrays[k+'_complex']=v.complex_envelope
    np.savez_compressed(a.out/'finest_arrays.npz',**arrays)
    summary={'status':'COMPLETED_THREE_GRID_CENTRE_REGIONAL_CHECK_NOT_ABSOLUTE_CONVERGENCE',
        'code_sha256':sha256(__file__),'verified_capsule_sha256':checks,
        'patch_results_sha256':sha256(a.patch_results/'summary.json'),
        'metrics':metrics,'early_crest_to_slope_change_L2_ratio':ratios,
        'scope':'Same geometry, physical PML and source; centre station only. Ratios are waveform-change sensitivity, not separated energy. No Richardson error budget, full-scan convergence, finite3D or field acceptance.'}
    save(a.out/'summary.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes=plt.subplots(3,1,figsize=(10,9),constrained_layout=True)
    for ax,k,label in zip(axes,('base','crest','slope'),('原界面减全覆盖层','局部抬高拱顶引起的变化','局部抬高坡面引起的变化')):
        for g,name in zip(grids,('2.5 cm','1/60 m','1.25 cm')):
            ax.plot(t[full],products[g+'_'+k].real_bandpass[full],label=name)
        ax.set_title(label); ax.set_xlabel('时间（ns）'); ax.set_ylabel('源归一化带符号响应'); ax.legend()
    fig.suptitle('15 m 航高中心道：三档原生网格的区域归因复核')
    fig.savefig(a.out/'three_grid_regional_check.png',dpi=150); plt.close(fig)
    print(json.dumps({'metrics':metrics,'ratios':ratios},indent=2))


if __name__=='__main__': main()
