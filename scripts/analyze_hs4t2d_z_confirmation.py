"""Check whether centre-trace changes shrink in a three-level Z sequence."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from analyze_hs4t2d_cause_controls import load_response, verify_manifest, weighted_comparison
from analyze_hs4t2d_source_controls import array_comparison
from hs_capsule_identity import sha256
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--first-results',type=Path,required=True)
    ap.add_argument('--confirmation',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new output required')
    verify_manifest(args.first_results); verify_manifest(args.confirmation)
    c=json.loads(args.contract.read_text('utf-8'))
    events=[json.loads(s) for s in (args.contract.parent/'execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('status')!='COMPLETED' or events[-1].get('traces')!=2 or events[0]['contract_sha256']!=sha256(args.contract):
        raise ValueError('complete frozen final pair required')
    rough,_=load_response(args.contract.parent/'z1p25_rough')
    half,_=load_response(args.contract.parent/'z1p25_halfspace')
    spectra={'Z1p25cm':rough.response-half.response}
    with np.load(args.first_results/'diagnostic_arrays.npz') as data:
        spectra['Z5cm']=data['narrow_halfspace_difference_frequency_complex'][:,6:7].copy()
        untapered=data['center_no_tail_taper_frequency_complex'].copy()
        truncated=data['center_record_300ns_no_taper_frequency_complex'].copy()
    with np.load(args.confirmation/'confirmation_arrays.npz') as data:
        spectra['Z2p5cm']=data['zfine13_difference_frequency_complex'][:,6:7].copy()
    products={name:reconstruct_time_response(replace(rough,response=value),window='hann',zero_pad_factor=64)
              for name,value in spectra.items()}
    t=products['Z5cm'].time*1e9
    window=(t>=160)&(t<=220)
    comparisons={}
    for a,b in [('Z5cm','Z2p5cm'),('Z2p5cm','Z1p25cm')]:
        comparisons[a+'_to_'+b]={'spectral':weighted_comparison(spectra[a],spectra[b]),
            'signed_window':array_comparison(products[a].real_bandpass[window],products[b].real_bandpass[window]),
            'envelope_window':array_comparison(2*np.abs(products[a].complex_envelope[window]),2*np.abs(products[b].complex_envelope[window]))}
    processing={}
    for name,value in [('no_tail_taper',untapered),('300ns_record_no_taper',truncated)]:
        p=reconstruct_time_response(replace(rough,response=value),window='hann',zero_pad_factor=64)
        processing[name]=array_comparison(products['Z5cm'].real_bandpass[window],p.real_bandpass[window])
    peaks={name:float(t[window][np.argmax(np.abs(p.complex_envelope[window,0]))]) for name,p in products.items()}
    args.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for name,p in products.items():
        axes[0].plot(t[window],p.real_bandpass[window,0],label=name)
        axes[1].plot(t[window],2*np.abs(p.complex_envelope[window,0]),label=name)
    for ax in axes:
        ax.legend(); ax.set_xlabel('重建时间（ns）')
    axes[0].set_ylabel('带符号场/电流响应'); axes[1].set_ylabel('2|复包络|')
    fig.suptitle('同一中心道 / 0.8 m界面 / 15 m航高 / 全覆盖层匹配差分\n'
                 'Z三档：5/2.5/1.25 cm；X/Y均保持5 cm；64倍补零仅时间插值，无增益')
    fig.savefig(args.out/'three_level_vertical_grid.png',dpi=145); plt.close(fig)
    data={'time_ns':t[window],'frequency_Hz':rough.frequency}
    for name,p in products.items():
        data[name+'_signed']=p.real_bandpass[window]
        data[name+'_complex_envelope']=p.complex_envelope[window]
        data[name+'_frequency_complex']=spectra[name]
    np.savez_compressed(args.out/'grid_arrays.npz',**data)
    report={'status':'COMPLETED','contract_sha256':sha256(args.contract),'new_traces':2,
            'comparisons':comparisons,'processing_only_center_window':processing,
            'window_envelope_max_time_ns':peaks,'source_all_bins_valid':bool(rough.source_valid.all() and half.source_valid.all()),
            'time_interpolation':'same official Hann/20-170MHz/501 bins; 64x padding diagnostic, not new bandwidth',
            'interpretation':'Z sensitivity sequence only. X remains coarse; no complete spatial convergence or absolute interface-arrival acceptance.',
            'code_sha256':sha256(__file__)}
    (args.out/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for g in c['groups']:
        for row in json.loads((args.contract.parent/g['id']/'audit.json').read_text('utf-8')):
            if sha256(args.contract.parent/g['id']/row['file'])!=row['sha256']:
                raise ValueError('raw output changed during analysis')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)}
        for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
