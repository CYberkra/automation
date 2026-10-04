"""Fixed-band paired diagnostics for portable HS4 stages and endpoint controls."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import h5py
import numpy as np

from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime
from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
    direct_frequency_response, reconstruct_time_response)
from hs4_cross_pc import ROOT, FINE, save, read


def response(path, component):
    s=load_source(path); r=load_receiver(path,receiver_path='/rxs/rx1',component=component)
    product=direct_frequency_response(s,r,FREQ,tail_taper_fraction=(round(200e-9/r.dt)-.25)/len(r.samples))
    if not product.source_valid.all(): raise ValueError('invalid source band bin')
    return product


def relative(a,b):
    return float(np.linalg.norm(b-a)/np.linalg.norm(a)) if np.linalg.norm(a)>0 else None


def reconstruct(reference, spectrum):
    return reconstruct_time_response(replace(reference,response=spectrum),window='hann',zero_pad_factor=8)


def input_path(capsule, group):
    p=Path(group['input'])
    # Historical contracts keep their original absolute bytes/hash. Resolve
    # their files within the cloned capsule, never in an old machine's folder.
    return capsule/group['id']/p.name if p.is_absolute() else capsule/p


def verify_raw(capsule):
    c=read(capsule/'execution_contract.json'); check=read(capsule/'completed_verification.json')
    if check['status']!='PASS' or check['contract_sha256']!=sha256(capsule/'execution_contract.json'):
        raise ValueError('verified complete capsule required')
    records={g['id']:g for g in check['groups']}
    for g in c['groups']:
        p=input_path(capsule,g)
        if sha256(p)!=g['input_sha256'] or sha256(p.with_suffix('.h5'))!=records[g['id']]['raw_sha256']:
            raise ValueError('frozen input/raw identity changed')
    return c


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--capsule',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); p.add_argument('--reference-capsule',type=Path)
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new result directory required')
    verify_official_runtime(); c=verify_raw(a.capsule)
    stage=c.get('stage',c.get('mode')); component='Ex' if stage.startswith('3d') else 'Ey'
    responses={}; native={}
    for g in c['groups']:
        raw=input_path(a.capsule,g).with_suffix('.h5')
        responses[g['id']]=response(raw,component)
        with h5py.File(raw) as h: native[g['id']]=h['rxs/rx1/'+component][:]
    ref=next(iter(responses.values())); spectra={}; comparisons={}; runtime=c.get('runtime',{})
    if stage in ('replay','centre1cm'):
        for name,r in responses.items():
            old=response(FINE/name/'profile.h5','Ey')
            comparisons[name]={'full501_complex_relative_L2':relative(old.response,r.response)}
            if stage=='replay':
                with h5py.File(FINE/name/'profile.h5') as h: original=h['rxs/rx1/Ey'][:]
                comparisons[name]['native_receiver_bitwise_equal']=bool(np.array_equal(original,native[name]))
        spectra['base']=responses['base'].response-responses['halfspace'].response
        for variant in ('crest','slope'):
            spectra[variant]=responses[variant].response-responses['base'].response
        for name in ('base','crest','slope'):
            oldbase=response(FINE/'base/profile.h5','Ey')
            oldspec=oldbase.response-response(FINE/'halfspace/profile.h5','Ey').response if name=='base' else response(FINE/name/'profile.h5','Ey').response-oldbase.response
            spectra['old_'+name]=oldspec
    elif stage=='endpoints':
        for station in ('left','right'):
            base=responses[station+'_base'].response
            spectra[station+'_base']=base-responses[station+'_halfspace'].response
            for name in ('crest','slope'):
                spectra[station+'_'+name]=responses[station+'_'+name].response-base
        old=ROOT/'artifacts/research_checks/2026-10-04_hs4_local_patch_results/patch_arrays.npz'
        with np.load(old) as h:
            for station,index in [('left',0),('right',12)]:
                spectra['old_'+station+'_base']=h['high_base_spectrum'][:,index]
                for name in ('crest','slope'):
                    spectra['old_'+station+'_'+name]=h['high_'+name+'_change_spectrum'][:,index]
    else:
        stations=sorted({name.rsplit('_',1)[0] for name in responses})
        for station in stations:
            rough=responses[station+'_rough08'].response
            for reference in ('flat','halfspace'):
                spectra[station+'_rough_minus_'+reference]=rough-responses[station+'_'+reference].response
            spectra[station+'_raw_rough']=rough
        if a.reference_capsule:
            other=verify_raw(a.reference_capsule)
            for g in other['groups']:
                name=g['id']; raw=input_path(a.reference_capsule,g).with_suffix('.h5')
                if name in responses:
                    old=response(raw,'Ex'); comparisons[name]={'full501_complex_relative_L2':relative(old.response,responses[name].response)}
    products={name:reconstruct(ref,s) for name,s in spectra.items()}
    t=next(iter(products.values())).time*1e9; metrics={}; arrays={'time_ns':t,'frequency_Hz':FREQ}
    for name,v in products.items():
        arrays[name+'_spectrum']=spectra[name]; arrays[name+'_signed']=v.real_bandpass; arrays[name+'_complex']=v.complex_envelope
    for label,(lo,hi) in {'early':(160,180),'later':(180,220),'full':(160,220)}.items():
        mask=(t>=lo)&(t<=hi); record={}
        for name,v in products.items():
            if 'old_'+name in products:
                old=products['old_'+name]
                record[name]={'signed_relative_L2':relative(old.real_bandpass[mask],v.real_bandpass[mask]),
                              'complex_relative_L2':relative(old.complex_envelope[mask],v.complex_envelope[mask]),
                              'magnitude_relative_L2':relative(abs(old.complex_envelope[mask]),abs(v.complex_envelope[mask]))}
            if name.endswith('crest') and not name.startswith('old_'):
                slope=name[:-5]+'slope'
                if slope in products:
                    den=np.linalg.norm(products[slope].complex_envelope[mask])
                    record[name+'_to_slope_L2_ratio']=float(np.linalg.norm(v.complex_envelope[mask])/den) if den>0 else None
        metrics[label]=record
    a.out.mkdir(parents=True); np.savez_compressed(a.out/'arrays.npz',**arrays)
    summary={'status':'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE','stage':stage,
             'code_sha256':sha256(__file__),'contract_sha256':sha256(a.capsule/'execution_contract.json'),
             'verification_sha256':sha256(a.capsule/'completed_verification.json'),
             'runtime_source_identities':runtime.get('source_identities',c.get('source_identities')),
             'historical_build_scope':runtime.get('comparison_scope','SAME_LOCAL_FROZEN_RUNTIME'),
             'comparisons':comparisons,'metrics':metrics,'arrays_sha256':sha256(a.out/'arrays.npz'),
             'source_normalization':'Native emitted electric-current source; line-source spatial scale for2D; not portS21.',
             'limitations':'Adjacent-grid differences are not absolute error bounds. Whole-cover difference includes desired interfaces and is not clean training truth.3D ideal dipole is not finite antenna. No migrated shape or field validity is asserted.'}
    save(a.out/'summary.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    shown=[name for name in products if not name.startswith('old_')]
    fig,axes=plt.subplots(len(shown),1,figsize=(11,max(3,2.4*len(shown))),squeeze=False,constrained_layout=True)
    mask=(t>=140)&(t<=240)
    maximum=max(np.max(abs(products[name].real_bandpass[mask])) for name in shown)
    for ax,name in zip(axes[:,0],shown):
        ax.plot(t[mask],products[name].real_bandpass[mask],label=name)
        if 'old_'+name in products: ax.plot(t[mask],products['old_'+name].real_bandpass[mask],ls='--',label='coarser archived grid')
        ax.set_ylim(-maximum*1.05,maximum*1.05); ax.set_xlabel('Time (ns)'); ax.set_ylabel('Source-normalized response'); ax.legend()
    fig.suptitle(stage+': shared amplitude scale; paired/change diagnostics')
    fig.savefig(a.out/'comparison.png',dpi=150); plt.close(fig)
    print(json.dumps({'stage':stage,'metrics':metrics,'comparisons':comparisons},indent=2))


if __name__=='__main__': main()
