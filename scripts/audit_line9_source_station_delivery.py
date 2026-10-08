"""Independent native, exact-transform, direct-convolution and causal-metric audit."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.private/'execution_contract.json');v=read(a.private/'completed_verification.json');r=read(a.public/'analysis.json');m=c['study_manifest']
    events=[json.loads(s) for s in (a.private/'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256']==v['contract_sha256']==r['contract_sha256']==sha(a.private/'execution_contract.json')
    assert events[-1]['status']=='COMPLETED' and events[-1]['verification_sha256']==r['verification_sha256']==sha(a.private/'completed_verification.json')
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e];done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    assert [e['group'] for e in starts]==[e['group'] for e in done]==[g['id'] for g in c['groups']]
    assert len(done)==c['max_runs']==events[-1]['traces']==9
    assert sha(a.prepared/'manifest.json')==c['file_identities'][c['package']+'\\manifest.json']
    natives={};sources={}
    for original,g,row,e in zip(m['groups'],c['groups'],v['groups'],done):
        assert sha(a.prepared/original['input'])==g['input_sha256'] and sha(a.prepared/g['geometry'])==g['geometry_sha256'] and sha(a.prepared/g['material'])==g['material_sha256']
        path=a.private/(g['id']+'.h5');assert sha(path)==row['native_sha256']==e['raw_sha256']
        with h5py.File(path) as h:
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];dt=float(h.attrs['dt']);ea=h['srcs/src1/excitation'].attrs
            assert str(h.attrs['gprMax'])=='4.0.0' and x.dtype==s.dtype==np.float64 and np.isfinite(x).all() and np.isfinite(s).all()
            assert x.shape==s.shape==(g['expected_samples'],) and h.attrs['Iterations']==g['expected_samples'] and dt==row['dt_s']
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[4400,1600,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-12,rtol=0)
            assert ea['WaveformType']==g['source_type'] and ea['WaveformFrequency']==g['source_frequency_Hz'] and ea['SpatialScale']==.025 and ea['WaveformAmplitude']==40
            assert abs((len(x)-1)*dt-g['time_window_ns']*1e-9)<=dt
            natives[g['id']]=x;sources[g['id']]=s
    for j,name in enumerate(['H1','H0']):
        assert sha(a.prior/(name+'.h5'))==r['native_sha256'][j]
        with h5py.File(a.prior/(name+'.h5')) as h:natives['c0008_'+name]=h['rxs/rx1/Ez'][:]
    assert r['script_sha256']==sha(Path(__file__).parent/'analyze_line9_source_station_controls.py') and r['numerical_sha256']==sha(a.private/'sfcw.h5')
    direct_convolution={}
    for name in ['ricker_800','ricker_1600']:
        n=len(natives[name]);take=np.unique(np.r_[np.arange(1,n,397),np.arange(round(250e-9/dt),round(450e-9/dt),137)])
        expected=np.array([np.dot(natives['impulse_1600'][:k+1][::-1],sources[name][:k+1])/40 for k in take])
        error=float(np.linalg.norm(expected-natives[name][take])/np.linalg.norm(natives[name][take]));assert error<1e-9
        direct_convolution[name]=error
    for short,long in [('c0008_H0','impulse_1600'),('ricker_800','ricker_1600')]:assert np.array_equal(natives[short],natives[long][:len(natives[short])])
    inverse_errors={}
    with h5py.File(a.private/'sfcw.h5') as h:
        f=h['frequency_Hz'][:];t=h['time_s'][:];response=h['response'][:];names=h.attrs['columns'].split(',');ix={s:j for j,s in enumerate(names)}
        np.testing.assert_array_equal(f,20e6+300000*np.arange(501));assert names==r['columns']
        for window in ['hann','blackman']:
            w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();z=h[window+'_complex_bandpass'][:];take=np.arange(0,len(t),113)
            manual=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*response)/501
            err=float(np.linalg.norm(z[take]-manual)/np.linalg.norm(manual));assert err<1e-9;inverse_errors[window]=err
            for gate,row in r['metrics'][window]['source_window'].items():
                lo,hi=row['gate_ns'];mask=(t*1e9>=lo)&(t*1e9<=hi);base=z[mask,1]
                for name in ['ricker_800','impulse_1600','ricker_1600']:
                    x=z[mask,ix[name]];norm=np.linalg.norm(x-base)/np.linalg.norm(base)
                    assert abs(norm-row['models'][name]['change_L2_over_impulse800'])<1e-12
            for station,row in r['metrics'][window]['stations'].items():
                lo,hi=row['gate_ns'];mask=(t*1e9>=lo)&(t*1e9<=hi);j,k=ix[station+'_H1'],ix[station+'_H0'];delta=z[mask,j]-z[mask,k]
                assert abs(np.linalg.norm(z[mask,k])/np.linalg.norm(delta)-row['H0_over_difference_L2'])<1e-12
                assert abs(np.angle(np.vdot(z[mask,k],delta),deg=True)-row['H0_difference_phase_deg'])<1e-10
                if station!='c0008':assert abs(np.linalg.norm(z[mask,ix[station+'_far']])/np.linalg.norm(z[mask,k])-row['far_remaining_L2_over_H0'])<1e-12
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:
        if (a.public/name).exists():assert sha(a.public/name)==sha(a.private/name)
        else:shutil.copyfile(a.private/name,a.public/name)
    result=dict(status='PASS_NATIVE_EVENTS_DIRECT_CONVOLUTION_AND_COMPLEX_METRICS_NOT_ENTIRE_LINE_VALIDATION',audit_script_sha256=sha(__file__),analysis_script_sha256=r['script_sha256'],contract_sha256=r['contract_sha256'],verification_sha256=r['verification_sha256'],
        solver_started=9,solver_completed=9,independent_direct_convolution_relative_L2=direct_convolution,independent_inverse_relative_L2=inverse_errors,
        completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30,native_sha256=e['raw_sha256']) for e in done],plot_sha256={p.name:sha(p) for p in a.public.glob('*.png')})
    (a.public/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],convolution=direct_convolution,inverse=inverse_errors)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','public','prior','prepared']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
