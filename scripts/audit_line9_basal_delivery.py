"""Independent delivery checks of native event identity and complex recomposition; no solver."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(private,analysis,wavefield,numerical=None):
    def read(p):return json.loads(p.read_text('utf-8'))
    original=private/'original_execution_contract.json'
    recovery=private/'recovery_execution_contract.json'
    paired=private/'paired_verification_contract.json'
    c=read(recovery);v=read(private/'completed_verification.json');r=read(analysis/'analysis.json');w=read(wavefield/'wavefield_analysis.json')
    old_events=[json.loads(x) for x in (private/'original_execution.jsonl').read_text('utf-8').splitlines()]
    new_events=[json.loads(x) for x in (private/'recovery_execution.jsonl').read_text('utf-8').splitlines()]
    assert old_events[-1]==dict(status='FAILED',error='SESSION_LEASE_EXPIRED')
    assert new_events[-1]['status']=='COMPLETED'
    assert new_events[-1]['verification_sha256']==sha(private/'completed_verification.json')
    assert v['contract_sha256']==sha(paired)
    completed=[e for e in old_events+new_events if e['status']=='COMPLETED' and 'group' in e]
    attempts=[e for e in old_events+new_events if e['status']=='STARTED' and 'group' in e]
    assert [e['group'] for e in completed]==['H1','H0'] and len(attempts)==3
    for name,row,event in zip(['H1','H0'],v['groups'],completed):
        assert row['raw_sha256']==event['raw_sha256']==sha(private/(name+'.h5'))
        assert row['dtype']=='float64' and row['snapshot_count']==256
    numerical=numerical or private/'paired_sfcw.h5'
    assert sha(numerical)==r['private_numerical_sha256']
    expected_script=sha(Path(__file__).parent/'analyze_line9_basal_pair.py')
    assert r['script_sha256']==w['script_sha256']==expected_script
    assert w['verification_sha256']==r['verification_sha256']==sha(private/'completed_verification.json')
    errors={}
    with h5py.File(numerical) as h:
        f=h['frequency_Hz'][:];t=h['time_s'][:];response=h['response'][:]
        np.testing.assert_array_equal(f,20e6+np.arange(501)*300000.)
        for window in ['hann','blackman']:
            weights=np.hanning(501) if window=='hann' else np.blackman(501)
            weights/=weights.mean()
            take=np.arange(0,len(t),101)
            expected=np.exp(2j*np.pi*t[take,None]*f)@(weights[:,None]*response)/501
            z=h[window+'_complex_bandpass'][:]
            err=float(np.linalg.norm(z[take]-expected)/np.linalg.norm(expected));assert err<1e-9
            gate=abs(t*1e9-r['windows'][window]['basal_sand']['primary_peak_ns'])<=12
            delta=z[:,0]-z[:,1]
            metrics=r['windows'][window]['basal_sand']
            norm=float(np.linalg.norm(delta[gate])/np.linalg.norm(z[gate,0]))
            phase=float(np.angle(np.vdot(z[gate,1],delta[gate]),deg=True))
            assert abs(norm-metrics['difference_over_H1_L2'])<1e-13 and abs(phase-metrics['H0_difference_inner_product_phase_deg'])<1e-10
            errors[window]=dict(independent_direct_inverse_relative_L2=err,complex_recomposition_relative_L2=float(np.linalg.norm(z[:,0]-(z[:,1]+delta))/np.linalg.norm(z[:,0])))
    result=dict(status='PASS_NATIVE_IDENTITY_AND_INDEPENDENT_ARRAY_CHECKS_NOT_FULL_PHYSICAL_VALIDATION',
        audit_script_sha256=sha(__file__),analysis_script_sha256=expected_script,
        original_contract_sha256=sha(original),recovery_contract_sha256=sha(recovery),paired_verification_contract_sha256=sha(paired),
        native_verification_sha256=sha(private/'completed_verification.json'),solver_attempts_started=len(attempts),
        completed_unique_models=2,interrupted_attempts_preserved=1,completed_solver_elapsed_s=[e['elapsed_s'] for e in completed],
        peak_owned_RSS_GiB=[e['peak_owned_RSS_bytes']/2**30 for e in completed],hardware=c['hardware_at_freeze'],
        original_session_lease_s=120,recovery_session_lease_s=c['max_lease_age_s'],native_dtype='float64',
        native_samples_per_trace=13569,dt_s=v['groups'][0]['dt_s'],snapshots_per_model=256,snapshot_components=['Ex','Ey','Ez','Hx','Hy','Hz'],
        checks=errors,raw_native_sha256=[row['raw_sha256'] for row in v['groups']],
        artifact_sha256={p.name:sha(p) for p in sorted(analysis.glob('*.png'))},
        wavefield_artifact_sha256={p.name:sha(p) for p in sorted(wavefield.iterdir()) if p.is_file()},
        private_raw_location='ROG E:/line9_basal_pair_20261008_r1; full snapshots retained there, not uploaded to Git',
        limitations='Single selected station, material contrast includes lower PML continuation; no independent PML amplitude ablation, no full-line or measured-data validation.')
    (analysis/'execution_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],checks=errors),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','analysis','wavefield']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--numerical',type=Path)
    a=p.parse_args();main(a.private,a.analysis,a.wavefield,a.numerical)
