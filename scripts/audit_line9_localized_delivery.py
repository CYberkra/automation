"""Independent completion/native/complex-metric audit of localized and height controls."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(private,public,prior):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(private/'execution_contract.json');v=read(private/'completed_verification.json');a=read(public/'analysis.json')
    events=[json.loads(s) for s in (private/'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256']==v['contract_sha256']==a['contract_sha256']==sha(private/'execution_contract.json')
    assert events[-1]['status']=='COMPLETED' and events[-1]['verification_sha256']==a['verification_sha256']==sha(private/'completed_verification.json')
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e];done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    assert [e['group'] for e in starts]==[e['group'] for e in done]==[g['id'] for g in c['groups']]
    assert len(done)==c['max_runs']==events[-1]['traces']==4
    for g,row,event in zip(c['groups'],v['groups'],done):
        path=private/(g['id']+'.h5');assert sha(path)==row['native_sha256']==event['raw_sha256']
        with h5py.File(path) as h:
            assert str(h.attrs['gprMax'])=='4.0.0' and h['rxs/rx1/Ez'].dtype==np.float64 and h['srcs/src1/excitation/samples'].dtype==np.float64
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            np.testing.assert_allclose(h['srcs/src1'].attrs['Position'],g['tx_m'],atol=1e-12,rtol=0)
            np.testing.assert_allclose(h['rxs/rx1'].attrs['Position'],g['rx_m'],atol=1e-12,rtol=0)
            assert h['rxs/rx1/Ez'].shape==(13569,) and h.attrs['dt']==row['dt_s']
    assert a['script_sha256']==sha(Path(__file__).parent/'analyze_line9_localized_height.py')
    assert sha(private/'sfcw.h5')==a['numerical_sha256']
    old=read(prior/'completed_verification.json')
    for name,row in zip(['H1','H0'],old['groups']):assert sha(prior/(name+'.h5'))==row['raw_sha256']
    errors={}
    with h5py.File(private/'sfcw.h5') as h:
        f=h['frequency_Hz'][:];t=h['time_s'][:];r=h['response'][:]
        np.testing.assert_array_equal(f,20e6+np.arange(501)*300000.)
        assert h.attrs['columns'].split(',')==a['columns']
        for window in ['hann','blackman']:
            w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();z=h[window+'_complex_bandpass'][:]
            take=np.arange(0,len(t),103);expected=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*r)/501
            e=float(np.linalg.norm(z[take]-expected)/np.linalg.norm(expected));assert e<1e-9;errors[window]=e
            lo,hi=c['study_manifest']['predefined_gates_ns']['high_basal'];mask=(t*1e9>=lo)&(t*1e9<=hi);base=z[mask,1];bn=np.linalg.norm(base)
            for j,name in [(2,'far_interbed_removed'),(3,'near_interbed_removed')]:
                row=a['metrics'][window]['localized'][name]
                assert abs(np.linalg.norm(z[mask,j])/bn-row['remaining_L2_over_high_H0'])<1e-12
                assert abs(np.linalg.norm(z[mask,j]-base)/bn-row['change_L2_over_high_H0'])<1e-12
            for height,j in [('high',0),('low',4)]:
                row=a['metrics'][window][height];lo,hi=row['gate_ns'];mask=(t*1e9>=lo)&(t*1e9<=hi);delta=z[:,j]-z[:,j+1]
                assert abs(np.linalg.norm(z[mask,j+1])/np.linalg.norm(delta[mask])-row['H0_over_difference_L2'])<1e-12
                assert abs(np.linalg.norm(delta[mask])/np.linalg.norm(z[mask,j])-row['difference_over_total_L2'])<1e-12
                phase=float(np.angle(np.vdot(z[mask,j+1],delta[mask]),deg=True));assert abs(phase-row['H0_difference_inner_product_phase_deg'])<1e-10
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:
        if (public/name).exists():assert sha(public/name)==sha(private/name)
        else:shutil.copyfile(private/name,public/name)
    result=dict(status='PASS_NATIVE_COMPLETION_AND_INDEPENDENT_COMPLEX_METRICS_NOT_WHOLE_LINE_VALIDATION',audit_script_sha256=sha(__file__),
        analysis_script_sha256=a['script_sha256'],contract_sha256=a['contract_sha256'],verification_sha256=a['verification_sha256'],solver_started=4,solver_completed=4,
        inverse_errors=errors,completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30,native_sha256=e['raw_sha256']) for e in done],
        plot_sha256={p.name:sha(p) for p in public.glob('*.png')})
    (public/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':errors}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['private','public','prior']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.private,a.public,a.prior)
