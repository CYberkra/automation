"""Independently recompute event identities, transforms and control metrics; no solver."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(private, public):
    read=lambda p: json.loads(p.read_text('utf-8'))
    c=read(private/'execution_contract.json');v=read(private/'completed_verification.json');a=read(public/'analysis.json')
    events=[json.loads(s) for s in (private/'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256']==v['contract_sha256']==a['contract_sha256']==sha(private/'execution_contract.json')
    assert events[-1]['status']=='COMPLETED' and events[-1]['verification_sha256']==a['verification_sha256']==sha(private/'completed_verification.json')
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e]
    done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    names=[g['id'] for g in c['groups']]
    assert [e['group'] for e in starts]==[e['group'] for e in done]==[g['id'] for g in v['groups']]==names
    assert len(done)==c['max_runs']==events[-1]['traces']
    for g,row,event in zip(c['groups'],v['groups'],done):
        p=private/(g['id']+'.h5')
        assert sha(p)==row['native_sha256']==event['raw_sha256']
        with h5py.File(p) as h:
            assert str(h.attrs['gprMax'])=='4.0.0' and h['rxs/rx1/Ez'].dtype==np.float64
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape'])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            np.testing.assert_allclose(h['srcs/src1'].attrs['Position'],g['tx_m'],rtol=0,atol=1e-12)
            np.testing.assert_allclose(h['rxs/rx1'].attrs['Position'],g['rx_m'],rtol=0,atol=1e-12)
            assert h['rxs/rx1/Ez'].shape==(row['samples'],) and h.attrs['dt']==row['dt_s']
    assert a['script_sha256']==sha(Path(__file__).parent/'analyze_line9_boundary_controls.py')
    numerical=private/'sfcw.h5';assert sha(numerical)==a['numerical_sha256']
    recomputed={}
    with h5py.File(numerical) as h:
        f=h['frequency_Hz'][:];t=h['time_s'][:];r=h['response'][:];columns=h.attrs['columns'].split(',')
        np.testing.assert_array_equal(f,20e6+300000*np.arange(501))
        assert columns==(['base']+names if a['reference_is_prior_H0'] else names)
        for window in ['hann','blackman']:
            w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean()
            z=h[window+'_complex_bandpass'][:];take=np.arange(0,len(t),89)
            expected=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*r)/501
            err=float(np.linalg.norm(z[take]-expected)/np.linalg.norm(expected));assert err<1e-9
            recomputed[window]={'inverse_relative_L2':err}
            for gate,(lo,hi) in a['gates_ns'].items():
                mask=(t*1e9>=lo)&(t*1e9<=hi);base=z[mask,0];bn=np.linalg.norm(base)
                for j,name in enumerate(columns):
                    x=z[mask,j];row=a['metrics'][window][gate][name]
                    vals=dict(L2_over_base=float(np.linalg.norm(x)/bn),change_L2_over_base=float(np.linalg.norm(x-base)/bn),
                              complex_correlation_with_base=float(abs(np.vdot(base,x))/(bn*np.linalg.norm(x))),
                              envelope_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(x))]]*1e9))
                    for key,val in vals.items():assert abs(val-row[key])<1e-10,(key,val,row[key])
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:
        if (public/name).exists():assert sha(public/name)==sha(private/name)
        else:shutil.copyfile(private/name,public/name)
    result=dict(status='PASS_NATIVE_EVENTS_AND_INDEPENDENT_TRANSFORM_METRICS_NOT_WHOLE_LINE_VALIDATION',audit_script_sha256=sha(__file__),
        contract_sha256=sha(private/'execution_contract.json'),verification_sha256=sha(private/'completed_verification.json'),
        solver_started=len(starts),solver_completed=len(done),checks=recomputed,
        completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30,native_sha256=e['raw_sha256']) for e in done],
        plot_sha256={p.name:sha(p) for p in public.glob('*.png')})
    (public/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'completed':len(done),'checks':recomputed}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--private',type=Path,required=True);p.add_argument('--public',type=Path,required=True)
    a=p.parse_args();main(a.private,a.public)
