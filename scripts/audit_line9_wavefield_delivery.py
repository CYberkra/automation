"""Independent receipt, source, event and complex-metric audit for the movie pair."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.private/'execution_contract.json');v=read(a.private/'completed_verification.json');r=read(a.public/'analysis.json');movie=read(a.public/'render_receipt.json')
    events=[json.loads(s) for s in (a.private/'execution.jsonl').read_text('utf-8').splitlines()]
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e];done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    assert [e['group'] for e in starts]==[e['group'] for e in done]==['H0','far_removed']
    assert events[-1]['status']=='COMPLETED' and events[-1]['verification_sha256']==sha(a.private/'completed_verification.json')==r['verification_sha256']
    assert events[0]['contract_sha256']==v['contract_sha256']==sha(a.private/'execution_contract.json')==r['contract_sha256']==movie['contract_sha256']
    assert sha(Path(__file__).parent/'analyze_line9_interbed_wavefield_delivery.py')==r['script_sha256']
    assert sha(Path(__file__).parent/'render_line9_interbed_wavefield.py')==movie['script_sha256']
    for path,digest in c['code_identities'].items():assert sha(Path(__file__).parent/Path(path).name)==digest
    assert sha(a.private/'sfcw.h5')==r['numerical_sha256']
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:assert sha(a.private/name)==sha(a.public/name)
    natives=[];raw_hashes=[]
    for g,row,e in zip(c['groups'],v['groups'],done):
        path=a.private/(g['id']+'.h5');digest=sha(path)
        assert digest==row['native_sha256']==e['raw_sha256'];raw_hashes.append(digest)
        with h5py.File(path) as h:
            dt=float(h.attrs['dt']);x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
            assert str(h.attrs['gprMax'])=='4.0.0' and x.dtype==s.dtype==np.float64 and x.shape==s.shape==(13569,)
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[4400,1600,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            for key,pos in [('srcs/src1',[78.6,35,0]),('rxs/rx1',[79.9,35,0])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-12,rtol=0)
            t=(np.arange(len(s))+.5)*dt;u=(np.pi*100e6*(t-np.sqrt(2)/100e6))**2
            np.testing.assert_allclose(s,40*(1-2*u)*np.exp(-u),atol=2e-12,rtol=1e-10)
            natives.append(x)
        assert len(row['snapshots'])==799
        np.testing.assert_array_equal([x['iteration'] for x in row['snapshots']],list(range(0,13569,17)))
    assert raw_hashes[0]==sha(a.reference)
    assert sha(a.public/'interbed_wavefield.gif')==movie['gif_sha256']
    receipt={(g['id'],x['iteration']):x['sha256'] for g in v['groups'] for x in g['snapshots']}
    for x in movie['rendered_snapshot_hashes']:assert x['sha256']==receipt[(x['group'],x['iteration'])]
    for x in movie['static_frames']:assert sha(a.public/x['file'])==x['sha256']
    errors={}
    with h5py.File(a.private/'sfcw.h5') as h:
        f=h['frequency_Hz'][:];t=h['time_s'][:];response=h['response'][:]
        np.testing.assert_array_equal(f,20e6+np.arange(501)*300000)
        np.testing.assert_array_equal(response[:,2],response[:,0]-response[:,1])
        gate=(t*1e9>=332.31137724550894)&(t*1e9<=356.31137724550894)
        for window in ['hann','blackman']:
            z=h[window+'_complex_bandpass'][:];w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();take=np.arange(3,len(t),131)
            manual=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*response)/501
            error=float(np.linalg.norm(z[take]-manual)/np.linalg.norm(manual));assert error<1e-9;errors[window]=error
            value=float(np.linalg.norm(z[gate,1])/np.linalg.norm(z[gate,0]))
            assert abs(value-r['metrics'][window]['far_remaining_L2_over_H0'])<1e-12
    result=dict(status='PASS_EVENTS_NATIVE_SOURCE_HASHES_MOVIE_RECEIPTS_AND_INDEPENDENT_INVERSE',audit_script_sha256=sha(__file__),independent_inverse_relative_L2=errors,
        solver_started=2,solver_completed=2,native_sha256=raw_hashes,total_six_component_snapshots=1598,rendered_frames=movie['rendered_frames'],
        full_snapshot_payload_audit='All six fields were checked for dtype/shape/finiteness and hashed on solver host. Local audit independently checks receipt identities, not a second download of9.3GB.',
        completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30) for e in done],
        public_media_sha256={p.name:sha(p) for p in a.public.iterdir() if p.suffix in ['.png','.gif']})
    (a.public/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],inverse=errors)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','public','reference']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
