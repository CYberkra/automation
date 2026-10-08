"""Independent native-metadata, direct-DFT and direct-inverse delivery audit.

No imports of the production SFCW, normalization, inverse, or analysis helpers.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text('utf-8'))


def native_response(p, scale, expected_hash):
    assert sha(p)==expected_hash
    with h5py.File(p) as h:
        assert h['rxs/rx1/Ez'].dtype==h['srcs/src1/excitation/samples'].dtype==np.float64
        x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
        assert x.shape==s.shape and len(x)==h.attrs['Iterations'] and np.isfinite(x).all()
        dt=float(h.attrs['dt']);ea=h['srcs/src1/excitation'].attrs;ra=h['rxs/rx1/Ez'].attrs
        assert ea['DrivingQuantity']=='electric_current' and ea['Polarisation']=='z'
        assert ea['WaveformAmplitude']==40 and ea['WaveformFrequency']==100e6
        assert ea['SampleInterval']==ra['SampleInterval']==dt
        ts=np.arange(len(s))*dt+ea['TimeSampleOffset'];tr=np.arange(len(x))*dt+ra['TimeSampleOffset']
    f=20e6+np.arange(501)*300000.;z=[]
    for k in range(0,501,13):
        freq=f[k:k+13,None]
        z.extend(((np.exp(-2j*np.pi*freq*tr)@x)/(np.exp(-2j*np.pi*freq*ts)@s))/scale)
    return np.array(z)


def main(a):
    assert not a.out.exists()
    pub=read(a.public/'analysis.json');ver=read(a.source/'completed_verification.json')
    assert ver['contract_sha256']==sha(a.source/'execution_contract.json')
    c=read(a.source/'execution_contract.json');f=20e6+np.arange(501)*300000.
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f)
        if a.mode=='version':z=h['v401_response'][:];base=h['v400_response'][:]
        else:z=h['response'][:]
    assert sha(a.numerical)==pub['numerical_sha256']
    identities=[];errors=[];native=[]
    for j,(g,v) in enumerate(zip(c['groups'],ver['groups'])):
        p=a.source/(g['id']+'.h5');q=native_response(p,.025,v['native_sha256']);native.append(q)
        error=float(np.linalg.norm(q-z[:,j])/np.linalg.norm(q));assert error<1e-9;errors.append(error)
        identities.append(dict(id=g['id'],native_sha256=sha(p)))
    if a.mode=='version':
        np.testing.assert_allclose(z[:,3],z[:,2]-z[:,1],rtol=0,atol=0)
        oldver=read(a.old/'completed_verification.json');old=[]
        for g in oldver['groups']:old.append(native_response(a.old/(g['id']+'.h5'),.025,g['native_sha256']))
        old=np.column_stack(old);old=np.column_stack([old,old[:,2]-old[:,1]])
        assert np.linalg.norm(old-base)/np.linalg.norm(old)<1e-9
        actual=np.column_stack(native);actual=np.column_stack([actual,actual[:,2]-actual[:,1]])
        reference=old
    else:
        actual=np.column_stack(native);actual=np.column_stack([actual,actual[:,1]-actual[:,0]])
        assert np.linalg.norm(actual-z)/np.linalg.norm(actual)<1e-9
        with h5py.File(a.reference) as h:reference=h['v401_response'][:,[1,2,3]]
    t=np.arange(4008)/(4008*300000.);lo,hi=c['study_manifest']['basal_gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi)
    metrics={}
    for name,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();e=np.exp(2j*np.pi*t[keep,None]*f);x=e@(actual*w[:,None])/501;r=e@(reference*w[:,None])/501
        norm=lambda value:np.linalg.norm(value)
        if a.mode=='version':
            changes=[float(norm(x[:,j]-r[:,j])/norm(r[:,j])) for j in range(4)]
            recorded=list(pub['metrics'][name]['relative_version_changes']['basal'].values())
            assert max(abs(np.array(changes)-recorded))<1e-9
            metrics[name]=dict(version_relative_changes=changes)
        else:
            changes=[float(norm(x[:,j]-r[:,j])/norm(r[:,j])) for j in [0,2]]
            recorded=pub['metrics'][name]
            assert abs(changes[0]-recorded['basal_H0_change_vs_wide'])<1e-9
            assert abs(changes[1]-recorded['basal_delta_change_vs_wide'])<1e-9
            metrics[name]=dict(H0_change_vs_wide=changes[0],delta_change_vs_wide=changes[1],
                               H0_over_delta=float(norm(x[:,0])/norm(x[:,2])),
                               H0_vs_delta_phase_deg=float(np.angle(np.vdot(x[:,2],x[:,0]),deg=True)))
    result=dict(status='PASS_INDEPENDENT_NATIVE_DFT_AND_DIRECT_INVERSE',audit_script_sha256=sha(__file__),
                analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),
                identities=identities,independent_DFT_relative_L2=errors,direct_inverse_metrics=metrics,
                limits='Audits data/processing and stated numerical comparisons; not field calibration or full numerical/physical certification.')
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['version','mesh'],required=True)
    for key in ['source','public','numerical','out','reference','old']:p.add_argument('--'+key,type=Path,required=key not in ['reference','old'])
    main(p.parse_args())
