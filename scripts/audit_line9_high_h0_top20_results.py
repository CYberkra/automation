"""Independent source-normalized DFT and direct inverse audit of one top-control."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from audit_line9_v401_delivery import native_response,sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    pub=read(a.public/'analysis.json');c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==pub['contract_sha256']==sha(a.source/'execution_contract.json')
    assert pub['verification_sha256']==sha(a.source/'completed_verification.json')
    assert sha(a.numerical)==pub['numerical_sha256']
    paths=[a.package/'baseline/native_H0.h5',a.source/'profile.h5']
    spectra=[native_response(p,.025,h) for p,h in zip(paths,pub['native_sha256'])]
    f=20e6+np.arange(501)*300000.
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);z=h['response'][:]
    assert z.shape==(501,3)
    errors=[float(np.linalg.norm(q-z[:,j])/np.linalg.norm(q)) for j,q in enumerate(spectra)]
    assert max(errors)<1e-9
    assert (z[:,1]-z[:,0]).tobytes()==z[:,2].tobytes()
    q=np.column_stack(spectra);t=np.arange(4008)/(4008*300000.);metrics={}
    for name,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();rows={}
        for gate,bounds in c['study_manifest']['gates_ns'].items():
            k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);x=np.exp(2j*np.pi*t[k,None]*f)@(q*w[:,None])/501
            value=float(np.linalg.norm(x[:,1]-x[:,0])/np.linalg.norm(x[:,0]))
            peaks=[float(t[k][abs(x[:,j]).argmax()]*1e9) for j in range(2)]
            r=pub['metrics'][name]['gates'][gate]
            assert abs(value-r['relative_complex_L2'])<1e-8
            np.testing.assert_allclose(peaks,r['peak_ns'],atol=1e-9,rtol=0)
            rows[gate]=dict(relative_complex_L2=value,peak_ns=peaks)
        metrics[name]=rows
    with h5py.File(paths[0]) as old,h5py.File(paths[1]) as new:
        s=old['srcs/src1/excitation/samples'][:];r=new['srcs/src1/excitation/samples'][:]
        assert s.dtype==r.dtype==np.float64 and s.tobytes()==r.tobytes()
        x,y=old['rxs/rx1/Ez'][:],new['rxs/rx1/Ez'][:]
        dt=c['study_manifest']['dt_s'];times=np.arange(len(x))*dt+old['rxs/rx1/Ez'].attrs['TimeSampleOffset']
        for gate,bounds in dict(early=[0,120],return_packet=[250,450],late=[450,1100]).items():
            k=(times*1e9>=bounds[0])&(times*1e9<=bounds[1]);value=float(np.linalg.norm(y[k]-x[k])/np.linalg.norm(x[k]))
            assert abs(value-pub['native_metrics'][gate]['relative_L2'])<1e-12
    frequency_change=float(np.linalg.norm(q[:,1]-q[:,0])/np.linalg.norm(q[:,0]))
    assert abs(frequency_change-pub['frequency_relative_complex_L2'])<1e-8
    result=dict(status='PASS_INDEPENDENT_ONE_TOP20_NATIVE_DFT_DIRECT_INVERSE_METRICS',audit_script_sha256=sha(__file__),
        analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),
        native_sha256=pub['native_sha256'],independent_DFT_relative_L2=errors,direct_inverse_metrics=metrics,
        frequency_relative_complex_L2=frequency_change,source_bitwise_equal=True,
        limits='Numerical processing agreement only. One top extension does not certify all PML or uniquely attribute cover paths.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','public','package','numerical','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
