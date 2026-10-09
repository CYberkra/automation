"""Describe complex change at the original H0 peak; not a new physical gate."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    r=read(a.public/'analysis.json');v=read(a.public/'independent_audit.json')
    assert v['analysis_sha256']==sha(a.public/'analysis.json') and v['status'].startswith('PASS')
    assert sha(a.numerical)==r['numerical_sha256']
    f=20e6+np.arange(501)*300000.;t=np.arange(4008)/(4008*300000.)
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);z=h['response'][:,:2]
    k=(t*1e9>=300)&(t*1e9<=450);metrics={}
    for name,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();x=np.exp(2j*np.pi*t[k,None]*f)@(z*w[:,None])/501
        j=abs(x[:,0]).argmax();p,q=x[j];peak=float(t[k][j]*1e9)
        assert peak==r['metrics'][name]['gates']['wide']['peak_ns'][0]
        metrics[name]=dict(original_peak_ns=peak,original_peak_amplitude=float(abs(p)),
            top20_amplitude_at_original_peak=float(abs(q)),
            amplitude_change_percent=float((abs(q)/abs(p)-1)*100),
            phase_change_at_original_peak_deg=float(np.angle(q/p,deg=True)),
            relative_complex_change_at_original_peak=float(abs(q-p)/abs(p)))
    result=dict(status='PASS_ORIGINAL_PEAK_COMPLEX_DIAGNOSTIC',script_sha256=sha(__file__),
        analysis_sha256=sha(a.public/'analysis.json'),independent_audit_sha256=sha(a.public/'independent_audit.json'),
        numerical_sha256=sha(a.numerical),metrics=metrics,
        limits='Peak-relative numerical diagnostic, no alignment/amplitude fit, energy attribution, sub-sample resolution or new acceptance threshold.')
    a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['public','numerical','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
