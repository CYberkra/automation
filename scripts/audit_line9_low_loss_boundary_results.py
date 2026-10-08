"""Independent low-loss expanded-domain native/DFT/inverse audit."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from audit_line9_v401_delivery import native_response,read,sha
from audit_line9_low_loss_boundary_inputs import audit as input_audit


def main(a):
    assert not a.out.exists()
    inputs=input_audit(a.package,a.parent);m=read(a.package/'manifest.json')
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');pub=read(a.public/'analysis.json')
    assert c['study_manifest']==m and v['completed']
    assert v['contract_sha256']==pub['contract_sha256']==sha(a.source/'execution_contract.json')
    assert pub['numerical_sha256']==sha(a.source/'sfcw.h5')
    paths={r['id']:(a.package/r['path'],r['native_sha256'],r['parent_group']) for r in m['reused']}
    groups={g['id']:g for g in m['groups']}
    for r in v['groups']:paths[r['id']]=(a.source/(r['id']+'.h5'),r['native_sha256'],groups[r['id']])
    names=['base','sides40','bottom20'];values=[];prefix=[];original={}
    for name in names:
        pair=[]
        for role in ['H0','H1']:
            label=name+'_'+role;p,digest,g=paths[label]
            with h5py.File(p) as h:
                assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==20352
                assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape'])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                    np.testing.assert_array_equal(h[key].attrs['GridPosition'],np.rint(np.array(pos)/.025))
                    np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
                source=h['srcs/src1/excitation/samples'][:];trace=h['rxs/rx1/Ez'][:]
                if not values and not pair:reference_source=source
                else:np.testing.assert_array_equal(source,reference_source)
                early=trace[np.arange(len(trace))*m['dt_s']<120e-9]
                if name=='base':original[role]=early
                else:prefix.append(dict(id=label,raw_0_120ns_bit_identical=bool(np.array_equal(early,original[role])),relative_L2=float(np.linalg.norm(early-original[role])/np.linalg.norm(original[role]))))
            pair.append(native_response(p,.025,digest))
        values.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(values);f=20e6+np.arange(501)*300000.;errors=[]
    with h5py.File(a.source/'sfcw.h5') as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);stored=h['response'][:]
        for j in range(9):
            err=float(np.linalg.norm(z[:,j]-stored[:,j])/np.linalg.norm(z[:,j]));assert err<1e-9;errors.append(err)
    t=np.arange(4008)/(4008*300000.);checks={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();checks[window]={}
        for gate,(lo,hi) in pub['gate_ns'].items():
            keep=(t*1e9>=lo)&(t*1e9<=hi)
            x=np.exp(2j*np.pi*t[keep,None]*f)@(z*w[:,None])/501;b=x[:,:3];worst=0.;worst_scaled=0.
            for j,name in enumerate(names):
                q=x[:,3*j:3*j+3];n=np.linalg.norm;want=pub['metrics'][window][gate][name]
                actual=dict(H0_over_delta=float(n(q[:,0])/n(q[:,2])),H1_vs_delta_correlation=float(abs(np.vdot(q[:,1],q[:,2]))/(n(q[:,1])*n(q[:,2]))),delta_norm_over_original=float(n(q[:,2])/n(b[:,2])),delta_peak_ns=float(t[keep][np.argmax(abs(q[:,2]))]*1e9))
                for key,value in actual.items():
                    err=abs(value-want[key]);scaled=err/max(1.,abs(want[key]))
                    # Early-window H0/delta can exceed 1e8: compare dimensionless
                    # metric error on its own scale, retaining absolute error.
                    assert scaled<1e-8,(gate,name,key,err,scaled)
                    worst=max(worst,err);worst_scaled=max(worst_scaled,scaled)
                for i,role in enumerate(['H0','H1','delta']):
                    change=n(q[:,i]-b[:,i])
                    for key,value in [('relative_complex_change',change/n(b[:,i])),('changes_over_original_delta',change/n(b[:,2]))]:
                        err=abs(float(value)-want[key][role]);scaled=err/max(1.,abs(want[key][role]))
                        assert scaled<1e-8,(gate,name,key,err,scaled)
                        worst=max(worst,err);worst_scaled=max(worst_scaled,scaled)
            checks[window][gate]=dict(max_absolute_metric_error=worst,max_scaled_metric_error=worst_scaled)
    result=dict(status='PASS_INDEPENDENT_LOW_LOSS_BOUNDARY_NATIVE_VOXEL_DFT_INVERSE',script_sha256=sha(__file__),
        input_audit=inputs,analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),
        independent_DFT_relative_L2=errors,direct_inverse_checks=checks,raw_early_prefix=prefix,
        limits='Processing and provenance audit, not global FDTD convergence or field/entire-line validation.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','parent','public','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
