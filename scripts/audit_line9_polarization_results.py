"""Independent native/DFT/direct inverse audit for the 2x2 constitutive controls."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from audit_line9_v401_delivery import native_response,read,sha
from audit_line9_polarization_inputs import audit as input_audit


def main(a):
    assert not a.out.exists()
    inputs=input_audit(a.package,a.parent)
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    m=read(a.package/'manifest.json');pub=read(a.public/'analysis.json')
    assert c['study_manifest']==m and v['completed']
    assert v['contract_sha256']==pub['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.source/'sfcw.h5')==pub['numerical_sha256']
    paths={r['id']:(a.package/r['path'],r['native_sha256']) for r in m['reused']}
    groups={g['id']:g for g in m['groups']}
    for r in v['groups']:paths[r['id']]=(a.source/(r['id']+'.h5'),r['native_sha256'])
    names=['span3_dc003','span3_dc0003','span03_dc003','span03_dc0003']
    vals=[]
    for name in names:
        pair=[]
        for role in ['H0','H1']:
            label=name+'_'+role;p,digest=paths[label]
            with h5py.File(p) as h:
                assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==20352
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
                source=h['srcs/src1/excitation/samples'][:]
                if not vals and not pair:reference_source=source
                else:np.testing.assert_array_equal(source,reference_source)
                for key,pos in [('srcs/src1',[177.95,39.375,0]),('rxs/rx1',[179.25,39.425,0])]:
                    np.testing.assert_array_equal(h[key].attrs['GridPosition'],np.rint(np.array(pos)/.025))
                    np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
            pair.append(native_response(p,.025,digest))
        vals.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(vals);f=20e6+np.arange(501)*300000.
    errors=[]
    with h5py.File(a.source/'sfcw.h5') as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f)
        stored=h['response'][:]
        for j in range(12):
            err=float(np.linalg.norm(stored[:,j]-z[:,j])/np.linalg.norm(z[:,j]));assert err<1e-9;errors.append(err)
    t=np.arange(4008)/(4008*300000.);checks={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        checks[window]={};w=w/w.mean()
        for gate in ['basal','deep']:
            lo,hi=m[gate+'_gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi)
            x=np.exp(2j*np.pi*t[keep,None]*f)@(z*w[:,None])/501
            expected=pub['metrics'][window][gate];norms=[];worst=0.
            for j,name in enumerate(names):
                h0,h1,d=x[:,3*j:3*j+3].T;n0=np.linalg.norm(h0);n1=np.linalg.norm(h1);nd=np.linalg.norm(d);norms.append(nd)
                q=dict(H0_norm=n0,H1_norm=n1,delta_norm=nd,H0_over_delta=n0/nd,
                    H1_vs_delta_correlation=abs(np.vdot(h1,d))/(n1*nd),
                    delta_vs_original_correlation=abs(np.vdot(x[:,2],d))/(np.linalg.norm(x[:,2])*nd),
                    delta_peak_time_ns=t[keep][np.argmax(abs(d))]*1e9)
                for k,value in q.items():
                    err=abs(value-expected['cases'][name][k]);assert err<1e-9;worst=max(worst,float(err))
            d=norms;q=dict(DC_reduction_at_span3=d[1]/d[0],DC_reduction_at_span03=d[3]/d[2],
                span_reduction_at_DC003=d[2]/d[0],span_reduction_at_DC0003=d[3]/d[1],
                both_reductions_over_original=d[3]/d[0],multiplicative_interaction=d[3]*d[0]/(d[2]*d[1]))
            for k,value in q.items():
                err=abs(value-expected['effects'][k]);assert err<1e-8;worst=max(worst,float(err))
            checks[window][gate]=dict(max_absolute_metric_error=worst,effects=q)
    # Scalar real-valued absorption check, independent from complex square-root analysis.
    for row in pub['constitutive_budget']:
        ei=row['epsilon_infinity'];de=row['delta_epsilon'];tau=row['tau_s'];dc=row['sigma_DC_S_m']
        f3=np.array([20,95,170])*1e6;u=2*np.pi*f3*tau
        er=ei+de/(1+u*u);loss=de*u/(1+u*u)+dc/(2*np.pi*f3*8.8541878128e-12)
        alpha=2*np.pi*f3/299792458*np.sqrt((np.hypot(er,loss)-er)/2)
        np.testing.assert_allclose(-20*alpha*14/np.log(10),row['absorption_7m_two_way_dB_20_95_170'],rtol=1e-10,atol=1e-10)
    result=dict(status='PASS_INDEPENDENT_FACTORIAL_NATIVE_MATERIAL_VOXEL_DFT_INVERSE',script_sha256=sha(__file__),
        input_audit=inputs,analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),
        independent_DFT_relative_L2=errors,direct_inverse_checks=checks,new_native=4,reused_native=4,
        limits='Audits provenance and numerical processing; not FDTD error-floor certification or field material validation.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','parent','public','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
