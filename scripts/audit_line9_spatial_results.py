"""Independent native, geometry, direct-DFT, Fresnel/Snell and inverse audit."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import brentq
from audit_line9_v401_delivery import native_response,read,sha
from audit_line9_spatial_inputs import audit as input_audit


def main(a):
    assert not a.out.exists()
    inputs=input_audit(a.package,a.parent);m=read(a.package/'manifest.json');pub=read(a.public/'analysis.json')
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==pub['contract_sha256']==sha(a.source/'execution_contract.json')
    assert sha(a.source/'sfcw.h5')==pub['numerical_sha256']
    assert sha(a.protocol)==pub['protocol_sha256']
    groups={g['id']:g for g in m['groups']};records={g['id']:g for g in v['groups']};values=[];identities=[]
    first_source=None
    for row in pub['stations']:
        pair=[]
        for role in ['H0','H1']:
            if row['id']=='anchor':
                entry=next(r for r in m['reused'] if r['id']=='anchor_'+role);p=a.package/entry['path'];digest=entry['native_sha256'];g=entry['parent_group']
            else:
                label=row['id']+'_'+role;p=a.source/(label+'.h5');digest=records[label]['native_sha256'];g=groups[label]
            with h5py.File(p) as h:
                assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==20352 and abs(h.attrs['dt']/m['dt_s']-1)<1e-14
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                    np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-11,rtol=0)
                    np.testing.assert_array_equal(h[key].attrs['GridPosition'],np.rint(np.array(pos)/.025))
                source=h['srcs/src1/excitation/samples'][:]
                if first_source is None:first_source=source
                else:np.testing.assert_array_equal(source,first_source)
            pair.append(native_response(p,.025,digest));identities.append(dict(id=row['id']+'_'+role,native_sha256=digest))
        values.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(values);f=20e6+np.arange(501)*300000.;omega=2*np.pi*f;db=read(a.package/m['groups'][0]['material'])['materials'];n=[]
    for mat in db.values():
        b=mat['base'];er=b['relative_permittivity']-1j*b['electric_conductivity_s_per_m']/(omega*8.8541878128e-12)
        for pole in mat.get('poles',[]):er+=pole['relative_permittivity_difference']/(1+1j*omega*pole['relaxation_time_s'])
        n.append(np.sqrt(er))
    n=np.column_stack(n);trans=lambda i,j:2*n[:,i]/(n[:,i]+n[:,j]);refl=lambda i,j:(n[:,i]-n[:,j])/(n[:,i]+n[:,j])
    coefficient=trans(0,1)*trans(1,0)*trans(1,2)*trans(2,1)*refl(2,3);templates=[]
    for row in pub['stations']:
        g=row['geometry'];d1=g['cover_base']['depth_m'];d2=g['basal_sand']['depth_m']-d1
        lengths=np.array([g['midpoint_agl_m'],d1,d2,0.]);active=lengths>0;nr=n[250,active].real;ll=lengths[active]
        sep=abs(row['rx_m'][0]-row['tx_m'][0])
        root=brentq(lambda q:2*np.sum(ll*q/np.sqrt(nr*nr-q*q))-sep,0,float(nr.min())*(1-1e-12),xtol=1e-14)
        tau=2*np.sum(ll*nr*nr/np.sqrt(nr*nr-root*root))/299792458.
        correction=tau-2*np.dot(lengths,n[250].real)/299792458.
        templates.append(coefficient*np.exp(-2j*omega*(n@lengths)/299792458.-1j*omega*correction))
    templates=np.column_stack(templates)
    with h5py.File(a.source/'sfcw.h5') as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);stored=h['response'][:];stored_templates=h['templates'][:]
        errors=[float(np.linalg.norm(z[:,j]-stored[:,j])/np.linalg.norm(z[:,j])) for j in range(12)];assert max(errors)<1e-9
        template_error=float(np.linalg.norm(templates-stored_templates)/np.linalg.norm(templates));assert template_error<1e-10
        grid=h['sparse_chainage_m'][:];mask=h['completed_mask'][:];positions=h['chainage_m'][:]
        assert len(grid)==262 and mask.sum()==4;np.testing.assert_allclose(sorted(grid[mask]),sorted(positions),atol=1e-10,rtol=0)
    t=np.arange(4008)/(4008*300000.);checks={};corr=lambda x,y:float(abs(np.vdot(x,y))/(np.linalg.norm(x)*np.linalg.norm(y)))
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();checks[window]=[];model=np.exp(2j*np.pi*t[:,None]*f)@(templates*w[:,None])/501
        common=(t*1e9>=180)&(t*1e9<=450);common_x=np.exp(2j*np.pi*t[common,None]*f)@(z*w[:,None])/501
        actual_peaks=[];model_peaks=[]
        for j,row in enumerate(pub['stations']):
            want=pub['metrics'][window][j];peak=float(t[np.argmax(abs(model[:,j]))]*1e9);assert abs(peak-want['template_peak_ns'])<1e-9
            worst=0.;actual_peak=None
            for label,bounds in [('basal',row['basal_gate_ns']),('wide',row['wide_gate_ns'])]:
                k=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);x=np.exp(2j*np.pi*t[k,None]*f)@(z[:,3*j:3*j+3]*w[:,None])/501
                h0,h1,delta=x.T;ticks=t[k]*1e9;imax=int(np.argmax(abs(delta)));hmax=int(np.argmax(abs(h1)))
                computed=dict(H0_over_delta_L2=float(np.linalg.norm(h0)/np.linalg.norm(delta)),H1_vs_delta_complex_correlation=corr(h1,delta),
                    delta_vs_local_primary_complex_correlation=corr(delta,model[k,j]),H1_vs_local_primary_complex_correlation=corr(h1,model[k,j]),
                    delta_peak_ns=float(ticks[imax]),H1_peak_ns=float(ticks[hmax]),actual_minus_template_peak_ns=float(ticks[imax]-peak),H1_minus_template_peak_ns=float(ticks[hmax]-peak),
                    delta_L2=float(np.linalg.norm(delta)),H0_L2=float(np.linalg.norm(h0)))
                for key,value in computed.items():
                    scaled=abs(value-want['gates'][label][key])/max(1.,abs(want['gates'][label][key]));assert scaled<1e-8,(window,row['id'],label,key,scaled);worst=max(worst,scaled)
                assert bool(imax in [0,len(ticks)-1])==want['gates'][label]['peak_on_gate_edge']
                if label=='basal':actual_peak=float(ticks[imax])
            for offset,key in [(0,'common_180_450_H0_peak_ns'),(1,'common_180_450_H1_peak_ns'),(2,'common_180_450_delta_peak_ns')]:
                actual=float(t[common][np.argmax(abs(common_x[:,3*j+offset]))]*1e9);assert abs(actual-want[key])<1e-9
            actual_peaks.append(actual_peak);model_peaks.append(peak);checks[window].append(dict(chainage_m=row['chainage_m'],max_scaled_metric_error=worst))
        for j,want in enumerate(pub['metrics'][window]):
            error=(actual_peaks[j]-actual_peaks[0])-(model_peaks[j]-model_peaks[0])
            assert abs(error-want['delta_anchor_relative_time_change_minus_model_time_change_ns'])<1e-9
    result=dict(status='PASS_INDEPENDENT_NATIVE_VOXELS_DIRECT_DFT_SNELL_FRESNEL_DIRECT_INVERSE_SPARSE_MASK',script_sha256=sha(__file__),
        input_audit=inputs,analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),
        independent_DFT_relative_L2=errors,independent_template_relative_L2=template_error,direct_inverse_checks=checks,native=identities,
        limits='Data, approximate template mathematics and processing audit; not independent physical applicability, full line, finite3D or field calibration.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','parent','public','protocol','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
