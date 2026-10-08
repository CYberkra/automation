"""Independent native DFT/inverse/availability audit; no production SFCW imports."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from audit_line9_v401_delivery import native_response,sha,read
from audit_line9_dense_loss_inputs import audit as audit_inputs

def main(a):
    assert not a.out.exists()
    inputs=audit_inputs(a.package);p=read(a.public/'analysis.json');m=read(a.package/'manifest.json');v=read(a.source/'snapshot_manifest.json');c=read(a.source/'execution_contract.json')
    assert inputs['bulk_execution_ready'] and c['study_manifest']==m
    assert p['contract_sha256']==v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert p['manifest_sha256']==sha(a.package/'manifest.json') and p['snapshot_sha256']==sha(a.source/'snapshot_manifest.json')
    assert p['numerical_sha256']==sha(a.numerical)
    records={r['id']:r for r in v['records']+m['reused']};items={r['id']:r for r in m['groups']+m['reused']}
    assert len(v['records'])==v['completed_new'] and len(records)==v['completed_new']+2
    if v['status']=='COMPLETE':
        complete=read(a.source/'completed_verification.json');assert v['completed_new']==46 and complete['completed'] and complete['contract_sha256']==v['contract_sha256']
        assert {r['id']:r['native_sha256'] for r in complete['groups']}=={r['id']:r['native_sha256'] for r in v['records']}
    f=20e6+np.arange(501)*300000.;t=np.arange(4008)/(4008*300000.);direct={};errors=[];source=None
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);stored=h['response'][:];ids=[x.decode() if isinstance(x,bytes) else x for x in h['ids'][:]]
    assert set(ids)==set(records) and len(ids)==len(set(ids))
    for j,sid in enumerate(ids):
        path=a.source/(sid+'.h5');g=items[sid]
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1';np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert h.attrs['Iterations']==20352 and abs(h.attrs['dt']/m['dt_s']-1)<1e-14
            for node,role in [('srcs/src1','tx'),('rxs/rx1','rx')]:
                np.testing.assert_allclose(h[node].attrs['Position'],g[role+'_m'],rtol=0,atol=1e-11)
                np.testing.assert_array_equal(h[node].attrs['GridPosition'],np.rint(np.array(g[role+'_m'])/.025))
            s=h['srcs/src1/excitation/samples'][:]
            if source is None:source=s
            else:np.testing.assert_array_equal(source,s)
        z=native_response(path,.025,records[sid]['native_sha256']);direct[sid]=z
        err=float(np.linalg.norm(stored[:,j]-z)/np.linalg.norm(z));assert err<1e-9;errors.append(err)
    availability={v:[f'{v}_x{round(x*100):05d}_H1' in direct for x in m['chainage_m']] for v in ['high','low']}
    assert availability==p['availability']
    checks={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();worst=0.;nmetrics=0
        for variant in ['high','low']:
            peakrecords={r['chainage_m']:r for r in p['H1_peaks'][window][variant]};anchors={r['chainage_m']:r for r in p['anchor_metrics'][window][variant]}
            for s in m['stations']:
                x=s['chainage_m'];prefix=f'{variant}_x{round(x*100):05d}_';sid=prefix+'H1';h0=prefix+'H0'
                if sid not in direct:assert x not in peakrecords and x not in anchors;continue
                lo,hi=s['templates'][variant]['basal_gate_ns'];k=(t*1e9>=lo)&(t*1e9<=hi);e=np.exp(2j*np.pi*t[k,None]*f)
                q=e@(direct[sid]*w)/501;peak=float(t[k][np.argmax(abs(q))]*1e9)
                assert abs(peak-peakrecords[x]['H1_peak_ns'])<1e-8
                if h0 not in direct:assert x not in anchors;continue
                b=e@(direct[h0]*w)/501;d=q-b;norm=np.linalg.norm
                values=dict(H0_over_delta=float(norm(b)/norm(d)),H1_vs_delta_complex_correlation=float(abs(np.vdot(q,d))/(norm(q)*norm(d))),delta_peak_ns=float(t[k][np.argmax(abs(d))]*1e9),H1_peak_ns=peak,delta_norm=float(norm(d)),H0_norm=float(norm(b)),H1_norm=float(norm(q)))
                for name,value in values.items():
                    err=abs(value-anchors[x][name])/max(1,abs(anchors[x][name]));assert err<1e-8,(window,variant,x,name,err);worst=max(worst,err);nmetrics+=1
        checks[window]=dict(anchor_scalar_checks=nmetrics,max_scaled_metric_error=worst)
    result=dict(status='PASS_INDEPENDENT_DENSE_NATIVE_DFT_INVERSE_MASKS',auditor_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),input_audit=inputs['status'],manifest_sha256=sha(a.package/'manifest.json'),contract_sha256=sha(a.source/'execution_contract.json'),completed_new=v['completed_new'],reused=2,independent_DFT_relative_L2=errors,checks=checks,availability=availability,limits='Identity/arithmetic/missing-column checks only; not complete FDTD error bound,continuous H0,unique propagation mechanism or site validation.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],completed_new=v['completed_new'],max_DFT_error=max(errors),checks=checks)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','public','numerical','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
