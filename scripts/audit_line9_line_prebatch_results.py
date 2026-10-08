"""Independent native/directDFT/direct-inverse arithmetic audit of220m controls."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from audit_line9_v401_delivery import native_response,sha,read
from audit_line9_line_prebatch_inputs import main as audit_inputs

def main(a):
    assert not a.out.exists();audit_inputs(argparse.Namespace(package=a.package,out=a.public/'independent_reaudit_inputs.json'))
    m=read(a.package/'manifest.json');c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');p=read(a.public/'analysis.json')
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==p['contract_sha256']==sha(a.source/'execution_contract.json')
    assert p['numerical_sha256']==sha(a.source/'sfcw.h5')
    groups={g['id']:g for g in m['groups']};records={r['id']:r for r in v['groups']};vals=[];names=['base','top20','right40'];source=None
    for name in names:
        pair=[]
        for role in ['H0','H1']:
            sid=name+'_'+role;path=a.source/(sid+'.h5');g=groups[sid]
            with h5py.File(path) as h:
                assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==20352
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
                for node,role2 in [('srcs/src1','tx'),('rxs/rx1','rx')]:
                    np.testing.assert_allclose(h[node].attrs['Position'],g[role2+'_m'],atol=1e-11,rtol=0)
                    np.testing.assert_array_equal(h[node].attrs['GridPosition'],np.rint(np.array(g[role2+'_m'])/.025))
                q=h['srcs/src1/excitation/samples'][:]
                if source is None:source=q
                else:np.testing.assert_array_equal(source,q)
            pair.append(native_response(path,.025,records[sid]['native_sha256']))
        vals.extend([*pair,pair[1]-pair[0]])
    z=np.column_stack(vals);f=20e6+np.arange(501)*300000.;errors=[]
    with h5py.File(a.source/'sfcw.h5') as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f);stored=h['response'][:]
        for j in range(9):
            err=float(np.linalg.norm(stored[:,j]-z[:,j])/np.linalg.norm(z[:,j]));assert err<1e-9;errors.append(err)
    t=np.arange(4008)/(4008*300000.);checks={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();checks[window]={}
        for gate,(lo,hi) in m['gates_ns'].items():
            k=(t*1e9>=lo)&(t*1e9<=hi);x=np.exp(2j*np.pi*t[k,None]*f)@(z*w[:,None])/501;b=x[:,:3];norm=np.linalg.norm;worst=0
            for j,name in enumerate(names):
                q=x[:,3*j:3*j+3];wanted=p['metrics'][window][gate][name]
                values=dict(H0_over_delta=float(norm(q[:,0])/norm(q[:,2])),H1_vs_delta_complex_correlation=float(abs(np.vdot(q[:,1],q[:,2]))/(norm(q[:,1])*norm(q[:,2]))),delta_peak_ns=float(t[k][np.argmax(abs(q[:,2]))]*1e9),H1_peak_ns=float(t[k][np.argmax(abs(q[:,1]))]*1e9))
                for key,value in values.items():
                    err=abs(value-wanted[key])/max(1,abs(wanted[key]));assert err<1e-8,(window,gate,name,key,err);worst=max(worst,err)
                for key,denom in [('relative_complex_change',None),('change_over_original_delta',norm(b[:,2]))]:
                    for i,r in enumerate(['H0','H1','delta']):
                        value=float(norm(q[:,i]-b[:,i])/(norm(b[:,i]) if denom is None else denom))
                        err=abs(value-wanted[key][r])/max(1,abs(wanted[key][r]));assert err<1e-8,(window,gate,name,key,r,err);worst=max(worst,err)
            checks[window][gate]=dict(max_scaled_metric_error=worst)
    result=dict(status='PASS_INDEPENDENT_NATIVE_DIRECT_DFT_DIRECT_INVERSE_INPUTS',script_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=sha(a.source/'execution_contract.json'),independent_DFT_relative_L2=errors,direct_inverse_checks=checks,limits='Arithmetic/source/input identity only, not complete FDTD error budget or field validation.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','public','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
