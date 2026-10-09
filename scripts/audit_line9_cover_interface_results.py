"""Independent native/direct-DFT/inverse and analytical-delay audit, no solver."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from audit_line9_v401_delivery import native_response,sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    pub=read(a.public/'analysis.json');c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json');m=c['study_manifest']
    assert v['completed'] and c['max_runs']==2 and c['no_retry'] and v['contract_sha256']==pub['contract_sha256']==sha(a.source/'execution_contract.json')
    assert pub['verification_sha256']==sha(a.source/'completed_verification.json') and sha(a.numerical)==pub['numerical_sha256']
    assert pub['native_sha256'][0]==m['baseline_sha256']['native_H0.h5'] and pub['native_sha256'][1:]==[r['native_sha256'] for r in v['groups']]
    paths=[a.package/'baseline/native_H0.h5']+[a.source/g['id']/'profile.h5' for g in c['groups']]
    spectra=[native_response(p,.025,h) for p,h in zip(paths,pub['native_sha256'])]
    f=20e6+np.arange(501)*300000.
    with h5py.File(a.numerical) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],f);z=h['response'][:]
    assert z.shape==(501,3)
    errors=[float(np.linalg.norm(s-z[:,j])/np.linalg.norm(s)) for j,s in enumerate(spectra)];assert max(errors)<1e-9
    q=np.column_stack(spectra);t=np.arange(4008)/(4008*300000.);metrics={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();rows={g['id']:{} for g in c['groups']}
        for gate,(lo,hi) in m['gates_ns'].items():
            k=(t*1e9>=lo)&(t*1e9<=hi);tt=t[k];x=np.exp(2j*np.pi*tt[:,None]*f)@(q*w[:,None])/501
            jo=abs(x[:,0]).argmax();original_index=int(np.flatnonzero(k)[jo])
            for j,g in enumerate(c['groups'],1):
                r=pub['metrics'][window][g['id']][gate];old,new=x[:,0],x[:,j]
                ratio=float(np.linalg.norm(new)/np.linalg.norm(old));change=float(np.linalg.norm(new-old)/np.linalg.norm(old))
                corr=np.vdot(old,new)/(np.linalg.norm(old)*np.linalg.norm(new));at=float(abs(new[jo])/abs(old[jo]))
                assert abs(ratio-r['new_over_original_complex_L2'])<1e-8 and abs(change-r['relative_complex_change_L2'])<1e-8 and abs(at-r['new_over_old_at_original_peak'])<1e-8
                assert abs(abs(corr)-r['complex_correlation_abs'])<1e-8 and abs(np.angle(corr,deg=True)-r['phase_deg'])<1e-6
                peaks=[float(tt[abs(old).argmax()]*1e9),float(tt[abs(new).argmax()]*1e9)]
                np.testing.assert_allclose(peaks,r['peak_ns'],atol=1e-9,rtol=0);assert original_index==r['original_peak_time_coordinate_index']
                rows[g['id']][gate]=dict(new_over_original_complex_L2=ratio,relative_complex_change_L2=change,peak_ns=peaks,new_over_old_at_original_peak=at)
        metrics[window]=rows
    raw=[]
    for path in paths:
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['dt']==m['dt_s'] and h.attrs['Iterations']==20352
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            for node,key in [('srcs/src1','tx_m'),('rxs/rx1','rx_m')]:np.testing.assert_allclose(h[node].attrs['Position'],c['groups'][0][key],rtol=0,atol=1e-12)
            s=h['srcs/src1/excitation/samples'][:]
            if raw:assert s.tobytes()==oldsource.tobytes()
            else:oldsource=s
            raw.append(h['rxs/rx1/Ez'][:]);assert raw[-1].dtype==s.dtype==np.float64
    times=np.arange(20352)*m['dt_s']
    for j,g in enumerate(c['groups'],1):
        for gate,(lo,hi) in dict(early=[0,120],return_packet=[250,450],late=[450,1100]).items():
            k=(times*1e9>=lo)&(times*1e9<=hi);old,new=raw[0][k],raw[j][k];r=pub['native_metrics'][g['id']][gate]
            assert abs(np.linalg.norm(new)/np.linalg.norm(old)-r['new_over_original_L2'])<1e-12
            assert abs(np.linalg.norm(new-old)/np.linalg.norm(old)-r['relative_change_L2'])<1e-12
            assert (new.tobytes()==old.tobytes())==r['bitwise_equal']
    cover=read(a.package/'baseline/materials.json')['materials']['material_001_cover'];b=cover['base'];p=cover['poles'][0]
    omega=2*np.pi*f;eps0=8.8541878128e-12
    eps=b['relative_permittivity']+p['relative_permittivity_difference']/(1+1j*omega*p['relaxation_time_s'])-1j*b['electric_conductivity_s_per_m']/(eps0*omega)
    derivative=-1j*p['relaxation_time_s']*p['relative_permittivity_difference']/(1+1j*omega*p['relaxation_time_s'])**2+1j*b['electric_conductivity_s_per_m']/(eps0*omega**2)
    n=np.sqrt(eps);d=(n+omega*derivative/(2*n)).real/299792458.*1e9
    pred=pub['nominal_delay_hypotheses']
    np.testing.assert_allclose([d.min(),d.max()],pred['one_bottom_roundtrip_ns_range'],rtol=1e-8,atol=0)
    np.testing.assert_allclose([2*d.min(),2*d.max()],pred['two_bottom_roundtrips_ns_range'],rtol=1e-8,atol=0)
    assert abs(d[250]-pred['one_bottom_roundtrip_at95MHz_ns'])<1e-7
    result=dict(status='PASS_INDEPENDENT_THREE_NATIVE_DFT_DIRECT_INVERSE_AND_DECLARED_DELAY',audit_script_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),contract_sha256=pub['contract_sha256'],native_sha256=pub['native_sha256'],independent_DFT_relative_L2=errors,direct_inverse_metrics=metrics,source_bitwise_equal=True,analytical_nominal_one_roundtrip_at95MHz_ns=float(d[250]),limits='Arithmetic and identity checks; global replacement modifies bulk propagation, and local vertical group delay cannot independently certify nonlocal bounce count or field validity.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],independent_DFT_relative_L2=errors)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','public','package','numerical','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
