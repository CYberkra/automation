"""Independent layered integrals and direct SFCW norms for material balance."""
import argparse
import json
import math
from pathlib import Path
import h5py
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import jv
from hs_capsule_identity import sha256 as sha

C=299792458.; EPS=8.8541878128e-12; MU=1/(EPS*C*C)


def eps(material, f):
    b=material['base']; real=b['relative_permittivity']; loss=b['electric_conductivity_s_per_m']/(2*np.pi*f*EPS)
    for p in material.get('poles',[]):
        x=2*np.pi*f*p['relaxation_time_s']; d=p['relative_permittivity_difference']
        real+=d/(1+x*x);loss+=d*x/(1+x*x)
    return real-1j*loss


def norm(z):
    return math.sqrt(math.fsum(float(v.real)**2+float(v.imag)**2 for v in np.ravel(z)))


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.public/'contract.json'); result=read(a.public/'analysis.json')
    assert sha(a.public/'contract.json')==result['contract_sha256'] and sha(a.arrays)==result['arrays_sha256']
    assert sha(Path(__file__).with_name('analyze_line9_material_return_balance.py'))==result['script_sha256']
    for name,digest in c['helper_sha256'].items():assert sha(Path(__file__).with_name(name))==digest
    with np.load(a.arrays) as h:f=h['frequency_Hz'];r=h['response'];epsilon=h['epsilon']
    assert np.array_equal(f,np.array(c['frequency_Hz']))
    constitutive=[]; adaptive=[]
    for j,rec in enumerate(c['materials']):
        materials=list(rec['materials'].values())
        e=np.array([[eps(m,frequency) for m in materials] for frequency in f])
        error=float(np.max(abs(e-epsilon[:,j])/abs(e)));assert error<1e-13;constitutive.append(error)
        for fi in [0,250,500]:
            freq=f[fi];k=2*np.pi*freq/C; e1=e[fi,1];e2=e[fi,2]
            def root(value):
                z=np.sqrt(value+0j);return -z if z.imag>0 else z
            def wave(u,evan=False):
                q=k*np.cosh(u) if evan else k*np.sin(u); ky=-1j*k*np.sinh(u) if evan else k*np.cos(u)
                k1=root(k*k*e1-q*q); k2=root(k*k*e2-q*q); phase=np.exp(-2j*k1*c['cover_thickness_m']);parts=[]
                for y1,y2 in [(k1,k2),(k1/e1,k2/e2)]:
                    top=(ky-y1)/(ky+y1);bottom=(y1-y2)/(y1+y2);first=(1-top*top)*bottom*phase;loop=-top*bottom*phase
                    parts.append(np.array([top,first,first*loop,first*loop*loop/(1-loop)]))
                air=np.exp(-1j*ky*c['air_height_sum_m']);beta=q*c['offset_m']
                line=-freq*MU*air*np.cos(beta)*(1j if evan else 1)*parts[0]
                point=-freq*MU/4*air*(1j*k*np.cosh(u) if evan else k*np.sin(u))*((jv(0,beta)-jv(2,beta))*parts[0]-(ky/k)**2*(jv(0,beta)+jv(2,beta))*parts[1])
                return np.stack([line,point])
            p,_=quad_vec(wave,0,np.pi/2,epsabs=1e-16,epsrel=1e-12)
            v,_=quad_vec(lambda u:wave(u,True),0,np.arcsinh(40/(k*c['air_height_sum_m'])),epsabs=1e-16,epsrel=1e-12)
            error=float(np.max(abs(r[fi,j]-(p+v))/abs(p+v)));assert error<1e-5
            adaptive.append({'case':rec['id'],'frequency_Hz':freq,'maximum_component_relative_error':error})
    old=read(Path(__file__).resolve().parents[1]/'artifacts/research_checks/2026-10-08_v401_polarization_dc_r1/analysis.json')
    assert sha(a.factorial_arrays)==old['numerical_sha256']
    assert c['source_analysis_sha256']==sha(Path(__file__).resolve().parents[1]/'artifacts/research_checks/2026-10-08_v401_polarization_dc_r1/analysis.json')
    with h5py.File(a.factorial_arrays) as h:
        assert np.array_equal(h['frequency_Hz'][:],f);z=h['response'][:]
    metric_errors=[];balance_errors=[];times=np.arange(4008)/(4008*300000.)
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();n=np.array([[norm(r[:,j,d,k]*w) for k in range(4)] for j in range(4) for d in range(2)]).reshape(4,2,4)
        ref=np.array(result['planar_component_relative_metrics'][window]['relative_to_original_by_case_dimension'])
        error=float(np.max(abs(n/n[0]-ref)/ref));assert error<1e-12;metric_errors.append(error)
        for gate in ['basal','deep']:
            lo,hi=old['gate_ns'][gate];t=times[(times*1e9>=lo)&(times*1e9<=hi)]
            values=np.exp(2j*np.pi*t[:,None]*f)@(z*w[:,None])/501
            hn=[norm(values[:,3*j]) for j in range(4)];dn=[norm(values[:,3*j+2]) for j in range(4)]
            for j,row in enumerate(result['actual_factorial_norm_balance'][window][gate]):
                expected=[hn[j]/hn[0],dn[j]/dn[0],dn[j]*hn[0]/(dn[0]*hn[j])]
                actual=[row[k] for k in ['H0_remaining_amplitude_norm_fraction','bottom_delta_amplitude_norm_gain','H0_over_delta_ratio_improvement']]
                error=float(np.max(abs(np.array(expected)-actual)/expected));assert error<1e-9;balance_errors.append(error)
                np.testing.assert_allclose(row['ratio_improvement_dB'],row['deep_gain_dB']+row['H0_reduction_dB'],rtol=1e-14,atol=1e-13)
    out={'status':'PASS_MATERIAL_BALANCE_12_ADAPTIVE_AND_16_NATIVE_NORM_ROWS','script_sha256':sha(__file__),
         'analysis_sha256':sha(a.public/'analysis.json'),'constitutive501_max_relative_error':max(constitutive),
         'adaptive_cases':adaptive,'adaptive_max_component_error':max(q['maximum_component_relative_error'] for q in adaptive),
         'planar_relative_metric_max_error':max(metric_errors),'native_norm_balance_rows':len(balance_errors),
         'native_direct_sum_max_relative_error':max(balance_errors),
         'limits':'Checks constitutive arithmetic and planar/reference balance, not site truth or nonflat path identity'}
    assert not a.out.exists();a.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k!='adaptive_cases'}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['public','arrays','factorial-arrays','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
