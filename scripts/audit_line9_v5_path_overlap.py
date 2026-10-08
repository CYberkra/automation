"""Independent transmission-product, Snell-root and direct-inverse path-time audit."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq,newton
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    r=json.loads(a.analysis.read_text('utf-8'));review=json.loads(a.review.read_text('utf-8'))
    assert sha(a.review)==r['review_sha256'] and sha(a.material)==r['material_sha256']
    db=json.loads(a.material.read_text('utf-8'))['materials'];f=20e6+np.arange(501)*300000.;omega=2*np.pi*f
    n=[]
    for mat in db.values():
        b=mat['base'];er=b['relative_permittivity']-1j*b['electric_conductivity_s_per_m']/(omega*8.8541878128e-12)
        for pole in mat.get('poles',[]):er+=pole['relative_permittivity_difference']/(1+1j*omega*pole['relaxation_time_s'])
        n.append(np.sqrt(er))
    n=np.column_stack(n);c0=299792458.;dt=review['dt_s'];kc=omega[:,None]*n/c0
    arg=n*.025/(c0*dt)*np.sin(np.pi*f*dt)[:,None];kg=2/.025*np.arcsin(arg);root_errors=[]
    for j in [0,250,500]:
        for medium in [0,1,2]:
            root=newton(lambda q:np.sin(q*.025/2)-arg[j,medium],kc[j,medium],fprime=lambda q:.025/2*np.cos(q*.025/2),tol=1e-12)
            error=abs(root-kg[j,medium])/abs(root);assert error<1e-12;root_errors.append(float(error))
    t=np.arange(4008)/(4008*300000.);matrix=np.exp(2j*np.pi*t[:,None]*f)/501
    trans=lambda a,b:2*n[:,a]/(n[:,a]+n[:,b])
    refl=lambda a,b:(n[:,a]-n[:,b])/(n[:,a]+n[:,b])
    basal_coef=trans(0,1)*trans(1,0)*trans(1,2)*trans(2,1)*refl(2,3)
    multiple_coef=trans(0,1)*trans(1,0)*refl(1,2)**2*refl(1,0)
    spectra={}
    for name,k in [('continuum',kc),('bulk_Yee',kg)]:
        primary=[];multiple=[]
        for row in review['records']:
            g=row['geometry'];h=g['midpoint_agl_m'];d1=g['cover_base']['depth_m'];d2=g['basal_sand']['depth_m']-d1
            # Native layouts audited in primary script are 1.3m horizontally separated.
            for lengths,coef,result in [(np.array([h,d1,d2,0]),basal_coef,primary),(np.array([h,2*d1,0,0]),multiple_coef,multiple)]:
                keep=lengths>0;v=lengths[keep];nr=n[250,keep].real
                p=brentq(lambda p:2*np.sum(v*p/np.sqrt(nr*nr-p*p))-1.3,0,float(nr.min())*(1-1e-12),xtol=1e-14)
                tau=2*np.sum(v*nr*nr/np.sqrt(nr*nr-p*p))/c0
                correction=tau-2*np.dot(lengths,n[250].real)/c0
                result.append(coef*np.exp(-2j*(k@lengths)-1j*omega*correction))
        spectra[name]=np.column_stack([*primary,*multiple])
    errors={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();errors[window]={}
        for name,z in spectra.items():
            x=matrix@(z*w[:,None]);peak=t[np.argmax(abs(x),axis=0)]*1e9
            expected=r['curves'][window][name]
            err=max(float(abs(peak[:194]-expected['basal_primary_ns']).max()),float(abs(peak[194:]-expected['cover_second_trip_ns']).max()))
            assert err<1e-9;errors[window][name]=err
    result=dict(status='PASS_INDEPENDENT_TRANSMISSION_SNELL_ROOT_BULK_ROOT_DIRECT_INVERSE',script_sha256=sha(__file__),
        analysis_sha256=sha(a.analysis),peak_time_max_absolute_error_ns=errors,bulk_sine_root_relative_error_max=max(root_errors),
        station_count=194,limits='Audits approximate path mathematics, not real event identity or validity of the local-plane approximation in nonflat geometry. No new FDTD.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['analysis','review','material','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
