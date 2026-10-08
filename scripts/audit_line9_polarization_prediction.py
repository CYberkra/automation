"""Independent input-admittance/adaptive integration of factorial planar predictions."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import hankel2
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    r=json.loads(a.analysis.read_text('utf-8'))
    assert sha(a.numerical)==r['numerical_sha256']
    with h5py.File(a.numerical) as h:f=h['frequency_Hz'][:];stored=h['response'][:]
    selected=[0,250,500];rows=[];c0=299792458.;e0=8.8541878128e-12;mu0=1/(e0*c0*c0)
    height=r['source_height_sum_m'];dx=r['horizontal_separation_m'];d1=r['cover_thickness_m'];d2=r['mud_thickness_m']
    for case,row in enumerate(r['cases']):
        errors=[]
        for j in selected:
            fj=f[j];omega=2*np.pi*fj;media=[]
            for mat in row['material_content']['materials'].values():
                db=mat['base'];er=db['relative_permittivity']-1j*db['electric_conductivity_s_per_m']/(omega*e0)
                for pole in mat.get('poles',[]):er+=pole['relative_permittivity_difference']/(1+1j*omega*pole['relaxation_time_s'])
                media.append(np.sqrt(er))
            k=omega/c0*np.array(media);k0=k[0].real
            def coefficient(q,ky0):
                ky=np.sqrt(k*k-q*q);ky=np.where(ky.imag>0,-ky,ky)
                answers=[]
                for bottom,last in [(ky[2],1),(ky[3],2)]:
                    y=bottom
                    for layer in range(last,0,-1):
                        tangent=np.tan(ky[layer]*[0,d1,d2][layer])
                        y=ky[layer]*(y+1j*ky[layer]*tangent)/(ky[layer]+1j*y*tangent)
                    answers.append((ky0-y)/(ky0+y))
                return np.array([*answers,answers[1]-answers[0]])
            def p(theta):
                q=k0*np.sin(theta);ky=k0*np.cos(theta)
                return np.cos(q*dx)*np.exp(-1j*ky*height)*coefficient(q,ky)
            def e(u):
                q=k0*np.cosh(u);ky=-1j*k0*np.sinh(u)
                return 1j*np.cos(q*dx)*np.exp(-1j*ky*height)*coefficient(q,ky)
            pv,_=quad_vec(p,0,np.pi/2,epsabs=1e-14,epsrel=1e-11)
            ev,_=quad_vec(e,0,np.arcsinh(40/(k0*height)),epsabs=1e-14,epsrel=1e-11)
            z=-fj*mu0/.025*(pv+ev)
            # Actual source heights differ by .05m; direct term only in H0/H1.
            direct=-omega*mu0/(4*.025)*hankel2(0,k0*np.hypot(dx,.05))
            z[:2]+=direct
            expected=stored[j,3*case:3*case+3]
            err=abs(z-expected)/abs(expected)
            assert max(err)<1e-5,(row['id'],fj,err)
            errors.append(err.tolist())
        rows.append(dict(id=row['id'],adaptive_vs_gauss_relative_error=errors))
    result=dict(status='PASS_INDEPENDENT_INPUT_ADMITTANCE_ADAPTIVE_QUADRATURE',script_sha256=sha(__file__),
        analysis_sha256=sha(a.analysis),selected_frequencies_Hz=f[selected].tolist(),cases=rows,
        maximum_relative_error=max(v for row in rows for errors in row['adaptive_vs_gauss_relative_error'] for v in errors),
        limits='Three audited frequencies per material condition; continuous planar local-column model only. No nonflat FDTD path/field validation.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['analysis','numerical','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
