"""Adaptive quadrature and independent input-admittance audit of planar2D diagnostic."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import quad_vec
from diagnose_line9_v5_planar_green import reflect_parts,root_outgoing,MU0
from review_line9_result_packages import C0,indices
from hs_capsule_identity import sha256 as sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    r=read(a.public/'analysis.json');m=read(a.prepared/'manifest.json');db=a.prepared/'base_H1/geometries/line9_research_materials_v1_smoothed.json'
    assert sha(db)==r['material_sha256']==m['material_sha256'] and sha(a.numerical)==r['numerical_sha256']
    assert sha(Path(__file__).parent/'diagnose_line9_v5_planar_green.py')==r['script_sha256']
    materials=read(db)['materials'];hs,hr=r['source_heights_m'];height=hs+hr;dx=r['horizontal_offset_m'];d1=r['cover_thickness_m'];d2=r['mud_thickness_m']
    with h5py.File(a.numerical) as h:f=h['frequency_Hz'][:];stored=h['response'][:]
    selected=[0,83,250,417,500];errors=[];admittance_errors=[]
    for j in selected:
        fj=f[j];k=2*np.pi*fj/C0*indices(materials,fj);k0=k[0].real
        def kernel(q,ky0):
            parts=reflect_parts(np.array([[q]]),np.array([[ky0]]),np.array([k[1]]),np.array([k[2]]),np.array([k[3]]),d1,d2)[0,0]
            ky=root_outgoing(k*k-q*q)
            # Independent input admittances; same-medium lower half-space forH0.
            for bottom,last in [(ky[2],1),(ky[3],2)]:
                y=bottom
                for layer in range(last,0,-1):
                    tangent=np.tan(ky[layer]*([0,d1,d2][layer]));y=ky[layer]*(y+1j*ky[layer]*tangent)/(ky[layer]+1j*y*tangent)
                reflected=(ky0-y)/(ky0+y)
                expected=parts[:4].sum() if last==1 else parts.sum()
                admittance_errors.append(float(abs(reflected-expected)))
            return parts
        def propagating(theta):
            q=k0*np.sin(theta);ky=k0*np.cos(theta)
            return np.cos(q*dx)*np.exp(-1j*ky*height)*kernel(q,ky)
        def evanescent(u):
            q=k0*np.cosh(u);ky=-1j*k0*np.sinh(u)
            return 1j*np.cos(q*dx)*np.exp(-1j*ky*height)*kernel(q,ky)
        p,_=quad_vec(propagating,0,np.pi/2,epsabs=1e-12,epsrel=1e-11)
        e,_=quad_vec(evanescent,0,np.arcsinh(40/(k0*height)),epsabs=1e-12,epsrel=1e-11)
        independent=-(2*np.pi*fj)*MU0/(2*np.pi*.025)*(p+e)
        err=abs(independent-stored[j,1:])/abs(stored[j,1:]);assert max(err)<1e-5;errors.append(err.tolist())
    assert max(admittance_errors)<1e-10
    result=dict(status='PASS_INDEPENDENT_ADAPTIVE_QUADRATURE_AND_ADMITTANCE_NOT_FDTD_AGREEMENT',audit_script_sha256=sha(__file__),selected_frequency_Hz=f[selected].tolist(),
        adaptive_vs_gauss_relative_error_by_component=errors,reflection_admittance_absolute_error_max=max(admittance_errors),finite_tail_exponent_main=36,finite_tail_exponent_independent=40,
        note='Only analytic planar model audited. Especially tiny higher-order component has its own relative quadrature error; no field-case or material validation implied.')
    (a.public/'independent_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],max_adaptive_relative=max(max(x) for x in errors),admittance=max(admittance_errors))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prepared','public','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
