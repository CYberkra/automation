"""Independent limits and angular controls for planar 3D spectral reference."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import jv
from line9_planar_point_layer_green import point_integral,layer_parts
from diagnose_line9_direct_dimensionality import point_dipole
from green_halfspace_hed_v0_1 import reflected_ex
from review_line9_result_packages import C0


def main():
    f=np.array([20e6,95e6,170e6]);k=2*np.pi*f/C0;checks={}
    for phi in [0.,np.pi/2]:
        def pec(q,ky):
            shape=(*q.shape,1);return -np.ones(shape),np.ones(shape)
        value=point_integral(f,16,1.3,phi,pec,order=512)[:,0]
        exact=-point_dipole(f,np.hypot(16,1.3*np.sin(phi)),np.full(3,1.3*np.cos(phi)))
        err=float(np.linalg.norm(value-exact)/np.linalg.norm(exact));assert err<1e-8
        checks['PEC_image_phi_'+str(phi)]=err
    # Direct numerical angular projection, independent of J0/J2 formula.
    for phi in [0.,.43,np.pi/2]:
        for beta in [.1,2.,5.]:
            s=quad(lambda a:np.sin(a)**2*np.cos(beta*np.cos(a-phi)),0,2*np.pi,epsabs=1e-12)[0]
            c=quad(lambda a:np.cos(a)**2*np.cos(beta*np.cos(a-phi)),0,2*np.pi,epsabs=1e-12)[0]
            expected=np.pi*np.array([jv(0,beta)+np.cos(2*phi)*jv(2,beta),jv(0,beta)-np.cos(2*phi)*jv(2,beta)])
            err=float(np.max(abs(np.array([s,c])-expected)));assert err<1e-11
            checks[f'angular_projection_phi{phi}_beta{beta}']=err
    q=k[:,None]*np.array([.1,.8,1.2]);ky=np.sqrt(k[:,None]**2-q*q+0j);ky=np.where(ky.imag>0,-ky,ky)
    te,tm=layer_parts(q,ky,k,np.ones(3),np.ones(3),6.575)
    assert np.max(abs(te))==0 and np.max(abs(tm))==0;checks['all_air_zero_reflection']=True
    eps1=np.full(3,11-.3j);eps2=np.full(3,13-.9j)
    # At zero thickness the complete stack must equal air/mud directly.
    te,tm=layer_parts(q,ky,k,eps1,eps2,0)
    kz=np.sqrt(k[:,None]**2*eps2[:,None]-q*q);kz=np.where(kz.imag>0,-kz,kz)
    expected=[(ky-kz)/(ky+kz),(eps2[:,None]*ky-kz)/(eps2[:,None]*ky+kz)]
    err=max(float(np.max(abs(z.sum(axis=-1)-e))) for z,e in zip([te,tm],expected));assert err<1e-12
    checks['zero_thickness_direct_halfspace']=err
    # With equal cover/mud the full response is a halfspace; compare legacy
    # independent midpoint quadrature at phi0, not the new Gauss machinery.
    def halfspace(q,ky):return layer_parts(q,ky,k,eps1,eps1,6.575)
    val=point_integral(f,16,1.3,0,halfspace)[:,:,].sum(axis=1)
    old=reflected_ex(f,1.3,16,eps1)
    err=float(np.linalg.norm(val-old)/np.linalg.norm(val));assert err<2e-6
    checks['legacy_midpoint_halfspace_comparison']=err
    for thickness in [-1,np.nan]:
        try:layer_parts(q,ky,k,eps1,eps2,thickness)
        except ValueError:checks['reject_thickness_'+str(thickness)]=True
        else:raise AssertionError('Invalid thickness accepted')
    print(json.dumps({'status':'PASS_PLANAR_POINT_LAYER_REFERENCE_CONTROLS','count':len(checks),'checks':checks,
                      'limits':'Analytic/source/angular/stack limits; not nonflat FDTD, antenna or site certification.'}))


if __name__=='__main__':main()
