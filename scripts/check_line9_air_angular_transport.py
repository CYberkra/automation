"""Periodic plane-wave controls for homogeneous air angular transport."""
import json
import numpy as np
from line9_air_angular_transport import project,sector_field,MU0,C0


def main():
    f=np.array([40e6,95e6,170e6]);dx=.5;n=512;dy=4.;x=np.arange(n)*dx
    checks={}
    for mode in [-20,-7,0,7,20]:
        q=2*np.pi*mode/(n*dx);k=2*np.pi*f/C0;ky=np.sqrt(k*k-q*q)
        wave=np.exp(1j*q*x)[None,:]*np.array([1+.2j,.4-.3j,-.1+.8j])[:,None]
        for kind,sign in [('up',1),('down',-1)]:
            hx=sign*ky[:,None]/(2*np.pi*f[:,None]*MU0)*wave
            result=project(wave,hx,f,dx,dy,.9)
            expected=wave*np.exp(-sign*1j*ky[:,None]*dy)
            error=np.linalg.norm(result['predicted_'+kind]-expected)/np.linalg.norm(expected)
            other='down' if kind=='up' else 'up'
            leakage=np.linalg.norm(result[other])/np.linalg.norm(wave)
            assert error<1e-12 and leakage<1e-12
            checks[kind+'_mode_'+str(mode)]={'propagation_relative_L2':float(error),'opposite_leakage':float(leakage)}
    # Counterpropagating superposition must retain both independent amplitudes.
    q=2*np.pi*7/(n*dx);ky=np.sqrt((2*np.pi*f/C0)**2-q*q)
    up=np.exp(1j*q*x)[None,:]*np.ones((3,1));down=.3*np.exp(-1j*q*x)[None,:]*np.ones((3,1))
    h=ky[:,None]/(2*np.pi*f[:,None]*MU0)*(up-down)
    result=project(up+down,h,f,dx,dy)
    err=max(np.linalg.norm(result['up']-up),np.linalg.norm(result['down']-down));assert err<1e-12
    checks['counterpropagating_superposition']=float(err)
    zero=project(np.zeros((3,n)),np.zeros((3,n)),f,dx,dy);assert all(np.max(abs(v))==0 for v in zero.values());checks['zero']=True
    narrow=sector_field(up,f,dx,.01);assert np.linalg.norm(narrow)<1e-12;checks['excluded_sector_zero']=True
    for name,args in [('negative_height',{'height_difference':-1}),('lightline',{'sector':1}),('zero_spacing',{'spacing':0})]:
        kw=dict(spacing=dx,height_difference=dy);kw.update(args)
        try:project(up,h,f,**kw)
        except ValueError:checks[name]=True
        else:raise AssertionError(name)
    print(json.dumps({'status':'PASS_AIR_ANGULAR_TRANSPORT_PERIODIC_CONTROLS','count':len(checks),'checks':checks,
                      'limits':'Periodic plane controls do not certify finite aperture, FDTD or unique ray identity'}))


if __name__=='__main__':main()
