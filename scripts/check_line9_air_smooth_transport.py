"""Known periodic modes check smooth-sector field amplitude and phase."""
import json
import numpy as np
from line9_air_angular_transport import MU0,C0
from line9_air_smooth_transport import smooth_transport,smooth_field


def main():
    f=np.array([40e6,95e6,170e6]);n=512;dx=.5;dy=4.;x=np.arange(n)*dx;checks={}
    for mode in [0,7,20,27,40]:
        q=2*np.pi*mode/(n*dx);k=2*np.pi*f/C0;prop=abs(q)<k
        if not prop.all():continue
        ky=np.sqrt(k*k-q*q);wave=np.exp(1j*q*x)[None,:]*np.ones((3,1))
        # Scalar explicit cosine transition, independent of FFT weight helper.
        ratio=abs(q)/k;u=np.clip((ratio-.65)/.25,0,1);weight=(1+np.cos(np.pi*u))/2
        for sign,name in [(1,'up'),(-1,'down')]:
            hx=sign*ky[:,None]/(2*np.pi*f[:,None]*MU0)*wave
            result=smooth_transport(wave,hx,f,dx,dy)
            expected=weight[:,None]*wave*np.exp(-sign*1j*ky[:,None]*dy)
            error=float(np.linalg.norm(result['predicted_'+name]-expected)/np.linalg.norm(wave));assert error<1e-12
            checks[f'{name}_mode{mode}']=error
        expected=weight[:,None]*wave;err=float(np.linalg.norm(smooth_field(wave,f,dx)-expected)/np.linalg.norm(wave));assert err<1e-12
        checks['field_mode'+str(mode)]=err
    for roll in [0,1,np.nan]:
        try:smooth_transport(wave,hx,f,dx,dy,rolloff=roll)
        except ValueError:checks['reject_roll_'+str(roll)]=True
        else:raise AssertionError(roll)
    print(json.dumps({'status':'PASS_SMOOTH_AIR_TRANSPORT_PERIODIC_CONTROLS','count':len(checks),'checks':checks,
                      'limits':'Known modes only, no finite-aperture/physical/energy certification'}))


if __name__=='__main__':main()
