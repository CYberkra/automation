"""CPU planar layered 3D horizontal electric point-dipole reference.

exp(+iwt), outgoing Im(ky)<=0. phi is the horizontal receiver bearing relative
to the horizontal dipole axis. The project's Ez perpendicular to the scan is
phi=pi/2 here; this is not a vertical dipole or a finite antenna port.

The angular projection gives TE[J0+cos(2phi)J2] and
-TM(ky0/k0)^2[J0-cos(2phi)J2]. Radial Weyl normalization is
-omega*mu0/(8pi); both propagating and evanescent terms are retained.
"""
import numpy as np
from scipy.special import jv, roots_legendre
from diagnose_line9_v5_planar_green import root_outgoing, MU0
from review_line9_result_packages import C0


def layer_parts(q,ky0,k0,eps1,eps2,thickness):
    """Air/cover/mud TE and TM: surface, one cover return, two, higher remainder."""
    if not np.isfinite(thickness) or thickness<0:
        raise ValueError('Nonnegative finite cover thickness required')
    ky1=root_outgoing(k0[:,None]**2*eps1[:,None]-q*q)
    ky2=root_outgoing(k0[:,None]**2*eps2[:,None]-q*q)
    te=[];tm=[]
    for polar,target in [('TE',te),('TM',tm)]:
        y0=ky0;y1=ky1;y2=ky2
        if polar=='TM':y1=ky1/eps1[:,None];y2=ky2/eps2[:,None]
        a=(y0-y1)/(y0+y1);b=(y1-y2)/(y1+y2)
        p=np.exp(-2j*ky1*thickness)
        full=(a+b*p)/(1+a*b*p)
        once=(1-a*a)*b*p;twice=-(1-a*a)*a*b*b*p*p
        target.extend([a,once,twice,full-a-once-twice])
    return np.stack(te,axis=-1),np.stack(tm,axis=-1)


def point_integral(freq,height_sum,offset,phi,kernel,order=512,tail_exponent=36.):
    f=np.asarray(freq,dtype=float)
    if f.ndim!=1 or not np.isfinite(f).all() or (f<=0).any() or height_sum<=0 or offset<0 or order<16 or tail_exponent<=0 or not np.isfinite(phi):
        raise ValueError('Positive frequencies/height/order/tail and finite bearing required')
    x,w=roots_legendre(order);k0=2*np.pi*f/C0
    theta=(x+1)*np.pi/4;wt=w*np.pi/4
    umax=np.arcsinh(tail_exponent/(k0*height_sum));u=(x[None,:]+1)*umax[:,None]/2;wu=w[None,:]*umax[:,None]/2
    parts=[]
    for q,ky,measure in [(k0[:,None]*np.sin(theta),k0[:,None]*np.cos(theta),k0[:,None]*np.sin(theta)*wt),
                          (k0[:,None]*np.cosh(u),-1j*k0[:,None]*np.sinh(u),1j*k0[:,None]*np.cosh(u)*wu)]:
        te,tm=kernel(q,ky);b=q*offset;j0=jv(0,b);j2=jv(2,b);c=np.cos(2*phi)
        angular_te=j0+c*j2;angular_tm=-(ky/k0[:,None])**2*(j0-c*j2)
        wave=np.exp(-1j*ky*height_sum)*measure
        parts.append(np.sum(wave[:,:,None]*(angular_te[:,:,None]*te+angular_tm[:,:,None]*tm),axis=1))
    return -2*np.pi*f[:,None]*MU0/(8*np.pi)*(parts[0]+parts[1])
