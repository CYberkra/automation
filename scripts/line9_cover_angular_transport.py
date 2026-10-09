"""Conditional TM split in passive Debye cover, restricted to air-radiating q.

exp(+i omega t); outgoing ky has positive real and nonpositive imaginary part.
This finite sampled aperture does not certify absence of aliased high-q fields.
"""
import numpy as np
from line9_air_angular_transport import C0,EPS0,MU0
from line9_air_smooth_transport import smooth_weight


def cover_permittivity(f):
    f=np.asarray(f,float)
    if f.ndim!=1 or not np.isfinite(f).all() or np.any(f<=0):
        raise ValueError('Finite positive frequency vector required')
    omega=2*np.pi*f
    return 11+.5/(1+1j*omega*6.4567e-9)-1j*.001/(omega*EPS0)


def transport(e,hx,f,spacing,dy,sector=.9,rolloff=.25,epsilon=None):
    e=np.asarray(e);hx=np.asarray(hx);f=np.asarray(f,float)
    if e.ndim!=2 or e.shape!=hx.shape or e.shape[0]!=len(f) or f.ndim!=1 or e.shape[1]<2:
        raise ValueError('Matching frequency-by-position E and H required')
    if not np.isfinite(e).all() or not np.isfinite(hx).all() or not np.isfinite(f).all() or np.any(f<=0):
        raise ValueError('Finite fields and positive frequencies required')
    if not np.isfinite(spacing) or spacing<=0 or not np.isfinite(dy):
        raise ValueError('Positive spacing and finite signed plane displacement required')
    eps=cover_permittivity(f) if epsilon is None else np.asarray(epsilon,complex)
    if eps.shape!=f.shape or not np.isfinite(eps).all() or np.any(eps.real<=0) or np.any(eps.imag>0):
        raise ValueError('Passive positive-real permittivity vector required')
    weight=smooth_weight(f,e.shape[1],spacing,sector,rolloff)
    q=2*np.pi*np.fft.fftfreq(e.shape[1],spacing)[None,:]
    omega=2*np.pi*f[:,None];k0=omega/C0
    valid=abs(q)<=sector*k0
    ky=np.sqrt(k0*k0*eps[:,None]-q*q+0j)
    # Passive positive epsilon selects this outgoing branch without clipping loss.
    assert np.all(ky[valid].real>0) and np.all(ky[valid].imag<=0)
    safe=np.where(valid,ky,1);imp=np.where(valid,omega*MU0/safe,0)
    ef=np.fft.fft(e,axis=1);hf=np.fft.fft(hx,axis=1)
    up=(ef+imp*hf)*weight/2;down=(ef-imp*hf)*weight/2
    phase=np.where(valid,ky,0)*dy
    result={'up':up,'down':down,'predicted_up':up*np.exp(-1j*phase),
            'predicted_down':down*np.exp(1j*phase),'total':ef*weight}
    return {name:np.fft.ifft(v,axis=1) for name,v in result.items()}
