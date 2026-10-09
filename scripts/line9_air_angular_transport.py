"""Finite-aperture propagating-sector angular transport in homogeneous air.

This is a diagnostic projection, not full-space ray/energy identification.
exp(+i omega t), Ez upgoing has Hx=+ky/(omega*mu) Ez.
Only |qx|<=sector*k0 is retained; grazing and evanescent fields are excluded.
"""
import numpy as np

C0=299792458.
EPS0=8.8541878128e-12
MU0=1/(EPS0*C0*C0)


def project(e, hx, frequencies, spacing, height_difference, sector=.9):
    e=np.asarray(e);hx=np.asarray(hx);f=np.asarray(frequencies,float)
    if e.ndim!=2 or e.shape!=hx.shape or e.shape[0]!=len(f) or f.ndim!=1:
        raise ValueError('Equal frequency by spatial-sample arrays required')
    if not np.isfinite(e).all() or not np.isfinite(hx).all() or not np.isfinite(f).all() or np.any(f<=0):
        raise ValueError('Finite fields and positive frequencies required')
    if not np.isfinite(spacing) or spacing<=0 or not np.isfinite(height_difference) or height_difference<0 or not 0<sector<1:
        raise ValueError('Positive spacing, nonnegative height and interior propagating sector required')
    q=2*np.pi*np.fft.fftfreq(e.shape[1],spacing)[None,:]
    omega=2*np.pi*f[:,None];k=omega/C0
    valid=abs(q)<=sector*k
    ky=np.sqrt(np.maximum(k*k-q*q,0))
    # No arbitrary epsilon at the light line; excluded modes use unity solely
    # as a computational placeholder before an exact zero mask is applied.
    impedance=np.where(valid,omega*MU0/np.where(valid,ky,1),0)
    ef=np.fft.fft(e,axis=1);hf=np.fft.fft(hx,axis=1)
    up=np.where(valid,(ef+impedance*hf)/2,0)
    down=np.where(valid,(ef-impedance*hf)/2,0)
    propagated_up=up*np.exp(-1j*ky*height_difference)
    propagated_down=down*np.exp(1j*ky*height_difference)
    return {name:np.fft.ifft(v,axis=1) for name,v in
            [('up',up),('down',down),('predicted_up',propagated_up),('predicted_down',propagated_down)]}


def sector_field(e, frequencies, spacing, sector):
    q=2*np.pi*np.fft.fftfreq(e.shape[1],spacing)[None,:]
    mask=abs(q)<=sector*(2*np.pi*np.asarray(frequencies)[:,None]/C0)
    return np.fft.ifft(np.fft.fft(e,axis=1)*mask,axis=1)
