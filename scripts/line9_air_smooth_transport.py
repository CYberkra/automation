"""Smooth propagating-sector diagnostic; preserves the original hard-cut helper."""
import numpy as np
from line9_air_angular_transport import C0,project


def smooth_weight(frequencies,samples,spacing,sector,rolloff):
    if not 0<rolloff<sector<1:
        raise ValueError('Smooth rolloff must lie inside an interior propagating sector')
    ratio=abs(2*np.pi*np.fft.fftfreq(samples,spacing))[None,:]/(2*np.pi*np.asarray(frequencies)[:,None]/C0)
    transition=np.clip((ratio-sector+rolloff)/rolloff,0,1)
    return .5*(1+np.cos(np.pi*transition))


def smooth_transport(e,hx,frequencies,spacing,height_difference,sector=.9,rolloff=.25):
    result=project(e,hx,frequencies,spacing,height_difference,sector)
    weight=smooth_weight(frequencies,e.shape[1],spacing,sector,rolloff)
    return {name:np.fft.ifft(np.fft.fft(v,axis=1)*weight,axis=1) for name,v in result.items()}


def smooth_field(e,frequencies,spacing,sector=.9,rolloff=.25):
    weight=smooth_weight(frequencies,e.shape[1],spacing,sector,rolloff)
    return np.fft.ifft(np.fft.fft(e,axis=1)*weight,axis=1)
