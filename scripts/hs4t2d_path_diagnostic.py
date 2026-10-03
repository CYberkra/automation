"""Model-informed planar-ground ray timings; not an electromagnetic inversion."""
import numpy as np
from scipy.optimize import brentq

C0=299792458.
EPS0=8.8541878128e-12


def cover_index(frequency):
    omega=2*np.pi*np.asarray(frequency)
    epsilon=18.017+7.878/(1+1j*omega*6.4567e-9)-1j*.003/(omega*EPS0)
    return np.sqrt(epsilon)


def ray_timings(profile_z, midpoints, height):
    """Minimum travel time to any interface point, versus vertical-under-midpoint.

    Ground is flat, cover index uses 95MHz. Rays cross the ground according to
    phase-index Snell stationary paths. Report group delay on those paths.
    Selecting the first ray is not selecting the strongest scattered echo.
    """
    frequency=95e6
    phase_index=float(cover_index(frequency).real)
    df=1e3
    group_index=phase_index+frequency*float((cover_index(frequency+df).real-cover_index(frequency-df).real)/(2*df))
    q=(np.arange(480)+.5)*.025
    depth=12-np.asarray(profile_z)[np.floor(q/.25).astype(int)]
    phase=[]
    group=[]
    reflection=[]
    vertical=[]
    for midpoint in midpoints:
        tx,rx=midpoint-.65,midpoint+.65
        phase_time=np.zeros(len(q))
        group_time=np.zeros(len(q))
        for station in (tx,rx):
            for j,(target,d) in enumerate(zip(q,depth)):
                if target==station:
                    entry=target
                else:
                    def stationary(entry):
                        return (entry-station)/np.hypot(height,entry-station)+phase_index*(entry-target)/np.hypot(d,entry-target)
                    entry=brentq(stationary,min(station,target),max(station,target))
                air=np.hypot(height,entry-station)
                earth=np.hypot(d,entry-target)
                phase_time[j]+=(air+phase_index*earth)/C0
                group_time[j]+=(air+group_index*earth)/C0
        k=int(np.argmin(phase_time))
        phase.append(phase_time[k]*1e9)
        group.append(group_time[k]*1e9)
        reflection.append(float(q[k]))
        d=12-profile_z[min(47,int(midpoint/.25))]
        vertical.append((2*height+2*group_index*d)/C0*1e9)
    return {'phase_index_95MHz':phase_index,'group_index_95MHz':group_index,
            'minimum_ray_phase_time_ns':phase,'minimum_ray_group_time_ns':group,
            'minimum_ray_reflection_x_m':reflection,'vertical_group_time_ns':vertical,
            'scope':'earliest stationary-path diagnostic using 95MHz index and flat ground; not envelope-peak truth'}
