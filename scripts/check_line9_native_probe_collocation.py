"""Analytic controls for spatial and native half-time-step TM collocation."""
import json
import numpy as np
from analyze_line9_native_probes import collocate


def main():
    checks=[];n=25;time=np.arange(n,dtype=float);field=1+.04*time
    for degrees in [0,45,90,180,270]:
        theta=np.deg2rad(degrees);h1=field*np.sin(theta)/377;h2=-field*np.cos(theta)/377
        # Linear histories sampled half a step earlier than E.
        half=1+.04*(time-.5)
        e,hx,hy=collocate(field,half*np.sin(theta)/377,half*np.sin(theta)/377,
                          -half*np.cos(theta)/377,-half*np.cos(theta)/377)
        np.testing.assert_allclose(e*hx,field[:-1]**2*np.sin(theta)/377,atol=1e-17)
        np.testing.assert_allclose(-e*hy,field[:-1]**2*np.cos(theta)/377,atol=1e-17)
        checks.append('progressive_linear_wave_'+str(degrees))
    # A field affine in space and time must reconstruct its value at Ez exactly.
    f=lambda x,y,t: 2+3*x+5*y+7*t
    e,hx,hy=collocate(field,f(4,8.5,time-.5),f(4,7.5,time-.5),
                      f(4.5,8,time-.5),f(3.5,8,time-.5))
    np.testing.assert_array_equal(hx,f(4,8,time[:-1]));np.testing.assert_array_equal(hy,f(4,8,time[:-1]))
    assert len(e)==n-1
    checks.append('affine_spacetime_exact_and_last_sample_discarded')
    phase=np.arange(4096)*2*np.pi/4096;up=np.cos(phase-.37);down=np.cos(phase+.37)
    flux=(up+down)*(up-down)/377
    assert abs(np.mean(flux))<1e-17 and np.max(abs(flux))>1e-4
    checks.append('standing_wave_instantaneous_sign_does_not_count_echoes')
    for name,args in [('unequal_length',[field,field[:-1],field,field,field]),
                      ('nonfinite',[field,field*np.nan,field,field,field]),
                      ('one_sample',[field[:1]]*5),('not_1d',[field[:,None]]*5)]:
        try:collocate(*args)
        except ValueError:checks.append('reject_'+name)
        else:raise AssertionError('Invalid histories accepted')
    print(json.dumps({'status':'PASS_NATIVE_COLLOCATION_ANALYTIC_CONTROLS','count':len(checks),
                      'checks':checks,'limits':'Arithmetic/sign controls only, not FDTD or physical path certification.'}))


if __name__=='__main__':main()
