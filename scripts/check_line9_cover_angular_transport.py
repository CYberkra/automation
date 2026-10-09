"""Analytic lossy plane-wave and independent sign/air-limit controls."""
import argparse,json
from pathlib import Path
import numpy as np
from line9_cover_angular_transport import transport,cover_permittivity,C0,EPS0,MU0
from line9_air_smooth_transport import smooth_transport
from hs_capsule_identity import sha256 as sha


def main(a):
    rows=[];f=np.array([20e6,95e6,170e6]);n=128;dx=.5;x=np.arange(n)*dx
    # Independently written material spectrum, including DC loss once only.
    w=2*np.pi*f;eps=11+.5*(1-1j*w*6.4567e-9)/(1+(w*6.4567e-9)**2)-1j*.001/(w*EPS0)
    np.testing.assert_allclose(cover_permittivity(f),eps,rtol=2e-15,atol=0)
    for mode in [0,1,2]:
        q=2*np.pi*mode/(n*dx);ky=np.sqrt((w/C0)**2*eps-q*q)
        sp=np.exp(-1j*q*x)[None,:];up=(1+.2j)*sp*np.ones((3,1));down=(.3-.1j)*sp*np.ones((3,1))
        e=up+down;h=ky[:,None]/(w[:,None]*MU0)*(up-down)
        ratio=q/(w/C0);weight=(1+np.cos(np.pi*np.clip((ratio-.9+.25)/.25,0,1)))/2
        for dy in [-3,-1.5,0,1.5,3]:
            got=transport(e,h,f,dx,dy)
            expected={'up':up,'down':down,'predicted_up':up*np.exp(-1j*ky[:,None]*dy),
                      'predicted_down':down*np.exp(1j*ky[:,None]*dy),'total':e}
            errs={}
            for key,v in expected.items():
                wanted=v*weight[:,None];err=float(np.linalg.norm(got[key]-wanted)/np.linalg.norm(wanted));assert err<2e-13;errs[key]=err
            rows.append({'mode':mode,'dy_m':dy,'errors':errs})
    e=np.tile(np.exp(-1j*2*np.pi*x/(n*dx)),(3,1));h=e*.002
    got=transport(e,h,f,dx,4,epsilon=np.ones(3));air=smooth_transport(e,h,f,dx,4)
    for key,v in air.items():np.testing.assert_allclose(got[key],v,rtol=1e-13,atol=1e-14)
    bad=[('active_medium',{'epsilon':np.array([11+.1j]*3)}),('negative_real',{'epsilon':np.array([-1-.1j]*3)}),
         ('epsilon_shape',{'epsilon':np.array([11-.1j])}),('zero_spacing',{'spacing':0}),('nan_displacement',{'dy':np.nan}),
         ('invalid_sector',{'sector':1}),('invalid_rolloff',{'rolloff':.95}),('nan_fields',{'e':e*np.nan}),
         ('different_H_shape',{'hx':h[:,:-1]}),('negative_frequency',{'f':-f})]
    rejected=[]
    for name,change in bad:
        kw=dict(e=e,hx=h,f=f,spacing=dx,dy=1.5);kw.update(change)
        try:transport(**kw)
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    v={'status':'PASS_15_LOSSY_PLANE_CASES_AIR_LIMIT_AND_10_REJECTIONS','plane_cases':rows,'rejected':rejected,
       'material_independent_formula':True,'air_limit':True,'code_sha256':{p:sha(Path(__file__).with_name(p)) for p in ['line9_cover_angular_transport.py','check_line9_cover_angular_transport.py']},
       'new_solver_runs':0,'limits':'Analytic formula controls; no finite-aperture geological or field validation.'}
    assert not a.out.exists();a.out.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');print(v['status'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
