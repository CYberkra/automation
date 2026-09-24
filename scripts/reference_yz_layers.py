"""Independent 2D electric line-current reflection in a nonmagnetic planar stack.

Engineering exp(+i*w*t) convention; no FDTD, field data, or fitted correction.
Air source/receiver heights are 15m each, transverse separation 1.3m.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.constants import c,mu_0,epsilon_0
from scipy.special import roots_legendre,hankel2


def vertical(w,ky,er,sigma):
    q=np.sqrt((w/c)**2*er-1j*w*mu_0*sigma-ky**2+0j)
    return np.where(q.imag>0,-q,q)


def coefficient(w,ky,q0,cover,rock,thickness):
    q1=vertical(w,ky,*cover);q2=vertical(w,ky,*rock)
    r01=(q0-q1)/(q0+q1);r12=(q1-q2)/(q1+q2)
    travel=np.exp(-2j*q1*thickness)
    return (r01+r12*travel)/(1+r01*r12*travel)


def reflected(frequency,cover=(16,.01),rock=(9,.001),thickness=3.,order=256,cutoff=40.,constant=None):
    w=2*np.pi*np.asarray(frequency)[:,None];k=w/c
    x,weights=roots_legendre(order);theta=(x+1)*np.pi/4
    ky=k*np.sin(theta);q0=k*np.cos(theta)
    r=constant if constant is not None else coefficient(w,ky,q0,cover,rock,thickness)
    propagating=np.sum(weights*np.pi/4*r*np.cos(ky*1.3)*np.exp(-1j*q0*30),axis=1)
    umax=np.arcsinh(cutoff/(k*30));u=(x+1)*umax/2
    ky=k*np.cosh(u);q0=-1j*k*np.sinh(u)
    r=constant if constant is not None else coefficient(w,ky,q0,cover,rock,thickness)
    evanescent=np.sum(weights*umax/2*1j*r*np.cos(ky*1.3)*np.exp(-1j*q0*30),axis=1)
    return -w[:,0]*mu_0/(2*np.pi)*(propagating+evanescent)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    f=np.linspace(20e6,170e6,501);scale=-2*np.pi*f*mu_0/4*hankel2(0,2*np.pi*f/c*np.hypot(30,1.3))
    pec=reflected(f,constant=-1);zero=reflected(f,cover=(1,0),rock=(1,0))
    bare=reflected(f,cover=(9,.001),rock=(9,.001));covered=reflected(f)
    checks={'PEC_image_relative':float(max(abs(pec+scale)/abs(scale))),
      'matched_air_relative':float(max(abs(zero)/abs(scale))),
      'zero_thickness_relative':float(max(abs(reflected(f,thickness=0)-bare)/abs(bare))),
      'same_medium_thickness_relative':float(max(abs(reflected(f,cover=(9,.001),thickness=7)-bare)/abs(bare)))}
    for n,kw in [('bare',{'cover':(9,.001)}),('covered',{})]:
        ref=bare if n=='bare' else covered
        for order in (128,512):checks[f'{n}_order_{order}_relative']=float(max(abs(reflected(f,order=order,**kw)-ref)/abs(ref)))
        checks[n+'_cutoff_50_relative']=float(max(abs(reflected(f,cutoff=50,**kw)-ref)/abs(ref)))
    assert all(v<1e-8 for v in checks.values()),checks
    result={'checks':checks,'threshold':'1e-8 numerical identity/convergence checks, not FDTD acceptance','materials_hypothetical':True,'cover':{'er':16,'sigma_S_m':.01,'thickness_m':3},'rock':{'er':9,'sigma_S_m':.001},'source_height_m':15,'receiver_height_m':15,'baseline_m':1.3,'solver_invoked':False,'units':'(V/m)/A','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(a.output/'spectra.npz',frequency_Hz=f,bare=bare,covered=covered,PEC=pec)
    np.savetxt(a.output/'spectra.csv',np.column_stack([f,bare.real,bare.imag,covered.real,covered.imag]),delimiter=',',header='frequency_Hz,bare_real,bare_imag,covered_real,covered_imag',comments='')
    print(json.dumps(result))


if __name__=='__main__':main()
