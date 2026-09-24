"""Full-wave halfspace reflection using the unmodified empymod1.10.6 kernel.

Branch-aware Gaussian quadrature is our integration wrapper, not an empymod
built-in transform. Only isotropic lossless upper vacuum, er>=1 lower halfspace,
x-directed electric point current and cross-track receiver are supported.
"""
import numpy as np
from scipy.constants import epsilon_0,mu_0,c
from scipy.special import j0,j1,roots_legendre
from empymod import kernel
from functools import lru_cache

@lru_cache(maxsize=16)
def nodes(n):return roots_legendre(n)

def transverse(f,r,length=.05):
 k=2*np.pi*np.asarray(f)/c
 return -1j*np.sqrt(mu_0/epsilon_0)*k*length/(4*np.pi*r)*(1+1/(1j*k*r)+1/(1j*k*r)**2)*np.exp(-1j*k*r)

def reflected(f,er=9.,height=15.,offset=1.3,length=.05,order=256,umax=2.,sigma_lower=0.):
 if er<1 or height<=0 or offset<=0:raise ValueError('Unsupported geometry/material')
 w=2*np.pi*float(f);k=w*np.sqrt(epsilon_0*mu_0)
 eta=np.array([[1j*w*epsilon_0, sigma_lower+1j*w*epsilon_0*er]])
 zeta=np.full((1,2),1j*w*mu_0);depth=np.array([-np.inf,0.])
 x,weights=nodes(order)
 def integrand(lam):
  p0,p1,p0b=kernel.wavenumber(-height,-height,0,0,depth,eta,eta,zeta,zeta,lam[None,:],11,True,False,False,False)
  # Angle is pi/2: cos(2*angle)=-1. ab11 has J2 recurrence factor1/r.
  return ((p0-p0b).ravel()*j0(lam*offset)-p1.ravel()*j1(lam*offset)/offset)
 # lambda=k*sin(theta), 0<theta<pi/2 removes upper-medium branch singularity.
 theta=(x+1)*np.pi/4;lam=k*np.sin(theta);jac=k*np.cos(theta)
 prop=np.dot(weights,integrand(lam)*jac)*np.pi/4
 # lambda=sqrt(k^2+u^2), upper gamma=u; split at lower branch if present.
 split=[0.,umax];uc=k*np.sqrt(er-1)
 if sigma_lower==0 and 0<uc<umax:split.insert(1,uc)
 evan=0j
 for lo,hi in zip(split[:-1],split[1:]):
  u=lo+(x+1)*(hi-lo)/2;lam=np.sqrt(k*k+u*u)
  evan+=np.dot(weights,integrand(lam)*u/lam)*(hi-lo)/2
 return complex((prop+evan)*length)
