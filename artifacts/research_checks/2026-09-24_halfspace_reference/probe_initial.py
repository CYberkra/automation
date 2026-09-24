from pathlib import Path
import time,json,numpy as np,empymod
from scipy.constants import c,epsilon_0,mu_0
f=np.array([20.,50.,95.,130.,170.])*1e6;k=2*np.pi*f/c;r=1.3;th=-1j*np.sqrt(mu_0/epsilon_0)*k*.05/(4*np.pi*r)*(1+1/(1j*k*r)+1/(1j*k*r)**2)*np.exp(-1j*k*r)
start=time.time();v=empymod.dipole(src=[0,0,-15],rec=[0,1.3,-15],depth=[],res=[np.inf],epermH=[1],freqtime=f,ab=11,xdirect=True,verb=1)*.05
print('fullspace relative',max(abs(v-th)/abs(th)),flush=True)
for meth,kw in [('dlf',{}),('qwe',{'nquad':101,'maxint':80,'rtol':1e-9}),('quad',{'a':1e-8,'b':20,'pts_per_dec':1000,'rtol':1e-8,'limit':1000})]:
 t=time.time();v=empymod.dipole(src=[0,0,-15],rec=[0,1.3,-15],depth=[0],res=[np.inf,np.inf],epermH=[1,9],freqtime=f,ab=11,xdirect=None,ht=meth,htarg=kw,verb=1)*.05
 print(meth,'s',time.time()-t,'real',v.real,'imag',v.imag,flush=True)
