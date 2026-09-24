"""Unfitted axial Yee plane-wave phase budget; not full dipole/PML theory."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.constants import c
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
f=np.linspace(20e6,170e6,501);dx=.05;dt=9.629166007732353e-11;courant=c*dt/dx;k=2*np.pi*f/c
kn=2/dx*np.arcsin(np.sin(np.pi*f*dt)/courant);phase=(kn-k)*30
result={'formula':'sin(omega*dt/2) = (c*dt/dx)*sin(k_num*dx/2), axial vacuum Yee wave','dx_m':dx,'dt_s':dt,'path_m':30,'assumptions':'Axial infinite uniform lattice, near-vertical ray approximation; no dipole amplitude, image interface or PML model; not fitted to output','frequency_Hz':f.tolist(),'excess_phase_rad':phase.tolist(),'predicted_response_phase_error_deg':(-phase*180/np.pi).tolist(),'source':'https://eecs.wsu.edu/~schneidj/ufdtd/chap7.pdf'}
with a.output.open('x') as out:json.dump(result,out,indent=2)
print('170MHz phase lag degrees',phase[-1]*180/np.pi)
