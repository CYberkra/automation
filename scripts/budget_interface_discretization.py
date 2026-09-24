"""Post-result, unfitted normal-incidence Yee interface and refinement budget."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.constants import c
p=argparse.ArgumentParser();p.add_argument('--comparison',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();v=np.genfromtxt(a.comparison,delimiter=',',names=True);f=v['frequency_Hz'];budgets=[]
for dz in (.05,.04,.025):
 dt=1/(c*np.sqrt(2/.05**2+1/dz**2));u=np.sin(np.pi*f*dt)/(c*dt/dz);q1=2*np.arcsin(u);q2=2*np.arcsin(3*u)
 gamma=(np.sin(q1)-np.sin(q2))/(np.sin(q1)+np.sin(q2));amp=20*np.log10(abs(gamma)/.5);phase=-(q1/dz-2*np.pi*f/c)*30*180/np.pi
 row=dict(dz_m=dz,dt_s=dt,cells=480*480*round(24/dz),relative_cells=.05/dz,phase_170MHz_deg=float(phase[-1]),interface_amplitude_170MHz_dB=float(amp[-1]))
 if dz==.05:row['max_observed_amplitude_minus_budget_abs_dB']=float(max(abs(v['amplitude_error_dB']-amp)))
 budgets.append(row)
r={'status':'post_result_mechanism_check_no_fitting_no_correction','model':'Normal incidence, er1=1 er2=9, mu=1, tangential-E interface arithmetic er=5. Gamma=(sin(q1)-sin(q2))/(sin(q1)+sin(q2)), qj=2asin(sqrt(er_j)*sin(omega dt/2)/(c dt/dz)). Equivalent to Schneider7.96. Not an exact 3D dipole reference.','source':'https://eecs.wsu.edu/~schneidj/ufdtd/chap7.pdf','budgets':budgets,'next_candidate':'dx=dy=.05 dz=.04,24m cube;20 lateral PML cells and25 vertical preserve1m thickness. New paired air and dielectric runs needed because dt changes. Estimated cells+25%; actual GPU/host/PML memory must be checked before launch. No new run authorized by this static file alone.','solver_invoked':False}
with a.output.open('x') as out:json.dump(r,out,indent=2)
print(json.dumps(r))
