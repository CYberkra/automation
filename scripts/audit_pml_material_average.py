"""Small official PML averaging probe; no grid creation or FDTD solver."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from scipy.constants import c
from gprMax.cython.pml_build import pml_average_er_mr
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
plane=np.zeros((480,480),dtype=np.uint32);plane[:,:80]=1
result={}
for name,ers in [('air',np.array([1.,1.])),('PEC',np.array([1.,1.])),('dielectric',np.array([1.,9.]))]:
 er,mr=pml_average_er_mr(480,480,8,plane,ers,np.ones(2));result[name]={'average_er':er,'average_mr':mr,'sigma_max_relative_to_air':1/np.sqrt(er*mr)}
result['side_boundary_geometric_paths_ns']={'ymin_inner_outer':[21.9/c*1e9,23.9/c*1e9],'ymax_inner_outer':[22.1/c*1e9,24.1/c*1e9],'x_inner_outer':[np.hypot(22,1.3)/c*1e9,np.hypot(24,1.3)/c*1e9]}
result['interpretation']='Source-derived mechanism candidate: changing lower 4m from air to er9 changes average er on full-height side PML, hence automatically selected sigma_max even in air portion. PEC nominal er remains1. Geometric side returns bracket observed77ns precursor. Not a causal isolation experiment; no solver invoked.'
source=Path('E:/gprMax-v.4.0.0/gprMax-v.4.0.0/gprMax');installed=Path('artifacts/local_checks/gprmax_v4_gpu_env/Lib/site-packages/gprMax')
result['source_sha256']={}
for name in ('grid/fdtd_grid.py','pml.py'):
 b=(source/name).read_bytes();assert hashlib.sha256(b).hexdigest()==hashlib.sha256((installed/name).read_bytes()).hexdigest();result['source_sha256'][name]=hashlib.sha256(b).hexdigest()
result['source_sha256']['cython/pml_build.pyx']=hashlib.sha256((source/'cython/pml_build.pyx').read_bytes()).hexdigest()
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result))
