"""Official small CFS profile equivalence probe; no FDTD grid or solve."""
import json,argparse
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from gprMax import config
from gprMax.pml import CFS
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
config.sim_config=SimpleNamespace(em_consts=config.SimulationConfig.em_consts,dtypes={'float_or_double':np.float64})
auto=CFS();auto.calculate_sigmamax(.05,1.,1.);text=repr(float(auto.sigma.max));fixed=CFS();fixed.sigma.max=float(text)
checks={}
for name in ('alpha','kappa','sigma'):
 av=auto.calculate_values(20,getattr(auto,name));fv=fixed.calculate_values(20,getattr(fixed,name))
 for field,x,y in zip(('E','H'),av,fv):checks[name+'_'+field]=bool(np.array_equal(x,y))
assert all(checks.values())
result={'sigma_max_air':float(auto.sigma.max),'literal':text,'z0':float(config.SimulationConfig.em_consts['z0']),'checks':checks,'solver_invoked':False,'scope':'21-sample scratch CFS profiles only, temporary constants/dtype namespace; does not instantiate SimulationConfig or FDTDGrid','reuse_air_reason':'Air er=mr=1 on all six faces and dx=dy=dz=.05; fixed numeric literal exactly equals automatic sigma and all six E/H profile arrays. Same dt, geometry, PML order, and thickness. Existing T800 is the matched air reference without another FDTD run.'}
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result))
