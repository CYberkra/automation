"""Official SFCW response versus air line-current and PEC image references."""
import argparse,hashlib,json
from pathlib import Path
from dataclasses import replace
import h5py
import numpy as np
from scipy.constants import c,mu_0
from scipy.special import hankel2
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response
from analyze_vertical_refinement import metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False,parents=True)
paths={n:Path(f'artifacts/research_checks/2026-09-25_YZ_{n}/YZ_{n}.h5') for n in ('AIR','PEC')}
sources={n:load_source(path) for n,path in paths.items()}
receivers={n:load_receiver(path,receiver_path='name:measurement',component='Ex') for n,path in paths.items()}
dt=receivers['AIR'].dt;f=np.linspace(20e6,170e6,501);omega=2*np.pi*f
def line(distance):return -omega*mu_0/4*hankel2(0,omega/c*distance)
theory={'air':line(1.3),'reflection':-line(np.hypot(1.3,30))}
checks={};meta={}
for n,path in paths.items():
 src=sources[n];rx=receivers[n]
 with h5py.File(path) as h:
  assert np.array_equal(h.attrs['nx_ny_nz'],[1,960,1760])
  assert h[rx.path].dtype==np.float64
  assert np.allclose(h[rx.path].parent.attrs['Position'][1:],[12.65,39],rtol=0,atol=1e-12)
 assert np.isclose(dt,1/(c*np.sqrt(2/.025**2)),rtol=1e-12,atol=0)
 assert rx.dt==src.dt==dt and rx.time_offset==0 and src.time_offset==dt/2
 assert np.array_equal(src.samples,sources['AIR'].samples)
 assert np.all(np.isfinite(rx.samples)) and np.any(rx.samples)
 meta[n]={'dt':dt,'samples':len(rx.samples),'source_units':src.units,'source_quantity':src.quantity,'source_spatial_scale':src.spatial_scale,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
 checks[n+'_grid_dtype_position_source_timing_finite']=True
nt=min(len(receivers['AIR'].samples),int(np.floor(800e-9/dt))+1)
arrays={'frequency_Hz':f,'theory_air':theory['air'],'theory_reflection':theory['reflection']};rows={}
for taper in (0,200,400):
 fraction=0 if taper==0 else (round(taper*1e-9/dt)-.25)/nt
 h={}
 for n in paths:
  result=direct_frequency_response(sources[n],replace(receivers[n],samples=receivers[n].samples[:nt]),f,tail_taper_fraction=fraction)
  assert np.all(result.source_valid) and np.all(np.isfinite(result.response));h[n]=result.response
 air=h['AIR'];reflection=h['PEC']-air
 arrays[f'air_{taper}']=air;arrays[f'reflection_{taper}']=reflection
 rows[str(taper)]={'air':metrics(air,theory['air']),'reflection':metrics(reflection,theory['reflection'])}
result={'checks':checks,'metadata':meta,'taper_ns':rows,'taper_200_vs_400':{n:metrics(arrays[n+'_200'],arrays[n+'_400']) for n in ('air','reflection')},'physical_acceptance_threshold':None,'solver_invoked':False,'reference':'Engineering convention: -omega*mu0/4 * Hankel2(0,k*r); PEC reflected image has opposite sign. No fitted factor or delay.'}
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**arrays)
np.savetxt(a.output/'spectra.csv',np.column_stack([f]+[v for key in ('air_200','reflection_200','theory_air','theory_reflection') for v in (arrays[key].real,arrays[key].imag)]),delimiter=',',header='frequency_Hz,air_re,air_im,reflection_re,reflection_im,theory_air_re,theory_air_im,theory_reflection_re,theory_reflection_im',comments='')
print(json.dumps(result))
