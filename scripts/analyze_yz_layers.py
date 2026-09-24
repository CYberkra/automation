"""Compare official native-time SFCW with independent layered line-current spectra."""
import argparse,json,hashlib
from dataclasses import replace
from pathlib import Path
import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response
from analyze_vertical_refinement import metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
reference=Path('artifacts/research_checks/2026-09-25_yz_layers_reference/spectra.npz')
with np.load(reference) as ref:arrays={k:ref[k] for k in ref.files}
f=arrays['frequency_Hz'];paths={k:Path(f'artifacts/research_checks/2026-09-25_YZ_{k}/YZ_{k}.h5') for k in ('AIR','ROCK','COVER')}
src={k:load_source(v) for k,v in paths.items()};rx={k:load_receiver(v,receiver_path='name:measurement',component='Ex') for k,v in paths.items()};dt=src['AIR'].dt
checks={}
for k,path in paths.items():
 with h5py.File(path) as h:
  assert np.array_equal(h.attrs['nx_ny_nz'],[1,960,1760]) and h[rx[k].path].dtype==np.float64
  assert np.allclose(h[rx[k].path].parent.attrs['Position'][1:],[12.65,39],rtol=0,atol=1e-12)
 assert src[k].dt==rx[k].dt==dt and src[k].time_offset==dt/2 and rx[k].time_offset==0
 assert np.array_equal(src[k].samples,src['AIR'].samples) and np.all(np.isfinite(rx[k].samples))
 checks[k+'_metadata_finite_source']=True
n=min(len(rx['AIR'].samples),int(np.floor(800e-9/dt))+1);rows={}
for taper in (200,400):
 h={}
 for k in paths:
  result=direct_frequency_response(src[k],replace(rx[k],samples=rx[k].samples[:n]),f,tail_taper_fraction=(round(taper*1e-9/dt)-.25)/n)
  assert np.all(result.source_valid) and np.all(np.isfinite(result.response));h[k]=result.response
 for k,reference_key in [('ROCK','bare'),('COVER','covered')]:
  reflection=h[k]-h['AIR'];arrays[f'{k}_{taper}']=reflection;rows[f'{k}_{taper}']=metrics(reflection,arrays[reference_key])
result={'checks':checks,'comparisons':rows,'taper_200_vs_400':{k:metrics(arrays[k+'_200'],arrays[k+'_400']) for k in ('ROCK','COVER')},'inputs':{str(v):hashlib.sha256(v.read_bytes()).hexdigest() for v in [*paths.values(),reference]},'physical_acceptance_threshold':None,'solver_invoked':False}
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**arrays)
np.savetxt(a.output/'spectra.csv',np.column_stack([f]+[v for k in ('bare','covered','ROCK_200','COVER_200') for v in (arrays[k].real,arrays[k].imag)]),delimiter=',',header='frequency_Hz,bare_re,bare_im,covered_re,covered_im,ROCK_re,ROCK_im,COVER_re,COVER_im',comments='')
print(json.dumps(result))
