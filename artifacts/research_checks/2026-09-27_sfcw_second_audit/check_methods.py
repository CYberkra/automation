"""Bounded CPU-only independent-method/mesh screens; no solver and no test families."""
from pathlib import Path
import hashlib,json,math,inspect
import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as api
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'START_HERE.md').exists())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for role,stem in [('BG','B2D-C3m-BG'),('TGT','B2D-C3m-D10m-W4m-T0.5m-E20-S0.02')]:
 path=ROOT/f'artifacts/research_checks/2026-09-26_{stem}-CO11-t05/{stem}-CO11-t05.h5'
 src=api.load_source(path);rx=api.load_receiver(path,receiver_path='name:measurement',component='Ex')
 for ns in [0,200]:
  frac=0 if ns==0 else (round(ns*1e-9/rx.dt)-.25)/len(rx.samples)
  f=np.array([20e6,95e6,170e6])
  d=api.direct_frequency_response(src,rx,f,tail_taper_fraction=frac)
  h=api.homodyne_frequency_response(src,rx,f,cycles=8,tail_taper_fraction=frac)
  rel=np.abs(d.response-h.response)/np.abs(d.response)
  assert np.all(rel<1e-8),rel
  rows.append({'input':str(path.relative_to(ROOT)),'sha256':sha(path),'tail_ns':ns,'frequencies_hz':f.tolist(),'relative_direct_homodyne_difference':rel.tolist()})
mesh=[]
for folder,name in [('batch2d_v1','B2D-C3m-D10m-W4m-T0.5m-E20-S0.02'),('batch2d_v1','B2D-C3m-D10m-W4m-T0.5m-E28-S0.05'),('a0_3d_v1','B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM')]:
 path=ROOT/f'configs/research/{folder}/{name}.in';text=path.read_text();lines=text.splitlines()
 d=max(map(float,next(x.split(':')[1] for x in lines if x.startswith('#dx_dy_dz:')).split()))
 for line in lines:
  if not line.startswith('#material:'):continue
  er,sigma,mur,_,material=line.split(':')[1].split();er=float(er);sigma=float(sigma);mu=float(mur)*1.25663706212e-6;eps=er*8.8541878128e-12;w=2*math.pi*170e6
  beta=w*math.sqrt(mu*eps/2)*math.sqrt(math.sqrt(1+(sigma/(w*eps))**2)+1)
  wavelength=2*math.pi/beta;b=2*math.pi/(10*d)
  cutoff=(2*b*b/math.sqrt(4*b*b*mu*eps+mu*mu*sigma*sigma))/(2*math.pi)
  mesh.append({'input':str(path.relative_to(ROOT)),'sha256':sha(path),'material':material,'dx_max_m':d,'epsilon_r':er,'conductivity_S_m':sigma,'wavelength_170MHz_m':wavelength,'cells_per_wavelength_170MHz':wavelength/d,'ten_cell_screen_cutoff_hz':cutoff,'passes_ten_cell_screen_170MHz':wavelength/d>=10})
tails=[]
for stem in ['B2D-C3m-BG','B2D-C3m-D10m-W4m-T0.5m-E20-S0.02']:
 path=ROOT/f'artifacts/research_checks/2026-09-26_{stem}-MT33/{stem}-MT33.h5'
 with h5py.File(path) as h:
  x=np.column_stack([h[f'rxs/rx{i}/Ex'][:] for i in range(1,34)])
  tail=20*np.log10(np.max(abs(x[-int(np.ceil(.05*len(x))):]),axis=0)/np.max(abs(x),axis=0))
  tails.append({'input':path.relative_to(ROOT).as_posix(),'receiver_names':[f'mt{i:02}' for i in range(1,34)],'per_trace_raw_tail_dB':tail.tolist(),'per_trace_range_dB':[float(tail.min()),float(tail.max())]})
defaults={'source_floor_db':inspect.signature(api.direct_frequency_response).parameters['source_floor_db'].default,'normalise_window':inspect.signature(api.reconstruct_time_response).parameters['normalise_window'].default,'scope':'omitted call arguments resolved against processing_sha256; original export manifest preserved'}
r={'MT_per_trace_tail_check':tails,'pinned_API_defaults':defaults,'scope':'algorithm consistency and homogeneous material wavelength screening only; no physical convergence certificate','solver_invoked':False,'processing_sha256':sha(Path(api.__file__)),'homodyne_comparisons':rows,'mesh_screens':mesh}
Path(__file__).with_name('results.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'max_relative_method_difference':max(max(x['relative_direct_homodyne_difference']) for x in rows),'A0_mesh_screen':[x for x in mesh if '3D5CM' in x['input']],'solver_invoked':False},indent=2))
