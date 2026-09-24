"""Compare completed CUDA HSG control; official transforms, native Yee time axes."""
import argparse,json,csv,hashlib
from pathlib import Path
from dataclasses import replace
import numpy as np
import h5py
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response,engineering_dft,apply_tail_taper
from audit_official_sfcw import transverse_dipole
from analyze_vertical_refinement import metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--hsg-path',type=Path);p.add_argument('--ratio',type=int,default=3);a=p.parse_args();a.output.mkdir(exist_ok=False)
paths={k:Path('artifacts/research_checks/2026-09-25_'+n)/(n+'.h5') for k,n in [('reference','HSG_reference'),('hsg','HSG_air')]}
if a.hsg_path:paths['hsg']=a.hsg_path
assert a.ratio in (1,3,5)
freq=np.linspace(20e6,170e6,501);sources={k:load_source(v) for k,v in paths.items()};responses={};checks={};meta={};raw={}
assert sources['reference'].dt==sources['hsg'].dt and np.array_equal(sources['reference'].samples,sources['hsg'].samples)
for k,path in paths.items():
 for name in ('measurement','inside'):
  rx=load_receiver(path,receiver_path='name:'+name,component='Ex');source=sources[k]
  with h5py.File(path) as f:
   assert f[rx.path].dtype==np.float64
   assert np.allclose(f[rx.path].parent.attrs['Position'],[3,3.65 if name=='measurement' else 3,3],rtol=0,atol=1e-12)
  assert np.all(np.isfinite(rx.samples)) and np.any(rx.samples!=0) and rx.time_offset==0
  assert source.time_offset==source.dt/2 and source.spatial_scale==.05
  expected=source.dt/a.ratio if k=='hsg' and name=='inside' else source.dt
  assert np.isclose(rx.dt,expected,rtol=1e-14,atol=0)
  n=min(len(rx.samples),int(np.floor(200e-9/rx.dt))+1);fraction=(round(100e-9/rx.dt)-.25)/n
  receiver=replace(rx,samples=rx.samples[:n]);srcspec=engineering_dft(source.samples,source.dt,freq,time_offset=source.time_offset)
  native=engineering_dft(apply_tail_taper(receiver.samples,fraction),rx.dt,freq,time_offset=rx.time_offset)/srcspec
  if rx.dt==source.dt:
   official=direct_frequency_response(source,receiver,freq,tail_taper_fraction=fraction)
   assert np.all(official.source_valid) and np.array_equal(native,official.response)
  assert np.all(np.isfinite(native))
  key=k+'_'+name;responses[key]=native;raw[key]=rx.samples
  meta[key]=dict(path=rx.path,dt_s=rx.dt,time_offset_s=rx.time_offset,samples=len(rx.samples),transformed_samples=n)
  checks[key+'_finite_double_timing_position']=True
rows={};arrays={'frequency_Hz':freq,**responses,**{'raw_'+k:v for k,v in raw.items()}}
for name,distance in [('measurement',1.3),('inside',.65)]:
 ref=responses['reference_'+name];hsg=responses['hsg_'+name];theory=transverse_dipole(freq,distance,.05)
 rows[name]=dict(hsg_vs_reference=metrics(hsg,ref),reference_vs_infinite_air=metrics(ref,theory),hsg_vs_infinite_air=metrics(hsg,theory))
 arrays['theory_'+name]=theory
checks['native_DFT_matches_official_direct_on_equal_dt']=True
checks['same_source_history']=True
result=dict(checks=checks,metadata=meta,comparisons=rows,raw_outer_difference_peak_over_reference_peak=float(max(abs(raw['hsg_measurement']-raw['reference_measurement']))/max(abs(raw['reference_measurement']))),
 inputs={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in paths.items()},
 inner_transform='Official engineering_dft and apply_tail_taper on each native dt, source division; project diagnostic adapter because direct_frequency_response rejects mixed dt. No resampling or phase fitting.',
 physical_acceptance_threshold=None,geological_run_approved_by_this_result=False,solver_invoked=False)
(a.output/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode());np.savez_compressed(a.output/'arrays.npz',**arrays)
with (a.output/'spectra.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['frequency_Hz']+[k+s for k in responses for s in ('_real','_imag')]);w.writerows(zip(freq,*[x for h in responses.values() for x in (h.real,h.imag)]))
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for name in ('measurement','inside'):
 ratio=responses['hsg_'+name]/responses['reference_'+name];ax[0].plot(freq/1e6,20*np.log10(abs(ratio)),label=name);ax[1].plot(freq/1e6,np.angle(ratio,deg=True),label=name)
for axis,label in zip(ax,('HSG/reference amplitude [dB]','HSG/reference phase [degrees]')):axis.set(xlabel='Frequency [MHz]',ylabel=label);axis.legend()
fig.savefig(a.output/'comparison.png',dpi=170);plt.close(fig);print(json.dumps(result))
