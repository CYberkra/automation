"""20m target differences under domain expansion and mesh refinement."""
import argparse,json,hashlib
from pathlib import Path
from dataclasses import replace
import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response
from analyze_vertical_refinement import metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--zfine',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
prior=Path('artifacts/research_checks/2026-09-25_yz_depth_analysis/arrays.npz')
with np.load(prior) as archive:arrays={'frequency_Hz':archive['frequency_Hz'],'baseline':archive['DEP_20_difference']}
f=arrays['frequency_Hz'];checks={};inputs={str(prior):hashlib.sha256(prior.read_bytes()).hexdigest()};tail={}
models=[('WIDE',[1,1920,2480],[24.65,51]),('FINE',[1,2560,4000],[16.65,45])]
if a.zfine:models.append(('ZFINE',[1,2560,8000],[16.65,45]))
for mode,shape,position in models:
 responses={};sources={}
 for kind in ('BG','20'):
  name=f'DEP_{kind}_{mode}';path=Path(f'artifacts/research_checks/2026-09-25_{name}/{name}.h5')
  src=load_source(path);rx=load_receiver(path,receiver_path='name:measurement',component='Ex');sources[kind]=src
  with h5py.File(path) as h:
   assert np.array_equal(h.attrs['nx_ny_nz'],shape) and h[rx.path].dtype==np.float64
   assert np.allclose(h[rx.path].parent.attrs['Position'][1:],position,rtol=0,atol=1e-12)
  assert src.dt==rx.dt and src.time_offset==src.dt/2 and rx.time_offset==0 and np.all(np.isfinite(rx.samples))
  n=min(len(rx.samples),int(np.floor(1200e-9/rx.dt))+1)
  for taper in (200,400):
   result=direct_frequency_response(src,replace(rx,samples=rx.samples[:n]),f,tail_taper_fraction=(round(taper*1e-9/rx.dt)-.25)/n)
   assert np.all(result.source_valid) and np.all(np.isfinite(result.response));responses[kind,taper]=result.response
  checks[name+'_metadata_source_finite']=True;inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 assert sources['BG'].dt==sources['20'].dt and np.array_equal(sources['BG'].samples,sources['20'].samples)
 for taper in (200,400):arrays[mode+'_'+str(taper)]=responses['20',taper]-responses['BG',taper]
 tail[mode]=float(np.linalg.norm(arrays[mode+'_200']-arrays[mode+'_400'])/np.linalg.norm(arrays[mode+'_200']))
rows={}
for mode,_,_ in models:
 response=arrays[mode+'_200'];baseline=arrays['baseline'];rows[mode]=metrics(response,baseline)
 for label,mask in [('all',f>0),('20_70MHz',f<70e6),('70_120MHz',(f>=70e6)&(f<120e6)),('120_170MHz',f>=120e6)]:
  rows[mode][label+'_relative_L2_difference']=float(np.linalg.norm((response-baseline)[mask])/np.linalg.norm(baseline[mask]))
 rows[mode]['spectral_RMS_ratio_dB']=float(20*np.log10(np.linalg.norm(response)/np.linalg.norm(baseline)))
t=np.arange(8192)/(8192*300e3);window=np.hanning(501);arrays['time_s']=t
for k in ['baseline']+[m+'_200' for m,_,_ in models]:arrays[k+'_envelope']=abs(np.fft.ifft(arrays[k]*window,n=8192))*8192/window.sum()
result={'checks':checks,'comparisons_to_original_2p5cm':rows,'tail_window_relative_L2':tail,'inputs':inputs,'solver_invoked':False,'physical_acceptance_threshold':None,'interpretation':'Differences relative to original grid/domain, not absolute errors versus exact target solution.'}
if a.zfine:
 result['zfine_vs_uniform_fine']=metrics(arrays['ZFINE_200'],arrays['FINE_200'])
 result['zfine_vs_uniform_fine']['relative_L2_difference']=float(np.linalg.norm(arrays['ZFINE_200']-arrays['FINE_200'])/np.linalg.norm(arrays['FINE_200']))
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**arrays)
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
curves=[('baseline','Original 2.5cm'),('WIDE_200','Larger domain 2.5cm'),('FINE_200','Original domain 1.25cm')]
if a.zfine:curves.append(('ZFINE_200','dy 1.25cm / dz 0.625cm'))
for k,label in curves:
 ax[0].plot(f/1e6,abs(arrays[k]),label=label);ax[1].plot(t*1e9,arrays[k+'_envelope'],label=label)
ax[0].set(xlabel='Frequency [MHz]',ylabel='Target difference [(V/m)/A]');ax[1].set(xlabel='Time [ns]',ylabel='Hann envelope [(V/m)/A]',xlim=(480,570))
for axis in ax:axis.legend(fontsize=8);axis.grid(alpha=.2)
fig.savefig(a.output/'comparison.png',dpi=160);plt.close(fig)
print(json.dumps(result))
