"""Matched versus deliberately mismatched background subtraction; constructed controls."""
import argparse,hashlib,json
from pathlib import Path
from dataclasses import replace
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response
from research_evaluation_contract import waveform_metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--zfine',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
cp=Path('artifacts/research_checks/2026-09-25_yz_depth_analysis/arrays.npz');fp=Path('artifacts/research_checks/2026-09-25_deep_controls_analysis/arrays.npz');bp=Path('artifacts/research_checks/2026-09-25_DEP_BG_FINE/DEP_BG_FINE.h5')
with np.load(cp) as d:f=d['frequency_Hz'];coarse=(d['background'],d['DEP_20_difference'])
with np.load(fp) as d:fine_target=d['FINE_200']
src=load_source(bp);rx=load_receiver(bp,receiver_path='name:measurement',component='Ex');n=min(len(rx.samples),int(np.floor(1200e-9/rx.dt))+1)
result=direct_frequency_response(src,replace(rx,samples=rx.samples[:n]),f,tail_taper_fraction=(round(200e-9/rx.dt)-.25)/n);assert np.all(result.source_valid)
models=[('coarse',coarse),('fine',(result.response,fine_target))];inputs=[cp,fp,bp]
if a.zfine:
 zp=Path('artifacts/research_checks/2026-09-25_deep_zfine_analysis/arrays.npz');zb=Path('artifacts/research_checks/2026-09-25_DEP_BG_ZFINE/DEP_BG_ZFINE.h5')
 with np.load(zp) as x:target=x['ZFINE_200']
 zs=load_source(zb);zr=load_receiver(zb,receiver_path='name:measurement',component='Ex');zn=min(len(zr.samples),int(np.floor(1200e-9/zr.dt))+1)
 zh=direct_frequency_response(zs,replace(zr,samples=zr.samples[:zn]),f,tail_taper_fraction=(round(200e-9/zr.dt)-.25)/zn);assert np.all(zh.source_valid)
 models.append(('zfine',(zh.response,target)));inputs.extend([zp,zb])
t=np.arange(8192)/(8192*300e3);window=np.hanning(501);mask=((t>=440e-9)&(t<=600e-9))[:,None]
def wave(s):return (2*np.real(np.fft.ifft(s*window,n=len(t))*len(t)/window.sum()*np.exp(2j*np.pi*f[0]*t)))[:,None]
conditions=[('amplitude',e,0) for e in (0,1e-5,1e-4,1e-3,1e-2)]+[('timing',0,tau) for tau in (.02,.1,.5)]+[('late_coherent',e,500) for e in (1e-5,1e-4,1e-3)]
rows=[];arrays={'time_s':t,'mask':mask};checks=[]
for grid,(b,s) in models:
 truth=wave(s);arrays[grid+'_truth']=truth
 for i,(kind,e,tau_ns) in enumerate(conditions):
  phase=np.exp(-2j*np.pi*f*tau_ns*1e-9)
  bhat=(1+e)*b if kind=='amplitude' else b*phase if kind=='timing' else b+e*phase*b
  negative=b-bhat;positive=s+negative;yn=wave(negative);yp=wave(positive)
  diagnostic=waveform_metrics(yp,truth,mask,reference_kind='paired_contrast',state='isolated',scope='contrast')
  # Only a constructed additive-reference diagnostic, never physical eligibility.
  difference_error=float(np.linalg.norm((yp-yn-truth)[mask])/np.linalg.norm(truth[mask]));assert difference_error<1e-9
  if kind=='amplitude' and e==0:assert np.count_nonzero(negative)==0 and diagnostic['metrics']['nrmse']<1e-12
  row={'grid':grid,'condition':kind,'relative_amplitude_error':e,'delay_ns':tau_ns,'constructed_target_window':diagnostic['metrics'],
   'negative_window_peak_over_true_target_peak':float(max(abs(yn[mask]))/max(abs(truth[mask]))),'fullband_residual_over_target_L2':float(np.linalg.norm(negative)/np.linalg.norm(s)),
   'paired_difference_relative_error':difference_error,'physical_label_eligible':False}
  rows.append(row);arrays[f'{grid}_{i}_positive']=yp;arrays[f'{grid}_{i}_negative']=yn
r={'rows':rows,'checks':{'matched_background_preserves_target_and_zero_negative':True,'all_paired_differences_preserved':True},'source_reference_state':'numerically_unresolved','scope':'Constructed perturbation diagnosis, no measured noise or false-alarm/detection-rate claim','training_eligible':False,'unique_label':None,'solver_invoked':False,'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
(a.output/'results.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**arrays)
print(json.dumps({'rows':len(rows),'fine_late_controls':[x for x in rows if x['grid']=='fine' and x['condition']=='late_coherent']}))
