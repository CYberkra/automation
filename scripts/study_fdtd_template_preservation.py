"""Constructed 32-column mechanism study, not simulated flight-line data or labels."""
import argparse,json,hashlib
from pathlib import Path
from dataclasses import replace
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response
from research_operator_contract import catalogue,apply_configuration,ConfigUnavailable
from research_evaluation_contract import waveform_metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--zfine',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
coarse_path=Path('artifacts/research_checks/2026-09-25_yz_depth_analysis/arrays.npz');fine_path=Path('artifacts/research_checks/2026-09-25_deep_controls_analysis/arrays.npz')
bg_path=Path('artifacts/research_checks/2026-09-25_DEP_BG_FINE/DEP_BG_FINE.h5')
with np.load(coarse_path) as x:f=x['frequency_Hz'];coarse=(x['background'],x['DEP_20_difference'])
with np.load(fine_path) as x:fine_target=x['FINE_200']
src=load_source(bg_path);rx=load_receiver(bg_path,receiver_path='name:measurement',component='Ex');n=min(len(rx.samples),int(np.floor(1200e-9/rx.dt))+1)
r=direct_frequency_response(src,replace(rx,samples=rx.samples[:n]),f,tail_taper_fraction=(round(200e-9/rx.dt)-.25)/n);assert np.all(r.source_valid)
models=[('coarse',coarse),('fine',(r.response,fine_target))];inputs=[coarse_path,fine_path,bg_path]
if a.zfine:
 zp=Path('artifacts/research_checks/2026-09-25_deep_zfine_analysis/arrays.npz');zb=Path('artifacts/research_checks/2026-09-25_DEP_BG_ZFINE/DEP_BG_ZFINE.h5')
 with np.load(zp) as x:target=x['ZFINE_200']
 zs=load_source(zb);zr=load_receiver(zb,receiver_path='name:measurement',component='Ex');zn=min(len(zr.samples),int(np.floor(1200e-9/zr.dt))+1)
 zh=direct_frequency_response(zs,replace(zr,samples=zr.samples[:zn]),f,tail_taper_fraction=(round(200e-9/zr.dt)-.25)/zn);assert np.all(zh.source_valid)
 models.append(('zfine',(zh.response,target)));inputs.extend([zp,zb])
t=np.arange(8192)/(8192*300e3);keep=t<=1200e-9;window=np.hanning(501)
def real_wave(s):return (2*np.real(np.fft.ifft(s*window,n=8192)*8192/window.sum()*np.exp(2j*np.pi*f[0]*t)))[keep]
profiles={'common':np.ones(32),'localized':np.exp(-.5*((np.arange(32)-15.5)/4)**2)}
mask=np.broadcast_to(((t[keep]>=440e-9)&(t[keep]<=600e-9))[:,None],(sum(keep),32)).copy()
rows=[];saved={'time_s':t[keep],'mask':mask,**profiles};configs=[c for c in catalogue() if c['end_gain']==1]
for grid,(background,target) in models:
 b=real_wave(background);s=real_wave(target);saved[grid+'_background']=b;saved[grid+'_target']=s
 for shape,weights in profiles.items():
  x0=np.broadcast_to(b[:,None],mask.shape).copy();truth=s[:,None]*weights; x1=x0+truth
  for config in configs:
   row={'grid':grid,'profile':shape,'config_id':config['id'],'physical_label_eligible':False}
   try:
    y0=apply_configuration(x0,config['id'])['output'];y1=apply_configuration(x1,config['id'])['output'];change=y1-y0
    diagnostic=waveform_metrics(change,truth,mask,reference_kind='paired_contrast',state='isolated',scope='contrast')
    assert diagnostic['available'] and not diagnostic['absolute_preservation_eligible']
    row.update(available=True,constructed_contrast_diagnostic=diagnostic['metrics'])
    saved[grid+'_'+shape+'_'+config['id']]=change
   except ConfigUnavailable as e:row.update(available=False,reason=str(e))
   rows.append(row)
# Meaningful algebraic controls: identity preserves the constructed difference,
# and global full-mean subtraction removes a column-common target.
for row in rows:
 if row['config_id']=='B0_G1_BG':assert row['constructed_contrast_diagnostic']['nrmse']<1e-9
 if row['profile']=='common' and row['config_id']=='B3_G1_BG':assert row['constructed_contrast_diagnostic']['nrmse']>0.999999
result={'rows':rows,'checks':{'identity_contrast_preserved':True,'common_target_deleted_by_full_mean':True},'design':'32 constructed columns, common or prescribed Gaussian target weights, no flight-coordinate mapping. q=1 only; no physical clutter-removal score.',
 'reference_state_physical':'numerically_unresolved','constructed_reference_state':'isolated by algebraic addition for mechanism diagnostic only',
 'unique_label':None,'training_eligible':False,'solver_invoked':False,'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**saved)
print(json.dumps({'rows':len(rows),'available':sum(r['available'] for r in rows),'checks':result['checks']}))
