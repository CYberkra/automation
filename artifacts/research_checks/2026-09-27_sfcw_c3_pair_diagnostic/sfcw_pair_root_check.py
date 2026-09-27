"""Independent complex pair-window recomputation using math.fsum norms."""
from pathlib import Path
import json,hashlib,math
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'START_HERE.md').exists())
C=ROOT/'configs/research/sfcw_c3_pair_diagnostic_v1.json';c=json.loads(C.read_text())
mp=ROOT/c['input_manifest'];assert hashlib.sha256(mp.read_bytes()).hexdigest()==c['input_manifest_sha256']
m=json.loads(mp.read_text());zp=mp.parent/m['array_archive']['path'];assert hashlib.sha256(zp.read_bytes()).hexdigest()==m['array_archive']['sha256']
norm=lambda a:math.sqrt(math.fsum(float(v.real)**2+float(v.imag)**2 for v in a.ravel()))
ratio=lambda a,b: None if b==0 else a/b
rows=[];sensitivity=[]
with np.load(zp,allow_pickle=False) as z:
 data={}
 for fam in ['CO','MT']:
  for role in ['BG','TGT']:
   for taper in ['no_taper','tail_200ns']:
    rr=[r for r in m['records'] if (r['family'],r['role'],r['taper'])==(fam,role,taper)]
    rr.sort(key=lambda r:r['receiver_positions_xyz_m'][0][1])
    cols=[]
    for r in rr:
     k=r['official_reconstruction']['output_keys'];t=z[k['envelope_time_s']];v=z[k['complex_envelope']];cols.append(v[:,None] if v.ndim==1 else v)
    data[(fam,role,taper)]=np.concatenate(cols,axis=1)
  for taper in ['no_taper','tail_200ns']:
   a=data[(fam,'BG',taper)];b=data[(fam,'TGT',taper)];d=b-a
   for w in c['windows']:
    ix=np.flatnonzero((t>=w['start_s'])&(t<w['stop_s']));na=norm(a[ix]);nd=norm(d[ix])
    for op in c['operators']:
     if op=='identity':u=a;v=b
     else:
      u=a-np.mean(a,axis=1,keepdims=True);v=b-np.mean(b,axis=1,keepdims=True)
     q=v-u
     rows.append({'geometry':fam,'taper':taper,'window':w['id'],'operator':op,'n_samples':len(ix),'first_time_s':float(t[ix[0]]),'last_time_s':float(t[ix[-1]]),'background_norm':na,'delta_norm':nd,'background_over_delta':ratio(na,nd),'background_residual_ratio':ratio(norm(u[ix]),na),'paired_change_error':ratio(norm(q[ix]-d[ix]),nd),'delta_retention':ratio(norm(q[ix]),nd)})
  a=data[(fam,'BG','no_taper')];b=data[(fam,'BG','tail_200ns')];d=data[(fam,'TGT','no_taper')]-a;e=data[(fam,'TGT','tail_200ns')]-b
  for w in c['windows']:
   ix=np.flatnonzero((t>=w['start_s'])&(t<w['stop_s']))
   sensitivity.append({'geometry':fam,'window':w['id'],'background_taper_change':ratio(norm(b[ix]-a[ix]),norm(a[ix])),'delta_taper_change':ratio(norm(e[ix]-d[ix]),norm(d[ix]))})
result={'contract_sha256':hashlib.sha256(C.read_bytes()).hexdigest(),'method':'direct NPZ loading; independent math.fsum complex norms','rows':rows,'taper_sensitivity':sensitivity}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'D10':[r for r in rows if r['window']=='D10_context' and r['operator']=='trace_mean_removal'],'D10_taper':[r for r in sensitivity if r['window']=='D10_context']},indent=2))
