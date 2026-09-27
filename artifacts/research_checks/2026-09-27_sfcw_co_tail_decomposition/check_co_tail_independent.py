"""Independent H5 read, raised-cosine weights and direct complex DFT; no gprMax import."""
from pathlib import Path
import hashlib,json,math
import h5py
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'START_HERE.md').exists())
cp=ROOT/'configs/research/sfcw_co_tail_decomposition_v1.json'
c=json.loads(cp.read_text(encoding='utf-8'))
mp=ROOT/c['manifest'];assert hashlib.sha256(mp.read_bytes()).hexdigest()==c['manifest_sha256']
m=json.loads(mp.read_text(encoding='utf-8'))
rr=sorted([r for r in m['records'] if r['family']=='CO' and r['role']=='BG' and r['taper']=='no_taper'],key=lambda r:r['receiver_positions_xyz_m'][0][1])
columns=[];source=None
for r in rr:
 p=ROOT/r['input_path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['input_sha256']
 with h5py.File(p) as h:
  columns.append(h['rxs/rx1/Ex'][:]);s=h['srcs/src1/excitation/samples'][:]
  dt=float(h.attrs['dt']);offset=float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
  assert float(h['rxs/rx1/Ex'].attrs['TimeSampleOffset'])==0
  if source is None:source=s
  else:assert np.array_equal(source,s)
x=np.column_stack(columns);n=len(x);t=np.arange(n)*dt
center=lambda a:a-a.mean(axis=1,keepdims=True)
norm=lambda a:math.sqrt(math.fsum(float(v.real)**2+float(v.imag)**2 for v in np.ravel(a)))
rat=lambda a,b:None if b==0 else a/b
count=round(200e-9/dt);w=np.ones(n);w[-count:]=.5*(1+np.cos(np.linspace(0,np.pi,count)))
u=(1-w[:,None])*x;r=center(x);wr=center(w[:,None]*x)
raw=[];total=norm(r)**2;coverage=np.zeros(n,dtype=int)
for win in c['raw_windows']:
 ix=(t>=win['start_s'])&(t<win['stop_s']);coverage+=ix;nr=norm(r[ix])
 raw.append({'window':win['id'],'n':int(ix.sum()),'raw_norm':norm(x[ix]),'centered_norm':nr,'centered_square_share':rat(nr**2,total)})
f=20e6+np.arange(501)*300e3;response=np.empty((501,11),complex)
for i in range(0,len(f),32):
 ff=f[i:i+32];kernel=np.exp(-2j*np.pi*ff[:,None]*t)
 denom=(kernel@source)*np.exp(-2j*np.pi*ff*offset)
 response[i:i+len(ff)]=(kernel@u)/denom[:,None]
removed=np.fft.ifft(response,n=2004,axis=0)*4
zp=mp.parent/m['array_archive']['path'];assert hashlib.sha256(zp.read_bytes()).hexdigest()==m['array_archive']['sha256']
with np.load(zp,allow_pickle=False) as z:
 env={}
 for taper in ['no_taper','tail_200ns']:
  rows=sorted([q for q in m['records'] if q['family']=='CO' and q['role']=='BG' and q['taper']==taper],key=lambda q:q['receiver_positions_xyz_m'][0][1])
  env[taper]=np.column_stack([z[q['official_reconstruction']['output_keys']['complex_envelope']] for q in rows])
  tt=z[rows[0]['official_reconstruction']['output_keys']['envelope_time_s']]
a=center(env['no_taper']);b=center(env['tail_200ns']);v=center(removed)
rec=[]
for win in c['reconstruction_windows']:
 ix=(tt>=win['start_s'])&(tt<win['stop_s']);na=norm(a[ix]);nb=norm(b[ix]);nv=norm(v[ix])
 rec.append({'window':win['id'],'n':int(ix.sum()),'before_norm':na,'after_norm':nb,'removed_norm':nv,'after_over_before':rat(nb,na),'closure_over_before':rat(norm(a[ix]-b[ix]-v[ix]),na)})
result={'contract_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'raw':raw,'reconstruction':rec,
 'raw_unassigned':{'count':int(np.sum(coverage==0)),'times_s':t[coverage==0].tolist(),'centered_square_share':rat(norm(r[coverage==0])**2,total)},
 'taper_sample_count':count,'first_changed_sample_s':float(t[np.flatnonzero(w!=1)[0]]),
 'raw_centered_square_share_in_taper_support':rat(norm(r[w!=1])**2,total),
 'raw_centered_norm':norm(r),'raw_centered_norm_after_taper':norm(wr),
 'transform_removed_vs_saved_change_relative':rat(norm(removed-(env['no_taper']-env['tail_200ns'])),norm(env['no_taper']-env['tail_200ns'])),
 'per_trace':[{'rx_y':rr[j]['receiver_positions_xyz_m'][0][1],'centered_norm':norm(r[:,j]),'last_sample':float(x[-1,j]),'raw_mean':float(x[:,j].mean())} for j in range(11)]}
out=Path(__file__).with_suffix('.json');out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
np.savez_compressed(Path(__file__).with_suffix('.npz'),removed_response=response,removed_envelope=removed,raw_centered=r,raw_time_s=t,weights=w)
print(json.dumps(result,indent=2))
