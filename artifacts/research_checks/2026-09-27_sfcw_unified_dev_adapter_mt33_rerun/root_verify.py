"""Independent H5 -> explicit DFT and Fourier-sum acceptance, no official transform calls."""
import hashlib,json
from importlib.util import find_spec
from pathlib import Path
import h5py
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'START_HERE.md').is_file())
OUT=ROOT/'artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun'
m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert m['inputs_n']==24 and m['records_n']==48
assert sha(Path(find_spec('gprMax.toolboxes.SFCW.processing').origin))==m['api']['processing_source_sha256']
assert sha(OUT/m['array_archive']['path'])==m['array_archive']['sha256']
assert sha(ROOT/m['api']['adapter_script_path'])==m['api']['adapter_script_sha256']
assert sha(ROOT/m['acquisition_contract']['path'])==m['acquisition_contract']['sha256']
f=np.arange(501)*300000.+20000000.
fi=np.array([0,1,2,3,100,250,400,499,500])
ti=np.array([0,1,17,60,180,420,720,1200,2003])
max_dft=max_inverse=0.; dft_values=0;tail=[]; hashes={}; anchors={};rows=[];gathers={}
with np.load(OUT/m['array_archive']['path'],allow_pickle=False) as z:
 for r in m['records']:
  path=ROOT/r['input_path'];assert 'C3m' in path.name and 'C5m' not in path.name and 'C8m' not in path.name
  hashes.setdefault(str(path),sha(path));assert hashes[str(path)]==r['input_sha256']
  with h5py.File(path) as h:
   sd=h['srcs/src1/excitation'];s=sd['samples'][:];dt=float(sd.attrs['SampleInterval']);so=float(sd.attrs['TimeSampleOffset'])
   names=['rx'+str(i) for i in range(1,34)] if r['family']=='MT' else ['rx1']
   x=np.stack([h['rxs/'+k+'/Ex'][:] for k in names],axis=1)
   ro=float(h['rxs/'+names[0]+'/Ex'].attrs['TimeSampleOffset'])
   assert x.shape==(20352,33 if r['family']=='MT' else 1)
   pos=np.array([h['rxs/'+k].attrs['Position'] for k in names]);srcpos=h['srcs/src1'].attrs['Position']
   if r['family']=='MT':
    np.testing.assert_allclose(pos[:,1],12.65+.25*np.arange(33),rtol=0,atol=1e-12)
    assert abs(srcpos[1]-15.35)<1e-12
   else: assert abs(pos[0,1]-srcpos[1]-1.3)<1e-12
  np.testing.assert_array_equal(np.asarray(r['receiver_positions_xyz_m']),pos)
  np.testing.assert_array_equal(np.asarray(r['source_position_xyz_m']),srcpos)
  assert r['source_dt_s']==dt and r['source_time_offset_s']==so and r['receiver_time_offset_s']==ro
  prefix=f"{r['family']}_{r['role']}_{path.stem}"
  np.testing.assert_array_equal(z[prefix+'_raw_source_time_s'],so+np.arange(len(s))*dt)
  np.testing.assert_array_equal(z[prefix+'_raw_receiver_time_s'],ro+np.arange(len(x))*dt)
  expected_tail=20*np.log10(np.max(abs(x[-int(np.ceil(.05*len(x))):]))/np.max(abs(x)))
  assert abs(expected_tail-r['tail_relative_db_before_taper'])<1e-10
  tail.append(float(expected_tail))
  if r['taper']=='tail_200ns':
   count=round(200e-9/dt)
   assert abs(r['tail_taper_fraction']-(count-.25)/len(x))<1e-15
   x[-count:]*=((1+np.cos(np.linspace(0,np.pi,count)))/2)[:,None]
  else: assert r['tail_taper_fraction']==0
  # Direct summation includes physical source/receiver sample origins.
  ss=dt*(np.exp(-2j*np.pi*f[fi,None]*(np.arange(len(s))*dt+so))@s)
  rr=dt*(np.exp(-2j*np.pi*f[fi,None]*(np.arange(len(x))*dt+ro))@x)
  expected=rr/ss[:,None]
  keys=r['official_reconstruction']['output_keys'];h=z[keys['response']];h=h[:,None] if h.ndim==1 else h
  assert h.shape==(501,len(names)) and np.isfinite(h).all()
  gathers.setdefault((r['family'],r['role'],r['taper']),[]).append(h)
  assert z[keys['source_valid']].shape==(501,) and z[keys['source_valid']].all()
  np.testing.assert_array_equal(z[keys['window_weights']],np.ones(501))
  np.testing.assert_allclose(z[keys['source_spectrum']][fi],ss,rtol=1e-12,atol=0)
  np.testing.assert_array_equal(z[keys['frequency_hz']],f)
  error=float(np.max(abs(h[fi]-expected))/np.max(abs(expected)))
  assert error<1e-8,error
  max_dft=max(max_dft,error);dft_values+=expected.size
  t=z[keys['envelope_time_s']];np.testing.assert_allclose(t,np.arange(2004)/(2004*300000.),rtol=1e-14,atol=1e-20)
  env=z[keys['complex_envelope']];env=env[:,None] if env.ndim==1 else env
  inv=(np.exp(2j*np.pi*np.outer(ti,np.arange(501))/2004)@h)/501
  ie=float(np.max(abs(env[ti]-inv))/np.max(abs(inv)));assert ie<1e-10,ie;max_inverse=max(max_inverse,ie)
  band=z[keys['complex_bandpass']];band=band[:,None] if band.ndim==1 else band
  np.testing.assert_allclose(band,env*np.exp(2j*np.pi*20e6*t[:,None]),rtol=1e-12,atol=0)
  real=z[keys['real_bandpass']];real=real[:,None] if real.ndim==1 else real
  np.testing.assert_allclose(real,2*(env*np.exp(2j*np.pi*20e6*t[:,None])).real,rtol=1e-12,atol=np.max(abs(real))*1e-12)
  if r['family']=='MT': anchors[(r['role'],r['taper'],'MT')]=h[:,16]
  elif '-t05' in path.stem: anchors[(r['role'],r['taper'],'CO')]=h[:,0]
  rows.append({'input':r['input_path'],'taper':r['taper'],'DFT_max_error_over_selected_peak':error,'inverse_max_error_over_selected_peak':ie})
 for role in ['BG','TGT']:
  for taper in ['no_taper','tail_200ns']:
   np.testing.assert_allclose(anchors[(role,taper,'CO')],anchors[(role,taper,'MT')],rtol=1e-12,atol=0)
sensitivity=[]
for family in ['MT','CO']:
 for role in ['BG','TGT']:
  a=np.column_stack(gathers[(family,role,'no_taper')]);b=np.column_stack(gathers[(family,role,'tail_200ns')])
  sensitivity.append({'family':family,'role':role,'relative_L2_200ns_change_from_no_taper':float(np.linalg.norm(b-a)/np.linalg.norm(a))})
result={'tail_taper_sensitivity':sensitivity,'passed':True,'source_H5_files':len(hashes),'records_checked':len(rows),'complex_DFT_values_checked':dft_values,'method':'explicit complex exponential sums from H5 with physical time origins; inverse Fourier sums at 9 samples for every trace','max_DFT_error_over_selected_peak':max_dft,'max_inverse_error_over_selected_peak':max_inverse,'CO_t05_MT_rx17_anchor_checks':4,'raw_tail_dB_range':[min(tail),max(tail)],'record_decay_certified':False,'solver_invoked':False,'rows':rows}
(OUT/'root_acceptance.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
