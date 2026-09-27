"""Accept the official tail decomposition against independent H5/DFT calculations."""
from pathlib import Path
import hashlib,json,math
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'START_HERE.md').exists())
OUT=ROOT/'artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition'
r=json.loads((OUT/'check_co_tail_independent.json').read_text(encoding='utf-8'))
a=json.loads((OUT/'results.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert a['contract_sha256']==r['contract_sha256']==sha(ROOT/a['contract_path'])
assert a['script_sha256']==sha(ROOT/a['script_path'])
assert a['diagnostics_npz']['sha256']==sha(OUT/a['diagnostics_npz']['path'])
checks=[]
def same(x,y,label):
 assert math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-20),(label,x,y)
 checks.append(label)
for raw,official in zip(r['raw'],a['fixed_raw_windows']):
 assert raw['window']==official['window_id']
 for x,y in [('n','sample_count'),('centered_norm','centered_l2')]:same(raw[x],official[y],raw['window']+x)
 same(raw['centered_square_share'],official['share_of_full_raw_centered_squared_norm']['value'],raw['window']+'share')
for independent,official in zip(r['reconstruction'],a['reconstruction']['windows']):
 assert independent['window']==official['window_id']
 for x,y in [('n','sample_count'),('before_norm','centered_before_norm'),('after_norm','centered_after_taper_norm')]:same(independent[x],official[y],independent['window']+x)
 same(independent['after_over_before'],official['after_over_before_norm']['value'],independent['window']+'after_ratio')
 # The independently computed DFT uses a different accumulation order.
 assert math.isclose(independent['removed_norm'],official['centered_removed_tail_norm'],rel_tol=1e-8,abs_tol=1e-20)
 checks.append(independent['window']+'removed_norm')
for independent,official in zip(r['per_trace'],a['per_trace']):
 for x,y in [('rx_y','receiver_y_m'),('centered_norm','centered_l2_full_record'),('last_sample','raw_last_sample_v_per_m'),('raw_mean','raw_mean_v_per_m')]:same(independent[x],official[y],str(independent['rx_y'])+x)
same(r['raw_centered_square_share_in_taper_support'],a['raw_global']['actual_taper_nonzero_support']['centered_squared_norm_share_of_full_raw']['value'],'support_share')
same(r['raw_unassigned']['count'],a['raw_global']['fixed_windows_unassigned_samples']['count'],'unassigned_count')
same(r['raw_unassigned']['centered_square_share'],a['raw_global']['fixed_windows_unassigned_samples']['centered_squared_norm_share_of_full_raw']['value'],'unassigned_share')
input_specs={}
for rel in a['input_h5_sha256']:
 p=(ROOT/rel).with_suffix('.in');lines=p.read_text(encoding='utf-8').splitlines()
 assert '#domain: inf 32 50' in lines and '#dx_dy_dz: 0.025 0.025 0.025' in lines
 assert '#pml_cells: 0 40 40 0 40 40' in lines
 invariant=[s for s in lines if not s.startswith(('#title:','#hertzian_dipole:','#rx:'))]
 if input_specs:assert invariant==baseline
 else:baseline=invariant
 input_specs[p.relative_to(ROOT).as_posix()]=sha(p)
arrays=[]
with np.load(OUT/'check_co_tail_independent.npz') as i,np.load(OUT/'diagnostics.npz') as o:
 for x,y in [('raw_time_s','raw_time_s'),('raw_centered','raw_centered_residual'),('weights','official_tail_weight_W')]:
  assert np.array_equal(i[x],o[y]),x
  arrays.append({'name':x,'count':int(i[x].size),'bit_identical':True})
 for x,y in [('removed_response','removed_tail_frequency_response'),('removed_envelope','removed_tail_complex_envelope')]:
  error=float(np.linalg.norm(i[x]-o[y])/np.linalg.norm(o[y]));assert error<1e-8,(x,error)
  arrays.append({'name':x,'complex_values':int(i[x].size),'relative_l2_difference':error})
result={'status':'passed_arithmetic_and_input_geometry_checks_only','scalar_checks':len(checks),'arrays':arrays,
 'scalar_rtol':1e-11,'scalar_atol':1e-20,'independent_dft_relative_tolerance':1e-8,
 'tolerance_scope':'numerical implementation comparison only, not solver/physical error certification',
 'input_specifications_sha256':input_specs,'results_sha256':sha(OUT/'results.json'),
 'independent_script_sha256':sha(OUT/'check_co_tail_independent.py'),
 'independent_json_sha256':sha(OUT/'check_co_tail_independent.json'),
 'independent_npz_sha256':sha(OUT/'check_co_tail_independent.npz'),
 'PML_physical_cause_certified':False,'solver_executed':False,'training_executed':False}
(OUT/'root_acceptance.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
