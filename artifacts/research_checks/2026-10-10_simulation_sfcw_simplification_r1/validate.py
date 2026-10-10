"""Reproduce the cleanup evidence using existing local raw/array files; no solve."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import gprMax
import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
relative=lambda x,y:float(np.linalg.norm(x-y)/np.linalg.norm(y))
old=ROOT/'artifacts/local_checks/2026-10-10_cover_plane_results_r1/arrays_r1.npz'
new=ROOT/'artifacts/local_checks/2026-10-10_cover_plane_results_r2/arrays_r2.npz'
raw=ROOT/'artifacts/local_checks/2026-10-10_cover_plane_results_r1/source/profile.h5'
cli=new.parent/'official_main_hann.h5'
with h5py.File(raw) as h:
    offset=float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
    dt=float(h.attrs['dt'])
with np.load(old) as a,np.load(new) as b,h5py.File(cli) as h:
    native={key:bool(np.array_equal(a[key],b[key])) for key in ['time_s','Ez','Hx','Hy','source','rx']}
    main={key:bool(np.array_equal(a[key],b[key])) for key in ['main_hann','main_blackman']}
    assert all(native.values()) and all(main.values())
    phase=np.exp(2j*np.pi*b['frequency_Hz']*offset)[:,None,None]
    phase_errors={key:relative(b['spectra_'+key],a['spectra_'+key]*phase)
        for key in ['full','native_40_230','native_180_400','native_220_380']}
    assert max(phase_errors.values())<1e-9
    official_checks={
        'frequency_equal':bool(np.array_equal(h['frequency'][:],b['frequency_Hz'])),
        'time_equal':bool(np.array_equal(h['time_response/time'][:],b['sfcw_time_s'])),
        'complex_bandpass_over_length_relative_error':relative(h['time_response/complex_bandpass'][:]/.025,b['main_hann']),
        'real_bandpass_over_2_length_relative_error':relative(h['time_response/real_bandpass'][:]/(2*.025),b['main_hann'].real),
        'all_source_tones_valid':bool(h['source_valid'][:].all()),
        'receiver_tail_relative_db':float(h.attrs['ReceiverTailRelativeDB'])}
    assert official_checks['frequency_equal'] and official_checks['time_equal'] and official_checks['all_source_tones_valid']
    assert official_checks['complex_bandpass_over_length_relative_error']<1e-12
    assert official_checks['real_bandpass_over_2_length_relative_error']<1e-12
r1=ROOT/'artifacts/research_checks/2026-10-10_cover_plane_closure_r1/analysis.json'
r2=ROOT/'artifacts/research_checks/2026-10-10_cover_plane_closure_r2/analysis.json'
left=json.loads(r1.read_text('utf-8'));right=json.loads(r2.read_text('utf-8'))
primary=json.loads(r2.with_name('contract.json').read_text('utf-8'))['primary']
rows=[];largest=0.
for old_row,row in zip(left['metrics'],right['metrics']):
    assert {k:v for k,v in old_row.items() if k!='gates'}=={k:v for k,v in row.items() if k!='gates'}
    for gate in row['gates']:
        for key,value in row['gates'][gate].items():
            largest=max(largest,abs(value-old_row['gates'][gate][key]))
    if row['variant']=='full' and all(row[k]==v for k,v in primary.items()):
        rows.append({'pair':row['pair'],'window':row['window'],
            'old_late_up_relative_error':old_row['gates']['late']['up_relative_error'],
            'new_late_up_relative_error':row['gates']['late']['up_relative_error'],
            'new_early_down_relative_error':row['gates']['early']['down_relative_error'],
            'new_late_sector_change':row['gates']['late']['sector_raw_relative_change']})
test=subprocess.run([sys.executable,str(ROOT/'scripts/check_line9_cover_sfcw.py')],cwd=ROOT,capture_output=True,text=True)
assert test.returncode==0,test.stderr
(OUT/'regression.log').write_text(test.stdout+test.stderr,encoding='utf-8')
value={'status':'PASS_CPU_OFFICIAL_SFCW_CLEANUP','new_solver_runs':0,'gprmax_version':gprMax.__version__,
    'source_time_offset_s':offset,'dt_s':dt,'omitted_phase_deg_at_20_170MHz':list(360*np.array([20e6,170e6])*offset),
    'unchanged_native_arrays':native,'unchanged_main_profiles':main,'corrected_spectra_phase_relation_errors':phase_errors,
    'official_CLI_checks':official_checks,'metric_rows_compared':len(rows),'all_configuration_rows_compared':len(right['metrics']),
    'max_absolute_change_over_all_gate_metrics':largest,'primary_rows':rows,'regression_tests_passed':3,
    'identities':{str(p.relative_to(ROOT)):sha(p) for p in [old,new,raw,cli,r1,r2,Path(__file__),ROOT/'scripts/check_line9_cover_sfcw.py']},
    'official_processing_sha256':sha(sf.__file__),
    'limits':'Array/official-library agreement only; no field calibration, full-dimensional geology or unique path certification.'}
assert not (OUT/'validation.json').exists()
(OUT/'validation.json').write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
lines=['|平面对|窗|更正前晚窗上行误差|更正后晚窗上行误差|更正后早窗下行误差|更正后扇区复场变化|',
       '|---|---|---:|---:|---:|---:|']
for r in rows:
    numbers=[r[k]*100 for k in ['old_late_up_relative_error','new_late_up_relative_error','new_early_down_relative_error','new_late_sector_change']]
    lines.append('|'+str(r['pair'])+'|'+r['window']+'|'+'|'.join(f'{v:.3f}%' for v in numbers)+'|')
(OUT/'primary_metrics.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({k:value[k] for k in ['status','official_CLI_checks','corrected_spectra_phase_relation_errors','primary_rows']}))
