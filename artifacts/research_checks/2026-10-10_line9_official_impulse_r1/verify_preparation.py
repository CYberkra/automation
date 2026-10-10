"""Independently check source-only cards and archive text evidence; no FDTD."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import gprMax
import h5py
from gprMax.hash_cmds_file import get_user_objects

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
from line9_v401_version_controls import audit_source

OUT=Path(__file__).resolve().parent
PARENT=ROOT/'artifacts/local_checks/2026-10-08_v401_dense_loss_prepared_r1'
PACKAGE=ROOT/'artifacts/local_checks/2026-10-10_line9_official_impulse_package_r1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
old=json.loads((PARENT/'manifest.json').read_text('utf-8'))
new=json.loads((PACKAGE/'manifest.json').read_text('utf-8'))
parents={g['id']:g for g in old['groups']}
assert new['execution_authorized'] is False and new['reused_solver_outputs']==[]
assert not (PACKAGE/'execution_contract.json').exists()
assert not list(PACKAGE.rglob('profile.h5'))
rows=[]
for g in new['groups']:
    original=parents[g['id']];source=PARENT/original['input'];card=PACKAGE/g['input']
    left=source.read_text('utf-8').splitlines();right=card.read_text('utf-8').splitlines()
    changed=[i for i,(a,b) in enumerate(zip(left,right)) if a!=b]
    assert len(left)==len(right) and len(changed)==2
    assert right[changed[0]]=='#waveform: impulse 1 1 impulse'
    assert left[changed[1]].split()[:-1]==right[changed[1]].split()[:-1]
    assert right[changed[1]].split()[-1]=='impulse'
    for key in ['geometry','material']:
        assert sha(PARENT/original[key])==sha(PACKAGE/g[key])==g[key+'_sha256']
    objects=get_user_objects(right,input_dir=card.parent)
    wave=next(o.kwargs for o in objects if type(o).__name__=='Waveform')
    dipole=next(o.kwargs for o in objects if type(o).__name__=='HertzianDipole')
    assert wave=={'wave_type':'impulse','amp':1.,'freq':1.,'id':'impulse'}
    assert dipole['waveform_id']=='impulse' and dipole['start'] is None
    dst=OUT/'inputs'/g['id']/'profile.in';dst.parent.mkdir(parents=True,exist_ok=True)
    assert not dst.exists();dst.write_bytes(card.read_bytes())
    rows.append(dict(id=g['id'],changed_line_numbers_1based=[i+1 for i in changed],
        input_sha256=sha(card),parent_input_sha256=sha(source),
        geometry_sha256_unchanged=g['geometry_sha256'],material_sha256_unchanged=g['material_sha256'],
        official_waveform=wave,official_dipole=dipole))
raw=ROOT/'artifacts/local_checks/2026-10-10_cover_plane_results_r1/source/profile.h5'
with h5py.File(raw) as h:
    audit_source(h['srcs/src1/excitation'],{'source_type':'ricker','source_frequency_Hz':100e6,'source_amplitude_A':40.},float(h.attrs['dt']))
test=subprocess.run([sys.executable,str(ROOT/'scripts/check_line9_official_impulse.py')],cwd=ROOT,capture_output=True,text=True)
assert test.returncode==0,test.stderr
(OUT/'regression.log').write_text(test.stdout+test.stderr,encoding='utf-8')
assert not (OUT/'manifest.json').exists()
(OUT/'manifest.json').write_bytes((PACKAGE/'manifest.json').read_bytes())
value=dict(status='PASS_SOURCE_ONLY_PREPARATION_OFFICIAL_PARSE_AND5_CPU_TESTS',new_solver_runs=0,
    gprmax_version=gprMax.__version__,prepared_inputs=rows,regression_tests_passed=5,
    actual_legacy_Ricker_source_audit_passed=True,legacy_raw_sha256=sha(raw),
    source_design_sha256=sha(ROOT/'configs/research/line9_source_v0_2.json'),
    preparation_manifest_sha256=sha(PACKAGE/'manifest.json'),generator_sha256=sha(ROOT/'scripts/prepare_line9_official_impulse.py'),
    native_auditor_sha256=sha(ROOT/'scripts/line9_v401_version_controls.py'),
    verification_script_sha256=sha(__file__),
    official_example_sha256=sha(Path(gprMax.__file__).parent/'toolboxes/SFCW/examples/cylinder_sfcw_2D.in'),
    official_waveforms_sha256=sha(Path(gprMax.__file__).parent/'waveforms.py'),
    execution_contract_created=False,impulse_receiver_output_available=False,
    limits='Input/CPU validation only; no FDTD response or physics/field validation.')
assert not (OUT/'validation.json').exists()
(OUT/'validation.json').write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(status=value['status'],prepared_inputs=len(rows),new_solver_runs=0)))
