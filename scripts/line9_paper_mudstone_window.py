"""Authorized 190 m paper-mudstone H0/H1, 2400 ns; preserve historical inputs."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from gprMax.material_database import load_material_spec
import line9_v401_version_controls as native

APPROVAL = '2026-10-11 用户明确批准：只改泥岩，190m同站位H0/H1两项，2400ns与前1200ns比较；其余条件及官方impulse/direct保持。'


def prepare(a):
    assert not a.out.exists()
    m = json.loads((a.package/'manifest.json').read_text('utf-8'))
    db = json.loads(a.material.read_text('utf-8'))
    assert db['database']['id'] == 'line9_materials_v0_2'
    mud = db['materials']['material_002_mudstone']
    assert mud['model'] == 'constant' and not mud.get('poles')
    assert mud['base'] == dict(relative_permittivity=18., electric_conductivity_s_per_m=.006,
                              relative_permeability=1., magnetic_conductivity_s_per_m=0.)
    assert [g['id'] for g in m['groups']] == ['high_x19000_H0','high_x19000_H1']
    a.out.mkdir(parents=True)
    checks = []
    for g in m['groups']:
        old = copy.deepcopy(g)
        for key in ['input','geometry','material']:
            assert native.sha(a.package/g[key]) == g[key+'_sha256']
            dst = a.out/g[key]; dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(a.package/g[key], dst)
        olddb = json.loads((a.package/g['material']).read_text('utf-8'))
        for key in ['material_000_air','material_001_cover','material_003_sandstone']:
            assert olddb['materials'][key] == db['materials'][key]
        material_path = a.out/g['material']
        g['material'] = str(Path(g['material']).with_name('line9_materials_v0_2.json')).replace('\\','/')
        material_path.rename(a.out/g['material'])
        native.save(a.out/g['material'], db)
        for key in db['materials']:
            load_material_spec('line9_materials_v0_2', key, search_directory=(a.out/g['geometry']).parent)
        with h5py.File(a.out/g['geometry'],'r+') as h:
            h.attrs['MaterialDatabase'] = db['database']['id']
        with h5py.File(a.package/g['geometry']) as oldh, h5py.File(a.out/g['geometry']) as newh:
            assert list(oldh) == list(newh)
            for key in oldh:
                np.testing.assert_array_equal(oldh[key][:],newh[key][:])
            assert set(oldh.attrs) == set(newh.attrs)
            for key in oldh.attrs:
                if key != 'MaterialDatabase':
                    np.testing.assert_array_equal(oldh.attrs[key],newh.attrs[key])
        p = a.out/g['input']; before = p.read_text('utf-8')
        after = before.replace('#time_window: 1.2e-06','#time_window: 2.4e-06').replace(
            'line9_research_materials_v1_smoothed y','line9_materials_v0_2 y')
        changed = [(x,y) for x,y in zip(before.splitlines(),after.splitlines()) if x != y]
        assert len(changed) == 2 and '#waveform: impulse 1 1 impulse' in after
        p.write_text(after,encoding='utf-8')
        # Remove old-material peak-derived labels; they are not new ground truth.
        for key in ['peak_ns','basal_gate_ns','wide_gate_ns']: g.pop(key,None)
        g.update(time_window_ns=2400.,expected_samples=int(np.ceil(2.4e-6/m['dt_s']))+1,
                 comparison_crop_samples=old['expected_samples'], parent_input_sha256=old['input_sha256'],
                 parent_geometry_sha256=old['geometry_sha256'], parent_material_sha256=old['material_sha256'])
        for key in ['input','geometry','material']: g[key+'_sha256'] = native.sha(a.out/g[key])
        checks.append(dict(id=g['id'],unchanged_voxels=True,unchanged_other_materials=True,changed_input_lines=changed))
    # A vertical normal-incidence delay guide, not a ray-traced B-scan truth.
    with h5py.File(a.out/m['groups'][1]['geometry']) as h:
        column=h['data'][6800,:,0]
        thickness={int(k):float(np.count_nonzero(column==k)*.025) for k in [1,2]}
    def eps(key,f):
        p=db['materials'][key]; b=p['base']
        return b['relative_permittivity']-1j*b['electric_conductivity_s_per_m']/(2*np.pi*f*8.8541878128e-12)+sum(
            q['relative_permittivity_difference']/(1+2j*np.pi*f*q['relaxation_time_s']) for q in p.get('poles',[]))
    f=95e6; delay=2*(((38.95+39.025)/2-31.)/299792458.)
    for k,key in [(1,'material_001_cover'),(2,'material_002_mudstone')]:
        n=np.sqrt(eps(key,f)); derivative=(np.sqrt(eps(key,f+1000))-np.sqrt(eps(key,f-1000)))/2000
        delay += 2*thickness[k]*(n+f*derivative).real/299792458.
    m.pop('source_only_checks',None)
    m.update(status='PREPARED_APPROVED_PAPER_MUDSTONE_2400NS',approval_basis=APPROVAL,execution_authorized=True,
             generator_sha256=native.sha(__file__),parent_manifest_sha256=native.sha(a.package/'manifest.json'),
             material_database_sha256=native.sha(a.material),preparation_checks=checks,
             predeclared_analysis_windows_ns={'early':[0,120],'deep_broad':[250,550],'historical_deep':[300,450]},
             vertical_95MHz_group_delay_guide_ns=float(delay*1e9),vertical_thickness_m=thickness,
             limits='Two configurations at one station; no spatial B-scan or unique path truth. Only mudstone constitutive law and record length change. H0/H1 contrast is connected bottom sandstone, not all-layer clean truth. Air reduction shelved. No extra observers/snapshots, fit, taper, AGC or new method.')
    native.save(a.out/'manifest.json',m)
    print(json.dumps(dict(groups=2,samples=m['groups'][0]['expected_samples'],delay_guide_ns=m['vertical_95MHz_group_delay_guide_ns'])))


def freeze(a):
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert m['approval_basis'] == APPROVAL and m['execution_authorized']
    assert len(m['groups']) == 2 and all(g['time_window_ns']==2400 for g in m['groups'])
    a.action='prepare'; native.prepare(a)
    p=a.out/'execution_contract.json'; c=json.loads(p.read_text('utf-8'))
    c['approval_basis']=APPROVAL; c['study_manifest']=m
    c['limits']=m['limits']; c['code_identities'][str(Path(__file__).resolve())]=native.sha(__file__)
    native.save(p,c); native.save(a.out/'preflight_verification.json',native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p),groups=2,version=c['expected_version'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['prepare','freeze','run','verify'])
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--package',type=Path)
    p.add_argument('--material',type=Path)
    a=p.parse_args(); a.out=a.out.resolve()
    if a.package: a.package=a.package.resolve()
    if a.action=='prepare': prepare(a)
    elif a.action=='freeze': freeze(a)
    elif a.action=='run':
        native.supervisor.audit=native.audit
        native.supervisor.run(a.out/'execution_contract.json')
    else: print(json.dumps(native.audit(a.out/'execution_contract.json',True)))
