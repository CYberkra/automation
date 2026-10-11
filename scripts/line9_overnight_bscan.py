"""Approved 190->170 m H1 line, two new H0 anchors; existing solver unchanged."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from gprMax.material_database import load_material_spec
import line9_v401_version_controls as native

APPROVAL = ('2026-10-11 用户明确批准上述批次，跑完验收并推送：190→170m，'
            '41道H1、0.5m间距，复用190m H0/H1，新增40项H1及180/170m两项H0；'
            '先1m骨架后半米补齐；保持现行材料/210×42.5m域/2.5cm/2400ns/'
            '约8m离地/4.0.1原生FP64/官方impulse/direct501/Hann及Blackman；总时限6小时。')


def guide(column, db, air_height):
    """Prospective local vertical group-delay guide, not a ray-tracing truth."""
    def index(key, f):
        p = db['materials'][key]; b = p['base']
        eps = b['relative_permittivity'] - 1j*b['electric_conductivity_s_per_m']/(2*np.pi*f*8.8541878128e-12)
        eps += sum(q['relative_permittivity_difference']/(1+2j*np.pi*f*q['relaxation_time_s']) for q in p.get('poles', []))
        return np.sqrt(eps)
    f = 95e6
    time = 2*air_height/299792458.
    result = {'surface': time*1e9}
    for label, material, key in [('cover_base', 1, 'material_001_cover'), ('basal_sandstone', 2, 'material_002_mudstone')]:
        n = index(key, f); derivative = (index(key, f+1000)-index(key, f-1000))/2000
        time += 2*np.count_nonzero(column == material)*.025*(n+f*derivative).real/299792458.
        result[label] = float(time*1e9)
    return result


def prepare(a):
    assert not a.out.exists()
    old = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert len(old['groups']) == 2 and old['SFCW']['method'] == 'direct'
    templates = {g['role']: g for g in old['groups']}
    db = json.loads((a.package/templates['H1']['material']).read_text('utf-8'))
    assert db['database']['id'] == 'line9_materials_v0_2'
    assert db['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m'] == .006
    # Reuse requires the completed native certificate, source files and every identity.
    prior_c = json.loads((a.prior/'execution_contract.json').read_text('utf-8'))
    prior_v = json.loads((a.prior/'completed_verification.json').read_text('utf-8'))
    assert prior_v['completed'] and prior_v['contract_sha256'] == native.sha(a.prior/'execution_contract.json')
    records = {r['id']: r for r in prior_v['groups']}
    reused = []
    for g in prior_c['groups']:
        p = a.prior/(g['id']+'.h5')
        assert native.sha(p) == records[g['id']]['native_sha256']
        assert g['geometry_sha256'] == templates[g['role']]['geometry_sha256']
        assert g['material_sha256'] == templates[g['role']]['material_sha256']
        with h5py.File(p) as h:
            assert h.attrs['gprMax'] == '4.0.1' and h['rxs/rx1/Ez'].dtype == np.float64
            native.audit_source(h['srcs/src1/excitation'], g, float(h.attrs['dt']))
        reused.append(dict(g, native_path=str(p.resolve()), native_sha256=native.sha(p),
                           parent_contract_sha256=prior_v['contract_sha256']))
    a.out.mkdir(parents=True)
    for role, g in templates.items():
        for key in ['input', 'geometry', 'material']:
            assert native.sha(a.package/g[key]) == g[key+'_sha256']
        geometry = a.out/'geometries'/f'{role}.h5'; geometry.parent.mkdir(exist_ok=True)
        shutil.copyfile(a.package/g['geometry'], geometry)
    material = a.out/'geometries/line9_materials_v0_2.json'
    shutil.copyfile(a.package/templates['H1']['material'], material)
    for key in db['materials']:
        load_material_spec(db['database']['id'], key, search_directory=material.parent)
    with h5py.File(a.out/'geometries/H1.h5') as h:
        data = h['data'][:, :, 0]
        np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025]*3)
        assert h.attrs['MaterialDatabase'] == db['database']['id']
    stations = []
    for x in np.arange(190., 169.99, -.5):
        positions = []
        for offset in [-.65, .65]:
            internal_x = float(x-20+offset)
            ground = (np.flatnonzero(data[round(internal_x/.025)] != 0)[-1]+1)*.025
            positions.append([internal_x, float(ground+8.), .0125])
        assert min(40.5-p[1] for p in positions) > 0
        col = data[round((x-20)/.025)]
        ground = (np.flatnonzero(col != 0)[-1]+1)*.025
        stations.append(dict(profile_x_m=float(x), distance_from_220m_m=float(220-x),
            input_tx_m=positions[0], input_rx_m=positions[1],
            top_PML_min_clearance_m=min(40.5-p[1] for p in positions),
            guide_ns=guide(col, db, (positions[0][1]+positions[1][1])/2-ground),
            interfaces_y_m={str(int(v)):float((np.flatnonzero(col == v)[-1]+1)*.025) for v in [1,2,3]}))
    station_map = {s['profile_x_m']: s for s in stations}
    for g in reused:
        s = station_map[190.]
        for key in ['tx', 'rx']:
            np.testing.assert_allclose(g[key+'_m'], s['input_'+key+'_m'][:2]+[0.], rtol=0, atol=1e-12)
    # Anchor pairs first, integer skeleton second, then half-metre infill.
    order = [(180., 'H1'), (180., 'H0'), (170., 'H1'), (170., 'H0')]
    order += [(float(x), 'H1') for x in range(189, 170, -1) if x != 180]
    order += [(float(x), 'H1') for x in np.arange(189.5, 170., -1)]
    assert len(order) == 42 and len(set(order)) == 42
    groups = []
    for j, (x, role) in enumerate(order):
        g = copy.deepcopy(templates[role]); s = station_map[x]
        for key in list(g):
            if key.startswith('parent_') or key in ['comparison_crop_samples', 'chainage_m']:
                g.pop(key)
        name = f'line_x{round(x*100):05d}_{role}'
        g.update(id=name, profile_x_m=x, distance_from_220m_m=220-x,
                 phase='anchor_pairs' if j < 4 else 'integer_skeleton' if j < 22 else 'half_metre_infill',
                 input=f'cases/{name}/profile.in', geometry=f'geometries/{role}.h5',
                 material='geometries/line9_materials_v0_2.json')
        lines = []
        for line in (a.package/templates[role]['input']).read_text('utf-8').splitlines():
            if line.startswith('#title:'):
                line = f'#title: {name}; approved overnight spatial line; no snapshots'
            elif line.startswith('#hertzian_dipole:'):
                line = '#hertzian_dipole: z '+' '.join(f'{v:.15g}' for v in s['input_tx_m'])+' impulse'
            elif line.startswith('#rx:'):
                line = '#rx: '+' '.join(f'{v:.15g}' for v in s['input_rx_m'])+f' {name}_rx1 Ez'
            elif line.startswith('#geometry_objects_read:'):
                line = f'#geometry_objects_read: 0 0 0 ../../geometries/{role}.h5 line9_materials_v0_2 y'
            lines.append(line)
        p = a.out/g['input']; p.parent.mkdir(parents=True)
        p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        for key in ['input','geometry','material']:
            g[key+'_sha256'] = native.sha(a.out/g[key])
        for key in ['tx','rx']:
            g['input_'+key+'_m'] = s['input_'+key+'_m']
            g[key+'_m'] = s['input_'+key+'_m'][:2]+[0.]
        g['spacing_m'] = [.025]*3
        groups.append(g)
    manifest = dict(status='PREPARED_APPROVED_OVERNIGHT_LINE', approval_basis=APPROVAL,
        execution_authorized=True, groups=groups, reused_solver_outputs=reused,
        parent_manifest_sha256=native.sha(a.package/'manifest.json'), generator_sha256=native.sha(__file__),
        dt_s=old['dt_s'], SFCW=old['SFCW'], stations=stations, no_retry=True,
        predeclared_analysis_windows_ns=old['predeclared_analysis_windows_ns'],
        limits='41 H1 stations and 3 matched H0 anchors only; anchor differences shown separately, no interpolated H0 or synthetic full-line delta. Vertical guides are prospective local estimates, not ray-tracing truth. No AGC, fit, new methods, snapshots, extra observers, material, mesh, PML or domain changes. No field validation or spatial convergence claim.')
    native.save(a.out/'manifest.json', manifest)
    print(json.dumps(dict(new_runs=len(groups), H1_stations=41, reused_runs=2,
                         top_PML_min_clearance_m=min(s['top_PML_min_clearance_m'] for s in stations))))


def freeze(a):
    manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert manifest['approval_basis'] == APPROVAL and manifest['execution_authorized']
    assert len(manifest['groups']) == 42
    a.action = 'prepare'; native.prepare(a)
    p = a.out/'execution_contract.json'; c = json.loads(p.read_text('utf-8'))
    c.update(approval_basis=APPROVAL, study_manifest=manifest, limits=manifest['limits'], max_batch_wall_s=21600)
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    native.save(p, c); native.save(a.out/'preflight_verification.json', native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p), groups=42, max_batch_wall_s=21600)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare','freeze','run','verify'])
    for key in ['out','package','prior']:
        p.add_argument('--'+key, type=Path, required=key == 'out')
    a = p.parse_args(); a.out = a.out.resolve()
    if a.action == 'prepare': prepare(a)
    elif a.action == 'freeze': freeze(a)
    elif a.action == 'run':
        native.supervisor.audit = native.audit
        native.supervisor.run(a.out/'execution_contract.json')
    else:
        print(json.dumps(native.audit(a.out/'execution_contract.json', True)))
