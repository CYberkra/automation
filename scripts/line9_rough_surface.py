"""Authorized surface micro-relief sensitivity pair at x190, AGL 5 m; never retried."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from gprMax.material_database import load_material_spec
import line9_v401_version_controls as native

ROOT = Path(__file__).resolve().parents[1]
DL = 0.025
DT = DL / (299792458. * np.sqrt(2))
APPROVAL = ('2026-10-11 用户显性批准：地表起伏敏感性试验（建议 S5）。官方 #fractal_box:+'
            '#add_surface_roughness: 只适用于平直 box 面，无法贴合 line9 真实起伏地形，'
            '故以同族分形算法（谱合成 fBm，H=0.5，种子 20261011，RMS 目标 0.1m，截断 ±0.2m）'
            '直接微扰导入体素地表；几何导入路径仍为官方 #geometry_objects_read:。'
            'x190 站位 H0/H1 粗糙对，对照为同批消融基线（5m 平滑）。')


def prepare(a):
    assert not a.out.exists()
    v0_2 = json.loads((a.model / 'geometries/line9_materials_v0_2.json').read_text('utf-8'))
    a.out.mkdir(parents=True)
    groups, checks = [], []
    for role in ['H0', 'H1']:
        cid = f'rough_x190_{role}'
        case = a.out / cid
        src_in = a.model / f'cases/x190_{role}/profile.in'
        lines = src_in.read_text('utf-8').splitlines(keepends=True)
        titles = [i for i, s in enumerate(lines) if s.startswith('#title:')]
        geos = [i for i, s in enumerate(lines) if s.startswith('#geometry_objects_read:')]
        assert len(titles) == 1 and len(geos) == 1
        lines[titles[0]] = f'#title: {cid}; approved surface micro-relief sensitivity AGL5\n'
        parts = lines[geos[0]].split()
        assert parts[4] == f'../../geometries/{role}.h5'
        parts[4] = f'../../geometries/{role}_rough.h5'
        end = '\r\n' if lines[geos[0]].endswith('\r\n') else '\n'
        lines[geos[0]] = ' '.join(parts) + end
        card = case / 'cases' / cid / 'profile.in'
        card.parent.mkdir(parents=True)
        card.write_text(''.join(lines), encoding='utf-8')
        geo = case / 'geometries' / f'{role}_rough.h5'
        geo.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(a.rough / f'{role}_rough.h5', geo)
        mat = case / 'geometries/line9_materials_v0_2.json'
        native.save(mat, v0_2)
        for key in v0_2['materials']:
            load_material_spec('line9_materials_v0_2', key, search_directory=geo.parent)
        with h5py.File(geo, 'r+') as h:
            h.attrs['MaterialDatabase'] = v0_2['database']['id']
        with h5py.File(a.model / f'geometries/{role}.h5') as h0, h5py.File(geo) as h1:
            diff = int(np.count_nonzero(h0['data'][:] != h1['data'][:]))
        d = next(s.split() for s in lines if s.startswith('#hertzian_dipole:'))
        r = next(s.split() for s in lines if s.startswith('#rx: '))
        groups.append(dict(id=cid, role=role, variant='rough', profile_x_m=190.0,
                           input=str(card.relative_to(a.out)).replace('\\', '/'),
                           geometry=str(geo.relative_to(a.out)).replace('\\', '/'),
                           material=str(mat.relative_to(a.out)).replace('\\', '/'),
                           tx_m=[float(d[2]), float(d[3]), 0.], rx_m=[float(r[1]), float(r[2]), 0.],
                           native_shape=[8400, 1700, 1],
                           expected_samples=int(np.ceil(2.4e-6 / DT)) + 1))
        for key in ['input', 'geometry', 'material']:
            groups[-1][key + '_sha256'] = native.sha(a.out / groups[-1][key])
        checks.append(dict(id=cid, changed_surface_voxels=diff,
                           rough_stats=json.loads((a.rough / 'rough_stats.json').read_text('utf-8'))[role]))
    m = dict(status='PREPARED_APPROVED_ROUGH_SURFACE_AGL5', approval_basis=APPROVAL,
             execution_authorized=True, generator_sha256=native.sha(__file__), dt_s=float(DT),
             groups=groups, preparation_checks=checks, no_retry=True,
             reference_smooth='line9_material_ablation_20261011_r1 的 abl_base_x190_H0/H1（同模型同航高平滑对照）',
             SFCW=dict(method='direct', frequency_start_Hz=20e6, frequency_stop_Hz=170e6,
                       frequency_step_Hz=300e3, tones=501, windows=['hann', 'blackman'],
                       zero_pad_factor=8, tail_taper_fraction=0, time_shift_s=0,
                       field_reference_length_m=0.025),
             limits='Sensitivity pair only; fractal parameters are a test choice, not a site claim. Same-station smooth reference comes from the ablation batch, not re-run here.')
    native.save(a.out / 'manifest.json', m)
    print(json.dumps(dict(status=m['status'], groups=len(groups))))


def freeze(a):
    m = json.loads((a.package / 'manifest.json').read_text('utf-8'))
    assert m['approval_basis'] == APPROVAL and m['execution_authorized']
    assert len(m['groups']) == 2 and all(g['expected_samples'] == 40703 for g in m['groups'])
    a.action = 'prepare'
    native.prepare(a)
    p = a.out / 'execution_contract.json'
    c = json.loads(p.read_text('utf-8'))
    c['approval_basis'] = APPROVAL
    c['study_manifest'] = m
    c['limits'] = m['limits']
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    native.save(p, c)
    native.save(a.out / 'preflight_verification.json', native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p), groups=2, version=c['expected_version'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--model', type=Path, default=ROOT / 'model_current')
    p.add_argument('--rough', type=Path, required=True)
    p.add_argument('--package', type=Path)
    a = p.parse_args()
    a.out = a.out.resolve()
    a.model = a.model.resolve()
    a.rough = a.rough.resolve()
    if a.package:
        a.package = a.package.resolve()
    if a.action == 'prepare':
        prepare(a)
    elif a.action == 'freeze':
        freeze(a)
    elif a.action == 'run':
        native.supervisor.audit = native.audit
        native.supervisor.run(a.out / 'execution_contract.json')
    else:
        print(json.dumps(native.audit(a.out / 'execution_contract.json', True)))
