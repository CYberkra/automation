"""Authorized material ablation at x190, AGL 5 m: base + 3 single-variable variants x H0/H1."""
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
VARIANTS = ROOT / 'configs/research/line9_material_ablation_v0_1.json'
DL = 0.025
DT = DL / (299792458. * np.sqrt(2))
APPROVAL = json.loads(VARIANTS.read_text('utf-8'))['approval_basis']


def build_database(v0_2, variant):
    db = copy.deepcopy(v0_2)
    if 'target' not in variant:
        return db
    key = variant['target']
    entry = db['materials'][key]
    kept = entry['metadata']
    entry.clear()
    entry['name'] = kept['original_id']
    entry['model'] = variant['replacement']['model']
    entry['base'] = variant['replacement']['base']
    if variant['replacement']['poles']:
        entry['poles'] = variant['replacement']['poles']
    entry['averagable'] = True
    entry['metadata'] = kept
    entry['metadata']['ablation'] = variant['rationale']
    entry['metadata']['ablation_evidence'] = variant['evidence']
    return db


def prepare(a):
    assert not a.out.exists()
    spec = json.loads(VARIANTS.read_text('utf-8'))
    assert spec['approval_basis'] == APPROVAL
    v0_2 = json.loads((a.model / 'geometries/line9_materials_v0_2.json').read_text('utf-8'))
    assert v0_2['database']['id'] == 'line9_materials_v0_2'
    a.out.mkdir(parents=True)
    groups, checks = [], []
    for vname, variant in spec['variants'].items():
        db = build_database(v0_2, variant)
        for role in ['H0', 'H1']:
            cid = f'abl_{vname}_x190_{role}'
            case = a.out / cid
            src_in = a.model / f'cases/x190_{role}/profile.in'
            lines = src_in.read_text('utf-8').splitlines(keepends=True)
            titles = [i for i, s in enumerate(lines) if s.startswith('#title:')]
            assert len(titles) == 1
            lines[titles[0]] = f'#title: {cid}; approved material ablation AGL5\n'
            assert sum(s.startswith(('#hertzian_dipole:', '#rx: ')) for s in lines) == 2
            card = case / 'cases' / cid / 'profile.in'
            card.parent.mkdir(parents=True)
            card.write_text(''.join(lines), encoding='utf-8')
            geo = case / 'geometries' / f'{role}.h5'
            geo.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(a.model / f'geometries/{role}.h5', geo)
            mat = case / 'geometries/line9_materials_v0_2.json'
            native.save(mat, db)
            for key in db['materials']:
                load_material_spec('line9_materials_v0_2', key, search_directory=geo.parent)
            with h5py.File(geo, 'r+') as h:
                h.attrs['MaterialDatabase'] = db['database']['id']
            d = next(s.split() for s in lines if s.startswith('#hertzian_dipole:'))
            r = next(s.split() for s in lines if s.startswith('#rx: '))
            if 'target' in variant:
                others = [k for k in db['materials'] if k != variant['target']]
                others_equal = all(db['materials'][k] == v0_2['materials'][k] for k in others)
                target_changed = db['materials'][variant['target']]['base'] == variant['replacement']['base']
            else:
                others_equal = db == v0_2
                target_changed = None
            groups.append(dict(id=cid, role=role, variant=vname, profile_x_m=190.0,
                               input=str(card.relative_to(a.out)).replace('\\', '/'),
                               geometry=str(geo.relative_to(a.out)).replace('\\', '/'),
                               material=str(mat.relative_to(a.out)).replace('\\', '/'),
                               tx_m=[float(d[2]), float(d[3]), 0.], rx_m=[float(r[1]), float(r[2]), 0.],
                               native_shape=[8400, 1700, 1],
                               expected_samples=int(np.ceil(2.4e-6 / DT)) + 1))
            for key in ['input', 'geometry', 'material']:
                groups[-1][key + '_sha256'] = native.sha(a.out / groups[-1][key])
            assert others_equal and (target_changed in (None, True))
            checks.append(dict(id=cid, variant=vname, role=role,
                               other_materials_byte_equal=others_equal,
                               target_material_replaced=target_changed,
                               geometry_sha256=groups[-1]['geometry_sha256']))
    m = dict(status='PREPARED_APPROVED_MATERIAL_ABLATION_AGL5', approval_basis=APPROVAL,
             execution_authorized=True, generator_sha256=native.sha(__file__),
             variant_spec_sha256=native.sha(VARIANTS), dt_s=float(DT), groups=groups,
             preparation_checks=checks, no_retry=True,
             SFCW=dict(method='direct', frequency_start_Hz=20e6, frequency_stop_Hz=170e6,
                       frequency_step_Hz=300e3, tones=501, windows=['hann', 'blackman'],
                       zero_pad_factor=8, tail_taper_fraction=0, time_shift_s=0,
                       field_reference_length_m=0.025),
             limits='Single-station ablation probes, not site calibration. Constant materials replace dispersive ones as probes only. H0/H1 contrast is connected bottom sandstone. 8 m historical results are not the reference for this batch (different air path).')
    native.save(a.out / 'manifest.json', m)
    print(json.dumps(dict(status=m['status'], groups=len(groups), samples=groups[0]['expected_samples'])))


def freeze(a):
    m = json.loads((a.package / 'manifest.json').read_text('utf-8'))
    assert m['approval_basis'] == APPROVAL and m['execution_authorized']
    assert len(m['groups']) == 8 and all(g['expected_samples'] == 40703 for g in m['groups'])
    a.action = 'prepare'
    native.prepare(a)
    p = a.out / 'execution_contract.json'
    c = json.loads(p.read_text('utf-8'))
    c['approval_basis'] = APPROVAL
    c['study_manifest'] = m
    c['limits'] = m['limits']
    c['max_batch_wall_s'] = 10800
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    native.save(p, c)
    native.save(a.out / 'preflight_verification.json', native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p), groups=8, version=c['expected_version'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--model', type=Path, default=ROOT / 'model_current')
    p.add_argument('--package', type=Path)
    a = p.parse_args()
    a.out = a.out.resolve()
    a.model = a.model.resolve()
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
