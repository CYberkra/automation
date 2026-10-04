"""Archive bounded HS4T2D controls with all raw traces and one grid per group.

Duplicate geometry exports stay on the execution machine. Their identities and
the original all-trace grid audit remain recorded, rather than being replaced
by a claim that omitted files were checked on a fresh clone.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256


def read_json(path):
    return json.loads(path.read_text('utf-8'))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def verify(study):
    selection = read_json(study / 'archive_selection.json')
    manifest = read_json(study / 'manifest.json')
    for row in manifest:
        path = study / row['file']
        if path.stat().st_size != row['bytes'] or sha256(path) != row['sha256']:
            raise ValueError('archive identity differs: ' + row['file'])
    contract = read_json(study / 'execution_contract.json')
    completed = read_json(study / selection['all_trace_grid_audit'])
    if completed['status'] != 'PASS' or completed['contract_sha256'] != sha256(study / 'execution_contract.json'):
        raise ValueError('historical all-trace grid audit identity differs')
    events = [json.loads(s) for s in (study / 'execution.jsonl').read_text('utf-8').splitlines()]
    if events[0]['contract_sha256'] != completed['contract_sha256'] or events[-1].get('status') != 'COMPLETED' or events[-1].get('traces') != contract['max_runs']:
        raise ValueError('execution is incomplete or its contract differs')
    raw_count = 0
    for group in contract['groups']:
        name = group['id']
        directory = study / name
        if sha256(directory / 'profile.in') != group['sha256']:
            raise ValueError('input identity differs')
        rows = read_json(directory / 'audit.json')
        if len(rows) != group['traces']:
            raise ValueError('raw trace count differs')
        for row in rows:
            path = directory / row['file']
            if sha256(path) != row['sha256']:
                raise ValueError('raw identity differs')
            with h5py.File(path, 'r') as h:
                values = h['rxs/rx1/Ey'][:]
                if str(h.attrs['gprMax']) != '4.0.0' or values.dtype != np.float64 or not np.isfinite(values).all():
                    raise ValueError('V4/FP64 raw validity differs')
            raw_count += 1
        representative = selection['representative_geometry'][name]
        with h5py.File(study / representative, 'r') as h:
            material_hash = hashlib.sha256(h['VTKHDF/CellData/Material'][:].tobytes()).hexdigest()
        record = completed['groups'][name]
        if material_hash != record['material_array_sha256'] or record['actual_material_grids_verified'] != group['traces']:
            raise ValueError('representative material array or historical count differs')
    if raw_count != contract['max_runs']:
        raise ValueError('total raw trace count differs')
    # On the execution machine verify omitted files too; a clone only retains
    # their original identities and the historical independent material audit.
    locally_rechecked = 0
    for row in selection['omitted_duplicate_geometry']:
        path = study / row['file']
        if path.exists():
            if sha256(path) != row['sha256']:
                raise ValueError('local omitted geometry changed')
            locally_rechecked += 1
    return {'status': 'PASS', 'raw_traces_verified': raw_count,
            'representative_material_arrays_verified': len(contract['groups']),
            'omitted_geometry_locally_rechecked': locally_rechecked,
            'archived_files_verified': len(manifest),
            'manifest_sha256': sha256(study / 'manifest.json')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    study = args.study
    if args.create:
        if (study / 'manifest.json').exists() or (study / 'archive_selection.json').exists():
            raise ValueError('archive already frozen')
        contract = read_json(study / 'execution_contract.json')
        audit_name = 'completed_verification_v0_2.json' if (study / 'completed_verification_v0_2.json').exists() else 'completed_verification.json'
        audit = read_json(study / audit_name)
        omitted = []
        representatives = {}
        for group in contract['groups']:
            name = group['id']
            geometry = audit['groups'][name]['geometry_exports']
            representatives[name] = name + '/' + geometry[0]['file']
            for row in geometry[1:]:
                path = study / name / row['file']
                if sha256(path) != row['sha256']:
                    raise ValueError('geometry identity differs before archive')
                omitted.append({'file': name + '/' + row['file'], 'bytes': path.stat().st_size,
                                'sha256': row['sha256'],
                                'material_array_sha256': audit['groups'][name]['material_array_sha256']})
        write_json(study / 'archive_selection.json', {
            'status': 'COMPLETED_ARCHIVE', 'all_raw_traces_archived': contract['max_runs'],
            'all_trace_grid_audit': audit_name,
            'representative_geometry': representatives,
            'omitted_duplicate_geometry': omitted,
            'scope': 'All raw H5 and one actual grid per group are portable. Other grids remain local; their independent execution-time verification is historical evidence.'})
        omitted_names = {row['file'] for row in omitted}
        write_json(study / 'manifest.json', [
            {'file': path.relative_to(study).as_posix(), 'bytes': path.stat().st_size, 'sha256': sha256(path)}
            for path in sorted(study.rglob('*'))
            if path.is_file() and path.name != 'manifest.json' and path.relative_to(study).as_posix() not in omitted_names])
    print(json.dumps(verify(study), indent=2))


if __name__ == '__main__':
    main()
