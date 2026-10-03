"""Verify PML-only input changes and same actual material grid as baseline."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256
from prepare_hs4t2d_cause_controls import create_input, ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--completed', action='store_true')
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new verification required')
    contract = json.loads(args.contract.read_text('utf-8'))
    records = {}
    for group in contract['groups']:
        directory = args.contract.parent / group['id']
        source = ROOT / group['source_input']
        path = directory / 'profile.in'
        if sha256(source) != group['source_sha256'] or sha256(path) != group['sha256']:
            raise ValueError('input identity differs')
        baseline = create_input(source, 36, 15, .05, group['halfspace'], group['traces']).splitlines()
        actual = path.read_text('utf-8').splitlines()
        pml_line = '#pml_cells: ' + ' '.join(map(str, group['pml_cells']))
        if pml_line not in actual or [s for s in baseline if not s.startswith('#pml_cells:')] != [s for s in actual if not s.startswith('#pml_cells:')]:
            raise ValueError('a non-PML factor changed')
        role = 'halfspace' if group['halfspace'] else 'rough'
        reference = ROOT / ('artifacts/research_checks/2026-10-03_hs4t2d_cause_controls/wide36_' + role)
        with h5py.File(reference / 'hs4t2d_geom1.vtkhdf') as h:
            material = h['VTKHDF/CellData/Material'][:]
        arrays = []
        if args.completed:
            audit = json.loads((directory / 'audit.json').read_text('utf-8'))
            if len(audit) != group['traces']:
                raise ValueError('incomplete group')
            for k, row in enumerate(audit, 1):
                raw = directory / row['file']
                if sha256(raw) != row['sha256']:
                    raise ValueError('raw changed')
                with h5py.File(raw) as h:
                    values = h['rxs/rx1/Ey'][:]
                    if str(h.attrs['gprMax']) != '4.0.0' or values.dtype != np.float64 or not np.isfinite(values).all():
                        raise ValueError('raw validity differs')
                geom = directory / (f'hs4t2d_geom{k}.vtkhdf' if group['traces'] > 1 else 'hs4t2d_geom.vtkhdf')
                with h5py.File(geom) as h:
                    if not np.array_equal(h['VTKHDF/CellData/Material'][:], material):
                        raise ValueError('material changed with PML')
                arrays.append({'file': geom.name, 'sha256': sha256(geom)})
        records[group['id']] = {'PML_only_input_change': 'PASS',
                               'material_array_sha256': hashlib.sha256(material.tobytes()).hexdigest(),
                               'actual_material_grids_verified': len(arrays), 'geometry_exports': arrays}
    if args.completed:
        events = [json.loads(s) for s in (args.contract.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
        if events[-1].get('status') != 'COMPLETED' or events[-1].get('traces') != contract['max_runs']:
            raise ValueError('execution incomplete')
    args.out.write_text(json.dumps({'status': 'PASS', 'contract_sha256': sha256(args.contract),
                                   'code_sha256': sha256(__file__), 'groups': records}, indent=2) + '\n', encoding='utf-8')
    print('PASS', len(records), 'PML-only groups')


if __name__ == '__main__':
    main()
