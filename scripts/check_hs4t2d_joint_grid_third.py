"""Independent factor/raw/grid checks for the third (1/60m) centre resolution."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]


def check(contract_path, completed=False):
    c = json.loads(contract_path.read_text('utf-8'))
    base = Path(c['baseline_directory'])
    if sha256(base / 'execution_contract.json') != c['baseline_contract_sha256']:
        raise ValueError('baseline contract identity differs')
    records = {}
    required = [61] if c['stage'] == 'centre' else list(range(1, 52, 10)) + list(range(71, 122, 10))
    for role in (False, True):
        actual = [i for g in c['groups'] if g['halfspace'] == role for i in g['old_station_indices']]
        if actual != required:
            raise ValueError('incomplete, duplicate or out-of-order station coverage')
    fixed = ('#domain:', '#time_window:', '#omp_threads:', '#title:', '#pml_formulation:',
             '#waveform:', '#material:', '#add_dispersion_debye:', '#box:')
    for g in c['groups']:
        folder = contract_path.parent / g['id']
        source = ROOT / g['source_input']
        path = folder / 'profile.in'
        if sha256(path) != g['sha256'] or sha256(source) != g['source_sha256']:
            raise ValueError('input identity differs')
        old, new = source.read_text('utf-8').splitlines(), path.read_text('utf-8').splitlines()
        if [s for s in old if s.startswith(fixed)] != [s for s in new if s.startswith(fixed)]:
            raise ValueError('geometry/material/source/time factor changed')
        if '#dx_dy_dz: 0.016666666666666666 0.05 0.016666666666666666' not in new or '#pml_cells: 120 0 60 120 0 60' not in new:
            raise ValueError('fine grid or physical boundary widths differ')
        spacing = np.array(g['spacing_m'])
        if not np.array_equal(spacing, [1/60, .05, 1/60]):
            raise ValueError('undeclared spacing')
        nx, nz = 2160, 1980
        expected = np.full((nz, 1, nx), 2, dtype=np.uint32)
        profile = np.genfromtxt(ROOT / 'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv', delimiter=',', names=True)
        x = (np.arange(nx) + .5)/60 - 12
        bottom = np.zeros(nx) if g['halfspace'] else profile['grid_z_m'][np.clip(np.floor(x/.25).astype(int), 0, 47)]
        expected[:720] = 3
        for j, z in enumerate(bottom):
            expected[round(z*60):720, 0, j] = 4
        with h5py.File(base / g['baseline_group'] / 'hs4t2d_geom1.vtkhdf') as h:
            repeated = np.repeat(np.repeat(h['VTKHDF/CellData/Material'][:], 3, axis=0), 3, axis=2)
        if not np.array_equal(repeated, expected):
            raise ValueError('physical staircase differs from independent baseline/CSV')
        parsed = np.full_like(expected, 2)
        for s in new:
            if s.startswith('#box:'):
                w = s.split()
                a = np.array(list(map(float, w[1:7]))).reshape(2, 3) / spacing
                if not np.allclose(a, np.rint(a), rtol=0, atol=1e-10):
                    raise ValueError('geometry not cell aligned')
                lo, hi = np.rint(a).astype(int)
                parsed[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]] = {'rock': 3, 'cover': 4}[w[-1]]
        if not np.array_equal(parsed, expected):
            raise ValueError('parsed input geometry differs')
        idx = g['old_station_indices']
        tx = next(s.split() for s in new if s.startswith('#hertzian_dipole:'))
        rx = next(s.split() for s in new if s.startswith('#rx:'))
        expected_x = 14.6 + .05*(idx[0]-1)
        if not np.allclose(list(map(float, tx[2:5])), [expected_x, .025, 27], rtol=0, atol=1e-12) or not np.allclose(list(map(float, rx[1:4])), [expected_x+1.3, .025, 27], rtol=0, atol=1e-12):
            raise ValueError('source/receiver position differs')
        if len(idx) != g['traces'] or len(idx) > 1 and np.diff(idx).tolist() != [10]*(len(idx)-1):
            raise ValueError('invalid station sequence')
        steps = [s for s in new if s.startswith(('#src_steps:', '#rx_steps:'))]
        if steps != (['#src_steps: 0.5 0 0', '#rx_steps: 0.5 0 0'] if len(idx)>1 else []):
            raise ValueError('trace movement differs')
        exports = []
        if completed:
            rows = json.loads((folder / 'audit.json').read_text('utf-8'))
            if len(rows) != g['traces']:
                raise ValueError('incomplete outputs')
            for k, row in enumerate(rows, 1):
                raw = folder / row['file']
                if sha256(raw) != row['sha256'] or row['old_station_index'] != idx[k-1]:
                    raise ValueError('raw identity or index differs')
                with h5py.File(raw) as h:
                    v = h['rxs/rx1/Ey'][:]
                    if str(h.attrs['gprMax']) != '4.0.0' or v.dtype != np.float64 or not np.isfinite(v).all() or len(v) != int(h.attrs['Iterations']):
                        raise ValueError('raw validity differs')
                    if not np.array_equal(h.attrs['dx_dy_dz'], spacing) or not np.isclose(h.attrs['dt'], (1/60)/(299792458*np.sqrt(2)), rtol=1e-12, atol=0):
                        raise ValueError('fine grid /2D timestep differs')
                    sx = 14.6 + .05*(idx[k-1]-1)
                    for name, px in [('srcs/src1', sx), ('rxs/rx1', sx+1.3)]:
                        if not np.array_equal(h[name].attrs['GridPosition'], [round(px*60), 0, 1620]):
                            raise ValueError('actual station differs')
                geom = folder / (f'hs4t2d_geom{k}.vtkhdf' if len(idx)>1 else 'hs4t2d_geom.vtkhdf')
                with h5py.File(geom) as h:
                    if not np.array_equal(h['VTKHDF/CellData/Material'][:], expected):
                        raise ValueError('actual material map differs')
                exports.append({'file': geom.name, 'sha256': sha256(geom)})
        records[g['id']] = {'factors_preserved': 'PASS', 'fine_equals_repeated_coarse_voxels': True,
                            'material_array_sha256': hashlib.sha256(expected.tobytes()).hexdigest(),
                            'actual_material_grids_verified': len(exports), 'geometry_exports': exports}
    if completed:
        events = [json.loads(s) for s in (contract_path.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
        if events[0]['contract_sha256'] != sha256(contract_path) or events[-1].get('status') != 'COMPLETED' or events[-1].get('traces') != c['max_runs']:
            raise ValueError('incomplete or changed execution')
    return {'status': 'PASS', 'contract_sha256': sha256(contract_path), 'code_sha256': sha256(__file__),
            'groups': records, 'completed_outputs_checked': completed,
            'scope': 'input invariants, voxel geometry, version/rawdtype/2D dt and actual positions; not full spatial convergence'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract', type=Path, required=True)
    p.add_argument('--completed', action='store_true')
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise ValueError('new report required')
    report = check(args.contract, args.completed)
    args.out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS', len(report['groups']), 'groups')


if __name__ == '__main__':
    main()
