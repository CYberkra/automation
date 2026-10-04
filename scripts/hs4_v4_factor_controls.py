"""Frozen 15 m controls separating Debye-edge averaging, dt and geometry.

Same native spatial grid, domain, boundaries, material and acquisition.
Reuses the owned-process/resource supervisor without historical attempt reuse.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import gprMax
import h5py
import numpy as np
import hs4_station_grid_controls as supervisor
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
CENTRE = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'
REMAINING = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_remaining'
LOCAL = ROOT/'artifacts/research_checks/2026-10-04_hs4_local_wavefield_validation_r1'
PROFILE = ROOT/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv'
STATIONS = {'left': 14.6, 'centre': 17.6, 'right': 20.6}
FACTORS = [('averaging_y', station, role) for station in STATIONS for role in ('rough', 'halfspace')]
FACTORS += [('dt_half', 'centre', role) for role in ('rough', 'halfspace')]
FACTORS += [('flat', station, 'rough') for station in STATIONS]


def reference_raw(station, role):
    if station == 'centre':
        return CENTRE/f'centre_{role}'/'profile.h5'
    return REMAINING/f'{station}_{role}'/('profile1.h5' if station == 'left' else 'profile6.h5')


def lines_for(factor, station, role):
    old = CENTRE/f'centre_{role}'/'profile.in'
    lines = []
    for line in old.read_text('utf-8').splitlines():
        if line.startswith('#hertzian_dipole:'):
            line = f'#hertzian_dipole: y {STATIONS[station]:.12g} 0.025 27 impulse'
        elif line.startswith('#rx:'):
            line = f'#rx: {STATIONS[station]+1.3:.12g} 0.025 27 t01 Ey'
        elif factor == 'flat' and line.startswith('#box:') and line.endswith(' cover'):
            continue
        lines.append(line)
    if factor == 'averaging_y':
        lines.append('#dispersive_averaging: y')
    elif factor == 'dt_half':
        lines.append('#time_step_stability_factor: 0.5')
    elif factor == 'flat':
        lines.append('#box: 0 0 9.05 36 0.05 12 cover')
    else:
        raise ValueError('unknown factor')
    return old, lines


def expected_material(factor, role):
    profile = np.genfromtxt(PROFILE, delimiter=',', names=True)
    cells = np.full((1320, 1, 1440), 2, dtype=np.uint32)
    cells[:480] = 3
    x = (np.arange(1440)+.5)*.025-12
    bottom = profile['grid_z_m'][np.clip(np.floor(x/.25).astype(int), 0, 47)]
    if role == 'halfspace':
        bottom[:] = 0
    elif factor == 'flat':
        bottom[:] = 9.05
    for ix, z in enumerate(bottom):
        cells[round(z/.025):480, 0, ix] = 4
    return cells


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    if [(g['factor'], g['station'], g['role']) for g in c['groups']] != FACTORS:
        raise ValueError('frozen 11-case factor matrix differs')
    if sha256(PROFILE) != c['profile_sha256']:
        raise ValueError('independent geometric profile changed')
    rows = []
    for g in c['groups']:
        p = Path(g['input'])
        old, expected = lines_for(g['factor'], g['station'], g['role'])
        if sha256(p) != g['input_sha256'] or sha256(old) != g['reference_input_sha256']:
            raise ValueError('frozen input changed')
        if p.read_text('utf-8').splitlines() != expected:
            raise ValueError('undeclared physical or numerical factor')
        _, shape, parsed = supervisor.raster(supervisor.commands(p))
        material = expected_material(g['factor'], g['role'])
        if not np.array_equal(parsed, material) or not np.array_equal(shape, [1440, 1, 1320]):
            raise ValueError('parsed geometry differs from independent CSV/cell map')
        raw_reference = Path(g['reference_raw'])
        if sha256(raw_reference) != g['reference_raw_sha256']:
            raise ValueError('reference receiver changed')
        row = {'id': g['id'], 'input_sha256': sha256(p), 'physical_factor_check': 'PASS'}
        if completed:
            raw, geometry = p.with_suffix('.h5'), p.parent/'hs4t2d_geom.vtkhdf'
            with h5py.File(raw) as h:
                samples = h['rxs/rx1/Ey'][:]
                dt = float(h.attrs['dt'])
                if (str(h.attrs['gprMax']) != '4.0.0' or samples.dtype != np.float64
                        or not np.isfinite(samples).all() or len(samples) != h.attrs['Iterations']
                        or not np.array_equal(h.attrs['nx_ny_nz'], [1440, 1, 1320])
                        or not np.array_equal(h.attrs['dx_dy_dz'], [.025, .05, .025])
                        or not np.isclose(dt, c['base_dt_s']*g['dt_factor'], rtol=1e-13, atol=0)):
                    raise ValueError('runtime/grid/clock/native precision differs')
                expected_n = int(np.ceil(600e-9/dt))+1
                if len(samples) != expected_n:
                    raise ValueError('native time window differs')
                for key, position in (('srcs/src1', g['tx_m']), ('rxs/rx1', g['rx_m'])):
                    actual = h[key].attrs['GridPosition']
                    if not np.array_equal(actual, np.rint(np.array(position)/[.025, .05, .025]).astype(int)):
                        raise ValueError('actual source/receiver differs')
            with h5py.File(geometry) as h:
                if not np.array_equal(h['VTKHDF/CellData/Material'][:], material):
                    raise ValueError('actual bulk material map differs')
            log = (p.parent/'stdout.log').read_text('utf-8')
            if g['factor'] == 'averaging_y' and 'Dispersive material averaging: enabled' not in log:
                raise ValueError('requested interface averaging absent from solver log')
            row.update({'raw_sha256': sha256(raw), 'geometry_sha256': sha256(geometry),
                        'dt_s': dt, 'iterations': len(samples), 'dtype': str(samples.dtype),
                        'bulk_material_matches_independent_csv': True})
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Controlled discrete numerical/geometry diagnostic; no continuum or 3D certification.'}


def freeze(out):
    if out.exists() or str(gprMax.__version__) != '4.0.0':
        raise ValueError('new capsule and V4.0.0 required')
    old = json.loads((LOCAL/'execution_contract.json').read_text('utf-8'))
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    if any(sha256(package/name) != digest for name, digest in old['source_identities'].items()):
        raise ValueError('audited local runtime changed')
    out.mkdir(parents=True)
    groups = []
    for factor, station, role in FACTORS:
        name = f'{factor}_{station}_{role}'
        folder = out/name; folder.mkdir()
        original, lines = lines_for(factor, station, role)
        p = folder/'profile.in'; p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        raw = reference_raw(station, role)
        groups.append({'id': name, 'factor': factor, 'station': station, 'role': role,
                       'input': str(p.resolve()), 'input_sha256': sha256(p),
                       'reference_input_sha256': sha256(original), 'reference_raw': str(raw),
                       'reference_raw_sha256': sha256(raw),
                       'dt_factor': .5 if factor == 'dt_half' else 1.,
                       'tx_m': [STATIONS[station], .025, 27], 'rx_m': [STATIONS[station]+1.3, .025, 27]})
    with h5py.File(CENTRE/'centre_rough/profile.h5') as h:
        dt = float(h.attrs['dt'])
    docs = Path('D:/gprMax-v.4.0.0/docs/source')
    manual = {str(docs/name): sha256(docs/name) for name in ('gprmodelling.rst', 'input_hash_cmds.rst', 'features.rst')}
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 读取gprMax手册V4，再进行仿真，找到因素。',
         'python': str(Path(sys.executable).resolve()), 'source_identities': old['source_identities'],
         'max_runs': len(groups), 'base_dt_s': dt, 'groups': groups,
         'profile_sha256': sha256(PROFILE), 'manual_reading_identities': manual,
         'code_identities': {str(ROOT/name): sha256(ROOT/name) for name in (
             'scripts/hs4_v4_factor_controls.py', 'scripts/run_hs4_v4_factor_controls.cmd',
             'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py')},
         'gpu_lock': old['gpu_lock'], 'min_available_RAM_GiB': 1.25, 'min_free_VRAM_GiB': 3.5,
         'max_owned_RSS_GiB': 1.75, 'min_system_available_during_run_GiB': .2,
         'max_group_wall_s': 600, 'max_batch_wall_s': 2400,
         'invariants': '36x.05x33m; dx/z.025; same15m positions; PML2m sides/1mZ; Debye bulk; impulse;600ns;nativefloat64.',
         'factors': 'Debye interface averaging only at3stations; dt only atcentre; planar geometry only at3stations.',
         'processing': 'Actual source-normalized complex50120-170MHz;200ns tail; Hann unitmean/8x; complex differences before envelope.',
         'scope': 'No material-change/source-change/height-change sweep; reuse archived controls; no3D/field/G4 claims.'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen {len(groups)} controls; spatial grid and bulk Debye unchanged.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('freeze', 'run', 'verify'))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--audit-output', type=Path)
    args = ap.parse_args()
    path = args.out/'execution_contract.json'
    if args.action == 'freeze':
        freeze(args.out)
    elif args.action == 'run':
        if not args.execute:
            raise ValueError('--execute required')
        supervisor.audit = audit
        supervisor.run(path)
    else:
        if args.audit_output is None or args.audit_output.exists():
            raise ValueError('new --audit-output required')
        supervisor.save(args.audit_output, audit(path, True))


if __name__ == '__main__':
    main()
