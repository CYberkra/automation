"""Ground-coupled GPR control: 13-station dense B-scan, antenna at the surface.

User (2026-10-04): 再跑个地面雷达的吧. Same HS4T2D 36 m model, grid, materials,
impulse source, PML and 13-station layout as the archived 15 m airborne study;
the ONLY change is antenna height: src/rx z 27 m -> 12.025 m (one grid step
above the flat surface at z=12 m). 26 single-trace groups (13 stations x
rough/halfspace), each paired with its archived 15 m raw for clock/source
identity checks. Inputs stay byte-identical to the archived style (no explicit
dispersive-averaging line; V4 default n). One fresh attempt; no retry.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np

from hs_capsule_identity import sha256
import hs4_station_grid_controls as supervisor
from hs4_large_domain_dense import (CENTRE, N_STATIONS, RX_OFFSET,
                                    reference_audit, reference_raw, station_x)

ROOT = Path(__file__).resolve().parents[1]
CHECKS = ROOT / 'artifacts/research_checks'
LOCAL = CHECKS / '2026-10-04_hs4_local_wavefield_validation_r1'
GROUND_Z = 12.025  # one 0.025 m grid step above the flat surface at z=12 m


def input_lines(k, role):
    source = CENTRE / f'centre_{role}' / 'profile.in'
    result = []
    for line in source.read_text('utf-8').splitlines():
        if line.startswith('#hertzian_dipole:'):
            line = f'#hertzian_dipole: y {station_x(k):.12g} 0.025 {GROUND_Z} impulse'
        elif line.startswith('#rx:'):
            line = f'#rx: {station_x(k) + RX_OFFSET:.12g} 0.025 {GROUND_Z} t01 Ey'
        result.append(line)
    return source, result


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    if len(c['groups']) != 2 * N_STATIONS:
        raise ValueError('ground 26-case matrix required')
    rows = []
    for g in c['groups']:
        p = Path(g['input'])
        source, expected = input_lines(g['station_index'], g['role'])
        if (sha256(p) != g['input_sha256'] or sha256(source) != g['source_input_sha256']
                or p.read_text('utf-8').splitlines() != expected
                or sha256(g['reference_raw']) != g['reference_raw_sha256']):
            raise ValueError('frozen input/reference or declared station differs')
        spacing, shape, material = supervisor.raster(supervisor.commands(p))
        _, _, original = supervisor.raster(supervisor.commands(source))
        if not np.array_equal(material, original):
            raise ValueError('geometry differs from the archived 36m model')
        if not np.array_equal(spacing, [.025, .05, .025]) or not np.array_equal(shape, [1440, 1, 1320]):
            raise ValueError('native grid changed')
        row = {'id': g['id'], 'station_index': g['station_index'], 'role': g['role'],
               'input_sha256': sha256(p), 'physical_factor_check': 'PASS',
               'full_material_array_sha256': __import__('hashlib').sha256(material.tobytes()).hexdigest(),
               'geometry_identical_to_36m_archive': True}
        del original, material
        if completed:
            raw = p.with_suffix('.h5'); geom = p.parent / 'hs4t2d_geom.vtkhdf'
            with h5py.File(raw) as h, h5py.File(g['reference_raw']) as old:
                y = h['rxs/rx1/Ey'][:]
                if (str(h.attrs['gprMax']) != '4.0.0' or y.dtype != np.float64
                        or not np.isfinite(y).all() or not np.array_equal(h.attrs['nx_ny_nz'], shape)
                        or not np.array_equal(h.attrs['dx_dy_dz'], spacing)
                        or h.attrs['dt'] != old.attrs['dt'] or h.attrs['Iterations'] != old.attrs['Iterations']
                        or len(y) != h.attrs['Iterations']):
                    raise ValueError('native version/grid/clock/precision differs')
                for key, pos in [('srcs/src1', g['tx_m']), ('rxs/rx1', g['rx_m'])]:
                    if not np.array_equal(h[key].attrs['GridPosition'], np.rint(np.array(pos) / spacing).astype(int)):
                        raise ValueError('acquisition position differs')
                for key in ('samples',):
                    if not np.array_equal(h['srcs/src1/excitation/' + key][:], old['srcs/src1/excitation/' + key][:]):
                        raise ValueError('native source samples changed')
                for key in ('TimeSampleOffset', 'SampleInterval', 'SpatialScale', 'Quantity', 'Units'):
                    if h['srcs/src1/excitation'].attrs.get(key) != old['srcs/src1/excitation'].attrs.get(key):
                        raise ValueError('native source normalization changed')
            spacing2, shape2, material = supervisor.raster(supervisor.commands(p))
            with h5py.File(geom) as h:
                if not np.array_equal(h['VTKHDF/CellData/Material'][:], material):
                    raise ValueError('actual whole-domain material map differs')
            row.update(raw_sha256=sha256(raw), geometry_sha256=sha256(geom), dtype=str(y.dtype),
                       dt_s=float(c['base_dt_s']), iterations=len(y))
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Antenna-height factor only (27 m -> 12.025 m); archived 36 m geometry/materials/grid unchanged; 2D line-source model.'}


def freeze(out):
    if out.exists():
        raise ValueError('new capsule required')
    old = json.loads((LOCAL / 'execution_contract.json').read_text('utf-8'))
    import gprMax
    if str(gprMax.__version__) != '4.0.0':
        raise ValueError('user-specified V4.0.0 required')
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    expected_python = {n: h for n, h in old['source_identities'].items() if n.endswith('.py')}
    actual_python = {p.relative_to(package).as_posix(): sha256(p) for p in package.rglob('*.py')}
    if actual_python != expected_python:
        raise ValueError('V4 Python source differs from audited reference; stop for review')
    actual_binary = {p.relative_to(package).as_posix(): sha256(p) for p in package.rglob('*.pyd')}
    if not actual_binary:
        raise ValueError('compiled Windows solver extensions absent')
    out.mkdir(parents=True)
    groups = []
    for k in range(N_STATIONS):
        for role in ('rough', 'halfspace'):
            source, lines = input_lines(k, role)
            folder = out / f'g{k:02d}_{role}'
            folder.mkdir()
            p = folder / 'profile.in'
            p.write_text('\n'.join(lines) + '\n', encoding='utf-8')
            raw = reference_raw(k, role)
            reference_audit(k, role)
            groups.append({'id': folder.name, 'station_index': k, 'role': role,
                           'src_x_m': station_x(k),
                           'input': str(p.resolve()), 'input_sha256': sha256(p),
                           'source_input_sha256': sha256(source),
                           'reference_raw': str(raw), 'reference_raw_sha256': sha256(raw),
                           'tx_m': [station_x(k), .025, GROUND_Z],
                           'rx_m': [station_x(k) + RX_OFFSET, .025, GROUND_Z]})
    with h5py.File(reference_raw(6, 'rough')) as h:
        dt = float(h.attrs['dt'])
    c = {'status': 'FROZEN_APPROVED',
         'approval_basis': 'User 2026-10-04: 再跑个地面雷达的吧 (ground-coupled control).',
         'python': str(Path(sys.executable).resolve()),
         'source_identities': {**actual_python, **actual_binary},
         'runtime_build_policy': 'Python source exactly equals audited V4.0.0; native local extensions freshly hashed; cross-build numerical equality requires result comparison.',
         'groups': groups, 'max_runs': len(groups), 'base_dt_s': dt, 'gpu_lock': old['gpu_lock'],
         'code_identities': {str(ROOT / name): sha256(ROOT / name) for name in (
             'scripts/hs4_ground_gpr_dense.py', 'scripts/run_hs4_ground_gpr_dense.cmd',
             'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py')},
         'min_available_RAM_GiB': 1.8, 'min_free_VRAM_GiB': 4.4,
         'max_owned_RSS_GiB': 2.35, 'min_system_available_during_run_GiB': .2,
         'max_group_wall_s': 900, 'max_batch_wall_s': 900 * len(groups), 'no_retry': True,
         'invariants': 'Archived 36m domain/geometry/materials/grid/PML/600ns impulse/native float64 unchanged; inputs byte-identical to archived style except src/rx lines.',
         'factor': 'Antenna height only: z 27 -> 12.025 m (one grid step above flat surface z=12 m).',
         'scope': 'Ground-vs-airborne factor control; 13-station segment; 2D line source; not a real ground-antenna coupling model; no field claims.',
         'processing': 'Same actual-source normalized 501 20-170MHz, 200ns tail, Hann unitmean/8x; complex H1-H0 first.',
         'capacity_basis': '1.9M cells/trace (1/3 of the 108m run); measured 108m traces 25.4-35.6s and <=1.46GiB RSS on this machine.'}
    supervisor.save(out / 'execution_contract.json', c)
    supervisor.save(out / 'preflight_verification.json', audit(out / 'execution_contract.json'))
    print('Frozen ground-GPR dense B-scan: 26 traces (13 stations x 2 roles), z=12.025 m.', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('freeze', 'run', 'verify'))
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--audit-output', type=Path)
    a = ap.parse_args()
    cp = a.out / 'execution_contract.json'
    if a.action == 'freeze':
        freeze(a.out)
    elif a.action == 'run':
        if not a.execute:
            raise ValueError('--execute required')
        supervisor.audit = audit
        supervisor.run(cp)
    else:
        if a.audit_output is None or a.audit_output.exists():
            raise ValueError('new audit output required')
        supervisor.save(a.audit_output, audit(cp, True))


if __name__ == '__main__':
    main()
