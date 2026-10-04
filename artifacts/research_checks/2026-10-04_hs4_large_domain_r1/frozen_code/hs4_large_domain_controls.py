"""Threefold X-domain control, with the original local geometry translated intact."""
import argparse
import json
from pathlib import Path
import sys
import h5py
import numpy as np
from hs_capsule_identity import sha256
import hs4_station_grid_controls as supervisor
from hs4_v4_factor_controls import ROOT, LOCAL, reference_raw, STATIONS

ORDER = [(s, r) for s in ('centre', 'left', 'right') for r in ('rough', 'halfspace')]
SHIFT, WIDTH, DX = 36., 108., .025


def input_lines(station, role):
    source = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'/f'centre_{role}'/'profile.in'
    result = []
    for line in source.read_text('utf-8').splitlines():
        w = line.split()
        if line.startswith('#domain:'):
            line = '#domain: 108 0.05 33'
        elif line.startswith('#hertzian_dipole:'):
            line = f'#hertzian_dipole: y {STATIONS[station]+SHIFT:.12g} 0.025 27 impulse'
        elif line.startswith('#rx:'):
            line = f'#rx: {STATIONS[station]+SHIFT+1.3:.12g} 0.025 27 t01 Ey'
        elif line.startswith('#box:'):
            lo, hi = float(w[1]), float(w[4])
            w[1] = f'{0 if lo == 0 else lo+SHIFT:.12g}'
            w[4] = f'{WIDTH if hi == 36 else hi+SHIFT:.12g}'
            line = ' '.join(w)
        elif line.startswith('#geometry_view:'):
            line = '#geometry_view: 0 0 0 108 0.05 33 0.025 0.05 0.025 hs4t2d_geom n'
        result.append(line)
    result.append('#dispersive_averaging: n')
    return source, result


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    if [(g['station'], g['role']) for g in c['groups']] != ORDER:
        raise ValueError('frozen six-case matrix differs')
    rows = []
    for g in c['groups']:
        p = Path(g['input']); source, expected = input_lines(g['station'], g['role'])
        if (sha256(p) != g['input_sha256'] or sha256(source) != g['source_input_sha256']
                or p.read_text('utf-8').splitlines() != expected
                or sha256(g['reference_raw']) != g['reference_raw_sha256']):
            raise ValueError('frozen input/reference or declared factor differs')
        spacing, shape, material = supervisor.raster(supervisor.commands(p))
        _, _, original = supervisor.raster(supervisor.commands(source))
        # Independent full-array construction: constant continuation of old edge columns.
        padded = np.pad(original, ((0, 0), (0, 0), (1440, 1440)), mode='edge')
        if not np.array_equal(material, padded) or not np.array_equal(shape, [4320, 1, 1320]):
            raise ValueError('geometry not an intact translation plus edge continuation')
        if not np.array_equal(spacing, [.025, .05, .025]):
            raise ValueError('native spacing changed')
        row = {'id': g['id'], 'input_sha256': sha256(p), 'physical_factor_check': 'PASS',
               'full_material_array_sha256': __import__('hashlib').sha256(padded.tobytes()).hexdigest(),
               'old_roi_exactly_preserved': True, 'new_cells': shape.tolist()}
        del original, padded
        if completed:
            raw = p.with_suffix('.h5'); geom = p.parent/'hs4t2d_geom.vtkhdf'
            with h5py.File(raw) as h, h5py.File(g['reference_raw']) as old:
                y = h['rxs/rx1/Ey'][:]
                if (str(h.attrs['gprMax']) != '4.0.0' or y.dtype != np.float64
                        or not np.isfinite(y).all() or not np.array_equal(h.attrs['nx_ny_nz'], shape)
                        or not np.array_equal(h.attrs['dx_dy_dz'], spacing)
                        or h.attrs['dt'] != old.attrs['dt'] or h.attrs['Iterations'] != old.attrs['Iterations']
                        or len(y) != h.attrs['Iterations']):
                    raise ValueError('native version/grid/clock/precision differs')
                for key, pos in [('srcs/src1', g['tx_m']), ('rxs/rx1', g['rx_m'])]:
                    if not np.array_equal(h[key].attrs['GridPosition'], np.rint(np.array(pos)/spacing).astype(int)):
                        raise ValueError('acquisition translation differs')
                for key in ('samples',):
                    if not np.array_equal(h['srcs/src1/excitation/'+key][:], old['srcs/src1/excitation/'+key][:]):
                        raise ValueError('native source samples changed')
                for key in ('TimeSampleOffset', 'SampleInterval', 'SpatialScale', 'Quantity', 'Units'):
                    if h['srcs/src1/excitation'].attrs.get(key) != old['srcs/src1/excitation'].attrs.get(key):
                        raise ValueError('native source normalization changed')
            with h5py.File(geom) as h:
                if not np.array_equal(h['VTKHDF/CellData/Material'][:], material):
                    raise ValueError('actual whole-domain material map differs')
            row.update(raw_sha256=sha256(raw), geometry_sha256=sha256(geom), dtype=str(y.dtype),
                       dt_s=float(c['base_dt_s']), iterations=len(y))
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'X-domain sensitivity only; unchanged top/bottom, grid and 2D line-source model.'}


def freeze(out):
    if out.exists():
        raise ValueError('new capsule required')
    old = json.loads((LOCAL/'execution_contract.json').read_text('utf-8'))
    out.mkdir(parents=True); groups = []
    for station, role in ORDER:
        source, lines = input_lines(station, role)
        folder = out/f'{station}_{role}'; folder.mkdir(); p = folder/'profile.in'
        p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        raw = reference_raw(station, role)
        groups.append({'id': folder.name, 'station': station, 'role': role, 'input': str(p.resolve()),
                       'input_sha256': sha256(p), 'source_input_sha256': sha256(source),
                       'reference_raw': str(raw), 'reference_raw_sha256': sha256(raw),
                       'tx_m': [STATIONS[station]+SHIFT, .025, 27],
                       'rx_m': [STATIONS[station]+SHIFT+1.3, .025, 27]})
    with h5py.File(reference_raw('centre', 'rough')) as h:
        dt = float(h.attrs['dt'])
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 能不能跑一个很大的域看看？',
         'python': str(Path(sys.executable).resolve()), 'source_identities': old['source_identities'],
         'groups': groups, 'max_runs': 6, 'base_dt_s': dt, 'gpu_lock': old['gpu_lock'],
         'code_identities': {str(ROOT/name): sha256(ROOT/name) for name in (
             'scripts/hs4_large_domain_controls.py', 'scripts/run_hs4_large_domain_controls.cmd',
             'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py')},
         'min_available_RAM_GiB': 1.8, 'min_free_VRAM_GiB': 4.4,
         'max_owned_RSS_GiB': 2.35, 'min_system_available_during_run_GiB': .2,
         'max_group_wall_s': 900, 'max_batch_wall_s': 5400, 'no_retry': True,
         'invariants': 'dx/z.025; y.05 invariant; Debye parameters; 600ns impulse; nativefloat64; 15m altitude;1.3m Tx/Rx; same PML2m sides/1mZ; dispersive averaging n.',
         'factor': 'X36->108m; rigid+36m translation; edge-column continuation; no local interface or height change.',
         'scope': 'Three stations only; not full B-scan/grid convergence/3D/all-boundary isolation.',
         'processing': 'Same actual-source normalized50120-170MHz,200ns tail,Hann unitmean/8x; complex H1-H0 first; fixed160-180/180-220/160-220ns.',
         'capacity_basis': '5.7024M cells,3x old; olddevice~1.5GB scales~4.5GB decimal; expected host below2GiB; live guarded; no precision reduction or termination of other apps.'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print('Frozen six 108m cases, centre first; original ROI exactly preserved.', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('freeze', 'run', 'verify'))
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--audit-output', type=Path)
    a = ap.parse_args(); cp = a.out/'execution_contract.json'
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
