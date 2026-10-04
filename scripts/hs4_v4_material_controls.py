"""Separate matched-carrier bulk spectrum and conductivity from edge averaging."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import gprMax
import h5py
import numpy as np
from scipy.constants import epsilon_0
import hs4_station_grid_controls as supervisor
from hs4_v4_factor_controls import ROOT, CENTRE, LOCAL, PROFILE, expected_material
from hs_capsule_identity import sha256

FACTORS = [('matched95', role) for role in ('rough', 'halfspace')]
FACTORS += [('low_sigma', role) for role in ('rough', 'halfspace')]


def parameters():
    omega = 2*np.pi*95e6
    ratio = omega*6.4567e-9
    er = 18.017+7.878/(1+ratio**2)
    sigma = .003+omega*epsilon_0*7.878*ratio/(1+ratio**2)
    return float(er), float(sigma)


def lines_for(factor, role):
    er, sigma = parameters()
    if factor == 'low_sigma':
        sigma = .003
    elif factor != 'matched95':
        raise ValueError('unknown material factor')
    path = CENTRE/f'centre_{role}'/'profile.in'
    lines = []
    for line in path.read_text('utf-8').splitlines():
        if line.startswith('#add_dispersion_debye:'):
            continue
        if line.startswith('#material:') and line.endswith(' cover'):
            line = f'#material: {er:.17g} {sigma:.17g} 1 0 cover'
        elif line.startswith('#box:') and line.endswith(' cover'):
            # Preserve the non-averaged Debye reference's component assignment.
            line += ' n'
        lines.append(line)
    return path, lines


def raster(path):
    commands = supervisor.commands(path)
    # The independent bulk-cell raster ignores the explicit averaging token.
    commands['#box:'] = [words[:-1] if words[-1] == 'n' else words for words in commands['#box:']]
    return supervisor.raster(commands)


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    if [(g['factor'], g['role']) for g in c['groups']] != FACTORS or sha256(PROFILE) != c['profile_sha256']:
        raise ValueError('frozen factor/profile identity differs')
    er, sigma = parameters()
    omega = 2*np.pi*95e6
    original = 18.017+7.878/(1+1j*omega*6.4567e-9)+.003/(1j*omega*epsilon_0)
    constant = er+sigma/(1j*omega*epsilon_0)
    if abs(original-constant) > 1e-13:
        raise ValueError('complex permittivity at declared carrier does not match')
    rows = []
    for g in c['groups']:
        p = Path(g['input']); old, lines = lines_for(g['factor'], g['role'])
        reference = old.with_suffix('.h5')
        if (sha256(p) != g['input_sha256'] or sha256(old) != g['reference_input_sha256']
                or p.read_text('utf-8').splitlines() != lines or sha256(reference) != g['reference_raw_sha256']):
            raise ValueError('input/reference or undeclared factor differs')
        _, shape, parsed = raster(p)
        expected = expected_material('baseline', g['role'])
        if not np.array_equal(parsed, expected) or not np.array_equal(shape, [1440, 1, 1320]):
            raise ValueError('bulk geometry differs')
        row = {'id': g['id'], 'physical_factor_check': 'PASS', 'input_sha256': sha256(p)}
        if completed:
            raw = p.with_suffix('.h5'); geometry = p.parent/'hs4t2d_geom.vtkhdf'
            with h5py.File(raw) as h:
                values = h['rxs/rx1/Ey'][:]
                if (str(h.attrs['gprMax']) != '4.0.0' or values.dtype != np.float64 or not np.isfinite(values).all()
                        or not np.array_equal(h.attrs['nx_ny_nz'], [1440, 1, 1320])
                        or not np.array_equal(h.attrs['dx_dy_dz'], [.025, .05, .025])
                        or float(h.attrs['dt']) != c['base_dt_s'] or len(values) != c['iterations']):
                    raise ValueError('runtime/grid/native precision/clock differs')
                for key, position in (('srcs/src1', [17.6, .025, 27]), ('rxs/rx1', [18.9, .025, 27])):
                    if not np.array_equal(h[key].attrs['GridPosition'], np.rint(np.array(position)/[.025, .05, .025]).astype(int)):
                        raise ValueError('actual acquisition position differs')
            with h5py.File(geometry) as h:
                if not np.array_equal(h['VTKHDF/CellData/Material'][:], expected):
                    raise ValueError('actual material-cell map differs')
            row.update({'raw_sha256': sha256(raw), 'geometry_sha256': sha256(geometry),
                        'dtype': str(values.dtype), 'dt_s': c['base_dt_s'], 'samples': len(values)})
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path), 'code_sha256': sha256(__file__),
            'carrier_complex_permittivity_matching_error': abs(original-constant),
            'groups': rows, 'scope': 'Matched-carrier numerical material controls, not field material truth.'}


def freeze(out):
    if out.exists() or str(gprMax.__version__) != '4.0.0':
        raise ValueError('new capsule and V4.0.0 required')
    old = json.loads((LOCAL/'execution_contract.json').read_text('utf-8'))
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    if any(sha256(package/name) != digest for name, digest in old['source_identities'].items()):
        raise ValueError('audited runtime changed')
    out.mkdir(parents=True); groups = []
    for factor, role in FACTORS:
        name = f'{factor}_centre_{role}'; folder = out/name; folder.mkdir()
        original, lines = lines_for(factor, role)
        p = folder/'profile.in'; p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        groups.append({'id': name, 'factor': factor, 'role': role, 'input': str(p.resolve()),
                       'input_sha256': sha256(p), 'reference_input_sha256': sha256(original),
                       'reference_raw_sha256': sha256(original.with_suffix('.h5'))})
    with h5py.File(CENTRE/'centre_rough/profile.h5') as h:
        dt, iterations = float(h.attrs['dt']), int(h.attrs['Iterations'])
    er, sigma = parameters()
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 读取V4手册并继续仿真找到因素；界面平均先导仅0.18dB，需同界面处理的材料对照。',
         'python': str(Path(sys.executable).resolve()), 'source_identities': old['source_identities'],
         'code_identities': {str(ROOT/name): sha256(ROOT/name) for name in (
             'scripts/hs4_v4_material_controls.py', 'scripts/run_hs4_v4_material_controls.cmd',
             'scripts/hs4_v4_factor_controls.py', 'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py')},
         'profile_sha256': sha256(PROFILE), 'groups': groups, 'max_runs': 4,
         'base_dt_s': dt, 'iterations': iterations, 'match_frequency_hz': 95e6,
         'epsilon0_F_per_m': epsilon_0, 'constant_eps_real': er, 'constant_sigma_equivalent95_S_per_m': sigma,
         'comparison_low_sigma_S_per_m': .003, 'averaging_policy': 'Explicit cover box n for BOTH nondispersive controls; reference Debye default non-averaged.',
         'gpu_lock': old['gpu_lock'], 'min_available_RAM_GiB': 1.25, 'min_free_VRAM_GiB': 3.5,
         'max_owned_RSS_GiB': 1.75, 'min_system_available_during_run_GiB': .2,
         'max_group_wall_s': 600, 'max_batch_wall_s': 1200,
         'invariants': 'Same36m domain/dx025/PML/acquisition15m/impulse600ns/bulk-cell geometry/nativefloat64/explicit unaveraged cover.',
         'factors': 'Debye vs nondispersive matched-complex95 isolates bulk spectrum frequency dependence. Two constant models differ only sigma.',
         'scope': 'One centre station; causal control of supplied material assumptions, not geological correction or full B-scan.'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen4 material controls; eps95={er:.9g}; sigma95={sigma:.9g}; explicit box n.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('freeze', 'run', 'verify'))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--audit-output', type=Path)
    args = ap.parse_args(); path = args.out/'execution_contract.json'
    if args.action == 'freeze':
        freeze(args.out)
    elif args.action == 'run':
        if not args.execute:
            raise ValueError('--execute required')
        supervisor.audit = audit; supervisor.run(path)
    else:
        if args.audit_output is None or args.audit_output.exists():
            raise ValueError('new --audit-output required')
        supervisor.save(args.audit_output, audit(path, True))


if __name__ == '__main__':
    main()
