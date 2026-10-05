"""Bounded local V4 half-space point-source validation; new immutable capsules only.

Two traces, identical domain/source/receiver: (a) lossy dispersive half-space
(archived cover material filling z=0..3 m), (b) free-space control. Source and
receiver fly 8 m above the surface. FDTD results are compared against the
independent Sommerfeld reference (scripts/green_halfspace_hed_v0_1.py).
"""
import argparse
import json
from pathlib import Path
import sys

import gprMax
import h5py
import numpy as np

import hs4_station_grid_controls as supervisor
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
DL = 0.025
DOMAIN = (12., 12., 13.)           # surface z=3; source z=11; 1 m clearance to top PML
SOURCE = np.array([6., 6., 11.])   # 8 m above surface at z=3
RECEIVER = np.array([7.3, 6., 11.])
SURFACE_Z = 3.
MATERIAL = ('#material: 18.017 0.003 1 0 cover',
            '#add_dispersion_debye: 1 7.878 6.4567e-09 cover')
GROUPS = ('halfspace', 'freespace')


def recipe(kind):
    lines = ['#title: Half-space x-HED validation v0.1; 8m altitude; archived cover material',
             f'#domain: {DOMAIN[0]:g} {DOMAIN[1]:g} {DOMAIN[2]:g}',
             f'#dx_dy_dz: {DL:g} {DL:g} {DL:g}', '#time_window: 200e-9',
             '#omp_threads: 2', '#pml_cells: 40', '#pml_formulation: HORIPML',
             f'#waveform: ricker {1/DL:.12g} 100e6 pulse',
             '#hertzian_dipole: x ' + ' '.join(f'{x:.12g}' for x in SOURCE) + ' pulse',
             '#rx: ' + ' '.join(f'{x:.12g}' for x in RECEIVER) + ' r1 Ex Ey Ez']
    if kind == 'halfspace':
        lines += list(MATERIAL)
        lines.append(f'#box: 0 0 0 {DOMAIN[0]:g} {DOMAIN[1]:g} {SURFACE_Z:g} cover')
    return lines


def audit(path, completed=False):
    contract = json.loads(path.read_text('utf-8'))
    if [g['id'] for g in contract['groups']] != list(GROUPS):
        raise ValueError('declared trace set required')
    rows = []
    for group in contract['groups']:
        p = Path(group['input'])
        lines = recipe(group['id'])
        if sha256(p) != group['input_sha256'] or p.read_text('utf-8').splitlines() != lines:
            raise ValueError('input recipe identity differs')
        for position in (SOURCE, RECEIVER):
            if not np.allclose(position/DL, np.rint(position/DL), atol=1e-10, rtol=0):
                raise ValueError('off-grid source/receiver')
        row = {'id': group['id'], 'input_sha256': sha256(p), 'dl_m': DL,
               'surface_z_m': SURFACE_Z if group['id'] == 'halfspace' else None,
               'source_current_moment_amplitude_Am': 1.,
               'material': 'archived_cover_eps18.017_s0.003_debye' if group['id'] == 'halfspace' else 'air'}
        if completed:
            raw = p.with_suffix('.h5')
            with h5py.File(raw) as h:
                dt = float(h.attrs['dt'])
                shape = [round(d/DL) for d in DOMAIN]
                if (str(h.attrs['gprMax']) != '4.0.0'
                        or not np.array_equal(h.attrs['nx_ny_nz'], shape)
                        or not np.allclose(h.attrs['dx_dy_dz'], [DL]*3, rtol=0, atol=1e-14)
                        or not np.isclose(dt, DL/(299792458*np.sqrt(3)), rtol=1e-12)):
                    raise ValueError('native version/grid/time differs')
                source = h['srcs/src1']
                if not np.array_equal(source.attrs['GridPosition'], np.rint(SOURCE/DL).astype(int)):
                    raise ValueError('source grid position differs')
                excitation = source['excitation/samples'][:]
                if excitation.dtype != np.float64 or not np.isfinite(excitation).all():
                    raise ValueError('native source is not finite double')
                receiver = h['rxs/rx1']
                if not np.array_equal(receiver.attrs['GridPosition'], np.rint(RECEIVER/DL).astype(int)):
                    raise ValueError('receiver position differs')
                for component in ('Ex', 'Ey', 'Ez'):
                    values = receiver[component][:]
                    if (values.dtype != np.float64 or not np.isfinite(values).all()
                            or len(values) != int(h.attrs['Iterations'])):
                        raise ValueError('native receiver validity differs')
                row.update(raw_sha256=sha256(raw), dt_s=dt,
                           iterations=int(h.attrs['Iterations']), dtype='float64')
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Native execution identity; half-space physics checked separately against Sommerfeld reference.'}


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and installed V4.0.0 required')
    out.mkdir(parents=True)
    groups = []
    for kind in GROUPS:
        folder = out/kind
        folder.mkdir()
        p = folder/'profile.in'
        p.write_text('\n'.join(recipe(kind))+'\n', encoding='utf-8')
        groups.append({'id': kind, 'input': str(p.resolve()), 'input_sha256': sha256(p),
                       'spacing_m': DL, 'polarisation': 'x',
                       'tx_m': SOURCE.tolist(), 'receivers_m': [RECEIVER.tolist()]})
    package = Path(gprMax.__file__).parent
    source_names = ['__init__.py', 'sources.py', 'fields_outputs.py', 'model.py', 'waveforms.py',
                    'cuda_opencl/knl_fields_updates.py', 'cuda_opencl/knl_source_updates.py',
                    'cuda_opencl/knl_pml_updates_electric_HORIPML.py',
                    'cuda_opencl/knl_pml_updates_magnetic_HORIPML.py',
                    'toolboxes/SFCW/processing.py']
    source_names += [str(p.relative_to(package)).replace('\\', '/') for p in (package/'cython').glob('*.pyd')]
    code_names = ['scripts/hs4_halfspace_point_v0_1.py', 'scripts/run_halfspace_point_4090.cmd',
                  'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py']
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 跑，你去跑完 (2026-10-05); prior autonomous authorization persists.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': {n: sha256(package/n) for n in source_names},
         'code_identities': {str(ROOT/n): sha256(ROOT/n) for n in code_names},
         'groups': groups, 'max_runs': len(groups), 'no_retry': True,
         'min_available_RAM_GiB': 4., 'min_free_VRAM_GiB': 11.,
         'max_owned_RSS_GiB': 24., 'min_system_available_during_run_GiB': 2.,
         'max_group_wall_s': 2400, 'max_batch_wall_s': 5400,
         'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'stage': 'round2_half_space_sommerfeld_validation',
         'geometry': '12x12x13m; surface z=3m; cover z 0..3m; tx/rx 8m altitude; baseline 1.3m inline-x; 40-cell HORIPML all faces; tx 1m clear of top PML',
         'source': '100MHz Ricker actual-source deconvolution; amp*dl=1Am; ideal x-Hertzian, not a UAV hardware model',
         'processing': 'complex501 tones20-170MHz, tail20ns, co-polar Ex/(I*dl), no empirical phase/amplitude fit',
         'reference': 'scripts/green_halfspace_hed_v0_1.py Sommerfeld integrals (selftest PASS 2026-10-05) + FDTD-validated free-space dipole_field',
         'diagnostic_tolerances': {'complex_band_relative_L2': .05, 'maximum_phase_error_deg': 5.},
         'tolerance_scope': 'Predeclared engineering diagnostic; not field or high-precision certification.'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen {len(groups)} half-space validation traces; native GPU double', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('freeze', 'run', 'verify'))
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    contract = args.out/'execution_contract.json'
    if args.action == 'freeze':
        freeze(args.out)
    elif args.action == 'run':
        supervisor.audit = audit
        supervisor.run(contract)
    else:
        print(json.dumps(audit(contract, True), indent=2))


if __name__ == '__main__':
    main()
