"""Bounded local V4 free-space dipole benchmark; new immutable capsules only."""
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
SOURCE = np.array([3.4, 3.4, 4.])
TARGETS = [SOURCE+[1.3, 0, 0], SOURCE+[0, 1.3, 0],
           SOURCE+[-1.3, 0, 0], SOURCE+[0, -1.3, 0]]
CASES = [(d, p, 'forward') for d in (.1, .05) for p in ('x', 'y')]
CASES += [(.05, 'y', 'reverse_x'), (.05, 'y', 'reverse_y')]


def recipe(dl, pol, direction):
    tx = SOURCE if direction == 'forward' else TARGETS[0 if direction == 'reverse_x' else 1]
    rx = TARGETS if direction == 'forward' else [SOURCE]
    lines = ['#title: Local ideal-dipole continuum benchmark; not a UAV hardware model',
             '#domain: 8 8 8', f'#dx_dy_dz: {dl} {dl} {dl}', '#time_window: 100e-9',
             '#omp_threads: 2', f'#pml_cells: {round(1/dl)}', '#pml_formulation: HORIPML',
             f'#waveform: ricker {1/dl:.12g} 100e6 pulse',
             f'#hertzian_dipole: {pol} '+ ' '.join(f'{x:.12g}' for x in tx)+' pulse']
    for i, position in enumerate(rx):
        lines.append('#rx: '+' '.join(f'{x:.12g}' for x in position)+f' r{i+1} Ex Ey Ez')
    return lines, tx, rx


def audit(path, completed=False):
    contract = json.loads(path.read_text('utf-8'))
    if len(contract['groups']) != len(CASES):
        raise ValueError('six-case matrix required')
    rows = []
    for group, (dl, pol, direction) in zip(contract['groups'], CASES):
        p = Path(group['input'])
        lines, tx, rx = recipe(dl, pol, direction)
        if sha256(p) != group['input_sha256'] or p.read_text('utf-8').splitlines() != lines:
            raise ValueError('input recipe identity differs')
        for position in [tx, *rx]:
            if not np.allclose(position/dl, np.rint(position/dl), atol=1e-10, rtol=0):
                raise ValueError('off-grid source/receiver')
        row = {'id': group['id'], 'input_sha256': sha256(p), 'declared_air_only': True,
               'source_current_moment_amplitude_Am': 1., 'dl_m': dl, 'polarisation': pol}
        if completed:
            raw = p.with_suffix('.h5')
            with h5py.File(raw) as h:
                dt = float(h.attrs['dt'])
                if (str(h.attrs['gprMax']) != '4.0.0'
                        or not np.array_equal(h.attrs['nx_ny_nz'], [round(8/dl)]*3)
                        or not np.allclose(h.attrs['dx_dy_dz'], [dl]*3, rtol=0, atol=1e-14)
                        or not np.isclose(dt, dl/(299792458*np.sqrt(3)), rtol=1e-12)):
                    raise ValueError('native version/grid/time differs')
                source = h['srcs/src1']
                if not np.array_equal(source.attrs['GridPosition'], np.rint(tx/dl).astype(int)):
                    raise ValueError('source grid position differs')
                excitation = source['excitation/samples'][:]
                if excitation.dtype != np.float64 or not np.isfinite(excitation).all():
                    raise ValueError('native source is not finite double')
                for i, position in enumerate(rx):
                    receiver = h[f'rxs/rx{i+1}']
                    if not np.array_equal(receiver.attrs['GridPosition'], np.rint(position/dl).astype(int)):
                        raise ValueError('receiver position differs')
                    for component in ('Ex', 'Ey', 'Ez'):
                        values = receiver[component][:]
                        if (values.dtype != np.float64 or not np.isfinite(values).all()
                                or len(values) != int(h.attrs['Iterations'])):
                            raise ValueError('native receiver validity differs')
                row.update(raw_sha256=sha256(raw), dt_s=dt,
                           iterations=int(h.attrs['Iterations']), dtype='float64', receivers=len(rx))
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Native execution identity; continuum agreement checked separately.'}


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and installed V4.0.0 required')
    out.mkdir(parents=True)
    groups = []
    for dl, pol, direction in CASES:
        name = f'd{dl:g}_{pol}_{direction}'
        folder = out/name; folder.mkdir()
        lines, tx, rx = recipe(dl, pol, direction)
        p = folder/'profile.in'
        p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        groups.append({'id': name, 'input': str(p.resolve()), 'input_sha256': sha256(p),
                       'spacing_m': dl, 'polarisation': pol, 'direction': direction,
                       'tx_m': tx.tolist(), 'receivers_m': [r.tolist() for r in rx]})
    package = Path(gprMax.__file__).parent
    source_names = ['__init__.py', 'sources.py', 'fields_outputs.py', 'model.py', 'waveforms.py',
                    'cuda_opencl/knl_fields_updates.py', 'cuda_opencl/knl_source_updates.py',
                    'cuda_opencl/knl_pml_updates_electric_HORIPML.py',
                    'cuda_opencl/knl_pml_updates_magnetic_HORIPML.py',
                    'toolboxes/SFCW/processing.py']
    source_names += [str(p.relative_to(package)).replace('\\', '/') for p in (package/'cython').glob('*.pyd')]
    code_names = ['scripts/uav_local_benchmark.py', 'scripts/run_uav_local_benchmark.cmd',
                  'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py']
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 若是本地方便快速跑完就本地做完吧 (2026-10-05)',
         'python': str(Path(sys.executable).resolve()), 'source_identities': {n: sha256(package/n) for n in source_names},
         'code_identities': {str(ROOT/n): sha256(ROOT/n) for n in code_names},
         'groups': groups, 'max_runs': len(groups), 'no_retry': True,
         'min_available_RAM_GiB': 1.6, 'min_free_VRAM_GiB': 1.5,
         'max_owned_RSS_GiB': 1.6, 'min_system_available_during_run_GiB': .25,
         'max_group_wall_s': 240, 'max_batch_wall_s': 1200,
         'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'stage': 'round1_free_space_analytic_dipole', 'source': '100MHz Ricker actual-source deconvolution; amp*dl=1Am',
         'processing': 'complex501 tones20-170MHz, tail20ns, co-polar continuum E/(I*dl), no empirical phase/amplitude fit',
         'diagnostic_tolerances': {'complex_band_relative_L2': .05, 'maximum_phase_error_deg': 5.},
         'tolerance_scope': 'Predeclared engineering diagnostic; not field or high-precision certification.'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print('Frozen six free-space cases; maximum 4.096M cells; native GPU double', flush=True)


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
