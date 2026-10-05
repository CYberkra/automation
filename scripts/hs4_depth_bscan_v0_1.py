"""Bounded local V4 2D depth B-scan batch; new immutable capsules only.

52 traces: 13 stations (tx 15.0..21.0 m step 0.5, rx = tx+1.3, centre 18.0
matching the depth-factor batch) x 4 cases (fullcover / flat interface at
3 / 10 / 18.5 m). Same shared domain, surface, air path, materials, grid,
PML and 650 ns window as 2026-10-05_depth_factor_2d_r1; only station and
interface depth change. 2D y-invariant TM line-source family; mechanism
screening, depths are mechanism points, not Line9 truth.
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
DX = 0.025
DOMAIN = (36., 0.05, 31.5)
SURFACE_Z = 21.5
SOURCE_Z = 29.5
STATIONS = [15.0 + 0.5*i for i in range(13)]
BASELINE = 1.3
TIME_WINDOW = '650e-9'
MATERIALS = ('#material: 9 0.001 1 0 rock',
             '#material: 18.017 0.003 1 0 cover',
             '#add_dispersion_debye: 1 7.878 6.4567e-09 cover')
DEPTHS = {'fullcover': None, 'depth3': 3., 'depth10': 10., 'depth18_5': 18.5}


def fmt(v):
    return f'{v:.12g}'


def recipe(case, tx_x):
    rx_x = tx_x + BASELINE
    lines = [f'#title: 2D depth B-scan v0.1; {case}; tx={fmt(tx_x)}; 8m altitude',
             f'#domain: {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(DOMAIN[2])}',
             f'#dx_dy_dz: {DX:g} 0.05 {DX:g}', f'#time_window: {TIME_WINDOW}',
             '#omp_threads: 2', '#pml_cells: 80 0 40 80 0 40', '#pml_formulation: HORIPML',
             f'#waveform: ricker {1/DX:.12g} 100e6 pulse',
             f'#hertzian_dipole: y {fmt(tx_x)} 0.025 {fmt(SOURCE_Z)} pulse',
             f'#rx: {fmt(rx_x)} 0.025 {fmt(SOURCE_Z)} r1 Ey']
    lines += MATERIALS
    depth = DEPTHS[case]
    if depth is None:
        lines.append(f'#box: 0 0 0 {fmt(DOMAIN[0])} 0.05 {fmt(SURFACE_Z)} cover')
    else:
        iface = SURFACE_Z - depth
        lines.append(f'#box: 0 0 0 {fmt(DOMAIN[0])} 0.05 {fmt(iface)} rock')
        lines.append(f'#box: 0 0 {fmt(iface)} {fmt(DOMAIN[0])} 0.05 {fmt(SURFACE_Z)} cover')
    return lines


def audit(path, completed=False):
    contract = json.loads(path.read_text('utf-8'))
    if len(contract['groups']) != 52:
        raise ValueError('declared 52-trace B-scan matrix required')
    rows = []
    for group in contract['groups']:
        p = Path(group['input'])
        lines = recipe(group['case'], group['tx_x_m'])
        if sha256(p) != group['input_sha256'] or p.read_text('utf-8').splitlines() != lines:
            raise ValueError('input recipe identity differs')
        tx = np.array([group['tx_x_m'], 0.025, SOURCE_Z])
        rx = np.array([group['tx_x_m']+BASELINE, 0.025, SOURCE_Z])
        for position in (tx, rx):
            xz = np.array([position[0], position[2]])
            if not np.allclose(xz/DX, np.rint(xz/DX), atol=1e-10, rtol=0):
                raise ValueError('off-grid source/receiver')
        row = {'id': group['id'], 'input_sha256': sha256(p), 'dx_m': DX,
               'case': group['case'], 'tx_x_m': group['tx_x_m'],
               'source_current_moment_amplitude_Am': 1.}
        if completed:
            raw = p.with_suffix('.h5')
            with h5py.File(raw) as h:
                dt = float(h.attrs['dt'])
                shape = [round(DOMAIN[0]/DX), round(DOMAIN[1]/0.05), round(DOMAIN[2]/DX)]
                if (str(h.attrs['gprMax']) != '4.0.0'
                        or not np.array_equal(h.attrs['nx_ny_nz'], shape)
                        or not np.allclose(h.attrs['dx_dy_dz'], [DX, 0.05, DX], rtol=0, atol=1e-14)
                        or not np.isclose(dt, DX/(299792458*np.sqrt(2)), rtol=1e-12)):
                    raise ValueError('native version/grid/time differs')
                if not np.array_equal(h['srcs/src1'].attrs['GridPosition'],
                                      np.rint(tx/[DX, 0.05, DX]).astype(int)):
                    raise ValueError('source grid position differs')
                if not np.array_equal(h['rxs/rx1'].attrs['GridPosition'],
                                      np.rint(rx/[DX, 0.05, DX]).astype(int)):
                    raise ValueError('receiver position differs')
                values = h['rxs/rx1/Ey'][:]
                if (values.dtype != np.float64 or not np.isfinite(values).all()
                        or len(values) != int(h.attrs['Iterations'])):
                    raise ValueError('native receiver validity differs')
                row.update(raw_sha256=sha256(raw), dt_s=dt,
                           iterations=int(h.attrs['Iterations']), dtype='float64')
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Native execution identity; B-scan physics analysed separately.'}


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and installed V4.0.0 required')
    out.mkdir(parents=True)
    groups = []
    for case in DEPTHS:
        for i, tx_x in enumerate(STATIONS):
            name = f'{case}_s{i:02d}'
            folder = out/name
            folder.mkdir()
            p = folder/'profile.in'
            p.write_text('\n'.join(recipe(case, tx_x))+'\n', encoding='utf-8')
            groups.append({'id': name, 'input': str(p.resolve()), 'input_sha256': sha256(p),
                           'case': case, 'station_index': i, 'tx_x_m': tx_x, 'spacing_m': DX,
                           'tx_m': [tx_x, 0.025, SOURCE_Z],
                           'receivers_m': [[tx_x+BASELINE, 0.025, SOURCE_Z]]})
    package = Path(gprMax.__file__).parent
    source_names = ['__init__.py', 'sources.py', 'fields_outputs.py', 'model.py', 'waveforms.py',
                    'cuda_opencl/knl_fields_updates.py', 'cuda_opencl/knl_source_updates.py',
                    'cuda_opencl/knl_pml_updates_electric_HORIPML.py',
                    'cuda_opencl/knl_pml_updates_magnetic_HORIPML.py',
                    'toolboxes/SFCW/processing.py']
    source_names += [str(p.relative_to(package)).replace('\\', '/') for p in (package/'cython').glob('*.pyd')]
    code_names = ['scripts/hs4_depth_bscan_v0_1.py', 'scripts/run_depth_bscan_4090.cmd',
                  'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py']
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: bscan呢？ (2026-10-05); B-scan visualization of the depth-factor findings.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': {n: sha256(package/n) for n in source_names},
         'code_identities': {str(ROOT/n): sha256(ROOT/n) for n in code_names},
         'groups': groups, 'max_runs': len(groups), 'no_retry': True,
         'min_available_RAM_GiB': 2., 'min_free_VRAM_GiB': 2.,
         'max_owned_RSS_GiB': 6., 'min_system_available_during_run_GiB': 1.,
         'max_group_wall_s': 120, 'max_batch_wall_s': 6240,
         'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'stage': 'round5_depth_bscan_2d',
         'geometry': {'domain_m': list(DOMAIN), 'surface_z_m': SURFACE_Z, 'source_z_m': SOURCE_Z,
                      'stations_tx_x_m': STATIONS, 'baseline_m': BASELINE, 'interface_depths_m': DEPTHS,
                      'parent_batch': '2026-10-05_depth_factor_2d_r1 (same domain/materials/grid/window; centre station 18.0 identical)'},
         'source': '100MHz Ricker actual-source deconvolution; amp*dl=1Am; 2D y-invariant TM line-source family, not finite 3D',
         'processing': 'native signed B-scans + per-station isolated interface (flat minus shared fullcover), declared display enhancement labelled',
         'depths_note': '3/10/18.5 m are mechanism points from the design doc, not Line9-inverted truth'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen {len(groups)} depth B-scan traces; native GPU double', flush=True)


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
