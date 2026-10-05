"""Bounded local V4 3D layout/polarization/scenario batch; new immutable capsules only.

12 center-station traces: 2 baselines (inline (1.3,0,0) / cross-track (0,1.3,0))
x 2 horizontal polarizations (x=flight, y=cross) x 3 subsurface scenarios
(full cover / flat interface / archived 0.8m peak-to-peak relief extruded
along y). tx/rx midpoint fixed; only baseline orientation, declared
polarization and scenario change. Ideal 3D Hertzian dipole, 8 m altitude,
archived cover/rock materials; not a UAV hardware model.
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
DOMAIN = (10., 10., 16.)           # x=flight, y=cross, z=up; surface z=6; top PML 15..16
MIDPOINT = np.array([5., 5., 14.])  # 8 m above surface
SURFACE_Z = 6.
FLAT_INTERFACE_Z = 3.05            # 122 cells; archived relief mean 3.06 m rounded to grid (1 cm documented)
RELIEF_SHIFT = -6.0                # archived bottoms 8.65..9.45 -> 2.65..3.45 m, all on 0.05 m grid
ARCHIVED_BOTTOMS = [8.85, 8.85, 8.85, 8.75, 8.75, 8.65, 8.65, 8.65, 8.65, 8.75,
                    8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.35, 9.45, 9.45,
                    9.45, 9.45, 9.45, 9.45, 9.35, 9.35, 9.25, 9.25, 9.15, 9.05,
                    9.05, 8.95, 8.95, 8.95, 8.95, 8.95, 8.95, 8.95, 9.05, 9.15]
MATERIALS = ('#material: 9 0.001 1 0 rock',
             '#material: 18.017 0.003 1 0 cover',
             '#add_dispersion_debye: 1 7.878 6.4567e-09 cover')
LAYOUTS = {'inline': np.array([1.3, 0., 0.]), 'cross': np.array([0., 1.3, 0.])}
SCENARIOS = ('fullcover', 'flat', 'rough')


def fmt(v):
    return f'{v:.12g}'


def relief_bottoms():
    """Per 0.25 m column (40 columns over x=0..10), interface z on grid."""
    return [b + RELIEF_SHIFT for b in ARCHIVED_BOTTOMS]


def recipe(scenario, layout, pol):
    base = LAYOUTS[layout]
    tx = MIDPOINT - base/2
    rx = MIDPOINT + base/2
    lines = [f'#title: 3D layout/pol control v0.1; {scenario} {layout} {pol}-pol; centre station; 8m altitude',
             f'#domain: {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(DOMAIN[2])}',
             f'#dx_dy_dz: {DL:g} {DL:g} {DL:g}', '#time_window: 200e-9',
             '#omp_threads: 2', '#pml_cells: 40', '#pml_formulation: HORIPML',
             f'#waveform: ricker {1/DL:.12g} 100e6 pulse',
             f'#hertzian_dipole: {pol} ' + ' '.join(fmt(x) for x in tx) + ' pulse',
             '#rx: ' + ' '.join(fmt(x) for x in rx) + ' r1 Ex Ey Ez']
    lines += MATERIALS
    if scenario == 'fullcover':
        lines.append(f'#box: 0 0 0 {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(SURFACE_Z)} cover')
    elif scenario == 'flat':
        lines.append(f'#box: 0 0 0 {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(FLAT_INTERFACE_Z)} rock')
        lines.append(f'#box: 0 0 {fmt(FLAT_INTERFACE_Z)} {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(SURFACE_Z)} cover')
    else:
        for i, bottom in enumerate(relief_bottoms()):
            x0, x1 = i*0.25, (i+1)*0.25
            lines.append(f'#box: {fmt(x0)} 0 0 {fmt(x1)} {fmt(DOMAIN[1])} {fmt(bottom)} rock')
            lines.append(f'#box: {fmt(x0)} 0 {fmt(bottom)} {fmt(x1)} {fmt(DOMAIN[1])} {fmt(SURFACE_Z)} cover')
    return lines, tx, rx


def audit(path, completed=False):
    contract = json.loads(path.read_text('utf-8'))
    if len(contract['groups']) != 12:
        raise ValueError('declared 12-trace matrix required')
    rows = []
    for group in contract['groups']:
        scenario, layout, pol = group['scenario'], group['layout'], group['polarisation']
        p = Path(group['input'])
        lines, tx, rx = recipe(scenario, layout, pol)
        if sha256(p) != group['input_sha256'] or p.read_text('utf-8').splitlines() != lines:
            raise ValueError('input recipe identity differs')
        for position in (tx, rx):
            if not np.allclose(position/DL, np.rint(position/DL), atol=1e-10, rtol=0):
                raise ValueError('off-grid source/receiver')
        row = {'id': group['id'], 'input_sha256': sha256(p), 'dl_m': DL,
               'scenario': scenario, 'layout': layout, 'polarisation': pol,
               'source_current_moment_amplitude_Am': 1.}
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
                if not np.array_equal(source.attrs['GridPosition'], np.rint(tx/DL).astype(int)):
                    raise ValueError('source grid position differs')
                excitation = source['excitation/samples'][:]
                if excitation.dtype != np.float64 or not np.isfinite(excitation).all():
                    raise ValueError('native source is not finite double')
                receiver = h['rxs/rx1']
                if not np.array_equal(receiver.attrs['GridPosition'], np.rint(rx/DL).astype(int)):
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
            'scope': 'Native execution identity; layout/polarization effects analysed separately.'}


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and installed V4.0.0 required')
    out.mkdir(parents=True)
    groups = []
    for scenario in SCENARIOS:
        for layout in ('inline', 'cross'):
            for pol in ('x', 'y'):
                name = f'{scenario}_{layout}_{pol}'
                folder = out/name
                folder.mkdir()
                lines, tx, rx = recipe(scenario, layout, pol)
                p = folder/'profile.in'
                p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
                groups.append({'id': name, 'input': str(p.resolve()), 'input_sha256': sha256(p),
                               'scenario': scenario, 'layout': layout, 'polarisation': pol,
                               'spacing_m': DL, 'tx_m': tx.tolist(), 'receivers_m': [rx.tolist()]})
    package = Path(gprMax.__file__).parent
    source_names = ['__init__.py', 'sources.py', 'fields_outputs.py', 'model.py', 'waveforms.py',
                    'cuda_opencl/knl_fields_updates.py', 'cuda_opencl/knl_source_updates.py',
                    'cuda_opencl/knl_pml_updates_electric_HORIPML.py',
                    'cuda_opencl/knl_pml_updates_magnetic_HORIPML.py',
                    'toolboxes/SFCW/processing.py']
    source_names += [str(p.relative_to(package)).replace('\\', '/') for p in (package/'cython').glob('*.pyd')]
    code_names = ['scripts/hs4_layout_pol_batch_v0_1.py', 'scripts/run_layout_pol_batch_4090.cmd',
                  'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py']
    relief = relief_bottoms()
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 跑，你去跑完 (2026-10-05); design docs/research/2026-10-05_sim2real_next_step_design.md section B stage 1.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': {n: sha256(package/n) for n in source_names},
         'code_identities': {str(ROOT/n): sha256(ROOT/n) for n in code_names},
         'groups': groups, 'max_runs': len(groups), 'no_retry': True,
         'min_available_RAM_GiB': 4., 'min_free_VRAM_GiB': 11.,
         'max_owned_RSS_GiB': 24., 'min_system_available_during_run_GiB': 2.,
         'max_group_wall_s': 900, 'max_batch_wall_s': 10800,
         'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'stage': 'round3_layout_polarization_scenario_3d',
         'geometry': {'domain_m': list(DOMAIN), 'surface_z_m': SURFACE_Z, 'tx_rx_midpoint_m': MIDPOINT.tolist(),
                      'baselines_m': {k: v.tolist() for k, v in LAYOUTS.items()},
                      'flat_interface_z_m': FLAT_INTERFACE_Z,
                      'relief': {'generator': 'archived digitized profile (material-scan profile.in) central x 12.5..22.5 m extruded along y, shifted by -6.0 m',
                                 'bottoms_m': relief, 'mean_z_m': float(np.mean(relief)),
                                 'peak_to_peak_m': float(max(relief)-min(relief)),
                                 'flat_deviation_note': 'flat interface 3.05 m = relief mean 3.0525 m rounded to nearest cell (0.1 cell, 2.5 mm, grid quantization of 0.05 m-step archived profile)'}},
         'grid_precheck': {'dl_m': DL, 'cells': [round(d/DL) for d in DOMAIN],
                           'baseline_cells': 52, 'altitude_cells': 320, 'relief_pp_cells': 32,
                           'cover_lambda_170MHz_m': 0.4136, 'cells_per_lambda_170MHz': 16.5,
                           'dl0.05_rejected': 'only 8.3 cells/lambda in cover at 170 MHz'},
         'source': '100MHz Ricker actual-source deconvolution; amp*dl=1Am; ideal Hertzian, not a UAV hardware model',
         'processing': 'native signed A-scans + complex501 tones20-170MHz, tail20ns; no empirical phase/amplitude fit',
         'reference': 'round2 half-space Sommerfeld validation PASS (2026-10-05_halfspace_point_3d_r2) grounds this geometry/material chain',
         'diagnostic_note': 'layout effect must exceed this-batch numerical floor (pre-first-arrival cross-scene difference and 1e-12 DFT identity); no preset dB threshold for "main cause"',
         'resource_estimate': {'cells': 102400000, 'vram_GiB_approx': 12, 'wall_per_trace_s_approx': 400}}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen {len(groups)} layout/polarization/scenario traces; native GPU double', flush=True)


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
