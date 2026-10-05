"""Bounded local V4 2D gentle-slope terrain B-scan batch, fine grid + 2 m AGL; new immutable capsules only.

Fine variant of hs4_slope_bscan_v0_1 (2026-10-05_slope_bscan_2d_r1): dl halved
to 0.0125/0.05/0.0125 (2880x1x2520 cells, PML cell counts doubled to keep
x+-2 m / z+-1 m physical thickness) and antennas terrain-follow at quantised
surface + 2 m instead of + 8 m. User direction 2026-10-05: skip single-factor
ablation, produce the good sample first - this is the single fine/2 m model;
coarse-8m comparison uses the existing r1 capsule, no c2/f8 intermediates.

26 traces: 13 stations (tx 15.0..21.0 m step 0.5, rx = tx+1.3) x 2 cases
(slope_rough: sloped surface + half-slope rough interface; slope_fullcover:
same surface, all-cover background for the ideal subtraction reference).
Surface: one-sided gentle slope, 2 m drop over x 12..24 m (~9.5 deg), flat
21.5 m (x<=12) and 19.5 m (x>=24). Interface: half-slope trend (1 m drop over
the same interval, 19.5->18.5 m) plus the archived 0.8 m digitised relief
(48 columns of 0.25 m from 2026-10-05_hs4_material_scan_r1, mean removed,
edge values extended outside 12..24). Materials: cover eps_inf 11 / sigma
0.001 / Debye(0.5, 6.4567ns); rock eps 9 / sigma 0.001 (user-frozen
2026-10-05). Ricker 100 MHz actual source, amp*dl = 1 Am. Centre station
(tx=18) of slope_rough carries 200 snapshots at 34-iteration (~1.0018 ns)
spacing for the wavefield GIF. 2D y-invariant TM line-source family;
mechanism screening, not Line9 truth.
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
DX = 0.0125
DY = 0.05
DOMAIN = (36., DY, 31.5)
SURFACE_LEFT = 21.5
SLOPE_DROP = 2.0
SLOPE_X0, SLOPE_X1 = 12.0, 24.0
IFACE_LEFT = 19.5
IFACE_DROP = 1.0
AGL = 2.0
BASELINE = 1.3
STATIONS = [15.0 + 0.5*i for i in range(13)]
TIME_WINDOW = '200e-9'
MATERIALS = ('#material: 9 0.001 1 0 rock',
             '#material: 11 0.001 1 0 cover',
             '#add_dispersion_debye: 1 0.5 6.4567e-09 cover')
ROUGH_SOURCE = ROOT/'artifacts/research_checks/2026-10-05_hs4_material_scan_r1/eps18.017_s0.003_debye_rough/profile.in'
# Cover-column bottoms (z, m) of the 48 archived relief columns over x 12..24 m,
# transcribed from ROUGH_SOURCE; freeze() re-verifies them against that file.
ROUGH_BOTTOMS = (8.85, 8.85, 8.85, 8.85, 8.85, 8.75, 8.75, 8.65, 8.65, 8.65,
                 8.65, 8.75, 8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.35,
                 9.45, 9.45, 9.45, 9.45, 9.45, 9.45, 9.35, 9.35, 9.25, 9.25,
                 9.15, 9.05, 9.05, 8.95, 8.95, 8.95, 8.95, 8.95, 8.95, 8.95,
                 9.05, 9.15, 9.15, 9.25, 9.25, 9.25, 9.25, 9.25)
CASES = ('slope_rough', 'slope_fullcover')
SNAP_GROUP = 'slope_rough_s06'
SNAP_ITER_STEP = 34
SNAP_COUNT = 200
SNAP_REGION = (8.0, 0.0, 14.0, 28.1, DY, 31.4)
SNAP_STEP = (0.15, DY, 0.15)


def fmt(v):
    return f'{v:.12g}'


def rough_zero_mean():
    b = np.array(ROUGH_BOTTOMS, dtype=float)
    return b - b.mean()


def surface_z(x):
    frac = np.clip((np.asarray(x, dtype=float)-SLOPE_X0)/(SLOPE_X1-SLOPE_X0), 0., 1.)
    return SURFACE_LEFT - SLOPE_DROP*frac


def interface_z(x):
    x = np.asarray(x, dtype=float)
    frac = np.clip((x-SLOPE_X0)/(SLOPE_X1-SLOPE_X0), 0., 1.)
    trend = IFACE_LEFT - IFACE_DROP*frac
    col = np.clip(np.floor((x-SLOPE_X0)/0.25), 0, len(ROUGH_BOTTOMS)-1).astype(int)
    return trend + rough_zero_mean()[col]


def grid_index(z):
    return int(np.rint(np.asarray(z, dtype=float)/DX))


def column_grids():
    """Per-column integer grid heights (surface, interface) at column centres."""
    nx = round(DOMAIN[0]/DX)
    centres = DX*(np.arange(nx)+0.5)
    ks = np.rint(surface_z(centres)/DX).astype(int)
    ki = np.rint(interface_z(centres)/DX).astype(int)
    if not (np.all(ki > 0) and np.all(ki < ks) and np.all(ks < round(DOMAIN[2]/DX))):
        raise ValueError('interface/surface ordering violated')
    return ks, ki


def antenna_z(tx_x):
    return grid_index(surface_z(tx_x))*DX + AGL


def merged_runs(keys):
    runs = []
    start = 0
    for i in range(1, len(keys)+1):
        if i == len(keys) or keys[i] != keys[start]:
            runs.append((start, i, keys[start]))
            start = i
    return runs


def recipe(case, tx_x):
    z_src = antenna_z(tx_x)
    lines = [f'#title: 2D gentle-slope B-scan v0.1 fine grid; {case}; tx={fmt(tx_x)}; terrain-follow 2m AGL',
             f'#domain: {fmt(DOMAIN[0])} {fmt(DOMAIN[1])} {fmt(DOMAIN[2])}',
             f'#dx_dy_dz: {DX:g} {DY:g} {DX:g}', f'#time_window: {TIME_WINDOW}',
             '#omp_threads: 2', '#pml_cells: 160 0 80 160 0 80', '#pml_formulation: HORIPML',
             f'#waveform: ricker {1/DX:.12g} 100e6 pulse',
             f'#hertzian_dipole: y {fmt(tx_x)} {DY/2:g} {fmt(z_src)} pulse',
             f'#rx: {fmt(tx_x+BASELINE)} {DY/2:g} {fmt(z_src)} r1 Ey']
    lines += MATERIALS
    ks, ki = column_grids()
    if case == 'slope_rough':
        for i0, i1, key in merged_runs(list(zip(ki, ks))):
            k_rock, k_surf = key
            lines.append(f'#box: {fmt(i0*DX)} 0 0 {fmt(i1*DX)} {DY:g} {fmt(k_rock*DX)} rock')
            lines.append(f'#box: {fmt(i0*DX)} 0 {fmt(k_rock*DX)} {fmt(i1*DX)} {DY:g} {fmt(k_surf*DX)} cover')
    elif case == 'slope_fullcover':
        for i0, i1, k_surf in merged_runs(list(ks)):
            lines.append(f'#box: {fmt(i0*DX)} 0 0 {fmt(i1*DX)} {DY:g} {fmt(k_surf*DX)} cover')
    else:
        raise ValueError('unknown case')
    if f'{case}_s{STATIONS.index(tx_x):02d}' == SNAP_GROUP:
        r = SNAP_REGION + SNAP_STEP
        for k in range(SNAP_COUNT):
            lines.append(f'#snapshot: {fmt(r[0])} {fmt(r[1])} {fmt(r[2])} {fmt(r[3])} {fmt(r[4])} {fmt(r[5])}'
                         f' {fmt(r[6])} {fmt(r[7])} {fmt(r[8])} {k*SNAP_ITER_STEP} snap{k:04d}.h5')
    return lines


def verify_rough_source():
    parsed = []
    for line in ROUGH_SOURCE.read_text('utf-8').splitlines():
        if not line.startswith('#box:') or not line.endswith(' cover'):
            continue
        p = [float(v) for v in line.split()[1:7]]
        if abs(p[5]-12.0) < 1e-12 and abs(p[3]-p[0]-0.25) < 1e-12 and SLOPE_X0-1e-12 <= p[0] < SLOPE_X1:
            parsed.append((p[0], p[2]))
    parsed.sort()
    bottoms = tuple(z for _, z in parsed)
    if len(bottoms) != len(ROUGH_BOTTOMS) or any(abs(a-b) > 1e-9 for a, b in zip(bottoms, ROUGH_BOTTOMS)):
        raise ValueError('archived relief columns differ from transcribed constants')


def audit(path, completed=False):
    contract = json.loads(path.read_text('utf-8'))
    if len(contract['groups']) != 26:
        raise ValueError('declared 26-trace slope B-scan matrix required')
    rows = []
    for group in contract['groups']:
        p = Path(group['input'])
        lines = recipe(group['case'], group['tx_x_m'])
        if sha256(p) != group['input_sha256'] or p.read_text('utf-8').splitlines() != lines:
            raise ValueError('input recipe identity differs')
        z_src = antenna_z(group['tx_x_m'])
        tx = np.array([group['tx_x_m'], DY/2, z_src])
        rx = np.array([group['tx_x_m']+BASELINE, DY/2, z_src])
        for position in (tx, rx):
            xz = np.array([position[0], position[2]])
            if not np.allclose(xz/DX, np.rint(xz/DX), atol=1e-10, rtol=0):
                raise ValueError('off-grid source/receiver')
        row = {'id': group['id'], 'input_sha256': sha256(p), 'dx_m': DX,
               'case': group['case'], 'tx_x_m': group['tx_x_m'], 'antenna_z_m': z_src,
               'source_current_moment_amplitude_Am': 1.}
        if completed:
            raw = p.with_suffix('.h5')
            with h5py.File(raw) as h:
                dt = float(h.attrs['dt'])
                shape = [round(DOMAIN[0]/DX), round(DOMAIN[1]/DY), round(DOMAIN[2]/DX)]
                if (str(h.attrs['gprMax']) != '4.0.0'
                        or not np.array_equal(h.attrs['nx_ny_nz'], shape)
                        or not np.allclose(h.attrs['dx_dy_dz'], [DX, DY, DX], rtol=0, atol=1e-14)
                        or not np.isclose(dt, DX/(299792458*np.sqrt(2)), rtol=1e-12)):
                    raise ValueError('native version/grid/time differs')
                if not np.array_equal(h['srcs/src1'].attrs['GridPosition'],
                                      np.rint(tx/[DX, DY, DX]).astype(int)):
                    raise ValueError('source grid position differs')
                if not np.array_equal(h['rxs/rx1'].attrs['GridPosition'],
                                      np.rint(rx/[DX, DY, DX]).astype(int)):
                    raise ValueError('receiver position differs')
                values = h['rxs/rx1/Ey'][:]
                if (values.dtype != np.float64 or not np.isfinite(values).all()
                        or len(values) != int(h.attrs['Iterations'])):
                    raise ValueError('native receiver validity differs')
                row.update(raw_sha256=sha256(raw), dt_s=dt,
                           iterations=int(h.attrs['Iterations']), dtype='float64')
            if group['id'] == SNAP_GROUP:
                snaps = sorted((p.parent/'profile_snaps').glob('snap*.h5'))
                expect_nx = round((SNAP_REGION[3]-SNAP_REGION[0])/SNAP_STEP[0])
                expect_nz = round((SNAP_REGION[5]-SNAP_REGION[2])/SNAP_STEP[2])
                if len(snaps) != SNAP_COUNT:
                    raise ValueError('snapshot count differs')
                for k, snap in enumerate(snaps):
                    with h5py.File(snap) as h:
                        if (int(h.attrs['iteration']) != k*SNAP_ITER_STEP
                                or h['Ey'].shape != (expect_nx, 1, expect_nz)
                                or not np.isfinite(h['Ey'][:]).all()):
                            raise ValueError('snapshot identity/validity differs')
                row['snapshots'] = SNAP_COUNT
        rows.append(row)
    return {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
            'code_sha256': sha256(__file__), 'groups': rows,
            'scope': 'Native execution identity; B-scan physics analysed separately.'}


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and installed V4.0.0 required')
    verify_rough_source()
    out.mkdir(parents=True)
    groups = []
    for case in CASES:
        for i, tx_x in enumerate(STATIONS):
            name = f'{case}_s{i:02d}'
            folder = out/name
            folder.mkdir()
            p = folder/'profile.in'
            p.write_text('\n'.join(recipe(case, tx_x))+'\n', encoding='utf-8')
            z_src = antenna_z(tx_x)
            groups.append({'id': name, 'input': str(p.resolve()), 'input_sha256': sha256(p),
                           'case': case, 'station_index': i, 'tx_x_m': tx_x, 'spacing_m': DX,
                           'tx_m': [tx_x, DY/2, z_src],
                           'receivers_m': [[tx_x+BASELINE, DY/2, z_src]]})
    package = Path(gprMax.__file__).parent
    source_names = ['__init__.py', 'sources.py', 'fields_outputs.py', 'model.py', 'waveforms.py',
                    'cuda_opencl/knl_fields_updates.py', 'cuda_opencl/knl_source_updates.py',
                    'cuda_opencl/knl_pml_updates_electric_HORIPML.py',
                    'cuda_opencl/knl_pml_updates_magnetic_HORIPML.py',
                    'toolboxes/SFCW/processing.py']
    source_names += [str(p.relative_to(package)).replace('\\', '/') for p in (package/'cython').glob('*.pyd')]
    code_names = ['scripts/hs4_slope_bscan_fine2m_v0_1.py', 'scripts/run_slope_bscan_fine2m_4090.cmd',
                  'scripts/hs4_station_grid_controls.py', 'scripts/gprmax_cached_cuda_entry.py']
    ks, ki = column_grids()
    c = {'status': 'FROZEN_APPROVED', 'approval_basis': 'User: 不必，我们现在不考虑这个单因素，我们需要的是先做出好的样本再消融。只跑一个 (2026-10-05); single fine-grid 2 m AGL gentle-slope terrain-following B-scan with ideal + catalogue-operator background removal and dense wavefield snapshots; no c2/f8 intermediates.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': {n: sha256(package/n) for n in source_names},
         'code_identities': {str(ROOT/n): sha256(ROOT/n) for n in code_names},
         'groups': groups, 'max_runs': len(groups), 'no_retry': True,
         'min_available_RAM_GiB': 4., 'min_free_VRAM_GiB': 2.,
         'max_owned_RSS_GiB': 8., 'min_system_available_during_run_GiB': 1.,
         'max_group_wall_s': 300, 'max_batch_wall_s': 7800,
         'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'stage': 'round7_slope_bscan_2d_fine2m',
         'parent_batch': '2026-10-05_slope_bscan_2d_r1 (coarse dl0.025 + 8 m AGL; reused as the coarse/tall comparison cell, no rerun)',
         'geometry': {'domain_m': list(DOMAIN), 'dx_dy_dz_m': [DX, DY, DX], 'pml_cells': [160, 0, 80, 160, 0, 80],
                      'surface': {'left_flat_z_m': SURFACE_LEFT, 'drop_m': SLOPE_DROP,
                      'slope_interval_m': [SLOPE_X0, SLOPE_X1], 'right_flat_z_m': SURFACE_LEFT-SLOPE_DROP},
                      'interface': {'left_trend_z_m': IFACE_LEFT, 'drop_m': IFACE_DROP,
                                    'trend_interval_m': [SLOPE_X0, SLOPE_X1],
                                    'relief': 'archived 48x0.25m columns, mean removed, edge values extended; quantized to 0.0125 m',
                                    'relief_source': str(ROUGH_SOURCE), 'relief_source_sha256': sha256(ROUGH_SOURCE),
                                    'relief_range_m': [float(min(ROUGH_BOTTOMS)-np.mean(ROUGH_BOTTOMS)),
                                                       float(max(ROUGH_BOTTOMS)-np.mean(ROUGH_BOTTOMS))]},
                      'cover_thickness_m': [float((ks-ki).min()*DX), float((ks-ki).max()*DX)],
                      'stations_tx_x_m': STATIONS, 'baseline_m': BASELINE, 'agl_m': AGL,
                      'antenna_z_m': [antenna_z(x) for x in STATIONS],
                      'snapshot': {'group': SNAP_GROUP, 'region_m': list(SNAP_REGION), 'step_m': list(SNAP_STEP),
                                   'iteration_step': SNAP_ITER_STEP, 'count': SNAP_COUNT}},
         'materials': {'cover': 'eps_inf 11, sigma 0.001, Debye dE 0.5 tau 6.4567ns (user-frozen default 2026-10-05)',
                       'rock': 'eps 9, sigma 0.001 (unchanged)'},
         'source': '100MHz Ricker actual-source deconvolution; amp*dl=1Am; 2D y-invariant TM line-source family, not finite 3D',
         'processing': 'native signed B-scans; ideal background = slope_rough minus slope_fullcover (declared pair); catalogue mean/SVD operators applied in analysis; t^2 declared analysis gain labelled non-amplitude',
         'scope_note': 'single synthetic gentle-slope mechanism model, fine grid + 2 m AGL sample production (ablation deferred per user), not Yingshan site adaptation'}
    supervisor.save(out/'execution_contract.json', c)
    supervisor.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(f'Frozen {len(groups)} fine-grid 2m-AGL slope B-scan traces; native GPU double', flush=True)


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
