"""Stage B: 3D line segment with bowtie antenna pair over cropped line9 terrain.

Approved scope (user, 2026-10-12, before sleep): after Stage A antenna
characterization, start a 3D line simulation with the new antenna so a line
result exists by morning; domain/hardware adjustments pre-authorized to fit
the ROG 4090 (16 GB). Goal ends once the batch is verified running normally.

Design (bounded):
- 3 anchor stations profile_x 170/180/190 m (170 anchors the window
  convergence set; 170-190 spans the slope segment studied earlier).
- Per station a local 3D crop of the H1 (with-landslide) geometry:
  x +/-6 m along track, y 6.0 m cross-track extrusion (cross-track invariant
  terrain: documented approximation, same spirit as earlier local3D work),
  z = local relief + 11.5 m ground (incl. 1.5 m PML) + 7.25 m air
  (5 m AGL antenna + 0.75 m clearance + 1.5 m PML).
- Antennas: Stage A unfocused bowtie approximation, arms along y (broadside),
  Tx/Rx baseline 1.3 m along x, each 5 m above its own local surface.
  Tx centre-fed 50 ohm voltage source with official impulse waveform;
  Rx bowtie passive PEC, gap field Ey recorded as the A-scan.
- Window 900 ns (covers bottom-sand band 360-400 ns and the surface multiples).
- FP64 GPU, PML 60 cells HORIPML, dispersive averaging on, same material
  database as the 2D line work (line9_materials_v0_2).
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np

import gprMax
import line9_v401_version_controls as native
import hs4_station_grid_controls as supervisor
from antenna_stage_a import DT, DL, PML, bowtie_lines, fmt, DB_ID
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]
PARENT_GEO = ROOT/'model_current/geometries/H1.h5'
DB_SOURCE = ROOT/'model_current/geometries/line9_materials_v0_2.json'

APPROVAL = ('2026-10-12 用户批准「做真实天线吧…一切妥当后启用3D带新设计的天线后的测线仿真」，'
            '并预授权为使 ROG 可运行而调整域大小等配套。3D 测线批：锚点站 170/180/190，'
            '横轨不变地形近似，领结对（未聚焦原天线近似），900 ns 窗。')

WINDOW_S = 9e-7
STEPS = int(math.ceil(WINDOW_S/DT))+1
STATIONS = [170.0, 180.0, 190.0]
X_HALF_M = 6.0
Y_SPAN_M = 6.0
GROUND_DEPTH_M = 11.5      # includes 1.5 m bottom PML
AIR_TOP_M = 7.25           # 5 m AGL + 0.75 m clearance + 1.5 m top PML
BASELINE_M = 1.3
AGL_M = 5.0
ELEV_OFFSET_M = 417.5
PROFILE_OFFSET_M = 20.0


def surface_elevations(d2, i0, i1):
    slab = d2[i0:i1, :]
    nonair = slab != 0
    assert nonair.any() and not nonair[:, -1].any(), 'crop touches top of parent geometry'
    jlast = slab.shape[1]-1-np.argmax(nonair[:, ::-1], axis=1)
    return ELEV_OFFSET_M+(jlast+1)*DL   # first air cell bottom = surface level


def make_crop(path, sx):
    with h5py.File(PARENT_GEO) as h:
        d2 = h['data'][:, :, 0]
        assert h.attrs['MaterialDatabase'] == DB_ID
        i0 = int(round((sx-X_HALF_M-PROFILE_OFFSET_M)/DL))
        i1 = i0+int(round(2*X_HALF_M/DL))
        surf = surface_elevations(d2, i0, i1)
        relief = float(surf.max()-surf.min())
        z_bot = float(surf.min())-GROUND_DEPTH_M
        z_top = float(surf.max())+AIR_TOP_M
        nz = int(round((z_top-z_bot)/DL))
        assert abs(z_bot+nz*DL-z_top) < 1e-12
        j0 = int(round((z_bot-ELEV_OFFSET_M)/DL))
        assert j0 >= 0 and j0+nz <= d2.shape[1]
        ny = int(round(Y_SPAN_M/DL))
        nx = i1-i0
        with h5py.File(path, 'x') as out:
            data = out.create_dataset('data', shape=(nx, ny, nz), dtype='int16',
                                      chunks=(min(nx, 8), min(ny, 64), min(nz, 64)),
                                      compression='gzip', compression_opts=1)
            for first in range(0, nx, 8):
                last = min(nx, first+8)
                plane = d2[first+i0:last+i0, j0:j0+nz]
                data[first:last] = np.broadcast_to(plane[:, None, :], (last-first, ny, nz))
            out.create_dataset('material_keys', data=np.asarray(
                ['material_000_air', 'material_001_cover', 'material_002_mudstone', 'material_003_sandstone'], dtype='S'))
            out.attrs['dx_dy_dz'] = [DL]*3
            out.attrs['origin_xyz'] = [0., 0., 0.]
            out.attrs['shape_nxyz'] = [nx, ny, nz]
            out.attrs['MaterialDatabase'] = DB_ID
            out.attrs['MaterialDatabaseSchemaVersion'] = 1
            out.attrs['parent_geometry'] = str(PARENT_GEO)
            out.attrs['parent_i0'] = i0
            out.attrs['parent_j0'] = j0
            out.attrs['station_profile_x_m'] = sx
            out.attrs['cross_track_approximation'] = 'terrain extruded constant across y'
        with h5py.File(path, 'r') as out:
            assert np.array_equal(out['data'][:, 0, :], d2[i0:i1, j0:j0+nz])
            assert np.array_equal(out['data'][:, -1, :], d2[i0:i1, j0:j0+nz])
    return dict(nx=nx, ny=ny, nz=nz, z_bot=z_bot, relief=relief, surf=surf,
                centre_i=int(round((sx-PROFILE_OFFSET_M)/DL))-i0)


def case_text(case_id, geo, geo_file):
    nx, ny, nz = geo['nx'], geo['ny'], geo['nz']
    return [f'#title: {case_id}; 3D bowtie pair over cropped line9 terrain; stage B',
            f'#domain: {fmt([nx*DL, ny*DL, nz*DL])}',
            f'#dx_dy_dz: {fmt([DL]*3)}',
            f'#time_window: {WINDOW_S:.12g}',
            '#omp_threads: 2',
            f'#pml_cells: {PML} {PML} {PML} {PML} {PML} {PML}',
            '#pml_formulation: HORIPML',
            '#dispersive_averaging: y',
            '#waveform: impulse 1 1 impulse']


def prepare(a):
    assert not a.out.exists()
    (a.out/'cases').mkdir(parents=True)
    (a.out/'geometries').mkdir()
    db = a.out/'geometries'/DB_SOURCE.name
    shutil.copyfile(DB_SOURCE, db)
    groups = []
    for sx in STATIONS:
        cid = f'ant3d_x{round(sx*100):05d}'
        geo = make_crop(a.out/'geometries'/f'{cid}.h5', sx)
        ny = geo['ny']
        cy = ny*DL/2
        feeds = {}
        for sign, label in [(-1, 'tx'), (1, 'rx')]:
            i = geo['centre_i']+int(round(sign*BASELINE_M/2/DL))
            surf_here = geo['surf'][i]
            z = surf_here-geo['z_bot']+AGL_M
            z = round(z/DL)*DL
            feeds[label] = (round(i*DL, 10), round(cy, 10), round(z, 10))
        extra = bowtie_lines(feeds['rx'], 'y', None, passive_id='rxgap')
        lines = case_text(cid, geo, f'{cid}.h5')
        lines += bowtie_lines(feeds['tx'], 'y', 'impulse', window_s=WINDOW_S) + extra
        lines.append(f'#geometry_objects_read: 0 0 0 ../../geometries/{cid}.h5 {DB_ID} n')
        lines.append('')
        case = a.out/'cases'/cid
        case.mkdir()
        p = case/'profile.in'
        p.write_text('\n'.join(lines), encoding='utf-8')
        for label, f in feeds.items():
            for axis, v in enumerate(f):
                assert abs(v/DL-round(v/DL)) < 1e-9, f'{cid} {label} off grid'
                lim = [geo['nx'], geo['ny'], geo['nz']][axis]*DL
                assert PML*DL < v < lim-PML*DL, f'{cid} {label} inside PML'
        groups.append(dict(id=cid, input=f'cases/{cid}/profile.in',
                           geometry=f'geometries/{cid}.h5', material='geometries/'+DB_SOURCE.name,
                           input_sha256=sha(p), geometry_sha256=sha(a.out/'geometries'/f'{cid}.h5'),
                           material_sha256=sha(db),
                           native_shape=[geo['nx'], geo['ny'], geo['nz']], spacing_m=[DL]*3,
                           expected_samples=STEPS, tx_m=list(feeds['tx']), rx_m=list(feeds['rx']),
                           station_profile_x_m=sx, terrain_relief_m=geo['relief'],
                           source_type='impulse', source_frequency_Hz=1., source_amplitude_A=1.))
    manifest = dict(status='PREPARED_APPROVED_ANTENNA_3D_B', approval_basis=APPROVAL,
                    execution_authorized=True, groups=groups, no_retry=True, dt_s=DT,
                    window_s=WINDOW_S, stations_m=STATIONS, baseline_m=BASELINE_M, agl_m=AGL_M,
                    antenna_model='unfocused_centre_fed_bowtie_approximation (stage A design)',
                    parent_geometry_sha256=sha(PARENT_GEO), material_db_sha256=sha(db),
                    cross_track_approximation='terrain extruded constant across y over 6.0 m',
                    generator_sha256=sha(Path(__file__)),
                    limits='3 anchor stations only; not a full 41-trace line. Cross-track invariant '
                           'terrain. Bowtie = stage A unfocused approximation. FP64 GPU. '
                           'PML 60 cells. Window 900 ns.')
    native.save(a.out/'manifest.json', manifest)
    print(json.dumps(dict(groups=len(groups), expected_samples=STEPS,
                          cells=sum(int(np.prod(g['native_shape'])) for g in groups))))


def audit(path, completed=False):
    c = json.loads(Path(path).read_text('utf-8'))
    assert gprMax.__version__ == c['expected_version']
    assert Path(sys.executable).resolve() == Path(c['python']).resolve()
    native_dir = Path(gprMax.__file__).parent
    for category in ['code_identities', 'file_identities', 'binary_identities']:
        for p, digest in c[category].items():
            assert sha(p) == digest, (category, p)
    for p, digest in c['source_identities'].items():
        assert sha(native_dir/p) == digest, p
    rows = []
    for g in c['groups']:
        assert sha(g['input']) == g['input_sha256'], g['id']
        row = dict(id=g['id'], input_sha256=g['input_sha256'])
        if completed:
            raw = Path(g['input']).with_suffix('.h5')
            with h5py.File(raw) as h:
                assert str(h.attrs['gprMax']) == c['expected_version']
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'], g['native_shape'])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'], g['spacing_m'])
                assert abs(float(h.attrs['dt'])/g['dt_s']-1) < 1e-13
                assert h.attrs['Iterations'] == g['expected_samples']
                v = h['rxs/rx1/Ey'][:]
                assert v.dtype == np.float64 and v.shape == (g['expected_samples'],)
                assert np.isfinite(v).all() and np.any(v)
                exc = h['srcs/src1/excitation']
                assert exc.attrs['WaveformType'] == 'impulse'
                np.testing.assert_array_equal(np.flatnonzero(exc['samples'][:]), [0])
                np.testing.assert_allclose(h['srcs/src1'].attrs['Position'], g['tx_m'], atol=1e-12, rtol=0)
            row.update(native_sha256=sha(raw), rxgap_finite=True)
        rows.append(row)
    return dict(status='PASS_NATIVE_IDENTITY_NOT_PHYSICAL_CERTIFICATION', completed=completed,
                contract_sha256=sha(path), version=c['expected_version'], groups=rows)


def freeze(a):
    import psutil
    manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert manifest['approval_basis'] == APPROVAL and manifest['execution_authorized']
    assert len(manifest['groups']) == len(STATIONS)
    prepared = a.out/'prepared'
    groups = []
    for g in manifest['groups']:
        for key in ['input', 'geometry', 'material']:
            src = a.package/g[key]
            dst = prepared/g[key]
            assert sha(src) == g[key+'_sha256']
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
        groups.append(dict(g, input=str((prepared/g['input']).resolve()), dt_s=manifest['dt_s']))
    native_dir = Path(gprMax.__file__).parent
    nvcc, cl = shutil.which('nvcc'), shutil.which('cl')
    assert nvcc and cl, 'CUDA/MSVC environment required'
    hw = supervisor.live_resources(dict(min_available_RAM_GiB=.5, min_free_VRAM_GiB=.5))
    cells = max(int(np.prod(g['native_shape'])) for g in groups)
    estimate = cells*a.vram_bytes_per_cell+2*2**30
    assert hw['free_VRAM_bytes'] > estimate, (hw['free_VRAM_bytes'], estimate)
    assert psutil.virtual_memory().available > 20*2**30
    c = dict(status='FROZEN_APPROVED', expected_version=gprMax.__version__, python=str(Path(sys.executable).resolve()),
             approval_basis=APPROVAL, study_manifest=manifest, groups=groups, max_runs=len(groups), no_retry=True,
             source_identities={p.relative_to(native_dir).as_posix(): sha(p) for p in native_dir.rglob('*')
                                if p.is_file() and p.suffix in ('.py', '.pyd', '.so', '.tmpl')},
             code_identities={str(ROOT/'scripts'/name): sha(ROOT/'scripts'/name) for name in
                              ['line9_v401_version_controls.py', 'hs4_station_grid_controls.py',
                               'hs_capsule_identity.py', 'gprmax_cached_cuda_entry.py',
                               'antenna_stage_a.py', 'line9_antenna_3d_b.py']},
             file_identities={str(p.resolve()): sha(p) for p in a.out.rglob('*') if p.is_file()},
             binary_identities={str(Path(p).resolve()): sha(p) for p in [sys.executable, nvcc, cl]},
             min_available_RAM_GiB=20, min_free_VRAM_GiB=estimate/2**30,
             max_owned_RSS_GiB=24, min_system_available_during_run_GiB=.5,
             max_group_wall_s=10800, max_batch_wall_s=43200,
             gpu_lock='artifacts/local_checks/hs4_gpu_exclusive.lock',
             cancel_file=str((a.out/'USER_STOP').resolve()), lease_file=str((a.out/'session_heartbeat').resolve()),
             max_lease_age_s=600, solver_entrypoint='gprmax_cached_cuda_entry.py', additional_solver_args=[],
             hardware_at_freeze=hw, conservative_device_estimate_bytes=estimate,
             vram_estimate_basis=dict(bytes_per_cell=a.vram_bytes_per_cell,
                                      source='stage A measured FP64+HORIPML usage x1.15 margin'),
             limits=manifest['limits'])
    (a.out/'session_heartbeat').write_text('Active authorized version verification\n')
    native.save(p := a.out/'execution_contract.json', c)
    native.save(a.out/'preflight_verification.json', audit(p))
    print(json.dumps(dict(contract_sha256=sha(p), groups=len(groups), estimate_GiB=estimate/2**30)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--package', type=Path)
    p.add_argument('--vram-bytes-per-cell', type=float, default=200.0)
    a = p.parse_args()
    a.out = a.out.resolve()
    if a.action == 'prepare':
        prepare(a)
    elif a.action == 'freeze':
        a.package = a.package.resolve()
        freeze(a)
    elif a.action == 'run':
        supervisor.audit = audit
        supervisor.run(a.out/'execution_contract.json')
    else:
        print(json.dumps(audit(a.out/'execution_contract.json', True)))
