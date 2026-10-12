"""Sigma sensitivity batch on the B-r3 3D bowtie station x190 (user-approved
2026-10-12: "按你的建议开始做吧", approving the proposed mudstone sigma
{0.002, 0.006} x cover sigma {0.0005, 0.001} sweep at a single station).

Goal: quantify the bottom-sand / direct-wave amplitude ratio vs cover and
mudstone conductivity, against the real-data reference of about -40 dB (raw,
user-provided 2026-10-12). Everything except the two sigma values is identical
to the B-r3 x190 case (contract 9d176513): same cropped terrain, bowtie pair,
1.3 m baseline, 5 m AGL, official impulse 1 A excitation, 900 ns window,
FP64 GPU, PML 60, dispersive averaging on, same Debye poles.

The base point (mud 0.006, cover 0.001) already exists as B-r3
ant3d_x19000 (raw sha256 91b361c2...) and is reused, so this batch solves
only the 3 new combinations. Per-variant material database JSON files sit
beside the per-variant geometry H5 (gprMax resolves geometry databases next
to the .h5); voxel content is identical across variants, only the
MaterialDatabase attribute differs.
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
from antenna_stage_a import DT, DL, PML, bowtie_lines, fmt
from line9_antenna_3d_b import (PARENT_GEO, DB_SOURCE, WINDOW_S, STEPS,
                                X_HALF_M, Y_SPAN_M, GROUND_DEPTH_M, AIR_TOP_M,
                                BASELINE_M, AGL_M, ELEV_OFFSET_M,
                                PROFILE_OFFSET_M, case_text, surface_elevations)
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]

APPROVAL = ('2026-10-12 用户批准 σ 敏感性批：泥岩 σ∈{0.002,0.006} × 覆盖层 σ∈{0.0005,0.001}，'
            'x190 单站；基准组合 (0.006,0.001) 复用 B-r3 ant3d_x19000 不重跑；'
            '激励保持官方 impulse 1 A，其余物理/采集配置与 B-r3 完全一致。'
            '目标：量化底砂/直达幅度比对 σ 的响应，对照实测约 −40 dB。')

# (variant tag, mudstone sigma, cover sigma) -- base (0.006, 0.001) excluded, it is B-r3 x190
VARIANTS = [
    ('m002_c001',  0.002, 0.001),
    ('m002_c0005', 0.002, 0.0005),
    ('m006_c0005', 0.006, 0.0005),
]
STATION = 190.0
BASE_REUSE = dict(case='ant3d_x19000', batch='line9_antenna_b_20261012_r3',
                  raw_sha256='91b361c2b20275e593afd6647463d2bc594298f9e20a8ebad08569cc8322d697')


def make_variant_db(path, dbid, mud_sigma, cover_sigma):
    doc = json.loads(DB_SOURCE.read_text('utf-8'))
    doc['database']['id'] = dbid
    doc['database']['name'] = f'Line9 sigma sensitivity variant {dbid}'
    doc['database']['description'] = (
        f'Sigma sweep variant of line9_materials_v0_2: mudstone sigma={mud_sigma}, '
        f'cover sigma={cover_sigma}; all other parameters identical. Sensitivity '
        'experiment values, not site calibration. Parent DB sha256='+sha(DB_SOURCE))
    doc['materials']['material_001_cover']['base']['electric_conductivity_s_per_m'] = cover_sigma
    doc['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m'] = mud_sigma
    for key in ['material_001_cover', 'material_002_mudstone']:
        doc['materials'][key]['metadata']['parameter_status'] = 'SIGMA_SENSITIVITY_VARIANT_NOT_SITE_CALIBRATION'
    path.write_text(json.dumps(doc, indent=1, ensure_ascii=False), 'utf-8')


def make_crop(path, sx, dbid):
    with h5py.File(PARENT_GEO) as h:
        d2 = h['data'][:, :, 0]
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
            out.attrs['MaterialDatabase'] = dbid
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


def prepare(a):
    assert not a.out.exists()
    (a.out/'cases').mkdir(parents=True)
    (a.out/'geometries').mkdir()
    groups = []
    for tag, mud_s, cov_s in VARIANTS:
        dbid = f'line9_sig_{tag}'
        cid = f'sig_x19000_{tag}'
        make_variant_db(a.out/'geometries'/f'{dbid}.json', dbid, mud_s, cov_s)
        geo = make_crop(a.out/'geometries'/f'{cid}.h5', STATION, dbid)
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
        # geometry read first: it would overwrite the PEC antenna cells if placed later
        lines.append(f'#geometry_objects_read: 0 0 0 ../../geometries/{cid}.h5 {dbid} n')
        lines += bowtie_lines(feeds['tx'], 'y', 'impulse', window_s=WINDOW_S) + extra
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
                           geometry=f'geometries/{cid}.h5', material=f'geometries/{dbid}.json',
                           input_sha256=sha(p), geometry_sha256=sha(a.out/'geometries'/f'{cid}.h5'),
                           material_sha256=sha(a.out/'geometries'/f'{dbid}.json'),
                           native_shape=[geo['nx'], geo['ny'], geo['nz']], spacing_m=[DL]*3,
                           expected_samples=STEPS, tx_m=list(feeds['tx']), rx_m=list(feeds['rx']),
                           station_profile_x_m=STATION, terrain_relief_m=geo['relief'],
                           mudstone_sigma=mud_s, cover_sigma=cov_s,
                           source_type='impulse', source_frequency_Hz=1., source_amplitude_A=1.))
    manifest = dict(status='PREPARED_APPROVED_SIGMA_SWEEP', approval_basis=APPROVAL,
                    execution_authorized=True, groups=groups, no_retry=True, dt_s=DT,
                    window_s=WINDOW_S, stations_m=[STATION], baseline_m=BASELINE_M, agl_m=AGL_M,
                    antenna_model='unfocused_centre_fed_bowtie_approximation (stage A design)',
                    base_point_reuse=BASE_REUSE,
                    parent_geometry_sha256=sha(PARENT_GEO),
                    parent_material_db_sha256=sha(DB_SOURCE),
                    cross_track_approximation='terrain extruded constant across y over 6.0 m',
                    generator_sha256=sha(Path(__file__)),
                    limits='Single station x190, 3 sigma variants; base point reused from B-r3. '
                           'Cross-track invariant terrain. Bowtie = stage A unfocused '
                           'approximation. FP64 GPU. PML 60 cells. Window 900 ns. '
                           'Sensitivity experiment, NOT site-calibrated materials.')
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
                # the solved model must carry this variant's material database id
                prov = {k: dict(h['material_database_provenance'][k].attrs)
                        for k in h['material_database_provenance']}
                dbids = {str(v['DatabaseID']) for v in prov.values()}
                expected_dbid = 'line9_sig_'+g['id'].split('sig_x19000_')[1]
                assert dbids == {expected_dbid}, (g['id'], dbids)
            row.update(native_sha256=sha(raw), rxgap_finite=True)
        rows.append(row)
    return dict(status='PASS_NATIVE_IDENTITY_NOT_PHYSICAL_CERTIFICATION', completed=completed,
                contract_sha256=sha(path), version=c['expected_version'], groups=rows)


def freeze(a):
    import psutil
    manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert manifest['approval_basis'] == APPROVAL and manifest['execution_authorized']
    assert len(manifest['groups']) == len(VARIANTS)
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
                               'antenna_stage_a.py', 'line9_antenna_3d_b.py',
                               'line9_antenna_3d_b_sigma.py']},
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
                                      source='B-r3 measured ~11.5 GiB for the same grid; 150 B/cell conservative'),
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
    p.add_argument('--vram-bytes-per-cell', type=float, default=150.0)
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
