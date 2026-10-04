"""Prepare/run/audit bounded cover-material (eps_r x sigma) scan at 8 m altitude.

Motivation (user 2026-10-05: "做，以及孔径上的也跑仿真看看情况" — approving the
proposed material sensitivity scan): the measured Line9(36) UAV data at 9 m
altitude shows late reflections only ~27 dB below the early peak while the
archived 8 m dispersive wet-clay model puts the interface echo at -65.9 dB.
This scan measures how echo/direct ratio, decay slope and relief-shape fidelity
vary with cover eps_r in {9, 12, 18.017} and sigma in {0.0005, 0.001, 0.003}
(non-dispersive cover), plus one dispersive control (eps 18.017 + Debye +
sigma 0.003) anchoring the grid to the archived ladder. Labelled site-adaptation
diagnostic, not blind-test performance.

Same 36 m domain, grid, PML, boxes, impulse source and 13-station layout
(src 14.6..20.6 m step 0.5 m, rx +1.3 m) as the frozen joint-grid study, with
the antenna pair lowered z 27 -> 20 m (8 m altitude, matching the measured
~9 m). 10 materials x 2 roles x 1 group of 13 traces = 20 groups, 260 traces.
One fresh attempt; no retry.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np
import psutil
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'
HEIGHTS_CONTRACT = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c/execution_contract.json'
WHEEL = ROOT/'artifacts/local_checks/wheelhouse/gprmax-4.0.0-cp312-cp312-win_amd64.whl'
ANT_Z = 20.0
STATION_INDICES = list(range(1, 122, 10))  # 13 stations, src x0=14.6 step 0.5 via src_steps
MATERIALS = [(f'eps{eps:g}_s{sig:g}', eps, sig, False)
             for eps in (9.0, 12.0, 18.017) for sig in (0.0005, 0.001, 0.003)]
MATERIALS.append(('eps18.017_s0.003_debye', 18.017, 0.003, True))
ROLES = ['rough', 'halfspace']


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def cover_line(eps, sig):
    return f'#material: {eps:g} {sig:g} 1 0 cover'


def variant_lines(role, eps, sig, dispersive):
    original = REFERENCE/f'centre_{role}'/'profile.in'
    lines = []
    for line in original.read_text('utf-8').splitlines():
        if line.startswith('#add_dispersion_debye:'):
            if dispersive:
                lines.append(line)
            continue
        if line.startswith('#title:'):
            line = (f'#title: HS4T2D cover-material scan v0.1 (cover eps_r={eps:g} sigma={sig:g}'
                    f'{" +Debye" if dispersive else ""}; 8m altitude; companion of joint-grid study)')
        elif line.startswith('#material:') and line.split()[-1] == 'cover':
            line = cover_line(eps, sig)
        elif line.startswith('#hertzian_dipole:'):
            line = '#hertzian_dipole: y 14.6 0.025 20 impulse'
        elif line.startswith('#rx:'):
            line = '#rx: 15.9 0.025 20 t01 Ey'
        lines.append(line)
    lines += ['#src_steps: 0.5 0 0', '#rx_steps: 0.5 0 0']
    return original, lines


def prepare(folder):
    if folder.exists():
        raise ValueError('new capsule required')
    old = json.loads(HEIGHTS_CONTRACT.read_text('utf-8'))
    with h5py.File(REFERENCE/'centre_rough/profile.h5') as h:
        n, dt = len(h['rxs/rx1/Ey']), float(h.attrs['dt'])
    folder.mkdir(parents=True)
    groups = []
    for tag, eps, sig, dispersive in MATERIALS:
        for role in ROLES:
            name = f'{tag}_{role}'
            original, lines = variant_lines(role, eps, sig, dispersive)
            directory = folder/name
            directory.mkdir()
            path = directory/'profile.in'
            path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
            refgeom = REFERENCE/f'centre_{role}'/'hs4t2d_geom.vtkhdf'
            groups.append({'id': name, 'input': str(path.resolve()), 'input_sha256': sha256(path),
                'material_tag': tag, 'cover_eps_r': eps, 'cover_sigma': sig, 'debye': dispersive,
                'role': role, 'halfspace': role == 'halfspace',
                'station_indices': STATION_INDICES, 'traces': len(STATION_INDICES),
                'station_x0_m': 14.6, 'antenna_z_m': ANT_Z, 'spacing_m': [.025, .05, .025],
                'reference_input': str(original.resolve()), 'reference_input_sha256': sha256(original),
                'reference_geometry': str(refgeom.resolve()), 'reference_geometry_sha256': sha256(refgeom)})
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    identities, py_checked = {}, 0
    for name, digest in old['source_identities'].items():
        actual = sha256(package/name)
        if name.endswith('.py'):
            if actual != digest:
                raise ValueError('old V4 Python source differs: '+name)
            py_checked += 1
        identities[name] = actual
    c = {'status': 'FROZEN_APPROVED',
         'approval_basis': 'User 2026-10-05: 做，以及孔径上的也跑仿真看看情况 (approving proposed cover-material sigma/eps scan against measured Line9 decay); bounded scan frozen before execution.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': identities,
         'source_identity_note': f'All .py files byte-identical to the frozen 8m wavefield contract ({py_checked} checked). The .pyd files are the fresh local build recorded there; wheel sha256 recorded.',
         'wheel': str(WHEEL.resolve()), 'wheel_sha256': sha256(WHEEL),
         'code_identities': {str(ROOT/name): sha256(ROOT/name) for name in
                             ['scripts/hs4_material_scan_v0_1.py', 'scripts/run_hs4_material_scan_v0_1.cmd',
                              'scripts/gprmax_cached_cuda_entry.py']},
         'max_runs': len(groups), 'max_traces': sum(g['traces'] for g in groups),
         'groups': groups, 'dt_s': dt, 'iterations': n,
         'min_available_RAM_GiB': 1.5, 'min_free_VRAM_GiB': 3.5,
         'max_group_wall_s': 1200, 'max_wall_s': 24000,
         'max_owned_RSS_GiB': 1.75, 'min_live_system_available_MiB': 100,
         'no_retry': True, 'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'processing': 'Matched saved-source-normalized 501 frequencies 20-170MHz; physical 200ns taper/Hann/8x; rough-halfspace complex subtraction; fixed relief gates. No AGC/SVD/per-trace scaling.',
         'scope': 'TM invariant Y; 8 m antenna height; isolates cover eps_r/sigma (non-dispersive) plus one dispersive control. Site-adaptation diagnostic vs measured Line9 decay, NOT blind-test performance. No 3D/field claim.'}
    save(folder/'execution_contract.json', c)
    check(folder/'execution_contract.json', False)
    print(f'Frozen {len(groups)} groups / {c["max_traces"]} traces; dt={dt:.6e} s; iterations={n}')


def check(path, completed):
    c = json.loads(path.read_text('utf-8'))
    rows = []
    skip = ('#title:', '#hertzian_dipole:', '#rx:', '#src_steps:', '#rx_steps:', '#add_dispersion_debye:')

    def strip(lines):
        return [s for s in lines if not s.startswith(skip)
                and not (s.startswith('#material:') and s.split()[-1] == 'cover')]

    for g in c['groups']:
        p = Path(g['input'])
        if sha256(p) != g['input_sha256'] or sha256(g['reference_input']) != g['reference_input_sha256'] \
                or sha256(g['reference_geometry']) != g['reference_geometry_sha256']:
            raise ValueError('input/reference identity differs')
        lines = p.read_text('utf-8').splitlines()
        ref = Path(g['reference_input']).read_text('utf-8').splitlines()
        if strip(lines) != strip(ref):
            raise ValueError('physical invariant differs')
        cover = [s for s in lines if s.startswith('#material:') and s.split()[-1] == 'cover']
        if cover != [cover_line(g['cover_eps_r'], g['cover_sigma'])]:
            raise ValueError('cover material line differs')
        debye = [s for s in lines if s.startswith('#add_dispersion_debye:')]
        if debye != (['#add_dispersion_debye: 1 7.878 6.4567e-09 cover'] if g['debye'] else []):
            raise ValueError('dispersion line differs')
        tx = next(s.split() for s in lines if s.startswith('#hertzian_dipole:'))
        rx = next(s.split() for s in lines if s.startswith('#rx:'))
        if not np.allclose(list(map(float, tx[2:5])), [14.6, .025, ANT_Z], rtol=0, atol=1e-12) \
                or not np.allclose(list(map(float, rx[1:4])), [15.9, .025, ANT_Z], rtol=0, atol=1e-12):
            raise ValueError('source/receiver position differs')
        steps = [s for s in lines if s.startswith(('#src_steps:', '#rx_steps:'))]
        if steps != ['#src_steps: 0.5 0 0', '#rx_steps: 0.5 0 0']:
            raise ValueError('trace movement differs')
        if not completed:
            continue
        audit = json.loads((p.parent/'audit.json').read_text('utf-8'))
        if len(audit) != g['traces']:
            raise ValueError('incomplete outputs')
        for k, row in enumerate(audit, 1):
            raw = p.parent/row['file']
            if sha256(raw) != row['sha256'] or row['station_index'] != g['station_indices'][k-1]:
                raise ValueError('raw identity or station index differs')
            with h5py.File(raw) as h:
                v = h['rxs/rx1/Ey'][:]
                if str(h.attrs['gprMax']) != '4.0.0' or float(h.attrs['dt']) != c['dt_s'] \
                        or int(h.attrs['Iterations']) != c['iterations'] \
                        or not np.array_equal(h.attrs['dx_dy_dz'], g['spacing_m']) \
                        or v.dtype != np.float64 or not np.isfinite(v).all() or len(v) != c['iterations']:
                    raise ValueError('raw version/grid/time/validity differs')
                sx = 14.6+.05*(g['station_indices'][k-1]-1)
                for name, px in [('srcs/src1', sx), ('rxs/rx1', sx+1.3)]:
                    if not np.array_equal(h[name].attrs['GridPosition'], [round(px/.025), 0, round(ANT_Z/.025)]):
                        raise ValueError('actual station grid differs')
        geom = p.parent/'hs4t2d_geom1.vtkhdf'
        with h5py.File(geom) as h, h5py.File(g['reference_geometry']) as r:
            if not np.array_equal(h['VTKHDF/CellData/Material'][:], r['VTKHDF/CellData/Material'][:]):
                raise ValueError('actual material geometry differs')
        rows.append({'id': g['id'], 'material_tag': g['material_tag'], 'role': g['role'],
                     'traces': g['traces'], 'geometry_sha256': sha256(geom)})
    if completed:
        events = [json.loads(s) for s in (path.parent/'execution.jsonl').read_text('utf-8').splitlines()]
        if events[0]['contract_sha256'] != sha256(path):
            raise ValueError('execution contract identity differs')
        done = {e['group']: e for e in events if e.get('status') == 'COMPLETED' and 'group' in e}
        if any(e.get('status') == 'FAILED' for e in events if 'group' in e):
            raise ValueError('a solver group failed')
        for g in c['groups']:
            if g['id'] not in done or done[g['id']].get('traces') != g['traces']:
                raise ValueError('incomplete or changed execution')
        if sum(e['traces'] for e in done.values()) != c['max_traces']:
            raise ValueError('total trace count differs')
    report = {'status': 'PASS', 'completed': completed, 'contract_sha256': sha256(path),
              'code_sha256': sha256(__file__), 'groups': rows,
              'scope': 'input invariants, cover-material substitution, antenna height z=20, station geometry, raw validity and realised material grids; not full spatial convergence'}
    save(path.parent/('completed_verification.json' if completed else 'preflight_verification.json'), report)
    return report


def run(path):
    import msvcrt
    c = json.loads(path.read_text('utf-8'))
    if c['status'] != 'FROZEN_APPROVED' or len(c['groups']) != c['max_runs'] \
            or Path(sys.executable).resolve() != Path(c['python']).resolve():
        raise ValueError('approved bounded V4 study required')
    for name, digest in c['code_identities'].items():
        if sha256(name) != digest:
            raise ValueError('frozen code differs')
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    for name, digest in c['source_identities'].items():
        if sha256(package/name) != digest:
            raise ValueError('V4 runtime differs')
    check(path, False)
    record = path.parent/'execution.jsonl'
    if record.exists() or any(list(Path(g['input']).parent.glob('profile*.h5')) for g in c['groups']):
        raise ValueError('attempt consumed; new contract required for new work')
    lock = ROOT/c['gpu_lock']
    lock.parent.mkdir(parents=True, exist_ok=True)
    with open(lock, 'a+b') as held:
        held.seek(0)
        msvcrt.locking(held.fileno(), msvcrt.LK_NBLCK, 1)

        def event(row):
            with record.open('a', encoding='utf-8') as f:
                f.write(json.dumps(row, allow_nan=False)+'\n')
        digest = sha256(path)
        start = time.monotonic()
        event({'status': 'STARTED', 'contract_sha256': digest, 'unix_s': time.time(), 'owned_runner_pid': os.getpid()})
        try:
            for g in c['groups']:
                p = Path(g['input'])
                free = psutil.virtual_memory().available
                vram = int(subprocess.check_output(['nvidia-smi', '--query-gpu=memory.free', '--format=csv,noheader,nounits'],
                                                   text=True).splitlines()[0])*2**20
                if free < c['min_available_RAM_GiB']*2**30 or vram < c['min_free_VRAM_GiB']*2**30 \
                        or sha256(path) != digest or time.monotonic()-start > c['max_wall_s']:
                    raise RuntimeError('live resource minimum, wall budget or contract identity differs')
                command = [sys.executable, str(ROOT/'scripts/gprmax_cached_cuda_entry.py'), str(p),
                           '-n', str(g['traces']), '--geometry-fixed', '-gpu', '0', '-gpu_precision', 'double',
                           '--hide-progress-bars']
                env = os.environ.copy()
                env['HS4_CUDA_CACHE_LOG'] = str(p.parent/'cuda_cache.jsonl')
                event({'status': 'STARTED', 'group': g['id'], 'command': command,
                       'RAM_available_bytes': free, 'VRAM_free_bytes': vram})
                with (p.parent/'stdout.log').open('xb') as out, (p.parent/'stderr.log').open('xb') as err:
                    process = subprocess.Popen(command, cwd=p.parent, stdout=out, stderr=err, env=env)
                    started = time.monotonic()
                    peak = 0
                    while process.poll() is None:
                        try:
                            owned = psutil.Process(process.pid)
                            rss = sum(q.memory_info().rss for q in [owned]+owned.children(recursive=True) if q.is_running())
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            rss = 0
                        peak = max(peak, rss)
                        if time.monotonic()-started > c['max_group_wall_s'] or rss > c['max_owned_RSS_GiB']*2**30 \
                                or psutil.virtual_memory().available < c['min_live_system_available_MiB']*2**20:
                            process.terminate()
                            process.wait(timeout=20)
                            raise RuntimeError('owned solver exceeded wall/RSS/system memory guard')
                        time.sleep(.5)
                    if process.returncode:
                        raise RuntimeError('solver failed; preserve logs, no retry')
                rows = []
                for k in range(1, g['traces']+1):
                    raw = p.parent/f'profile{k}.h5'
                    with h5py.File(raw) as h:
                        rows.append({'file': raw.name, 'sha256': sha256(raw), 'trace': k,
                                     'station_index': g['station_indices'][k-1],
                                     'dt_s': float(h.attrs['dt']), 'samples': len(h['rxs/rx1/Ey'][:]),
                                     'source_m': h['srcs/src1'].attrs['Position'].tolist(),
                                     'receiver_m': h['rxs/rx1'].attrs['Position'].tolist()})
                (p.parent/'audit.json').write_text(json.dumps(rows, indent=2)+'\n', encoding='utf-8')
                event({'status': 'COMPLETED', 'group': g['id'], 'traces': len(rows),
                       'elapsed_s': time.monotonic()-started, 'peak_owned_RSS_bytes': peak})
                print('Completed', g['id'], flush=True)
            result = check(path, True)
            event({'status': 'COMPLETED', 'traces': c['max_traces'],
                   'elapsed_s': time.monotonic()-start,
                   'verification_sha256': sha256(path.parent/'completed_verification.json')})
        except Exception as exc:
            event({'status': 'FAILED', 'error': str(exc)})
            raise
        finally:
            held.seek(0)
            msvcrt.locking(held.fileno(), msvcrt.LK_UNLCK, 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'check'])
    parser.add_argument('--out', type=Path)
    parser.add_argument('--contract', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args.out)
    elif args.action == 'run':
        if not args.execute:
            raise ValueError('--execute required')
        run(args.contract)
    else:
        print(json.dumps(check(args.contract, True), ensure_ascii=False)[:500])
