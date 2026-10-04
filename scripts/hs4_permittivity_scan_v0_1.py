"""Prepare/run/audit bounded cover-permittivity scan on the 13-station 15 m HS4T2D grid.

Motivation (user 2026-10-04): the 8 m dense-wavefield unit showed the
cover->air critical angle (13.6 deg at eps_r=18.017) strips the high-angle
interface-elevation spectrum before it reaches the airborne antenna, while the
Kirchhoff diagnostic showed the large-scale (>~1.6 m) relief is still
recoverable. The lateral resolution bound lambda_cover/(2 sin theta_c) = c/(2f)
is permittivity-independent, but the transmitted in-band power fraction and
echo amplitude grow as eps_r drops. This scan measures how migration-based
relief recovery changes with cover eps_r in {6, 9, 18.017-non-dispersive}.

Same domain, grid, PML, boxes, impulse source, 13 stations (src 14.6..20.6 m
step 0.5 m, rx +1.3 m, z=27 m) as the frozen 15 m joint-grid study; only the
cover material eps_r changes and the Debye dispersion line is removed (the
eps18nd group also isolates the dispersion contribution against the archived
dispersive baseline). 18 groups = 3 eps x 2 roles x 3 segments, 78 traces.
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
EPSILONS = [('eps6', '6'), ('eps9', '9'), ('eps18nd', '18.017')]
ROLES = ['rough', 'halfspace']
SEGMENTS = [('centre', [61]), ('left', [1, 11, 21, 31, 41, 51]), ('right', [71, 81, 91, 101, 111, 121])]


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def variant_lines(role, eps, indices):
    original = REFERENCE/f'centre_{role}'/'profile.in'
    x0 = 14.6+.05*(indices[0]-1)
    lines = []
    for line in original.read_text('utf-8').splitlines():
        if line.startswith('#add_dispersion_debye:'):
            continue
        if line.startswith('#title:'):
            line = f'#title: HS4T2D cover-permittivity scan v0.1 (eps_r={eps}, non-dispersive cover; companion of joint-grid 15 m study)'
        elif line.startswith('#material:') and line.split()[-1] == 'cover':
            line = f'#material: {eps} 0.003 1 0 cover'
        elif line.startswith('#hertzian_dipole:'):
            line = f'#hertzian_dipole: y {x0:.12g} 0.025 27 impulse'
        elif line.startswith('#rx:'):
            line = f'#rx: {x0+1.3:.12g} 0.025 27 t01 Ey'
        lines.append(line)
    if len(indices) > 1:
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
    for tag, eps in EPSILONS:
        for role in ROLES:
            for segment, indices in SEGMENTS:
                name = f'{tag}_{segment}_{role}'
                original, lines = variant_lines(role, eps, indices)
                directory = folder/name
                directory.mkdir()
                path = directory/'profile.in'
                path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
                refgeom = REFERENCE/f'centre_{role}'/'hs4t2d_geom.vtkhdf'
                groups.append({'id': name, 'input': str(path.resolve()), 'input_sha256': sha256(path),
                    'eps_tag': tag, 'cover_eps_r': eps, 'role': role, 'halfspace': role == 'halfspace',
                    'segment': segment, 'station_indices': indices, 'traces': len(indices),
                    'station_x0_m': 14.6+.05*(indices[0]-1), 'spacing_m': [.025, .05, .025],
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
         'approval_basis': 'User2026-10-04: 大尺度形状仍可成像、临界角随覆盖层介电常数变 能否从这些角度上尝试做好我们的模型; bounded cover-permittivity scan frozen before execution.',
         'python': str(Path(sys.executable).resolve()), 'source_identities': identities,
         'source_identity_note': f'All .py files byte-identical to the frozen 8m wavefield contract ({py_checked} checked). The .pyd files are the fresh local build recorded there; wheel sha256 recorded.',
         'wheel': str(WHEEL.resolve()), 'wheel_sha256': sha256(WHEEL),
         'code_identities': {str(ROOT/name): sha256(ROOT/name) for name in
                             ['scripts/hs4_permittivity_scan_v0_1.py', 'scripts/run_hs4_permittivity_scan_v0_1.cmd',
                              'scripts/gprmax_cached_cuda_entry.py']},
         'max_runs': len(groups), 'max_traces': sum(g['traces'] for g in groups),
         'groups': groups, 'dt_s': dt, 'iterations': n,
         'min_available_RAM_GiB': 1.5, 'min_free_VRAM_GiB': 3.5,
         'max_group_wall_s': 3600, 'max_wall_s': 20000,
         'max_owned_RSS_GiB': 1.75, 'min_live_system_available_MiB': 100,
         'no_retry': True, 'gpu_lock': 'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'processing': 'Matched saved-source-normalized 501 frequencies 20-170MHz; physical 200ns taper/Hann/8x; rough-halfspace complex subtraction; Kirchhoff migration per eps with v=c/sqrt(eps_r). No AGC/SVD/per-trace scaling.',
         'scope': 'TM invariant Y; 15 m antenna height; isolates cover eps_r (non-dispersive) effect on airborne recoverability of the buried relief. eps18nd additionally separates the Debye dispersion contribution against the archived dispersive baseline. No 3D/field claim.'}
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
        if cover != [f'#material: {g["cover_eps_r"]} 0.003 1 0 cover'] or any(s.startswith('#add_dispersion_debye:') for s in lines):
            raise ValueError('cover material or dispersion line differs')
        x0 = g['station_x0_m']
        tx = next(s.split() for s in lines if s.startswith('#hertzian_dipole:'))
        rx = next(s.split() for s in lines if s.startswith('#rx:'))
        if not np.allclose(list(map(float, tx[2:5])), [x0, .025, 27], rtol=0, atol=1e-12) \
                or not np.allclose(list(map(float, rx[1:4])), [x0+1.3, .025, 27], rtol=0, atol=1e-12):
            raise ValueError('source/receiver position differs')
        steps = [s for s in lines if s.startswith(('#src_steps:', '#rx_steps:'))]
        if steps != (['#src_steps: 0.5 0 0', '#rx_steps: 0.5 0 0'] if g['traces'] > 1 else []):
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
                    if not np.array_equal(h[name].attrs['GridPosition'], [round(px/.025), 0, 1080]):
                        raise ValueError('actual station grid differs')
        geom = p.parent/('hs4t2d_geom1.vtkhdf' if g['traces'] > 1 else 'hs4t2d_geom.vtkhdf')
        with h5py.File(geom) as h, h5py.File(g['reference_geometry']) as r:
            if not np.array_equal(h['VTKHDF/CellData/Material'][:], r['VTKHDF/CellData/Material'][:]):
                raise ValueError('actual material geometry differs')
        rows.append({'id': g['id'], 'eps_tag': g['eps_tag'], 'role': g['role'], 'traces': g['traces'],
                     'geometry_sha256': sha256(geom)})
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
              'scope': 'input invariants, cover-material substitution, station geometry, raw validity and realised material grids; not full spatial convergence'}
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
                    raw = p.parent/('profile.h5' if g['traces'] == 1 else f'profile{k}.h5')
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
