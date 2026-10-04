"""Prepare/run/audit bounded 15 m per-permittivity wavefield snapshots (V4, no solver edits).

Companion to the cover-permittivity scan capsule: same HS4T2D model, 15 m
antenna height (z=27), centre station, ricker 95 MHz excitation and the same
snapshot ROI/hop (17 iterations ~ 1.0024 ns, 599 frames) as the verified 8 m
wavefield capsule — but only for the three NON-dispersive scan variants
(eps6, eps9, eps18nd) x (rough, halfspace) = 6 groups. The dispersive 15 m
case was already consumed by the heights-continuation capsule on the previous
machine and is not re-run; its archived receiver H5 plus the 8 m dispersive
GIF remain the dispersive references. Probe-based snapshot/native closure and
material-geometry equality with the tracked reference vtkhdf are enforced.
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
CELLS = [(i, 6) for i in (10,20,30,40,50)] + [(i,26) for i in (10,30,50)] + [(35,39),(35,126),(20,60),(40,60)]
FIELDS = ('Ex','Ey','Ez','Hx','Hy','Hz')
HOP = 17
EPSILONS = [('eps6', '6'), ('eps9', '9'), ('eps18nd', '18.017')]
ROLES = ['rough', 'halfspace']


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


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
            name = f'w15_{tag}_{role}'
            original = REFERENCE/f'centre_{role}'/'profile.in'
            lines = []
            for line in original.read_text('utf-8').splitlines():
                if line.startswith('#add_dispersion_debye:'):
                    continue
                if line.startswith('#title:'):
                    line = f'#title: HS4T2D 15m per-eps wavefield v0.1 (cover eps_r={eps}, non-dispersive; ricker95; centre station)'
                elif line.startswith('#waveform:'):
                    line = '#waveform: ricker 1 95e6 ricker95'
                elif line.startswith('#material:') and line.split()[-1] == 'cover':
                    line = f'#material: {eps} 0.003 1 0 cover'
                elif line.startswith('#hertzian_dipole:'):
                    line = '#hertzian_dipole: y 17.6 0.025 27 ricker95'
                lines.append(line)
            probes = []
            for k,(ix,iz) in enumerate(CELLS):
                x,z = 13.5+.15*ix+.075, 8+.15*iz+.075
                for suffix,px,pz,components in [('c',x,z,'Ey Hx Hz'),('z',x,z-.025,'Hx'),('x',x-.025,z,'Hz')]:
                    lines.append(f'#rx: {px:.12g} 0.025 {pz:.12g} p{k:02d}{suffix} {components}')
                probes.append({'id':f'p{k:02d}', 'snapshot_ix':ix, 'snapshot_iz':iz,'centre_xz_m':[x,z]})
            lines.extend(f'#snapshot: 13.5 0 8 22.5 0.05 28.1 0.15 0.05 0.15 {j} snap{j:05d}.h5' for j in range(0,n,HOP))
            directory = folder/name
            directory.mkdir()
            path = directory/'profile.in'
            path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            refgeom = REFERENCE/f'centre_{role}'/'hs4t2d_geom.vtkhdf'
            groups.append({'id':name,'input':str(path.resolve()),'input_sha256':sha256(path),
                'eps_tag':tag,'cover_eps_r':eps,'role':role,'height_m':15,'snapshot':True,'probes':probes,
                'reference_input':str(original.resolve()),'reference_input_sha256':sha256(original),
                'reference_geometry':str(refgeom.resolve()),'reference_geometry_sha256':sha256(refgeom)})
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    identities, py_checked = {}, 0
    for name,digest in old['source_identities'].items():
        actual = sha256(package/name)
        if name.endswith('.py'):
            if actual != digest:
                raise ValueError('old V4 Python source differs: '+name)
            py_checked += 1
        identities[name] = actual
    c = {'status':'FROZEN_APPROVED',
         'approval_basis':'User2026-10-04: 我需要看的是各个bscan、波场快照; bounded per-eps wavefield study frozen before execution. Dispersive 15m not re-run (consumed by heights-continuation capsule).',
         'python':str(Path(sys.executable).resolve()),'source_identities':identities,
         'source_identity_note':'All .py files byte-identical to the frozen 8m wavefield contract (%d checked). The .pyd files are the fresh local build recorded there; wheel sha256 recorded.'%py_checked,
         'wheel':str(WHEEL.resolve()),'wheel_sha256':sha256(WHEEL),
         'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in ['scripts/hs4_permittivity_wavefield_v0_1.py','scripts/run_hs4_permittivity_wavefield_v0_1.cmd','scripts/gprmax_cached_cuda_entry.py']},
         'max_runs':6,'groups':groups,'dt_s':dt,'iterations':n,'snapshot_hop_iterations':HOP,
         'snapshot_interval_s':HOP*dt,'snapshot_iterations':list(range(0,n,HOP)),
         'snapshot_shape_xyz':[60,1,134],'snapshot_extent_m':[13.5,0,8,22.5,.05,28.1],'snapshot_spacing_m':[.15,.05,.15],
         'snapshot_six_field_payload_GiB':60*134*6*8*len(range(0,n,HOP))/2**30,
         'min_available_RAM_GiB':1.25,'min_free_VRAM_GiB':3.5,
         'max_group_wall_s':1800,'max_wall_s':10000,'max_owned_RSS_GiB':1.75,'min_live_system_available_MiB':100,
         'no_retry':True,'gpu_lock':'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'processing':'Visualisation-only wavefields (ricker 95MHz); no SFCW chain claim. Probe/snapshot closure and geometry equality enforced as in the 8m capsule.',
         'scope':'TM invariant Y; 15 m antenna height; non-dispersive cover eps_r in {6,9,18.017}. Passive no-snapshot receiver-identity control not repeated (already proven on identical solver/build in the 8m capsule). No 3D/field claim.'}
    save(folder/'execution_contract.json',c)
    check(folder/'execution_contract.json',False)
    print('Frozen6runs;',c['snapshot_six_field_payload_GiB'],'GiB snapshot payload per trace;',len(c['snapshot_iterations']),'frames')


def receivers(h):
    return {str(rx.attrs['Name']):rx for rx in h['rxs'].values()}


def check(path,completed):
    c = json.loads(path.read_text('utf-8'))
    rows = []
    skip = ('#title:','#waveform:','#hertzian_dipole:','#rx:','#snapshot:','#add_dispersion_debye:')

    def strip(lines):
        return [s for s in lines if not s.startswith(skip)
                and not (s.startswith('#material:') and s.split()[-1] == 'cover')]

    for g in c['groups']:
        p = Path(g['input'])
        if sha256(p)!=g['input_sha256'] or sha256(g['reference_input'])!=g['reference_input_sha256'] \
                or sha256(g['reference_geometry'])!=g['reference_geometry_sha256']:
            raise ValueError('input/reference identity differs')
        lines = p.read_text('utf-8').splitlines()
        ref = Path(g['reference_input']).read_text('utf-8').splitlines()
        if strip(lines) != strip(ref):
            raise ValueError('physical invariant differs')
        cover = [s for s in lines if s.startswith('#material:') and s.split()[-1] == 'cover']
        if cover != [f'#material: {g["cover_eps_r"]} 0.003 1 0 cover'] or any(s.startswith('#add_dispersion_debye:') for s in lines):
            raise ValueError('cover material or dispersion line differs')
        if '#waveform: ricker 1 95e6 ricker95' not in lines \
                or '#hertzian_dipole: y 17.6 0.025 27 ricker95' not in lines \
                or '#rx: 18.9 0.025 27 t01 Ey' not in lines:
            raise ValueError('source/receiver lines differ')
        expected_n = len(c['snapshot_iterations'])
        if len([s for s in lines if s.startswith('#snapshot:')])!=expected_n:
            raise ValueError('snapshot input count differs')
        if any(len(s.split())!=12 for s in lines if s.startswith('#snapshot:')):
            raise ValueError('hash snapshot requires eleven arguments including filename suffix')
        if not completed:
            continue
        raw = p.with_suffix('.h5')
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax'])!='4.0.0' or float(h.attrs['dt'])!=c['dt_s'] or int(h.attrs['Iterations'])!=c['iterations'] or not np.array_equal(h.attrs['dx_dy_dz'],[.025,.05,.025]):
                raise ValueError('actual V4/grid/time differs')
            rx = receivers(h)
            if not np.array_equal(rx['t01'].attrs['GridPosition'],[756,0,1080]) or not np.array_equal(h['srcs/src1'].attrs['GridPosition'],[704,0,1080]):
                raise ValueError('station differs')
            for r in rx.values():
                for v in r.values():
                    if isinstance(v,h5py.Dataset) and (v.dtype!=np.float64 or v.shape!=(c['iterations'],) or not np.isfinite(v[:]).all()):
                        raise ValueError('receiver raw validity differs')
            snap_paths = sorted(p.parent.glob('profile_snaps/*.h5'))
            if len(snap_paths)!=expected_n:
                raise ValueError('snapshot output count differs')
            maximum = {'Ey':0.,'Hx':0.,'Hz':0.}
            snaps = []
            probe_values = {k:{f:r[f][:] for f in r} for k,r in rx.items() if k.startswith('p')}
            for j,s in zip(c['snapshot_iterations'],snap_paths):
                with h5py.File(s) as sh:
                    if int(sh.attrs['iteration'])!=j or float(sh.attrs['time'])!=j*c['dt_s'] or float(sh.attrs['magnetic_time'])!=(j-.5)*c['dt_s'] or not np.array_equal(sh.attrs['origin'],c['snapshot_extent_m'][:3]):
                        raise ValueError('snapshot time/origin differs')
                    arrays = {f:sh[f][:] for f in FIELDS}
                    if any(a.dtype!=np.float64 or list(a.shape)!=c['snapshot_shape_xyz'] or not np.isfinite(a).all() for a in arrays.values()):
                        raise ValueError('snapshot shape/dtype/finite differs')
                    for probe in g['probes']:
                        k,ix,iz=probe['id'],probe['snapshot_ix'],probe['snapshot_iz']
                        expected = {'Ey':probe_values[k+'c']['Ey'][j],
                            'Hx':.5*(probe_values[k+'c']['Hx'][j]+probe_values[k+'z']['Hx'][j]),
                            'Hz':.5*(probe_values[k+'c']['Hz'][j]+probe_values[k+'x']['Hz'][j])}
                        for f,value in expected.items():
                            actual = arrays[f][ix,0,iz]
                            maximum[f]=max(maximum[f],abs(actual-value))
                            if not np.isclose(actual,value,rtol=2e-13,atol=2e-15):
                                raise ValueError('native/snapshot spatial closure failed')
                snaps.append({'file':s.relative_to(p.parent).as_posix(),'sha256':sha256(s),'iteration':j})
        with h5py.File(p.parent/'hs4t2d_geom.vtkhdf') as h,h5py.File(g['reference_geometry']) as ref:
            if not np.array_equal(h['VTKHDF/CellData/Material'][:],ref['VTKHDF/CellData/Material'][:]):
                raise ValueError('actual material geometry differs')
        rows.append({'id':g['id'],'eps_tag':g['eps_tag'],'role':g['role'],'raw_sha256':sha256(raw),
                     'geometry_sha256':sha256(p.parent/'hs4t2d_geom.vtkhdf'),
                     'snapshot_probe_max_abs_difference':maximum,'snapshots':snaps})
    if completed:
        events = [json.loads(s) for s in (path.parent/'execution.jsonl').read_text('utf-8').splitlines()]
        if events[0]['contract_sha256'] != sha256(path):
            raise ValueError('execution contract identity differs')
        done = {e['group']: e for e in events if e.get('status') == 'COMPLETED' and 'group' in e}
        if any(e.get('status') == 'FAILED' for e in events if 'group' in e):
            raise ValueError('a solver group failed')
        if sorted(done) != sorted(g['id'] for g in c['groups']):
            raise ValueError('incomplete or changed execution')
    report={'status':'PASS','completed':completed,'contract_sha256':sha256(path),'code_sha256':sha256(__file__),'groups':rows}
    save(path.parent/('completed_verification.json' if completed else 'preflight_verification.json'),report)
    return report


def run(path):
    import msvcrt
    c=json.loads(path.read_text('utf-8'))
    if c['status']!='FROZEN_APPROVED' or len(c['groups'])!=c['max_runs'] or Path(sys.executable).resolve()!=Path(c['python']).resolve():
        raise ValueError('approved bounded V4 study required')
    for name,digest in c['code_identities'].items():
        if sha256(name)!=digest:
            raise ValueError('frozen code differs')
    package=Path(importlib.util.find_spec('gprMax').origin).parent
    for name,digest in c['source_identities'].items():
        if sha256(package/name)!=digest:
            raise ValueError('V4 runtime differs')
    check(path,False)
    record=path.parent/'execution.jsonl'
    if record.exists() or any(Path(g['input']).with_suffix('.h5').exists() for g in c['groups']):
        raise ValueError('attempt consumed; new contract required for new work')
    lock=ROOT/c['gpu_lock']; lock.parent.mkdir(parents=True,exist_ok=True)
    with open(lock,'a+b') as held:
        held.seek(0); msvcrt.locking(held.fileno(),msvcrt.LK_NBLCK,1)
        def event(row):
            with record.open('a',encoding='utf-8') as f:
                f.write(json.dumps(row,allow_nan=False)+'\n')
        digest=sha256(path)
        start=time.monotonic()
        event({'status':'STARTED','contract_sha256':digest,'unix_s':time.time(),'owned_runner_pid':os.getpid()})
        try:
            for g in c['groups']:
                p=Path(g['input'])
                free=psutil.virtual_memory().available
                vram=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).splitlines()[0])*2**20
                if free<c['min_available_RAM_GiB']*2**30 or vram<c['min_free_VRAM_GiB']*2**30 or sha256(path)!=digest or time.monotonic()-start>c['max_wall_s']:
                    raise RuntimeError('live resource minimum, wall budget or contract identity differs')
                command=[sys.executable,str(ROOT/'scripts/gprmax_cached_cuda_entry.py'),str(p),'-gpu','0','-gpu_precision','double','--hide-progress-bars']
                env=os.environ.copy(); env['HS4_CUDA_CACHE_LOG']=str(p.parent/'cuda_cache.jsonl')
                event({'status':'STARTED','group':g['id'],'command':command,'RAM_available_bytes':free,'VRAM_free_bytes':vram})
                with (p.parent/'stdout.log').open('xb') as out,(p.parent/'stderr.log').open('xb') as err:
                    process=subprocess.Popen(command,cwd=p.parent,stdout=out,stderr=err,env=env)
                    started=time.monotonic(); peak=0
                    while process.poll() is None:
                        try:
                            owned=psutil.Process(process.pid)
                            rss=sum(q.memory_info().rss for q in [owned]+owned.children(recursive=True) if q.is_running())
                        except (psutil.NoSuchProcess,psutil.AccessDenied):
                            rss=0
                        peak=max(peak,rss)
                        if time.monotonic()-started>c['max_group_wall_s'] or rss>c['max_owned_RSS_GiB']*2**30 or psutil.virtual_memory().available<c['min_live_system_available_MiB']*2**20:
                            process.terminate(); process.wait(timeout=20)
                            raise RuntimeError('owned solver exceeded wall/RSS/system memory guard')
                        time.sleep(.5)
                    if process.returncode:
                        raise RuntimeError('solver failed; preserve logs, no retry')
                event({'status':'COMPLETED','group':g['id'],'elapsed_s':time.monotonic()-started,'peak_owned_RSS_bytes':peak,'raw_sha256':sha256(p.with_suffix('.h5'))})
                print('Completed',g['id'],flush=True)
            result=check(path,True)
            event({'status':'COMPLETED','traces':c['max_runs'],'elapsed_s':time.monotonic()-start,'verification_sha256':sha256(path.parent/'completed_verification.json')})
        except Exception as exc:
            event({'status':'FAILED','error':str(exc)})
            raise
        finally:
            held.seek(0); msvcrt.locking(held.fileno(),msvcrt.LK_UNLCK,1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run','check'])
    parser.add_argument('--out',type=Path)
    parser.add_argument('--contract',type=Path)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.action=='prepare': prepare(args.out)
    elif args.action=='run':
        if not args.execute: raise ValueError('--execute required')
        run(args.contract)
    else: print(json.dumps(check(args.contract,True),ensure_ascii=False)[:500])
