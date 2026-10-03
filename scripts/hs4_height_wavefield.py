"""Prepare/run/audit bounded matched-height V4 wavefields without solver edits."""
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
CELLS = [(i, 6) for i in (20,30,40,50,60)] + [(i,26) for i in (20,40,60)] + [(46,39),(46,126),(30,60),(50,60)]
FIELDS = ('Ex','Ey','Ez','Hx','Hy','Hz')


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def prepare(folder):
    if folder.exists():
        raise ValueError('new capsule required')
    old = json.loads((REFERENCE/'execution_contract.json').read_text('utf-8'))
    with h5py.File(REFERENCE/'centre_rough/profile.h5') as h:
        n, dt = len(h['rxs/rx1/Ey']), float(h.attrs['dt'])
    folder.mkdir(parents=True)
    groups = []
    for name,height,role,snaps in [('passive_reference',15,'rough',False),
            ('high_rough',15,'rough',True),('high_halfspace',15,'halfspace',True),
            ('low_rough',2,'rough',True),('low_halfspace',2,'halfspace',True)]:
        original = REFERENCE/('centre_'+role)/'profile.in'
        lines = []
        for line in original.read_text('utf-8').splitlines():
            w = line.split()
            if line.startswith('#waveform:'):
                line = '#waveform: ricker 1 95e6 ricker95'
            elif line.startswith('#hertzian_dipole:'):
                w[4],w[-1] = str(12+height),'ricker95'
                line = ' '.join(w)
            elif line.startswith('#rx:'):
                w[3] = str(12+height)
                line = ' '.join(w)
            lines.append(line)
        probes = []
        for k,(ix,iz) in enumerate(CELLS):
            # Output step is six native cells: Ey is exactly at centre node;
            # Hx/Hz interpolate two half-cell samples along Z/X respectively.
            x,z = 12+.15*ix+.075, 8+.15*iz+.075
            for suffix,px,pz,components in [('c',x,z,'Ey Hx Hz'),('z',x,z-.025,'Hx'),('x',x-.025,z,'Hz')]:
                lines.append(f'#rx: {px:.12g} 0.025 {pz:.12g} p{k:02d}{suffix} {components}')
            probes.append({'id':f'p{k:02d}', 'snapshot_ix':ix, 'snapshot_iz':iz,'centre_xz_m':[x,z]})
        if snaps:
            lines.extend(f'#snapshot: 12 0 8 24 0.05 28.1 0.15 0.05 0.15 {j} snap{j:05d} .h5' for j in range(0,n,10))
        directory = folder/name
        directory.mkdir()
        path = directory/'profile.in'
        path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        groups.append({'id':name,'input':str(path.resolve()),'input_sha256':sha256(path),
            'height_m':height,'role':role,'snapshot':snaps,'probes':probes,
            'reference_input':str(original),'reference_input_sha256':sha256(original),
            'reference_geometry':str(REFERENCE/('centre_'+role)/'hs4t2d_geom.vtkhdf'),
            'reference_geometry_sha256':sha256(REFERENCE/('centre_'+role)/'hs4t2d_geom.vtkhdf')})
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    for name,digest in old['source_identities'].items():
        if sha256(package/name)!=digest:
            raise ValueError('old V4 source differs: '+name)
    c = {'status':'FROZEN_APPROVED','approval_basis':'User2026-10-04: 自主推进直到彻底解决问题; bounded wavefield pilot frozen before execution.',
         'python':old['python'],'source_identities':old['source_identities'],
         'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in ['scripts/hs4_height_wavefield.py','scripts/run_hs4_height_wavefield.cmd','scripts/gprmax_cached_cuda_entry.py']},
         'max_runs':5,'groups':groups,'dt_s':dt,'iterations':n,'snapshot_iterations':list(range(0,n,10)),
         'snapshot_shape_xyz':[80,1,134],'snapshot_extent_m':[12,0,8,24,.05,28.1],'snapshot_spacing_m':[.15,.05,.15],
         'snapshot_six_field_payload_GiB':80*134*6*8*len(range(0,n,10))/2**30,
         'output_variant_reason':'Available RAM fell below old design budget. Only output spacing changes .1->.15m; fine2.5cm FDTD, geometry, source/height comparisons and native interpolation closure preserved. One ROI; actual auto GPU2CPU mode logged, no implementation edits.',
         'min_available_RAM_GiB_snapshot':1.25,'min_available_RAM_GiB_plain':.8,'min_free_VRAM_GiB':3.5,
         'max_group_wall_s':600,'max_owned_RSS_GiB':1.75,'min_live_system_available_MiB':100,
         'no_retry':True,'gpu_lock':'artifacts/local_checks/hs4_gpu_exclusive.lock',
         'processing':'Matched saved-source-normalized 501 frequencies20-170MHz; physical200ns taper/Hann/8x; native probes, E/H half-step and source equivalence required. No AGC/SVD/perframe scaling.',
         'scope':'TM invariant Y; two heights, matched cover background, passive replay; no unique spatial attribution or 3D/field claim from centre alone.'}
    save(folder/'execution_contract.json',c)
    check(folder/'execution_contract.json',False)
    print('Frozen5runs;',c['snapshot_six_field_payload_GiB'],'GiB snapshot payload per trace')


def receivers(h):
    return {str(rx.attrs['Name']):rx for rx in h['rxs'].values()}


def check(path,completed):
    c = json.loads(path.read_text('utf-8'))
    rows = []
    for g in c['groups']:
        p = Path(g['input'])
        if sha256(p)!=g['input_sha256'] or sha256(g['reference_input'])!=g['reference_input_sha256'] or sha256(g['reference_geometry'])!=g['reference_geometry_sha256']:
            raise ValueError('input/reference identity differs')
        lines = p.read_text('utf-8').splitlines()
        ref = Path(g['reference_input']).read_text('utf-8').splitlines()
        preserved = ('#domain:','#dx_dy_dz:','#time_window:','#pml_cells:','#pml_formulation:','#material:','#add_dispersion_debye:','#box:','#geometry_view:')
        if [s for s in lines if s.startswith(preserved)]!=[s for s in ref if s.startswith(preserved)]:
            raise ValueError('physical invariant differs')
        expected_n = len(c['snapshot_iterations']) if g['snapshot'] else 0
        if len([s for s in lines if s.startswith('#snapshot:')])!=expected_n:
            raise ValueError('snapshot input count differs')
        if not completed:
            continue
        raw = p.with_suffix('.h5')
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax'])!='4.0.0' or float(h.attrs['dt'])!=c['dt_s'] or int(h.attrs['Iterations'])!=c['iterations'] or not np.array_equal(h.attrs['dx_dy_dz'],[.025,.05,.025]):
                raise ValueError('actual V4/grid/time differs')
            rx = receivers(h)
            if not np.array_equal(rx['t01'].attrs['GridPosition'],[756,0,round((12+g['height_m'])/.025)]) or not np.array_equal(h['srcs/src1'].attrs['GridPosition'],[704,0,round((12+g['height_m'])/.025)]):
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
        rows.append({'id':g['id'],'raw_sha256':sha256(raw),'geometry_sha256':sha256(p.parent/'hs4t2d_geom.vtkhdf'),
                     'snapshot_probe_max_abs_difference':maximum,'snapshots':snaps})
    if completed:
        a,b = path.parent/'passive_reference/profile.h5',path.parent/'high_rough/profile.h5'
        with h5py.File(a) as h,h5py.File(b) as r:
            for name,v in receivers(h).items():
                for field in v:
                    if not np.array_equal(v[field][:],receivers(r)[name][field][:]):
                        raise ValueError('snapshot observer changed native receiver')
    report={'status':'PASS','completed':completed,'contract_sha256':sha256(path),'code_sha256':sha256(__file__),'groups':rows,'passive_receiver_bit_identical':completed}
    if not completed:
        save(path.parent/'preflight_verification.json',report)
    else:
        save(path.parent/'completed_verification.json',report)
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
        event({'status':'STARTED','contract_sha256':digest,'unix_s':time.time(),'owned_runner_pid':os.getpid()})
        try:
            for g in c['groups']:
                p=Path(g['input'])
                free=psutil.virtual_memory().available
                vram=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).splitlines()[0])*2**20
                minimum=c['min_available_RAM_GiB_snapshot' if g['snapshot'] else 'min_available_RAM_GiB_plain']
                if free<minimum*2**30 or vram<c['min_free_VRAM_GiB']*2**30 or sha256(path)!=digest:
                    raise RuntimeError('live resource minimum or contract identity differs')
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
            event({'status':'COMPLETED','traces':c['max_runs'],'verification_sha256':sha256(path.parent/'completed_verification.json')})
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
