"""Bounded native-grid controls; frozen historical inputs are read-only."""
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

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'artifacts/research_checks/2026-10-04_hs4_patch_finest'


def save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def commands(path):
    result={}
    for line in path.read_text('utf-8').splitlines():
        if line.startswith('#'):
            key,*values=line.split(); result.setdefault(key,[]).append(values)
    return result


def raster(c):
    spacing=np.array(list(map(float,c['#dx_dy_dz:'][0])))
    shape=np.rint(np.array(list(map(float,c['#domain:'][0])))/spacing).astype(int)
    material=np.full(tuple(shape[::-1]),2,dtype=np.uint32)
    for words in c['#box:']:
        limits=np.array(list(map(float,words[:6]))).reshape(2,3)/spacing
        if not np.allclose(limits,np.rint(limits),rtol=0,atol=1e-9): raise ValueError('off-grid box')
        lo,hi=np.rint(limits).astype(int)
        material[lo[2]:hi[2],lo[1]:hi[1],lo[0]:hi[0]]={'rock':3,'cover':4}[words[-1]]
    return spacing,shape,material


def prepare(out,mode):
    if out.exists(): raise ValueError('new capsule required')
    old=json.loads((REFERENCE/'execution_contract.json').read_text('utf-8'))
    dx=.0125 if mode=='endpoints' else .01
    stations=[('left',1),('right',121)] if mode=='endpoints' else [('centre',61)]
    out.mkdir(parents=True); groups=[]
    for station,index in stations:
        for variant in ('base','halfspace','crest','slope'):
            source=REFERENCE/variant/'profile.in'
            tx=14.6+.05*(index-1); rx=tx+1.3
            lines=[]
            for line in source.read_text('utf-8').splitlines():
                w=line.split()
                if line.startswith('#dx_dy_dz:'): line=f'#dx_dy_dz: {dx} 0.05 {dx}'
                elif line.startswith('#pml_cells:'): line=f'#pml_cells: {round(2/dx)} 0 {round(1/dx)} {round(2/dx)} 0 {round(1/dx)}'
                elif line.startswith('#geometry_view:'): line=f'#geometry_view: 0 0 0 36 0.05 33 {dx} 0.05 {dx} hs4t2d_geom n'
                elif line.startswith('#hertzian_dipole:'): w[2]=f'{tx:.12g}'; line=' '.join(w)
                elif line.startswith('#rx:'): w[1]=f'{rx:.12g}'; line=' '.join(w)
                lines.append(line)
            folder=out/(station+'_'+variant); folder.mkdir(); p=folder/'profile.in'
            p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            groups.append({'id':folder.name,'station':station,'old_station_index':index,'variant':variant,
                'input':str(p.resolve()),'input_sha256':sha256(p),'source_input':str(source),'source_input_sha256':sha256(source),
                'tx_m':[tx,0.,27.],'rx_m':[rx,0.,27.],'spacing_m':[dx,.05,dx],
                'reference_geometry':str(REFERENCE/variant/'hs4t2d_geom.vtkhdf'),
                'reference_geometry_sha256':sha256(REFERENCE/variant/'hs4t2d_geom.vtkhdf')})
    c={'status':'FROZEN_APPROVED','approval_basis':'User2026-10-04 says do the next grid/station and3D reference controls; prior autonomous authorization persists.',
        'mode':mode,'python':old['python'],'source_identities':old['source_identities'],'groups':groups,'max_runs':len(groups),
        'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in ['scripts/hs4_station_grid_controls.py','scripts/run_hs4_station_grid_controls.cmd','scripts/gprmax_cached_cuda_entry.py']},
        'min_available_RAM_GiB':2.25 if mode=='endpoints' else 3.75,'min_free_VRAM_GiB':3.5 if mode=='endpoints' else 4.3,
        'max_owned_RSS_GiB':3. if mode=='endpoints' else 4.3,'min_system_available_during_run_GiB':.2,
        'max_group_wall_s':900,'max_batch_wall_s':900*len(groups),'no_retry':True,
        'gpu_lock':'artifacts/local_checks/hs4_gpu_exclusive.lock',
        'processing':'Native source-normalized501complex20-170MHz tones,200ns tail taper,Hann/8x; fixed160-180/180-220/160-220ns; no peak-picked windows or normalization.',
        'invariants':'36x.05x33m; same physical box/material geometry and impulse source;15m height;2m sides1m top/bottom PML. Only declared native X/Z spacing and actual station change. Invariant-Y TM line source, not finite3D.'}
    save(out/'execution_contract.json',c)
    save(out/'preflight_verification.json',audit(out/'execution_contract.json',False))
    print('Frozen',len(groups),'traces',mode,flush=True)


def audit(path,completed):
    c=json.loads(path.read_text('utf-8')); rows=[]
    invariant=('#domain:','#time_window:','#material:','#add_dispersion_debye:','#waveform:','#box:')
    for g in c['groups']:
        p=Path(g['input']); source=Path(g['source_input'])
        if sha256(p)!=g['input_sha256'] or sha256(source)!=g['source_input_sha256']: raise ValueError('frozen input changed')
        a,b=commands(p),commands(source)
        if any(a.get(k)!=b.get(k) for k in invariant): raise ValueError('physical/source factor changed')
        spacing,shape,material=raster(a)
        if not np.array_equal(spacing,g['spacing_m']): raise ValueError('spacing differs')
        pml=np.array(list(map(int,a['#pml_cells:'][0])))
        if not np.allclose(pml[[0,2,3,5]]*spacing[[0,2,0,2]],[2,1,2,1],rtol=0,atol=1e-12): raise ValueError('PML widths differ')
        for position in (g['tx_m'],g['rx_m']):
            if not np.allclose(np.array(position)/spacing,np.rint(np.array(position)/spacing),rtol=0,atol=1e-9): raise ValueError('off-grid station')
        record={'id':g['id'],'input_sha256':sha256(p),'native_grid_cells':shape.tolist(),'physical_factor_check':'PASS'}
        if completed:
            raw=p.with_suffix('.h5'); geom=p.parent/'hs4t2d_geom.vtkhdf'
            with h5py.File(raw) as h:
                v=h['rxs/rx1/Ey'][:]; dt=float(h.attrs['dt']); iterations=int(h.attrs['Iterations'])
                if str(h.attrs['gprMax'])!='4.0.0' or not np.array_equal(h.attrs['dx_dy_dz'],spacing) or not np.array_equal(h.attrs['nx_ny_nz'],shape) or not np.isclose(dt,spacing[0]/(299792458*np.sqrt(2)),rtol=1e-12,atol=0): raise ValueError('runtime/native grid/time differs')
                if v.dtype!=np.float64 or not np.isfinite(v).all() or len(v)!=iterations: raise ValueError('raw precision/validity differs')
                for key,position in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                    if not np.array_equal(h[key].attrs['GridPosition'],np.rint(np.array(position)/spacing).astype(int)): raise ValueError('actual station differs')
            with h5py.File(geom) as h:
                if not np.array_equal(h['VTKHDF/CellData/Material'][:],material): raise ValueError('actual material map differs')
            if c['mode']=='endpoints':
                if sha256(g['reference_geometry'])!=g['reference_geometry_sha256']: raise ValueError('reference grid changed')
                with h5py.File(g['reference_geometry']) as h:
                    if not np.array_equal(h['VTKHDF/CellData/Material'][:],material): raise ValueError('station changed material geometry')
            record.update({'raw_sha256':sha256(raw),'geometry_sha256':sha256(geom),'dt_s':dt,'iterations':iterations,'dtype':str(v.dtype)})
        rows.append(record)
    return {'status':'PASS','completed':completed,'code_sha256':sha256(__file__),'contract_sha256':sha256(path),'groups':rows}


def live_resources(c):
    ram=psutil.virtual_memory().available
    vram=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).splitlines()[0])*2**20
    if ram<c['min_available_RAM_GiB']*2**30 or vram<c['min_free_VRAM_GiB']*2**30: raise RuntimeError('live RAM/VRAM preflight: no solve authorized by capacity')
    return {'available_RAM_bytes':ram,'free_VRAM_bytes':vram}


def terminate_owned_tree(process):
    try: children=psutil.Process(process.pid).children(recursive=True)
    except psutil.NoSuchProcess: children=[]
    for child in reversed(children):
        try: child.terminate()
        except psutil.NoSuchProcess: pass
    process.terminate() if process.poll() is None else None
    _,alive=psutil.wait_procs(children,timeout=5)
    for child in alive:
        try: child.kill()
        except psutil.NoSuchProcess: pass
    process.wait(timeout=10)


def run(path):
    import msvcrt
    c=json.loads(path.read_text('utf-8')); package=Path(importlib.util.find_spec('gprMax').origin).parent
    if c['status']!='FROZEN_APPROVED' or len(c['groups'])!=c['max_runs'] or Path(sys.executable).resolve()!=Path(c['python']).resolve(): raise ValueError('frozen bounded runtime required')
    for name,digest in c['code_identities'].items():
        if sha256(name)!=digest: raise ValueError('execution code changed')
    for name,digest in c['source_identities'].items():
        if sha256(package/name)!=digest: raise ValueError('runtime changed')
    audit(path,False)
    log=path.parent/'execution.jsonl'
    if log.exists() or any(Path(g['input']).with_suffix('.h5').exists() for g in c['groups']): raise ValueError('attempt consumed; no retry')
    first_resources=live_resources(c)  # No attempt consumed when preflight alone rejects.
    with (ROOT/c['gpu_lock']).open('a+b') as lock:
        lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        def event(value):
            with log.open('a',encoding='utf-8') as f: f.write(json.dumps(value,allow_nan=False)+'\n')
        start_batch=time.monotonic(); process=None
        event({'status':'STARTED','contract_sha256':sha256(path),'unix_s':time.time(),'runner_pid':os.getpid(),**first_resources})
        try:
            for g in c['groups']:
                resources=live_resources(c); p=Path(g['input'])
                extra=c.get('additional_solver_args',[])
                if extra != []:
                    raise ValueError('Unsupported additional solver arguments')
                entry=c.get('solver_entrypoint','gprmax_cached_cuda_entry.py')
                if entry not in ('gprmax_cached_cuda_entry.py','gprmax_snapshot_cuda_entry.py'):
                    raise ValueError('Unsupported solver entrypoint')
                command=[sys.executable,str(ROOT/'scripts'/entry),str(p),'-gpu','0','-gpu_precision','double','--hide-progress-bars',*extra]
                env=os.environ.copy(); env['HS4_CUDA_CACHE_LOG']=str(p.parent/'cuda_cache.jsonl')
                event({'status':'STARTED','group':g['id'],'command':command,**resources})
                with (p.parent/'stdout.log').open('xb') as out,(p.parent/'stderr.log').open('xb') as err:
                    process=subprocess.Popen(command,cwd=p.parent,env=env,stdout=out,stderr=err)
                    start=time.monotonic(); peak=0
                    while process.poll() is None:
                        try:
                            parent=psutil.Process(process.pid); rss=0
                            for owned in [parent]+parent.children(recursive=True):
                                try: rss+=owned.memory_info().rss
                                except psutil.NoSuchProcess: pass
                            peak=max(peak,rss)
                        except psutil.NoSuchProcess: pass
                        if time.monotonic()-start>c['max_group_wall_s'] or time.monotonic()-start_batch>c['max_batch_wall_s'] or peak>c['max_owned_RSS_GiB']*2**30 or psutil.virtual_memory().available<c['min_system_available_during_run_GiB']*2**30:
                            event({'status':'GUARD_STOP','group':g['id'],'elapsed_s':time.monotonic()-start,
                                'peak_owned_RSS_bytes':peak,'available_RAM_bytes':psutil.virtual_memory().available,
                                'max_owned_RSS_GiB':c['max_owned_RSS_GiB']})
                            terminate_owned_tree(process); raise RuntimeError('owned solver resource/wall guard')
                        time.sleep(.5)
                    if process.returncode: raise RuntimeError('solver failed; preserve attempt')
                event({'status':'COMPLETED','group':g['id'],'elapsed_s':time.monotonic()-start,'peak_owned_RSS_bytes':peak,'raw_sha256':sha256(p.with_suffix('.h5'))})
                print('Completed',g['id'],flush=True)
                process=None
            verification=audit(path,True)
            save(path.parent/'completed_verification.json',verification)
            event({'status':'COMPLETED','traces':len(c['groups']),'verification_sha256':sha256(path.parent/'completed_verification.json')})
        except BaseException as exc:
            if process is not None and process.poll() is None: terminate_owned_tree(process)
            event({'status':'FAILED','error':str(exc)}); raise
        finally:
            lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=['prepare','run','verify'])
    p.add_argument('--out',type=Path); p.add_argument('--mode',choices=['endpoints','centre1cm'])
    p.add_argument('--contract',type=Path); p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if a.action=='prepare': prepare(a.out,a.mode)
    elif a.action=='verify':
        if a.out.exists(): raise ValueError('new audit output required')
        save(a.out,audit(a.contract,True))
    elif a.execute: run(a.contract)
    else: raise ValueError('--execute required')
