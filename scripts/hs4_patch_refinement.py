"""Frozen localized interface controls; preserve all unperturbed geometry."""
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
BASE=ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def prepare(out):
    if out.exists(): raise ValueError('new capsule required')
    old=json.loads((BASE/'execution_contract.json').read_text('utf-8'))
    out.mkdir(parents=True); groups=[]
    third=ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_third/centre_rough'
    for name,patch,sign,reference,spacing in [('fine_crest','crest',1,third,1/60),
            ('fine_slope','slope',1,third,1/60),('reverse_crest','crest',-1,BASE/'centre_rough',.025)]:
        height,count,step=15,1,0
        source=reference/'profile.in'; text=source.read_text('utf-8')
        lines=[]; changes=[]
        for line in text.splitlines():
            w=line.split()
            if line.startswith('#hertzian_dipole:'):
                w[2],w[4]='17.6',str(12+height); line=' '.join(w)
            elif line.startswith('#rx:'):
                w[1],w[3]='18.9',str(12+height); line=' '.join(w)
            elif line.startswith('#box:') and w[-1]=='cover' and abs(float(w[4])-float(w[1])-.25)<1e-12:
                index=round((float(w[1])-12)/.25)
                start={'crest':20,'slope':32}.get(patch,-100)
                if start<=index<start+6:
                    delta=sign*(.05 if index in (start,start+5) else .1)
                    before=float(w[3]); w[3]=f'{before+delta:.12g}'; line=' '.join(w)
                    changes.append({'original_bin':index,'x0_m':float(w[1]),'x1_m':float(w[4]),'original_z_m':before,'new_z_m':float(w[3]),'delta_z_m':delta})
            lines.append(line)
        folder=out/name; folder.mkdir(); p=folder/'profile.in'
        p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        groups.append({'id':name,'input':str(p.resolve()),'sha256':sha256(p),'height_m':height,'patch':patch,
            'traces':count,'step_m':step,'old_station_indices':[61],'spacing_m':[spacing,.05,spacing],
            'source_input':str(source),'source_input_sha256':sha256(source),
            'reference_geometry':str(reference/'hs4t2d_geom.vtkhdf'),'reference_geometry_sha256':sha256(reference/'hs4t2d_geom.vtkhdf'),
            'modified_boxes':changes,'changed_cross_section_m2':sum(r['delta_z_m']*(r['x1_m']-r['x0_m']) for r in changes)})
    contract={'status':'FROZEN_APPROVED','approval_basis':'User2026-10-04 explicitly authorized autonomous investigation until attribution resolved. Local crest/slope controls frozen before execution.',
        'python':old['python'],'source_identities':old['source_identities'],
        'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in ['scripts/hs4_patch_refinement.py','scripts/run_hs4_patch_refinement.cmd','scripts/gprmax_cached_cuda_entry.py']},
        'max_runs':3,'groups':groups,
        'min_available_RAM_GiB':1.5,'min_free_VRAM_GiB':3.5,'max_group_wall_s':600,'max_owned_RSS_GiB':2.5,
        'gpu_lock':'artifacts/local_checks/hs4_gpu_exclusive.lock','no_retry':True,
        'invariants':'36x.05x33m; .025/.05/.025m native grid; source impulse1 and invariant line-source scale; same materials, offset1.3m, fixed physicalPML2m sides1m top/bottom. Height only15/2m. Geometry perturbation limited to6interior .25m boxes; adjacent endpoint continuation/PML unchanged.',
        'interpretation':'Centre61 crest/right-slope equal geometry perturbation at1/60m; reverse crest at2.5cm. Compare to already archived same-grid centre reference; no baseline replay. Physical PML widths/staircase preserved. Supervisor includes actual Python child RSS; historical35run launcher-only RSS is not an actual solver peak.',
        'processing':'Source-normalized complex501tones20-170MHz/Hann/8x/200ns taper; signed and complex differences; high fixed160-220ns, low same window minus fixed air double delay. NoSVD/AGC/pertrace normalization. Numerical grid robustness remains separate.'}
    save(out/'execution_contract.json',contract)
    audit(out/'execution_contract.json',False)
    print('Frozen3traces: finer-centre crest/slope and reverse-crest control')


def expected_geometry(c,g):
    dx,dy,dz=g['spacing_m']
    with h5py.File(g['reference_geometry']) as h:
        expected=h['VTKHDF/CellData/Material'][:]
    for r in g['modified_boxes']:
        a,b=round(r['x0_m']/dx),round(r['x1_m']/dx)
        lo,hi=round(r['original_z_m']/dz),round(r['new_z_m']/dz)
        expected[min(lo,hi):max(lo,hi),0,a:b]=3 if hi>lo else 4
    return expected


def audit(path,completed):
    c=json.loads(path.read_text('utf-8')); rows=[]
    fixed=('#domain:','#dx_dy_dz:','#time_window:','#pml_cells:','#pml_formulation:','#material:','#add_dispersion_debye:','#waveform:')
    for g in c['groups']:
        if sha256(g['source_input'])!=g['source_input_sha256'] or sha256(g['reference_geometry'])!=g['reference_geometry_sha256']:
            raise ValueError('reference identity changed')
        original=Path(g['source_input']).read_text('utf-8').splitlines()
        dx,dy,dz=g['spacing_m']
        p=Path(g['input']); lines=p.read_text('utf-8').splitlines()
        if sha256(p)!=g['sha256'] or [s for s in lines if s.startswith(fixed)]!=[s for s in original if s.startswith(fixed)]:
            raise ValueError('frozen factor changed')
        expected=expected_geometry(c,g); parsed=np.full_like(expected,2)
        for s in lines:
            if s.startswith('#box:'):
                w=s.split(); limits=np.array(list(map(float,w[1:7]))).reshape(2,3)/[dx,dy,dz]
                if not np.allclose(limits,np.rint(limits),atol=1e-10,rtol=0):
                    raise ValueError('off-grid physical box')
                a,b=np.rint(limits).astype(int)
                parsed[a[2]:b[2],a[1]:b[1],a[0]:b[0]]={'rock':3,'cover':4}[w[-1]]
        if not np.array_equal(parsed,expected) or len(g['modified_boxes'])!=(0 if g['patch']=='none' else 6):
            raise ValueError('nonlocal/unexpected geometry change')
        outputs=[]
        if completed:
            for j,index in enumerate(g['old_station_indices'],1):
                raw=p.parent/'profile.h5'; geom=p.parent/'hs4t2d_geom.vtkhdf'
                with h5py.File(raw) as h:
                    v=h['rxs/rx1/Ey'][:]; x=14.6+.05*(index-1); z=12+g['height_m']
                    if str(h.attrs['gprMax'])!='4.0.0' or v.dtype!=np.float64 or not np.isfinite(v).all() or len(v)!=int(h.attrs['Iterations']) or not np.isclose(float(h.attrs['dt']),dx/(299792458*np.sqrt(2)),rtol=1e-12,atol=0):
                        raise ValueError('raw output validity/time differs')
                    for name,xx in [('srcs/src1',x),('rxs/rx1',x+1.3)]:
                        if not np.array_equal(h[name].attrs['GridPosition'],[round(xx/dx),0,round(z/dz)]):
                            raise ValueError('actual station differs')
                with h5py.File(geom) as h:
                    if not np.array_equal(h['VTKHDF/CellData/Material'][:],expected):
                        raise ValueError('actual perturbed material map differs')
                outputs.append({'file':raw.name,'sha256':sha256(raw),'old_station_index':index,'geometry_file':geom.name,'geometry_sha256':sha256(geom)})
        rows.append({'group':g['id'],'geometry_factor_check':'PASS','outputs':outputs})
    report={'status':'PASS','contract_sha256':sha256(path),'code_sha256':sha256(__file__),'completed':completed,'groups':rows}
    save(path.parent/('completed_verification.json' if completed else 'preflight_verification.json'),report)
    return report


def run(path):
    import msvcrt
    c=json.loads(path.read_text('utf-8'))
    if c['status']!='FROZEN_APPROVED' or sum(g['traces'] for g in c['groups'])!=c['max_runs'] or Path(sys.executable).resolve()!=Path(c['python']).resolve():
        raise ValueError('approved bounded runtime required')
    package=Path(importlib.util.find_spec('gprMax').origin).parent
    for name,digest in c['code_identities'].items():
        if sha256(name)!=digest: raise ValueError('code changed')
    for name,digest in c['source_identities'].items():
        if sha256(package/name)!=digest: raise ValueError('runtime changed')
    audit(path,False); record=path.parent/'execution.jsonl'
    if record.exists() or any(list(Path(g['input']).parent.glob('profile*.h5')) for g in c['groups']):
        raise ValueError('attempt consumed')
    with (ROOT/c['gpu_lock']).open('a+b') as lock:
        lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        def event(r):
            with record.open('a',encoding='utf-8') as f: f.write(json.dumps(r,allow_nan=False)+'\n')
        event({'status':'STARTED','contract_sha256':sha256(path),'unix_s':time.time(),'owned_pid':os.getpid()})
        try:
            for g in c['groups']:
                p=Path(g['input']); free=psutil.virtual_memory().available
                vram=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).splitlines()[0])*2**20
                if free<c['min_available_RAM_GiB']*2**30 or vram<c['min_free_VRAM_GiB']*2**30: raise RuntimeError('live resource guard')
                command=[sys.executable,str(ROOT/'scripts/gprmax_cached_cuda_entry.py'),str(p),'-n',str(g['traces']),'--geometry-fixed','-gpu','0','-gpu_precision','double','--hide-progress-bars']
                env=os.environ.copy(); env['HS4_CUDA_CACHE_LOG']=str(p.parent/'cuda_cache.jsonl')
                event({'status':'STARTED','group':g['id'],'command':command,'available_RAM_bytes':free,'free_VRAM_bytes':vram})
                with (p.parent/'stdout.log').open('xb') as out,(p.parent/'stderr.log').open('xb') as err:
                    process=subprocess.Popen(command,cwd=p.parent,env=env,stdout=out,stderr=err); start=time.monotonic(); peak=0
                    while process.poll() is None:
                        try:
                            owned=psutil.Process(process.pid)
                            peak=max(peak,sum(q.memory_info().rss for q in [owned]+owned.children(recursive=True) if q.is_running()))
                        except psutil.NoSuchProcess: pass
                        if time.monotonic()-start>c['max_group_wall_s'] or peak>c['max_owned_RSS_GiB']*2**30 or psutil.virtual_memory().available<100*2**20:
                            process.terminate(); process.wait(timeout=20); raise RuntimeError('owned solver wall/RSS/system memory guard')
                        time.sleep(.5)
                    if process.returncode: raise RuntimeError('solver failed; no retry')
                event({'status':'COMPLETED','group':g['id'],'traces':g['traces'],'elapsed_s':time.monotonic()-start,'peak_owned_RSS_bytes':peak})
                print('Completed',g['id'],flush=True)
            audit(path,True); event({'status':'COMPLETED','traces':c['max_runs'],'verification_sha256':sha256(path.parent/'completed_verification.json')})
        except Exception as exc:
            event({'status':'FAILED','error':str(exc)}); raise
        finally:
            lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=['prepare','run','check'])
    p.add_argument('--out',type=Path); p.add_argument('--contract',type=Path); p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if a.action=='prepare': prepare(a.out)
    elif a.action=='check': audit(a.contract,True)
    elif a.execute: run(a.contract)
    else: raise ValueError('--execute required')
