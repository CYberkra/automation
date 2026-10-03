"""Run the approved fresh three-group CUDA-double comparison without retries."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

from hs_capsule_identity import sha256


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--execute',action='store_true')
    args = ap.parse_args()
    c = json.loads(args.contract.read_text('utf-8'))
    if not args.execute or c['status'] != 'FROZEN_APPROVED' or c['backend'] != 'CUDA':
        raise ValueError('approved CUDA study and --execute required')
    if Path(c['python']).resolve() != Path(sys.executable).resolve():
        raise ValueError('wrong Python')
    if 'runner_sha256' in c and sha256(__file__) != c['runner_sha256']:
        raise ValueError('runner identity differs')
    root = Path(__file__).resolve().parents[1]
    run = Path(c['run_directory'])
    record = run/'series_execution.jsonl'
    if record.exists():
        raise ValueError('attempt consumed; no automatic rerun')
    import importlib.util
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    for name,digest in c['source_identities'].items():
        if sha256(package/name) != digest:
            raise ValueError('solver source identity differs')
    digest = sha256(args.contract)
    cache_entry = root/'scripts/gprmax_cached_cuda_entry.py'
    if 'cuda_cache_entry_sha256' in c and sha256(cache_entry) != c['cuda_cache_entry_sha256']:
        raise ValueError('cache entry identity differs')
    def event(row):
        with record.open('a',encoding='utf-8') as f:
            f.write(json.dumps(row,allow_nan=False)+'\n')
    start=time.monotonic()
    event({'status':'STARTED','contract_sha256':digest,'start_unix_s':time.time()})
    for group in c['groups']:
        path = Path(group['input'])
        if sha256(path) != group['sha256'] or sha256(root/group['source_input']) != group['source_sha256']:
            raise RuntimeError('input changed')
        if list(path.parent.glob('profile*.h5')):
            raise RuntimeError('outputs already exist')
        memory_wait_start = time.monotonic()
        free=psutil.virtual_memory().available
        while free < c['memory_min_available_GiB']*2**30:
            if time.monotonic()-memory_wait_start >= 60:
                event({'group':group['id'],'status':'PREFLIGHT_STOPPED','available_RAM_bytes':free})
                raise RuntimeError('available memory below contract for 60 seconds')
            time.sleep(2)
            free=psutil.virtual_memory().available
        remaining = c['max_wall_s']-(time.monotonic()-start)
        if remaining <= 0 or sha256(args.contract) != digest:
            raise RuntimeError('wall budget or contract identity violation')
        command=([sys.executable,str(cache_entry)] if 'cuda_cache_entry_sha256' in c else [sys.executable,'-m','gprMax']) + [str(path),'-n',str(group['traces']),
                 '--geometry-fixed','-gpu',str(c['gpu_device']),'-gpu_precision','double','--hide-progress-bars']
        event({'group':group['id'],'status':'STARTED','command':command,'available_RAM_bytes':free})
        with (path.parent/'stdout.log').open('xb') as stdout,(path.parent/'stderr.log').open('xb') as stderr:
            import os
            environment = os.environ.copy()
            environment['HS4_CUDA_CACHE_LOG'] = str(path.parent/'cuda_cache.jsonl')
            result=subprocess.run(command,cwd=path.parent,stdout=stdout,stderr=stderr,timeout=remaining,env=environment)
        if result.returncode != 0:
            event({'group':group['id'],'status':'FAILED','exit_code':result.returncode})
            raise RuntimeError('solver failed; see saved logs')
        import h5py
        import numpy as np
        rows=[]
        for k in range(1,group['traces']+1):
            output=path.parent/f'profile{k}.h5'
            original=root/f'artifacts/research_checks/2026-10-02_hs4t2d_transect/hs4t2d_t{k:02d}.h5'
            with h5py.File(output,'r') as h,h5py.File(original,'r') as old:
                raw=h['rxs/rx1/Ey'][:]
                if raw.dtype != np.float64 or raw.shape != (5089,) or not np.isfinite(raw).all():
                    raise RuntimeError('raw dtype/shape/finite audit failed')
                if str(h.attrs['gprMax']) != '4.0.0' or h.attrs['dt'] != old.attrs['dt']:
                    raise RuntimeError('version/dt differs')
                for name in ['srcs/src1','rxs/rx1']:
                    if not np.array_equal(h[name].attrs['Position'],old[name].attrs['Position']):
                        raise RuntimeError('station differs from archived input')
                rows.append({'trace':k,'file':output.name,'sha256':sha256(output),'bytes':output.stat().st_size,
                             'dt':float(h.attrs['dt']),'dtype':str(raw.dtype),'solver':'4.0.0',
                             'source_m':h['srcs/src1'].attrs['Position'].tolist(),
                             'receiver_m':h['rxs/rx1'].attrs['Position'].tolist(),
                             'archived_raw_max_abs_difference':float(np.max(np.abs(raw-old['rxs/rx1/Ey'][:]))) if group['id']=='baseline' else None})
        (path.parent/'audit.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        event({'group':group['id'],'status':'COMPLETED','traces':len(rows),'elapsed_s':time.monotonic()-start})
        print('Completed',group['id'],len(rows),'traces',flush=True)
    event({'status':'COMPLETED','traces':sum(g['traces'] for g in c['groups']),'elapsed_s':time.monotonic()-start})


if __name__=='__main__':
    main()
