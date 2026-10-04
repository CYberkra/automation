"""Execute frozen joint-grid groups; preserve every attempt and raw output."""
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


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args()
    c=json.loads(args.contract.read_text('utf-8'))
    if not args.execute or c['status']!='FROZEN_APPROVED' or sum(g['traces'] for g in c['groups'])!=c['max_runs']:
        raise ValueError('approved bounded study and --execute required')
    if Path(c['python']).resolve()!=Path(sys.executable).resolve():
        raise ValueError('Python differs')
    root=Path(__file__).resolve().parents[1]
    package=Path(importlib.util.find_spec('gprMax').origin).parent
    for name,digest in c['code_identities'].items():
        if sha256(name)!=digest:
            raise ValueError('runner identity differs')
    for name,digest in c['source_identities'].items():
        if sha256(package/name)!=digest:
            raise ValueError('V4 source identity differs')
    from check_hs4t2d_joint_grid_third import check
    check(args.contract, completed=False)
    prerequisite=c.get('prerequisite')
    if prerequisite:
        pilot=Path(prerequisite['directory'])
        if sha256(pilot/'execution_contract.json')!=prerequisite['contract_sha256'] or sha256(pilot/'completed_verification.json')!=prerequisite['completed_verification_sha256']:
            raise ValueError('pilot identity differs')
    record=args.contract.parent/'execution.jsonl'
    if record.exists():
        raise ValueError('attempt consumed; no automatic retry')
    digest=sha256(args.contract)
    def event(row):
        with record.open('a',encoding='utf-8') as f:
            f.write(json.dumps(row,allow_nan=False)+'\n')
    start=time.monotonic()
    event({'status':'STARTED','contract_sha256':digest,'start_unix_s':time.time()})
    try:
        for g in c['groups']:
            path=Path(g['input'])
            if sha256(path)!=g['sha256'] or sha256(root/g['source_input'])!=g['source_sha256']:
                raise ValueError('input changed')
            if list(path.parent.glob('profile*.h5')):
                raise ValueError('outputs already exist')
            check_start=time.monotonic()
            while True:
                free=psutil.virtual_memory().available
                vram=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).splitlines()[0])*2**20
                if free>=c['min_available_RAM_GiB']*2**30 and vram>=g['minimum_free_VRAM_GiB']*2**30:
                    break
                if time.monotonic()-check_start>60:
                    raise RuntimeError('RAM/VRAM below frozen preflight minimum')
                time.sleep(2)
            remaining=c['max_wall_s']-(time.monotonic()-start)
            if remaining<=0 or sha256(args.contract)!=digest:
                raise ValueError('wall budget or frozen identity differs')
            command=[sys.executable,str(root/'scripts/gprmax_cached_cuda_entry.py'),str(path),
                     '-n',str(g['traces']),'--geometry-fixed','-gpu','0','-gpu_precision','double','--hide-progress-bars']
            event({'group':g['id'],'status':'STARTED','command':command,'RAM_available_bytes':free,'VRAM_free_bytes':vram})
            environment=os.environ.copy()
            environment['HS4_CUDA_CACHE_LOG']=str(path.parent/'cuda_cache.jsonl')
            with (path.parent/'stdout.log').open('xb') as out,(path.parent/'stderr.log').open('xb') as err:
                solved=subprocess.run(command,cwd=path.parent,stdout=out,stderr=err,timeout=remaining,env=environment)
            if solved.returncode:
                raise RuntimeError('solver failed; saved stdout/stderr')
            rows=[]
            for k in range(1,g['traces']+1):
                p=path.parent/('profile.h5' if g['traces']==1 else f'profile{k}.h5')
                with h5py.File(p,'r') as h:
                    raw=h['rxs/rx1/Ey'][:]
                    expected_tx=g['station_x_m'][k-1]
                    dx,dy,dz=g['spacing_m']
                    if (str(h.attrs['gprMax'])!='4.0.0' or raw.dtype!=np.float64 or not np.isfinite(raw).all()
                            or len(raw)!=int(h.attrs['Iterations'])):
                        raise ValueError('raw version/dtype/shape/finite audit failed')
                    for name,x in [('srcs/src1',expected_tx),('rxs/rx1',expected_tx+1.3)]:
                        expected=[round(x/dx),0,round(27/dz)]
                        if not np.array_equal(h[name].attrs['GridPosition'],expected):
                            raise ValueError('actual station grid differs')
                    if not np.array_equal(h.attrs['dx_dy_dz'],g['spacing_m']):
                        raise ValueError('actual grid spacing differs')
                    rows.append({'file':p.name,'sha256':sha256(p),'trace':k,'old_station_index':g['old_station_indices'][k-1],
                                 'dt_s':float(h.attrs['dt']),'samples':len(raw),'dtype':str(raw.dtype),
                                 'source_m':h['srcs/src1'].attrs['Position'].tolist(),
                                 'receiver_m':h['rxs/rx1'].attrs['Position'].tolist()})
            (path.parent/'audit.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
            event({'group':g['id'],'status':'COMPLETED','traces':len(rows),'elapsed_s':time.monotonic()-start})
            print('Completed',g['id'],len(rows),flush=True)
        event({'status':'COMPLETED','traces':c['max_runs'],'elapsed_s':time.monotonic()-start})
    except Exception as exc:
        event({'status':'FAILED','group':g['id'],'error':str(exc),'elapsed_s':time.monotonic()-start})
        raise


if __name__=='__main__':
    main()
