"""Execute only the separately frozen HS4T2D comparison, serially in FP64."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--contract', type=Path, required=True)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    contract = json.loads(args.contract.read_text('utf-8'))
    if not args.execute or contract['status'] != 'FROZEN_APPROVED':
        raise ValueError('explicit --execute and approved frozen experiment required')
    digest = sha256(args.contract)
    run = Path(contract['run_directory'])
    if (run / 'execution.jsonl').exists():
        raise ValueError('execution record exists; refuse automatic retry/overwrite')
    if Path(sys.executable).resolve() != Path(contract['python']).resolve():
        raise ValueError('Python identity differs')
    start = time.monotonic()
    cases = contract['cases']
    for case in cases:
        if sha256(case['input']) != case['sha256'] or Path(case['input']).with_suffix('.h5').exists():
            raise ValueError('input differs or output already exists')
    completed = []
    def event(value):
        with (run / 'execution.jsonl').open('a', encoding='utf-8') as f:
            f.write(json.dumps(value, allow_nan=False) + '\n')
    event({'status':'STARTED','contract_sha256':digest,'started_unix_s':time.time(),'cases':len(cases)})
    try:
        for case in cases:
            if time.monotonic() - start >= contract['batch_wall_limit_s']:
                raise RuntimeError('batch wall budget reached')
            available = psutil.virtual_memory().available
            if available < contract['memory_min_available_GiB'] * 2**30:
                raise RuntimeError(f'insufficient available RAM: {available/2**30:.3f} GiB')
            path = Path(case['input'])
            if sha256(args.contract) != digest or sha256(path) != case['sha256']:
                raise RuntimeError('contract/input changed during execution')
            command = [sys.executable,'-m','gprMax',str(path),'-cpu_precision','double','--hide-progress-bars']
            begin = time.monotonic()
            event({'id':case['id'],'status':'STARTED','command':command,'available_RAM_bytes':available})
            with path.with_suffix('.stdout.log').open('xb') as stdout, path.with_suffix('.stderr.log').open('xb') as stderr:
                result = subprocess.run(command,cwd=path.parent,stdout=stdout,stderr=stderr,
                                        timeout=contract['per_trace_timeout_s'])
            if result.returncode != 0:
                raise RuntimeError(f"{case['id']} exit {result.returncode}; see saved logs")
            import h5py
            import numpy as np
            output = path.with_suffix('.h5')
            with h5py.File(output,'r') as h:
                raw = h['rxs/rx1/Ey'][:]
                if raw.dtype != np.float64 or raw.shape != (5089,) or not np.isfinite(raw).all():
                    raise RuntimeError('raw receiver dtype/shape/finite audit failed')
                metadata = {'dt':float(h.attrs['dt']),'dtype':str(raw.dtype),'samples':len(raw),
                            'solver_version':str(h.attrs['gprMax']),
                            'tx_m':h['srcs/src1'].attrs['Position'].tolist(),
                            'rx_m':h['rxs/rx1'].attrs['Position'].tolist()}
                if metadata['solver_version'] != '4.0.0' or metadata['dt'] != 1.1793271683748422e-10:
                    raise RuntimeError('solver version/time-step differs')
            if case['id'] == 'baseline_replay/t61':
                original = ROOT/'artifacts/research_checks/2026-10-02_hs4t2d_transect/hs4t2d_t61.h5'
                with h5py.File(original,'r') as old:
                    error = float(np.max(np.abs(raw-old['rxs/rx1/Ey'][:])))
                    metadata['archived_raw_max_abs_error'] = error
                if error != 0:
                    raise RuntimeError(f'archived center trace is not reproduced exactly: {error}')
            row = {'id':case['id'],'status':'COMPLETED','wall_s':time.monotonic()-begin,
                   'input_sha256':sha256(path),'output_sha256':sha256(output),**metadata}
            event(row)
            completed.append(row)
            print(f"{len(completed)}/{len(cases)} {case['id']} {row['wall_s']:.2f}s",flush=True)
        status = 'COMPLETED'
    except Exception as exc:
        status = 'STOPPED'
        event({'status':status,'reason':str(exc),'completed':len(completed)})
        raise
    finally:
        (run/'execution_summary.json').write_text(json.dumps({'status':locals().get('status','STOPPED'),
            'completed':len(completed),'expected':len(cases),'contract_sha256':digest,
            'wall_s':time.monotonic()-start,'records':completed},indent=2)+'\n',encoding='utf-8')
    event({'status':status,'completed':len(completed)})


if __name__ == '__main__':
    main()
