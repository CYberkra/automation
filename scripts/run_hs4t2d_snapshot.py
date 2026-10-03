"""Execute one frozen V4 CUDA FP64 trace with passive field snapshots."""
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract', type=Path, required=True)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    c = json.loads(args.contract.read_text('utf-8'))
    if not args.execute or c['status'] != 'FROZEN_APPROVED' or c['max_runs'] != 1:
        raise ValueError('approved one-run contract and --execute required')
    if Path(sys.executable).resolve() != Path(c['python']).resolve():
        raise ValueError('Python differs')
    for name, digest in c['code_identities'].items():
        if sha256(name) != digest:
            raise ValueError('code differs')
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    for name, digest in c['source_identities'].items():
        if sha256(package/name) != digest:
            raise ValueError('V4 source differs')
    source, reference = Path(c['input']), Path(c['reference_raw'])
    if sha256(source) != c['input_sha256'] or sha256(reference) != c['reference_raw_sha256']:
        raise ValueError('input/reference differs')
    event = source.parent/'execution.json'
    if event.exists() or list(source.parent.glob('*.h5')):
        raise ValueError('attempt already consumed')
    available = psutil.virtual_memory().available
    if available < c['min_available_RAM_GiB']*2**30:
        raise ValueError('RAM below frozen minimum')
    command = [sys.executable, c['cache_entry'], str(source), '-gpu', '0',
               '-gpu_precision', 'double', '--hide-progress-bars']
    digest = sha256(args.contract)
    record = {'status':'STARTED', 'contract_sha256':digest, 'command':command,
              'available_RAM_bytes':available, 'started_unix_s':time.time()}
    event.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    environment = os.environ.copy()
    environment['HS4_CUDA_CACHE_LOG'] = str(source.parent/'cuda_cache.jsonl')
    started = time.monotonic()
    try:
        with (source.parent/'stdout.log').open('xb') as out, (source.parent/'stderr.log').open('xb') as err:
            solved = subprocess.run(command, cwd=source.parent, stdout=out, stderr=err,
                                    timeout=c['max_wall_s'], env=environment)
        if solved.returncode:
            raise RuntimeError('solver failed; see saved logs')
        output = source.with_suffix('.h5')
        with h5py.File(output,'r') as h, h5py.File(reference,'r') as old:
            raw = h['rxs/rx1/Ey'][:]
            if (str(h.attrs['gprMax']) != '4.0.0' or raw.dtype != np.float64
                    or not np.isfinite(raw).all() or not np.array_equal(raw,old['rxs/rx1/Ey'][:])):
                raise ValueError('raw V4/dtype/passive-observer replay audit failed')
            dt = float(h.attrs['dt'])
        snapshots = sorted(source.parent.glob('*_snaps/*.h5'))
        if len(snapshots) != len(c['snapshot_times_ns']):
            raise ValueError('snapshot count differs')
        rows = []
        for p in snapshots:
            with h5py.File(p,'r') as h:
                ey = h['Ey'][:]
                if str(h.attrs['gprMax']) != '4.0.0' or ey.dtype != np.float64 or not np.isfinite(ey).all():
                    raise ValueError('snapshot version/dtype/finite audit failed')
                requested = float(p.stem.split('_')[-1])
                actual_ns = float(h.attrs['time'])*1e9
                if abs(actual_ns-requested) > dt*1e9/2+1e-9:
                    raise ValueError('snapshot time quantization differs')
                rows.append({'file':p.relative_to(source.parent).as_posix(),'sha256':sha256(p),
                             'requested_time_ns':requested,'actual_time_ns':actual_ns,
                             'iteration':int(h.attrs['iteration']),'Ey_shape':list(ey.shape),
                             'Ey_dtype':str(ey.dtype),'Ey_peak_abs':float(np.max(np.abs(ey)))})
        if sha256(source) != c['input_sha256'] or sha256(args.contract) != digest:
            raise ValueError('frozen identities changed')
        record.update(status='COMPLETED',elapsed_s=time.monotonic()-started,raw_sha256=sha256(output),
                      raw_bit_identical_to_non_snapshot_trace=True,snapshots=rows)
    except Exception as exc:
        record.update(status='FAILED',elapsed_s=time.monotonic()-started,error=str(exc))
        raise
    finally:
        event.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
