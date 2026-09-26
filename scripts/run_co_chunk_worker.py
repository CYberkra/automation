"""Chunk worker for batch2d_v1_co: run N cases sequentially inside ONE process.

Launched (supervised) by run_approved_batch2d_co.py. Reads a job JSON:
  {device_id, precision, cases: [{run_id, input_path, input_sha256, run_directory}]}
For each case: hash-lock verify input, copy into run directory, gprMax.run() with
output into the run directory, verify output exists. Any failure exits 1 immediately;
the launcher aborts later chunks.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import traceback


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def main():
    job_path = Path(sys.argv[1])
    job = json.loads(job_path.read_text(encoding='utf-8'))
    device_id = str(job['device_id'])
    precision = job['precision']
    import gprMax  # imported once, reused across cases in this process

    for case in job['cases']:
        run_id = case['run_id']
        src = Path(case['input_path'])
        run_dir = Path(case['run_directory'])
        run_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(run_dir)
        data = src.read_bytes()
        assert sha256_bytes(data) == case['input_sha256'], f'input hash mismatch {run_id}'
        local_in = run_dir / f'{run_id}.in'
        local_in.write_bytes(data)
        out_h5 = run_dir / f'{run_id}.h5'
        try:
            sys.argv = ['gprMax']
            gprMax.run(inputfile=str(local_in), outputfile=str(out_h5),
                       gpu=device_id, gpu_precision=precision,
                       hide_progress_bars=True)
        except Exception:
            print('CHUNK_CASE_FAILED ' + json.dumps({'run_id': run_id,
                                                     'error': traceback.format_exc()}),
                  flush=True)
            sys.exit(1)
        if not out_h5.exists():
            print('CHUNK_CASE_FAILED ' + json.dumps(
                {'run_id': run_id, 'error': 'output h5 missing'}), flush=True)
            sys.exit(1)
        print('CHUNK_CASE_DONE ' + json.dumps({'run_id': run_id,
                                               'output_bytes': out_h5.stat().st_size}),
              flush=True)
    print('CHUNK_ALL_DONE', flush=True)


if __name__ == '__main__':
    main()
