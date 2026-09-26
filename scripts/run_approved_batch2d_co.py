"""Run the frozen batch2d_v1_co batch (chunked in-process execution) through the supervisor.

Default is preflight only. --execute consumes each case's one-attempt record immediately
before its chunk launches. Chunks run serially in frozen order (BG all, then A0); within a
chunk a single supervised process runs all cases via gprMax.run(). Any chunk failure aborts
the batch: later chunks untouched, unfinished cases of the failed chunk marked in their
attempt records. No retries.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import shutil
from datetime import datetime, timezone

from bounded_windows_process import supervise

ROOT = Path(__file__).resolve().parents[1]
BATCH = 'batch2d_v1_co'
CO_DIR = ROOT / 'configs/research/batch2d_v1_co'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--chunks', default=None,
                   help='comma-separated chunk_ids to run (default: all chunks in frozen order)')
    a = p.parse_args()
    gate = json.loads((ROOT / 'configs/research/gprmax_v4_execution_gate.json').read_text(encoding='utf-8'))
    contracts = gate.get('approved_execution_contracts')
    if gate.get('approved_to_simulate') is not True or not isinstance(contracts, list) or not contracts:
        raise SystemExit('No pending approved batch; consumed attempts cannot be reused')
    if gate.get('batch_id') != BATCH:
        raise SystemExit('Gate batch_id mismatch')
    budget = gate['approved_compute_budget']
    if gate['approved_run_ids'] != [c['run_id'] for c in contracts]:
        raise SystemExit('approved_run_ids must list the batch contracts in execution order')
    if budget['max_fdtd_runs'] < len(contracts):
        raise SystemExit('Batch exceeds max_fdtd_runs')
    for c in contracts:
        run_id = c['run_id']
        if not run_id.replace('_', '').replace('-', '').replace('.', '').isalnum():
            raise SystemExit('Invalid run identifier')
        assert sha(Path(__file__)) == c['launcher_sha256'], 'Launcher changed'
        assert sha(ROOT / 'scripts/run_co_chunk_worker.py') == c['worker_sha256'], 'Worker changed'
        assert sha(ROOT / 'scripts/bounded_windows_process.py') == c['supervisor_sha256'], 'Supervisor changed'
        assert sha(ROOT / c['input_path']) == c['input_sha256'], 'Input changed'
    if gate.get('execution_policy', {}).get('required_solver_backend') != 'cuda':
        raise SystemExit('This launcher requires the recorded CUDA-only policy')
    if budget['backend'] != 'CUDA' or budget['precision'] != 'double':
        raise SystemExit('Only the reviewed CUDA double calibration is supported')
    if not shutil.which('cl.exe') or not shutil.which('nvcc.exe'):
        raise SystemExit('Use the documented MSVC x64 compiler environment; no CPU fallback')
    first = contracts[0]
    identity_path = ROOT / first['runtime_identity_path']
    assert sha(identity_path) == first['runtime_identity_sha256']
    identity = json.loads(identity_path.read_text(encoding='utf-8'))
    assert Path(sys.executable).resolve() == Path(identity['executable']).resolve(), 'Wrong Python'
    assert importlib.metadata.version('gprMax') == '4.0.0', 'Wrong distribution'
    cuda_identity_path = ROOT / first['cuda_identity_path']
    assert sha(cuda_identity_path) == first['cuda_identity_sha256']
    gpu_identity = json.loads(cuda_identity_path.read_text(encoding='utf-8'))
    assert gpu_identity['passed'] and gpu_identity['kernel_float64']
    assert sha(Path(gpu_identity['native_module_path'])) == gpu_identity['native_module_sha256']
    import psutil
    import pycuda.driver as cuda
    cuda.init()
    device = cuda.Device(budget['device_id'])
    assert device.name() == gpu_identity['device_name'], 'Unexpected CUDA device'

    chunks_doc = json.loads((CO_DIR / 'chunks.json').read_text(encoding='utf-8'))
    chunks = chunks_doc['chunks']
    contract_by_id = {c['run_id']: c for c in contracts}
    ordered_ids = [rid for ch in chunks for rid in ch['run_ids']]
    assert ordered_ids == gate['approved_run_ids'], 'chunks.json order must match approved_run_ids'
    if a.chunks:
        wanted = [s.strip() for s in a.chunks.split(',') if s.strip()]
        known = {ch['chunk_id'] for ch in chunks}
        assert all(w in known for w in wanted), f'unknown chunk id: {wanted}'
        chunks = [ch for ch in chunks if ch['chunk_id'] in wanted]

    report = []
    for ch in chunks:
        chunk_id = ch['chunk_id']
        chunk_contracts = [contract_by_id[rid] for rid in ch['run_ids']]
        pending = [c for c in chunk_contracts
                   if not (ROOT / c['attempt_record']).exists()]
        if not pending:
            print(json.dumps({'chunk': chunk_id, 'skipped': 'attempts already consumed'}),
                  flush=True)
            continue
        if len(pending) != len(chunk_contracts):
            raise SystemExit(f'{chunk_id} partially consumed; refusing to run')
        context = device.make_context()
        try:
            gpu_free, gpu_total = cuda.mem_get_info()
        finally:
            context.pop()
            context.detach()
        free = psutil.virtual_memory().available
        assert gpu_free >= budget['minimum_available_VRAM_GiB'] * 2**30, 'Insufficient available VRAM'
        assert free >= budget['minimum_available_RAM_GiB'] * 2**30, 'Insufficient available RAM'
        chunk_dir = (ROOT / 'artifacts/simulations' / f'2026-09-26_{chunk_id}').resolve()
        assert chunk_dir.is_relative_to((ROOT / 'artifacts/simulations').resolve())
        if chunk_dir.exists():
            raise SystemExit(f'{chunk_dir} already exists')
        job = {'device_id': budget['device_id'], 'precision': budget['precision'],
               'cases': []}
        for c in pending:
            run_dir = (ROOT / c['run_directory']).resolve()
            assert run_dir.is_relative_to((ROOT / 'artifacts/simulations').resolve())
            if run_dir.exists():
                raise SystemExit(f'{run_dir} already exists')
            preflight = dict(run_id=c['run_id'], batch=BATCH, chunk_id=chunk_id,
                             input_sha256=sha(ROOT / c['input_path']),
                             available_RAM_bytes=free, available_VRAM_bytes=gpu_free,
                             total_VRAM_bytes=gpu_total, backend='CUDA',
                             device_id=budget['device_id'], precision='double',
                             cpu_fallback=False, approved=True, execute=a.execute,
                             utc=datetime.now(timezone.utc).isoformat())
            print(json.dumps({'chunk': chunk_id, 'preflight': preflight}), flush=True)
            if not a.execute:
                continue
            attempt_path = ROOT / c['attempt_record']
            attempt_path.parent.mkdir(parents=True, exist_ok=True)
            with attempt_path.open('x', encoding='utf-8') as f:
                json.dump(preflight, f, indent=2)
            job['cases'].append({
                'run_id': c['run_id'], 'input_path': str(ROOT / c['input_path']),
                'input_sha256': c['input_sha256'],
                'run_directory': str(ROOT / c['run_directory'])})
        if not a.execute:
            continue
        chunk_dir.mkdir(parents=True)
        job_path = chunk_dir / 'chunk_job.json'
        job_path.write_text(json.dumps(job, indent=1), encoding='utf-8')
        wall_s = 300 + 150 * len(job['cases'])
        result = supervise([sys.executable, '-u',
                            str(ROOT / 'scripts/run_co_chunk_worker.py'), str(job_path)],
                           chunk_dir, wall_s=wall_s,
                           memory_bytes=budget['job_commit_GiB'] * 2**30,
                           output_bytes=budget['output_GiB'] * 2**30, poll_s=.25)
        outcome = dict(chunk=chunk_id, **result)
        print(json.dumps(outcome), flush=True)
        if result['reason'] != 'completed':
            for c in pending:
                attempt_path = ROOT / c['attempt_record']
                if attempt_path.exists() and 'outcome' not in json.loads(
                        attempt_path.read_text(encoding='utf-8')):
                    attempt_path.write_text(json.dumps(
                        dict(json.loads(attempt_path.read_text(encoding='utf-8')),
                             outcome=result, note='chunk aborted before this case ran'),
                        indent=2), encoding='utf-8')
            raise SystemExit(f'{chunk_id} did not complete ({result["reason"]}); '
                             'batch aborted, later chunks untouched')
        for c in pending:
            run_dir = ROOT / c['run_directory']
            h5 = run_dir / f"{c['run_id']}.h5"
            if not h5.exists():
                raise SystemExit(f"{c['run_id']} output missing after chunk completion")
            report.append(dict(run_id=c['run_id'], chunk=chunk_id, h5_bytes=h5.stat().st_size))
            archive = ROOT / 'artifacts/research_checks' / f"2026-09-26_{c['run_id']}"
            archive.mkdir(parents=True, exist_ok=True)
            for item in run_dir.iterdir():
                if item.is_file():
                    shutil.copy2(item, archive / item.name)
    print(json.dumps({'batch': BATCH, 'runs_completed': len(report)}), flush=True)


if __name__ == '__main__':
    main()
