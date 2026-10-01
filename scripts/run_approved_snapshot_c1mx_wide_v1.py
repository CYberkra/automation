"""Run the snapshot_c1mx_v1 diagnostic case (no VRAM sampler) through the Windows supervisor.

Default is preflight only. --execute consumes each contract's one-attempt
record immediately before that run's launch. Contracts run in gate order
(2D_BG, 2D_TGT, 3D_BG, 3D_TGT); any failed run aborts the batch, so 3D is
never touched when a 2D smoke case fails. No retries, no extra cases.
Already-consumed attempts are reported and skipped, never rerun.
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
BATCH = 'snapshot_c1mx_wide_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solver_code(source, run_id, device_id):
    """Copy a hash-locked text input into the supervised run directory."""
    if source.suffix != '.in':
        raise ValueError('Unsupported input type')
    name = run_id + '.in'
    prefix = ("import hashlib,pathlib,runpy,sys; "
              f"p=pathlib.Path({str(source)!r}); b=p.read_bytes(); "
              f"assert hashlib.sha256(b).hexdigest()=={sha(source)!r}; "
              f"pathlib.Path({name!r}).write_bytes(b); ")
    return prefix + (f"sys.argv=['gprMax',{name!r},'-gpu',{str(device_id)!r},'-gpu_precision','double']; "
                     "runpy.run_module('gprMax',run_name='__main__')")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute', action='store_true')
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
        assert sha(ROOT / 'scripts/bounded_windows_process.py') == c['supervisor_sha256'], 'Supervisor changed'
        assert sha(ROOT / c['input_path']) == c['input_sha256'], 'Input changed'
    if gate.get('execution_policy', {}).get('required_solver_backend') != 'cuda':
        raise SystemExit('This launcher now requires the recorded CUDA-only policy')
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
    for item in identity['compiled_modules_imported'].values():
        assert sha(Path(item['path'])) == item['sha256'], 'Native module changed'
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
    report = []
    for c in contracts:
        run_id = c['run_id']
        attempt = ROOT / c['attempt_record']
        if attempt.exists():
            report.append(dict(run_id=run_id, skipped='attempt already consumed; no rerun'))
            print(json.dumps(report[-1]), flush=True)
            continue
        context = device.make_context()
        try:
            gpu_free, gpu_total = cuda.mem_get_info()
        finally:
            context.pop()
            context.detach()
        free = psutil.virtual_memory().available
        assert gpu_free >= budget['minimum_available_VRAM_GiB'] * 2**30, 'Insufficient available VRAM'
        assert free >= budget['minimum_available_RAM_GiB'] * 2**30, 'Insufficient available RAM'
        run_dir = (ROOT / c['run_directory']).resolve()
        assert run_dir.is_relative_to((ROOT / 'artifacts/simulations').resolve())
        if run_dir.exists():
            raise SystemExit('Run directory already exists')
        preflight = dict(run_id=run_id, batch=BATCH, input_sha256=sha(ROOT / c['input_path']),
                         available_RAM_bytes=free, available_VRAM_bytes=gpu_free,
                         total_VRAM_bytes=gpu_total, backend='CUDA', device_id=budget['device_id'],
                         precision='double', cpu_fallback=False, approved=True, execute=a.execute,
                         utc=datetime.now(timezone.utc).isoformat())
        print(json.dumps(preflight), flush=True)
        if not a.execute:
            continue
        attempt.parent.mkdir(parents=True, exist_ok=True)
        with attempt.open('x', encoding='utf-8') as f:
            json.dump(preflight, f, indent=2)
        result = supervise([sys.executable, '-u', '-c', solver_code(ROOT / c['input_path'], run_id, budget['device_id'])],
                           run_dir, wall_s=budget['wall_minutes'] * 60,
                           memory_bytes=budget['job_commit_GiB'] * 2**30,
                           output_bytes=budget['output_GiB'] * 2**30, poll_s=.25)
        outcome = dict(run_id=run_id, **result)
        report.append(outcome)
        print(json.dumps(outcome), flush=True)
        if result['reason'] != 'completed':
            (ROOT / c['attempt_record']).write_text(json.dumps(dict(preflight, outcome=result), indent=2), encoding='utf-8')
            raise SystemExit(f'{run_id} did not complete ({result["reason"]}); batch aborted, later attempts untouched')
    print(json.dumps({'batch': BATCH, 'runs_reported': len(report)}), flush=True)


if __name__ == '__main__':
    main()
