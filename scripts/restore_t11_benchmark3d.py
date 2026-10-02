"""Restore B3D5CM-C3mR2-BG-CO13-t11 (benchmark3d_r2_co).

Basis (user sign-off): attempt t11 was consumed 2026-10-02 ~02:56 but the
supervised solver was killed by an external process-terminate at ~03:00
(DBG_TERMINATE_PROCESS during agent-runtime recovery), leaving no h5.
User instruction 2026-10-02 09:25 "恢复t11契约" explicitly authorises this
out-of-contract restore. Zero-retry discipline otherwise unchanged; this
script only touches the t11 contract entry.

Reuses the frozen gate contract and performs the same hash/identity checks
as run_approved_benchmark3d_r2.py before consuming the fresh attempt record.
"""
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import sys

from bounded_windows_process import supervise

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
from run_approved_benchmark3d_r2 import solver_code  # noqa: E402

RUN_ID = 'B3D5CM-C3mR2-BG-CO13-t11'
BATCH = 'benchmark3d_r2_co'
RESTORE_BASIS = ('user sign-off 2026-10-02 09:25 (restore t11 contract); '
                 'original attempt consumed but solver killed externally '
                 '~03:00 before producing h5; no rerun of other contracts')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    gate = json.loads((ROOT / 'configs/research/gprmax_v4_execution_gate.json').read_text(encoding='utf-8'))
    contracts = gate.get('approved_execution_contracts')
    contract = next((c for c in contracts if c['run_id'] == RUN_ID), None)
    if contract is None:
        raise SystemExit('t11 contract not found in gate')
    if gate.get('batch_id') != BATCH or gate.get('approved_to_simulate') is not True:
        raise SystemExit('gate state mismatch')
    budget = gate['approved_compute_budget']

    # same integrity assertions as the batch launcher
    assert sha(ROOT / 'scripts/run_approved_benchmark3d_r2.py') == contract['launcher_sha256'], 'Launcher changed'
    assert sha(ROOT / 'scripts/bounded_windows_process.py') == contract['supervisor_sha256'], 'Supervisor changed'
    assert sha(ROOT / contract['input_path']) == contract['input_sha256'], 'Input changed'
    if gate.get('execution_policy', {}).get('required_solver_backend') != 'cuda':
        raise SystemExit('CUDA-only policy required')
    identity_path = ROOT / contract['runtime_identity_path']
    assert sha(identity_path) == contract['runtime_identity_sha256']
    identity = json.loads(identity_path.read_text(encoding='utf-8'))
    assert Path(sys.executable).resolve() == Path(identity['executable']).resolve(), 'Wrong Python'
    assert importlib.metadata.version('gprMax') == '4.0.0'
    for item in identity['compiled_modules_imported'].values():
        assert sha(Path(item['path'])) == item['sha256']
    cuda_identity_path = ROOT / contract['cuda_identity_path']
    assert sha(cuda_identity_path) == contract['cuda_identity_sha256']
    gpu_identity = json.loads(cuda_identity_path.read_text(encoding='utf-8'))
    assert gpu_identity['passed'] and gpu_identity['kernel_float64']
    assert sha(Path(gpu_identity['native_module_path'])) == gpu_identity['native_module_sha256']
    if not shutil.which('cl.exe') or not shutil.which('nvcc.exe'):
        raise SystemExit('MSVC+nvcc environment required; run via run_benchmark3d_r2_gpu.cmd')

    import psutil
    import pycuda.driver as cuda
    cuda.init()
    device = cuda.Device(budget['device_id'])
    assert device.name() == gpu_identity['device_name'], 'Unexpected CUDA device'
    context = device.make_context()
    try:
        gpu_free, gpu_total = cuda.mem_get_info()
    finally:
        context.pop()
        context.detach()
    free = psutil.virtual_memory().available
    assert gpu_free >= budget['minimum_available_VRAM_GiB'] * 2**30, 'Insufficient VRAM'
    assert free >= budget['minimum_available_RAM_GiB'] * 2**30, 'Insufficient RAM'

    run_dir = (ROOT / contract['run_directory']).resolve()
    assert run_dir.is_relative_to((ROOT / 'artifacts/simulations').resolve())
    if run_dir.exists():
        shutil.rmtree(run_dir)  # burned partial dir from the killed run; input hash re-verified above
    preflight = dict(run_id=RUN_ID, batch=BATCH, input_sha256=sha(ROOT / contract['input_path']),
                     available_RAM_bytes=free, available_VRAM_bytes=gpu_free,
                     total_VRAM_bytes=gpu_total, backend='CUDA', device_id=budget['device_id'],
                     precision='double', cpu_fallback=False, approved=True, execute=True,
                     restore_basis=RESTORE_BASIS,
                     utc=datetime.now(timezone.utc).isoformat())
    print(json.dumps(preflight), flush=True)
    attempt = ROOT / contract['attempt_record']
    attempt.parent.mkdir(parents=True, exist_ok=True)
    if attempt.exists():
        burned = json.loads(attempt.read_text(encoding='utf-8'))
        assert 'outcome' not in burned, 'attempt already has outcome; refuse to restore'
        archived = attempt.with_name(attempt.stem + '_burned_20261002_0300kill' + attempt.suffix)
        attempt.rename(archived)
        print(json.dumps({'archived_burned_attempt': str(archived)}), flush=True)
    with attempt.open('x', encoding='utf-8') as f:
        json.dump(preflight, f, indent=2)
    result = supervise([sys.executable, '-u', '-c', solver_code(ROOT / contract['input_path'], RUN_ID, budget['device_id'])],
                       run_dir, wall_s=budget['wall_minutes'] * 60,
                       memory_bytes=budget['job_commit_GiB'] * 2**30,
                       output_bytes=budget['output_GiB'] * 2**30, poll_s=.25)
    outcome = dict(run_id=RUN_ID, **result)
    print(json.dumps(outcome), flush=True)
    if result['reason'] != 'completed':
        (ROOT / contract['attempt_record']).write_text(json.dumps(dict(preflight, outcome=result), indent=2), encoding='utf-8')
        raise SystemExit(f'{RUN_ID} did not complete ({result["reason"]})')


if __name__ == '__main__':
    main()
