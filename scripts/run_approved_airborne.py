"""Launch the single explicitly approved 3D run through the Windows supervisor.

Default is preflight only. --execute consumes a one-attempt record before launch.
No retries or subsequent models are dispatched by this script.
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

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def solver_code(source, run_id, device_id):
    """Copy a hash-locked text input or reviewed standalone official API model."""
    suffix=source.suffix
    if suffix not in ('.in','.py'):raise ValueError('Unsupported input type')
    name=run_id+suffix
    prefix=("import hashlib,pathlib,runpy,sys; "
            f"p=pathlib.Path({str(source)!r}); b=p.read_bytes(); "
            f"assert hashlib.sha256(b).hexdigest()=={sha(source)!r}; "
            f"pathlib.Path({name!r}).write_bytes(b); ")
    if suffix=='.py':
        return prefix+f"sys.argv=[{name!r}]; runpy.run_path({name!r},run_name='__main__')"
    return prefix+(f"sys.argv=['gprMax',{name!r},'-gpu',{str(device_id)!r},'-gpu_precision','double']; "
                   "runpy.run_module('gprMax',run_name='__main__')")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    gate=json.loads((ROOT/'configs/research/gprmax_v4_execution_gate.json').read_text(encoding='utf-8'))
    contract=gate['approved_execution_contract']
    run_id=contract.get('run_id','M00_x_3d')
    if not run_id.replace('_','').isalnum():raise SystemExit('Invalid run identifier')
    if not gate['approved_to_simulate'] or gate['approved_run_ids'] != [run_id]:
        raise SystemExit('No recorded authorization for this single run')
    source=ROOT/contract['input_path']
    assert sha(source)==contract['input_sha256'],'Input changed'
    assert sha(Path(__file__))==contract['launcher_sha256'],'Launcher changed'
    assert sha(ROOT/'scripts/bounded_windows_process.py')==contract['supervisor_sha256'],'Supervisor changed'
    identity=json.loads((ROOT/contract['runtime_identity_path']).read_text(encoding='utf-8'))
    assert sha(ROOT/contract['runtime_identity_path'])==contract['runtime_identity_sha256']
    assert Path(sys.executable).resolve()==Path(identity['executable']).resolve(),'Wrong Python'
    assert importlib.metadata.version('gprMax')=='4.0.0','Wrong distribution'
    for item in identity['compiled_modules_imported'].values():
        assert sha(Path(item['path']))==item['sha256'],'Native module changed'
    import psutil
    free=psutil.virtual_memory().available
    budget=gate['approved_compute_budget']
    if gate.get('execution_policy',{}).get('required_solver_backend')!='cuda':
        raise SystemExit('This launcher now requires the recorded CUDA-only policy')
    if budget['backend']!='CUDA' or budget['precision']!='double':
        raise SystemExit('Only the reviewed CUDA double calibration is supported')
    if not shutil.which('cl.exe') or not shutil.which('nvcc.exe'):
        raise SystemExit('Use the documented MSVC x64 compiler environment; no CPU fallback')
    cuda_identity_path=ROOT/contract['cuda_identity_path']
    assert sha(cuda_identity_path)==contract['cuda_identity_sha256']
    gpu_identity=json.loads(cuda_identity_path.read_text(encoding='utf-8'))
    assert gpu_identity['passed'] and gpu_identity['kernel_float64']
    assert sha(Path(gpu_identity['native_module_path']))==gpu_identity['native_module_sha256']
    import pycuda.driver as cuda
    cuda.init()
    device=cuda.Device(budget['device_id'])
    assert device.name()==gpu_identity['device_name'],'Unexpected CUDA device'
    context=device.make_context()
    try:
        gpu_free,gpu_total=cuda.mem_get_info()
    finally:
        context.pop();context.detach()
    assert gpu_free>=budget['minimum_available_VRAM_GiB']*2**30,'Insufficient available VRAM'
    assert free>=budget['minimum_available_RAM_GiB']*2**30,'Insufficient available RAM'
    attempt=ROOT/contract['attempt_record']
    if attempt.exists():raise SystemExit('Attempt already consumed; no automatic retry')
    run_dir=(ROOT/contract['run_directory']).resolve()
    assert run_dir.is_relative_to((ROOT/'artifacts/simulations').resolve())
    if run_dir.exists():raise SystemExit('Run directory already exists')
    preflight=dict(run_id=run_id,input_sha256=sha(source),available_RAM_bytes=free,
                   available_VRAM_bytes=gpu_free,total_VRAM_bytes=gpu_total,
                   backend='CUDA',device_id=budget['device_id'],precision='double',cpu_fallback=False,
                   approved=True,execute=a.execute,utc=datetime.now(timezone.utc).isoformat())
    print(json.dumps(preflight),flush=True)
    if not a.execute:return
    attempt.parent.mkdir(parents=True,exist_ok=True)
    with attempt.open('x',encoding='utf-8') as f:json.dump(preflight,f,indent=2)
    # Input is copied into the fresh supervised directory, so all official outputs
    # remain there. No shell interpretation, solver code edits, or private API.
    code=solver_code(source,run_id,budget['device_id'])
    result=supervise([sys.executable,'-u','-c',code],run_dir,
        wall_s=budget['wall_minutes']*60,memory_bytes=budget['job_commit_GiB']*2**30,
        output_bytes=budget['output_GiB']*2**30,poll_s=.25)
    print(json.dumps(result),flush=True)
    if result['reason']!='completed':raise SystemExit(1)


if __name__=='__main__':main()
