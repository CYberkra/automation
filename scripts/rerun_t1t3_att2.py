"""att2 rerun for the 40 batch2d_slope_t1t3_co cases killed mid-solver.

Accident record: during chunked execution (300 s tool-call windows) the last case
started inside a chunk window was killed by the chunk timeout WHILE THE SOLVER WAS
RUNNING; the attempt had already been consumed at launch, so the batch runner skips
these runs forever and their h5 was never written. Diagnosis per run dir: stdout.log
progress bars stop at 5-13%, no completion line, no h5. 40/396 cases affected, spread
across all 12 mothers (each chunk boundary claimed ~1 case).

Recovery follows the batch2d_slope_t2 NC att2 precedent (decision_log 2026-09-27):
the SAME hash-locked input is rerun under the SAME authorization (user "按建议开始",
gate batch2d_slope_t1t3_co) via a NEW run directory suffix _att2 and a NEW attempt
record suffix _attempt_att2; the frozen batch runner/launcher is NOT modified. This
script is operational tooling, not part of the frozen contract; its sha is recorded
in each att2 attempt record.

Discipline: input sha256 verified against cases.json before every run; run dir must
not exist; one attempt per run (att2), no retries; same caps as the gate budget
(wall 20 min, job commit 6 GiB, output 2 GiB); CUDA double device 0; serial in
cases.json order; abort on first failure. Default is preflight; --execute consumes
the att2 attempt immediately before launch. Already-consumed att2 attempts are
reported and skipped, never rerun.
"""
import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BATCH = 'batch2d_slope_t1t3_co'
CASES = ROOT / 'configs/research/batch2d_slope_t1t3_co/cases.json'
DATE = '2026-09-29'
BUDGET = {'wall_s': 1200, 'memory_bytes': 6 * 2 ** 30, 'output_bytes': 2 * 2 ** 30,
          'device_id': 0}

# Environment fix (same as run_batch2d_slope_t1t3_co_gpu.cmd): without the CUDA
# v13.3 bin directories on PATH, pycuda cannot invoke nvcc and the solver exits
# within ~2 s (nonzero_exit). Two att2 attempts (t12, t23) were consumed by this
# launch-environment failure before the fix; they are rerun under a further suffix.
_CUDA = Path(r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3')
os.environ['PATH'] = (str(_CUDA / 'bin') + os.pathsep + str(_CUDA / 'bin' / 'x64')
                      + os.pathsep + os.environ['PATH'])

sys.path.insert(0, str(ROOT / 'scripts'))
from bounded_windows_process import supervise  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solver_code(source, run_id, device_id):
    if source.suffix != '.in':
        raise ValueError('Unsupported input type')
    name = run_id + '.in'
    prefix = ("import hashlib,pathlib,runpy,sys; "
              f"p=pathlib.Path({str(source)!r}); b=p.read_bytes(); "
              f"assert hashlib.sha256(b).hexdigest()=={sha(source)!r}; "
              f"pathlib.Path({name!r}).write_bytes(b); ")
    return prefix + (f"sys.argv=['gprMax',{name!r},'-gpu',{str(device_id)!r},'-gpu_precision','double']; "
                     "runpy.run_module('gprMax',run_name='__main__')")


def missing_runs():
    cases = json.loads(CASES.read_text(encoding='utf-8'))['cases']
    out = []
    for c in cases:
        rid = c['run_id']
        ok = any((ROOT / 'artifacts/simulations' / f'{DATE}_{rid}{sfx}' / f'{rid}.h5').exists()
                 for sfx in ('', '_att2', '_att3', '_att4', '_att5'))
        if not ok:
            out.append(c)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--suffix', default='att2',
                   help='attempt/run-dir suffix (att3 for the two env-failure runs)')
    a = p.parse_args()
    sfx = a.suffix
    todo = missing_runs()
    print(json.dumps({'att2_pending': len(todo), 'execute': a.execute, 'suffix': sfx,
                      'launcher_sha256': sha(Path(__file__))}), flush=True)
    report = []
    for c in todo:
        run_id = c['run_id']
        att2 = ROOT / 'artifacts/research_checks' / f'{DATE}_{run_id}_attempt_{sfx}.json'
        run_dir = ROOT / 'artifacts/simulations' / f'{DATE}_{run_id}_{sfx}'
        if att2.exists():
            report.append(dict(run_id=run_id, skipped=f'{sfx} attempt already consumed'))
            print(json.dumps(report[-1]), flush=True)
            continue
        if sfx != 'att2':
            prior = ROOT / 'artifacts/research_checks' / f'{DATE}_{run_id}_attempt_att2.json'
            assert prior.exists(), f'{run_id}: {sfx} requires a consumed att2 attempt'
        src = ROOT / 'configs/research/batch2d_slope_t1t3_co' / c['file']
        assert sha(src) == c['input_sha256'], f'input hash mismatch {run_id}'
        assert not run_dir.exists(), f'{sfx} run directory already exists {run_dir}'
        preflight = dict(run_id=run_id, batch=BATCH, attempt=sfx,
                         input_sha256=c['input_sha256'],
                         basis='chunk-timeout mid-solver kill; same hash-locked input, '
                               'same authorization (gate note 2026-09-29)',
                         launcher_sha256=sha(Path(__file__)),
                         backend='CUDA', device_id=BUDGET['device_id'], precision='double',
                         approved=True, execute=a.execute,
                         utc=datetime.now(timezone.utc).isoformat())
        print(json.dumps(preflight), flush=True)
        if not a.execute:
            continue
        att2.parent.mkdir(parents=True, exist_ok=True)
        with att2.open('x', encoding='utf-8') as f:
            json.dump(preflight, f, indent=2)
        result = supervise([sys.executable, '-u', '-c', solver_code(src, run_id, BUDGET['device_id'])],
                           run_dir, wall_s=BUDGET['wall_s'],
                           memory_bytes=BUDGET['memory_bytes'],
                           output_bytes=BUDGET['output_bytes'], poll_s=.25)
        outcome = dict(run_id=run_id, **result)
        report.append(outcome)
        print(json.dumps(outcome), flush=True)
        if result['reason'] != 'completed':
            att2.write_text(json.dumps(dict(preflight, outcome=result), indent=2), encoding='utf-8')
            raise SystemExit(f'{run_id} {sfx} did not complete ({result["reason"]}); aborted')
    print(json.dumps({'batch': f'{BATCH}_{sfx}', 'runs_reported': len(report)}), flush=True)


if __name__ == '__main__':
    main()
