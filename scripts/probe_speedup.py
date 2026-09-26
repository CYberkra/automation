"""Timing probes for gprMax acceleration options (local checks, not gated runs).

Modes:
  single_double <in> <tag>   one model, CUDA double, one process
  seq_double <in> <tag>      two models sequentially in ONE process (option 3)
  single_single <in> <tag>   one model, CUDA single precision (option 1)
  concurrent <in> <tag>      two processes at once, each one model double (option 2)

Prints one JSON line per phase with wall seconds; h5 written under
artifacts/local_checks/2026-09-26_speedup_probes/<tag>/.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / 'artifacts/local_checks/2026-09-26_speedup_probes'


def run_once(in_path, precision, outdir, stem):
    import gprMax
    sys.argv = ['gprMax']
    gprMax.run(inputfile=str(in_path), outputfile=str(outdir / f'{stem}.h5'),
               gpu='0', gpu_precision=precision,
               hide_progress_bars=True)


def main():
    mode, in_path, tag = sys.argv[1], Path(sys.argv[2]).resolve(), sys.argv[3]
    outdir = SCRATCH / tag
    outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(outdir)
    rec = dict(mode=mode, tag=tag, in_file=in_path.name)

    if mode == 'single_double':
        t0 = time.perf_counter()
        run_once(in_path, 'double', outdir, tag)
        rec['wall_s'] = time.perf_counter() - t0
    elif mode == 'seq_double':
        walls = []
        for i in (1, 2):
            t0 = time.perf_counter()
            run_once(in_path, 'double', outdir, f'{tag}_{i}')
            walls.append(time.perf_counter() - t0)
        rec['wall_each_s'] = walls
        rec['wall_total_s'] = sum(walls)
    elif mode == 'single_single':
        t0 = time.perf_counter()
        run_once(in_path, 'single', outdir, tag)
        rec['wall_s'] = time.perf_counter() - t0
    elif mode == 'concurrent':
        t0 = time.perf_counter()
        procs = [subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                                   'single_double', str(in_path), f'{tag}_w{i}'],
                                  cwd=outdir)
                 for i in (1, 2)]
        codes = [p.wait() for p in procs]
        rec['wall_s'] = time.perf_counter() - t0
        rec['exit_codes'] = codes
    else:
        raise SystemExit(f'unknown mode {mode}')
    print('PROBE_RESULT ' + json.dumps(rec), flush=True)


if __name__ == '__main__':
    main()
