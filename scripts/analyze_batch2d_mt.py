"""Post-run CPU audit of the batch2d_v1_mt multi-trace batch (22 cases x 33 rx).

Verifies, for every approved contract:
  * 33 receiver groups named mt01..mt33, y = 12.65 + 0.25*k m, z = 45 m,
    GridPosition (0, 506+10*k, 1800), Ex-only, 20352 float64 samples.
  * Anchor trace mt17 (y 16.65 m, grid y 666) is bit-for-bit identical to the
    archived single-trace mother model h5 (max_abs_difference must be 0.0).
  * NC (no-cavity) role: all 33 traces bit-for-bit equal to the same-family
    BG case, i.e. the NC-BG difference is exactly zero on every trace.
  * Supervision outcome per case (reason/exit/wall/peak job commit).

Deterministic output: no timestamps, fixed key order; running twice must give
byte-identical results.json (r1/r2 reproducibility check).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
N_TRACES = 33
N_SAMPLES = 20352
ANCHOR_K = 17  # 1-based trace index of the anchor (mt17)
ARCHIVE_DIR = 'artifacts/research_checks/2026-09-26_{mother}/{mother}.h5'


def check_layout(path, run_id):
    """Validate receiver array layout of one MT output file."""
    issues = []
    with h5py.File(path, 'r') as f:
        if int(f.attrs['nrx']) != N_TRACES:
            issues.append(f'nrx={f.attrs["nrx"]} != {N_TRACES}')
        for k in range(1, N_TRACES + 1):
            gname = f'rxs/rx{k}'
            if gname not in f:
                issues.append(f'missing {gname}')
                continue
            g = f[gname]
            name = g.attrs.get('Name')
            if name != f'mt{k:02d}':
                issues.append(f'{gname} Name={name!r}')
            pos = np.asarray(g.attrs.get('Position'), dtype=float)
            y_expected = 12.65 + 0.25 * (k - 1)
            if pos.shape != (3,) or abs(pos[0]) > 1e-12 \
                    or abs(pos[1] - y_expected) > 1e-9 or abs(pos[2] - 45.0) > 1e-12:
                issues.append(f'{gname} Position={pos.tolist()} expected y={y_expected}')
            grid = np.asarray(g.attrs.get('GridPosition'), dtype=int)
            grid_expected = [0, 506 + 10 * (k - 1), 1800]
            if grid.tolist() != grid_expected:
                issues.append(f'{gname} GridPosition={grid.tolist()} expected {grid_expected}')
            keys = sorted(g.keys())
            if keys != ['Ex']:
                issues.append(f'{gname} datasets={keys}')
                continue
            dset = g['Ex']
            if dset.shape != (N_SAMPLES,) or dset.dtype != np.float64:
                issues.append(f'{gname}/Ex shape={dset.shape} dtype={dset.dtype}')
    return issues


def anchor_diff(mt_path, mother):
    """max |mt17 - archived single-trace| (0.0 iff bit-for-bit equal)."""
    arch_path = ROOT / ARCHIVE_DIR.format(mother=mother)
    with h5py.File(mt_path, 'r') as f:
        a = f[f'rxs/rx{ANCHOR_K}/Ex'][:]
    with h5py.File(arch_path, 'r') as f:
        b = f['rxs/rx1/Ex'][:]
    if a.shape != b.shape:
        return None, f'shape mismatch mt{a.shape} vs archive{b.shape}'
    if not np.array_equal(a, b):
        return float(np.max(np.abs(a - b))), 'not bit-for-bit equal'
    return 0.0, None


def bg_pair_zero_check(nc_path, family):
    """NC role: traces must equal the same-family BG case bit-for-bit
    (both are the no-cavity background; NC-BG difference is exactly zero).
    Returns (max_abs_diff, mismatch_count)."""
    bg_path = ROOT / 'artifacts/simulations' / f'2026-09-26_B2D-{family}-BG-MT33' / f'B2D-{family}-BG-MT33.h5'
    worst = 0.0
    mismatched = 0
    with h5py.File(nc_path, 'r') as fa, h5py.File(bg_path, 'r') as fb:
        for k in range(1, N_TRACES + 1):
            a = fa[f'rxs/rx{k}/Ex'][:]
            b = fb[f'rxs/rx{k}/Ex'][:]
            if not np.array_equal(a, b):
                mismatched += 1
                worst = max(worst, float(np.max(np.abs(a - b))))
    return worst, mismatched


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path,
                   default=ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_mt_analysis/results.json')
    a = p.parse_args()

    gate = json.loads((ROOT / 'configs/research/gprmax_v4_execution_gate.json').read_text(encoding='utf-8'))
    contracts = gate['approved_execution_contracts']
    assert gate['batch_id'] == 'batch2d_v1_mt'

    cases = []
    for c in contracts:
        run_id = c['run_id']
        mother = run_id[:-5]  # strip '-MT33'
        run_dir = ROOT / c['run_directory']
        h5_path = run_dir / f'{run_id}.h5'
        supervision = json.loads((run_dir / 'supervision.json').read_text(encoding='utf-8'))

        case = {
            'run_id': run_id,
            'mother': mother,
            'supervision_reason': supervision['reason'],
            'exit_code': supervision['exit_code'],
            'wall_s': round(supervision['wall_s'], 3),
            'peak_job_commit_bytes': supervision['peak_job_commit_bytes'],
        }
        issues = check_layout(h5_path, run_id)
        diff, anchor_issue = anchor_diff(h5_path, mother)
        case['anchor_max_abs_difference'] = diff
        case['anchor_bit_for_bit'] = (diff == 0.0)
        if anchor_issue:
            issues.append(anchor_issue)
        if '-NC-' in run_id:
            family = mother.split('-')[1]  # B2D-C1m-D10m-NC -> C1m
            worst, mismatched = bg_pair_zero_check(h5_path, family)
            case['nc_minus_bg_max_abs'] = worst
            case['nc_minus_bg_mismatched_traces'] = mismatched
            if mismatched:
                issues.append(f'NC differs from BG on {mismatched} traces (max {worst})')
        case['layout_ok'] = not issues
        case['issues'] = issues
        cases.append(case)

    result = {
        'batch': 'batch2d_v1_mt',
        'n_cases': len(cases),
        'n_traces_per_case': N_TRACES,
        'n_samples_per_trace': N_SAMPLES,
        'anchor_trace': f'mt{ANCHOR_K:02d}',
        'all_cases_pass': all(c['layout_ok'] and c['anchor_bit_for_bit']
                              and c['supervision_reason'] == 'completed' and c['exit_code'] == 0
                              for c in cases),
        'cases': cases,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'all_cases_pass': result['all_cases_pass'], 'n_cases': result['n_cases'],
                      'output': str(a.output)}, sort_keys=True))
    if not result['all_cases_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
