"""Post-run acceptance analysis for batch2d_slope_t1t3_co (396 runs).

Read-only over artifacts/simulations. Checks:
1. All 396 output h5 files exist with the expected 20352-sample Ex trace.
2. NC - BG zero-difference per trace (4 tier/TZ pairs x 33 traces): the nullcontrast
   target box is numerically identical to rock, so NC must equal its BG twin
   bit-for-bit (t2/v1 lineage property, re-verified on the new geometries).
3. TGT - BG sanity: difference energy concentrated in the target-reflection band,
   reported as relative L2 (mechanism evidence, not detectability).
4. Wall/memory summary from the runner chunk log + supervision records.
5. tz_pair: S?TZX-BG vs S?X-BG are DIFFERENT materials (tzone present) - reported
   as expected-nonzero, only recorded.

Deterministic: no timestamps; r1/r2 byte-identical.
"""
import json
import re
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIMS = ROOT / 'artifacts' / 'simulations'
CHUNK_LOG = ROOT / 'artifacts' / 'local_checks' / 't1t3_chunk.log'
OUT = ROOT / 'artifacts' / 'research_checks' / 'batch2d_slope_t1t3_co_analysis'
TIERS = ['S3X', 'S3TZX', 'S1X', 'S1TZX']
SCENES = {'BG': 'BG', 'NC': 'NC', 'TGT': 'D10m-W4m-T0.5m-E20-S0.02'}
N_TRACES = 33
DATE = '2026-09-29'


def load_ex(run_id):
    for sfx in ('', '_att2', '_att3', '_att4', '_att5'):
        p = SIMS / f'{DATE}_{run_id}{sfx}' / f'{run_id}.h5'
        if p.exists():
            with h5py.File(p, 'r') as f:
                path = 'rxs/rx1/Ex'
                if path not in f:
                    keys = list(f.keys())
                    raise KeyError(f'{path} not in {run_id}: {keys}')
                return np.asarray(f[path]), sfx
    return None, None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {'schema': 'batch2d_slope_t1t3_co_analysis/1', 'date': DATE,
           'batch': 'batch2d_slope_t1t3_co', 'n_expected': 396,
           'h5_presence': [], 'nc_bg_zero_diff': [], 'tgt_bg': [],
           'tz_material_pair': []}
    n_ok = 0
    for tier in TIERS:
        for k in range(1, N_TRACES + 1):
            tid = f't{k:02d}'
            bg_id = f'B2D-C3m{tier}-{SCENES["BG"]}-CO33-{tid}'
            nc_id = f'B2D-C3m{tier}-{SCENES["NC"]}-CO33-{tid}'
            tg_id = f'B2D-C3m{tier}-{SCENES["TGT"]}-CO33-{tid}'
            bg, _s1 = load_ex(bg_id); nc, _s2 = load_ex(nc_id); tg, _s3 = load_ex(tg_id)
            for rid, arr, sfx in ((bg_id, bg, _s1), (nc_id, nc, _s2), (tg_id, tg, _s3)):
                present = arr is not None and arr.shape[0] == 20352
                res['h5_presence'].append({'run_id': rid, 'present': bool(present),
                                           'source': sfx or 'primary'})
                n_ok += bool(present)
            if bg is not None and nc is not None:
                d = float(np.max(np.abs(nc - bg)))
                res['nc_bg_zero_diff'].append(
                    {'tier': tier, 'trace': tid, 'max_abs_diff': d, 'zero': d == 0.0})
            if bg is not None and tg is not None:
                diff = tg - bg
                rel = float(np.linalg.norm(diff) / max(np.linalg.norm(bg), 1e-30))
                res['tgt_bg'].append({'tier': tier, 'trace': tid,
                                      'rel_L2_tgt_minus_bg': rel})
    # tz pairs: tzone variant BG must differ from non-tz BG (materials differ)
    for base, tzv in (('S3X', 'S3TZX'), ('S1X', 'S1TZX')):
        diffs = []
        for k in range(1, N_TRACES + 1):
            tid = f't{k:02d}'
            a, _ = load_ex(f'B2D-C3m{base}-BG-CO33-{tid}')
            b, _ = load_ex(f'B2D-C3m{tzv}-BG-CO33-{tid}')
            if a is not None and b is not None:
                diffs.append(float(np.max(np.abs(b - a))))
        res['tz_material_pair'].append(
            {'pair': f'{tzv} vs {base}', 'max_abs_diff_range':
             [min(diffs), max(diffs)] if diffs else None})

    text = CHUNK_LOG.read_text(encoding='utf-8', errors='replace')
    walls = [float(m) for m in re.findall(r'"wall_s": ([0-9.]+)', text)]
    peaks = [int(m) for m in re.findall(r'"peak_job_commit_bytes": (\d+)', text)]
    fails = re.findall(r'"reason": "(?!completed")([^"]+)"', text)
    res['execution_summary'] = {
        'h5_present': n_ok, 'n_completed_records': len(walls),
        'wall_s_min': min(walls) if walls else None,
        'wall_s_max': max(walls) if walls else None,
        'job_commit_GiB_max': max(peaks) / 2 ** 30 if peaks else None,
        'noncompleted_reasons': sorted(set(fails)),
        'note': 'records cover cases whose supervision finalized; cases in flight '
                'at a chunk timeout carry h5 + stdout completion lines instead '
                '(t01 explicitly: solver completed 26.5 s, supervision record '
                'lost to the first chunk timeout)'}
    (OUT / 'results.json').write_text(
        json.dumps(res, ensure_ascii=False, indent=1, sort_keys=True) + '\n',
        encoding='utf-8')
    nc_bad = [r for r in res['nc_bg_zero_diff'] if not r['zero']]
    print(f'h5 present {n_ok}/396 | NC-BG zero-diff {len(res["nc_bg_zero_diff"]) - len(nc_bad)}/'
          f'{len(res["nc_bg_zero_diff"])} | walls {res["execution_summary"]["wall_s_min"]:.1f}-'
          f'{res["execution_summary"]["wall_s_max"]:.1f} s | commit max '
          f'{res["execution_summary"]["job_commit_GiB_max"]:.2f} GiB')
    print('NC violations:', nc_bad[:5])
    print('noncompleted:', res['execution_summary']['noncompleted_reasons'])


if __name__ == '__main__':
    main()
