"""Empirical sign-symmetry study of sample_shift damage levels (evidence for G4
mission tolerance v0.2 determination, user 2026-09-26 22:23).

Read-only aggregation over the r1 damage-ladder chunks: for every
candidate x geometry x dev-family x |shift| group, compare the median D of
positive-shift records against negative-shift records. Prints the count and
the distribution of |median(D+) - median(D-)|.

Finding (frozen into configs/research/g4_mission_tolerance_v0.2.json
v02_rationale.evidence_sign_symmetry): 484/486 groups differ by <= 0.01 in
median D; the worst is 0.012, below the value-jitter epsilon (~0.03)
calibrated from the reference-uncertainty budget. Operator response to shift
sign is a nuisance direction at the reference-uncertainty scale.
"""
import json
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHUNK_GLOB = 'artifacts/research_checks/2026-09-26_damage_ladder_r1_c{:02d}'
N_CHUNKS = 8
DEV_FAMILIES = ('c1', 'c3')
THRESHOLD = 0.01


def main():
    groups = defaultdict(lambda: defaultdict(list))
    n_records = 0
    for k in range(N_CHUNKS):
        d = json.loads((ROOT / CHUNK_GLOB.format(k) / 'records.json')
                       .read_bytes().decode('utf-8'))
        for r in d['records']:
            n_records += 1
            if (r['family'] not in DEV_FAMILIES or r['availability'] != 'ran'
                    or r['row_type'] != 'event'
                    or r['damage']['type'] != 'sample_shift'):
                continue
            lvl = r['damage']['level']
            if lvl == 0:
                continue
            sign = '+' if lvl > 0 else '-'
            key = (r['candidate_id'], r['geometry'], r['family'], abs(lvl))
            groups[key][sign].append(r['metrics']['waveform']['metrics']['nrmse'])

    rows = []
    for key in sorted(groups):
        pos, neg = groups[key].get('+', []), groups[key].get('-', [])
        if not pos or not neg:
            continue
        rows.append((key, statistics.median(pos), statistics.median(neg),
                     len(pos), len(neg)))
    diffs = [abs(a - b) for _key, a, b, _np, _nn in rows]
    within = sum(1 for x in diffs if x <= THRESHOLD)

    print(f'ladder records scanned: {n_records}')
    print(f'comparable groups (candidate x geometry x family x |shift|): {len(rows)}')
    print(f'|median(D+) - median(D-)|: min={min(diffs):.4f} '
          f'median={statistics.median(diffs):.4f} max={max(diffs):.4f}')
    print(f'within {THRESHOLD}: {within}/{len(diffs)} = {within / len(diffs) * 100:.1f}%')
    print('\nworst 10 asymmetries:')
    for key, a, b, np_, nn in sorted(rows, key=lambda r: -abs(r[1] - r[2]))[:10]:
        print(f'  {key[0]:18s} {key[1]} {key[2]} k={key[3]:2d} '
              f'med(+)={a:.4f} med(-)={b:.4f} (n={np_}/{nn})')


if __name__ == '__main__':
    main()
