"""Acceptance for the damage ladder: merge r1/r2 chunks in chunk order, strip
per-record resource fields, verify byte-identity across the two runs; sanity
checks (identity anchor D==0, availability totals, manifest consistency)."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
N_CHUNKS = 8


def merged(run, pattern='artifacts/research_checks/2026-09-26_damage_ladder_{run}_c{k:02d}'):
    recs = []
    manifests = []
    for k in range(N_CHUNKS):
        d = ROOT / pattern.format(run=run, k=k)
        rs = json.loads((d / 'records.json').read_text(encoding='utf-8'))['records']
        for r in rs:
            r.pop('resource', None)
        recs.extend(rs)
        m = json.loads((d / 'run_manifest.json').read_text(encoding='utf-8'))
        for key in ('run_tag', 'total_wall_s'):
            m.pop(key, None)
        manifests.append(m)
    return recs, manifests


def comparison_bytes(records, *, cross_version=False):
    # Only the runner identity differs legitimately across runner revisions.
    # Preserve every numerical result, input hash and other provenance field.
    rows = copy.deepcopy(records)
    for row in rows:
        row.pop('resource', None)
        if cross_version:
            row['provenance'] = {k: v for k, v in row['provenance'].items()
                                 if k != 'runner_script_sha256'}
    return json.dumps(rows, sort_keys=True, ensure_ascii=False).encode('utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left-pattern', default='artifacts/research_checks/2026-09-26_damage_ladder_{run}_c{k:02d}')
    p.add_argument('--right-pattern', default='artifacts/research_checks/2026-09-26_damage_ladder_{run}_c{k:02d}')
    p.add_argument('--cross-version', action='store_true',
                   help='Numerical regression only: exclude runner_script_sha256, preserve archived identity')
    p.add_argument('--split', choices=['dev', 'test'], default='dev')
    a = p.parse_args()
    r1, m1 = merged('r1', a.left_pattern)
    r2, m2 = merged('r2', a.right_pattern)
    for rows in (r1, r2):
        hashes = {r['provenance']['runner_script_sha256'] for r in rows}
        assert len(hashes) == 1 and all(len(h) == 64 and all(c in '0123456789abcdef' for c in h) for h in hashes)
        print('recorded runner SHA256:', sorted(hashes))
    b1 = comparison_bytes(r1, cross_version=a.cross_version)
    b2 = comparison_bytes(r2, cross_version=a.cross_version)
    assert b1 == b2, 'record regression differs (outside permitted provenance/resource exclusions)'
    assert m1 == m2, 'manifest differs outside run_tag/total_wall_s'
    from run_g4_mission_capability import validate_records, DEV_FAMILIES, TEST_FAMILIES
    families = DEV_FAMILIES if a.split == 'dev' else TEST_FAMILIES
    validate_records(r1, families)
    validate_records(r2, families)
    print('comparison mode:', 'cross-version numerical regression' if a.cross_version else 'same-version reproducibility')
    print('records:', len(r1), 'sha256:', hashlib.sha256(b1).hexdigest())
    print('records/manifests/coverage/identity checks passed')


if __name__ == '__main__':
    main()
