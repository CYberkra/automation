"""Acceptance for the damage ladder: merge r1/r2 chunks in chunk order, strip
per-record resource fields, verify byte-identity across the two runs; sanity
checks (identity anchor D==0, availability totals, manifest consistency)."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
N_CHUNKS = 8


def merged(run):
    recs = []
    manifests = []
    for k in range(N_CHUNKS):
        d = ROOT / f'artifacts/research_checks/2026-09-26_damage_ladder_{run}_c{k:02d}'
        rs = json.loads((d / 'records.json').read_text(encoding='utf-8'))['records']
        for r in rs:
            r.pop('resource', None)
        recs.extend(rs)
        m = json.loads((d / 'run_manifest.json').read_text(encoding='utf-8'))
        for key in ('run_tag', 'total_wall_s'):
            m.pop(key, None)
        manifests.append(m)
    return recs, manifests


r1, m1 = merged('r1')
r2, m2 = merged('r2')
b1 = json.dumps(r1, sort_keys=True, ensure_ascii=False).encode('utf-8')
b2 = json.dumps(r2, sort_keys=True, ensure_ascii=False).encode('utf-8')
print('merged records stripped byte-identical:', b1 == b2, f'({len(r1)} records)')
print('stripped sha256:', hashlib.sha256(b1).hexdigest())

for i, (a, b) in enumerate(zip(m1, m2)):
    assert json.dumps(a, sort_keys=True, ensure_ascii=False) == \
           json.dumps(b, sort_keys=True, ensure_ascii=False), f'manifest chunk {i} differs'
print('per-chunk manifests identical except run_tag/total_wall_s: True')

ran = [r for r in r1 if r['availability'] == 'ran']
idents = [r for r in ran if r.get('candidate_id') == 'B0_G1_BG']
ok = all(r['metrics']['waveform']['metrics'] is not None
         and r['metrics']['waveform']['metrics']['nrmse'] == 0.0 for r in idents)
print('identity anchor D==0 across all ran records:', ok, f'({len(idents)} records)')
tot = {'ran': 0, 'unavailable': 0}
for r in r1:
    tot[r['availability']] += 1
print('merged availability:', tot)
refs = {json.dumps(r['damage']['reference']) for r in r1}
print('damage reference declarations:', refs)
contrast = sum(1 for r in ran
               if r['metrics']['waveform']['available']
               and r['metrics']['waveform']['reference_provenance']
               == 'identity_of_damaged_input_constructed')
print('ran records with constructed-reference D/A/H:', contrast, f'/{len(ran)}')
