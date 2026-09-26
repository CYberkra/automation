"""Acceptance check for the CO evaluation runner: r1/r2 byte-identity after
stripping per-record resource fields; manifest identical except
run_tag/total_wall_s; key metric sanity vs the frozen event table."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-26_eval_batch2d_co_r1'
R2 = ROOT / 'artifacts/research_checks/2026-09-26_eval_batch2d_co_r2'


def stripped_records(run_dir):
    recs = json.loads((run_dir / 'records.json').read_text(encoding='utf-8'))['records']
    for r in recs:
        r.pop('resource', None)
    return recs


def stripped_manifest(run_dir):
    m = json.loads((run_dir / 'run_manifest.json').read_text(encoding='utf-8'))
    m.pop('run_tag', None)
    m.pop('total_wall_s', None)
    return m


r1, r2 = stripped_records(R1), stripped_records(R2)
b1 = json.dumps(r1, sort_keys=True, ensure_ascii=False).encode('utf-8')
b2 = json.dumps(r2, sort_keys=True, ensure_ascii=False).encode('utf-8')
print('records stripped byte-identical:', b1 == b2, f'({len(r1)} records)')

m1, m2 = stripped_manifest(R1), stripped_manifest(R2)
mb1 = json.dumps(m1, sort_keys=True, ensure_ascii=False).encode('utf-8')
mb2 = json.dumps(m2, sort_keys=True, ensure_ascii=False).encode('utf-8')
print('manifest stripped byte-identical:', mb1 == mb2)

# --- key-number sanity on r1 ---
m = json.loads((R1 / 'run_manifest.json').read_text(encoding='utf-8'))
print('event_rows:', m['event_rows'], 'n_records:', m['n_records'],
      'availability:', m['availability'])
print('gates:', m['gates'])
ran = [r for r in json.loads((R1 / 'records.json').read_text(encoding='utf-8'))['records']
       if r['availability'] == 'ran']
print('ran rows (event,case):', sorted({(r["event_id"], r["case"]) for r in ran}))
nb = {}
for r in ran:
    if r['row_type'] == 'negative_control':
        nb.setdefault(r['candidate_id'], []).append(r['metrics']['N_b']['energy_ratio'])
print('n gain records:', sum(1 for r in ran if r['gain_diagnostics'] is not None))
print('clipped/overflow:', sum(1 for r in ran
                               if r['gain_diagnostics'] and
                               (r['gain_diagnostics']['clipped_samples'] or r['gain_diagnostics']['overflow'])))
print('nc rows ran:', sorted({r['event_id'] for r in ran if r['row_type'] == 'negative_control'}))
import hashlib
h = hashlib.sha256(b1).hexdigest()
print('stripped records sha256:', h)
