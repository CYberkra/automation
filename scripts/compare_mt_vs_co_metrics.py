"""Same-metric comparison between the MT (variable-offset, 33-trace) and CO
(common-offset, 11-trace) evaluations, restricted to rows that exist in both:
  * target event EV-C3m-D10m-W4m-T0.5m-E20-S0.02 (identical frozen window, applied
    unshifted in both runners);
  * c3 BG-NEG rows D5m/D10m/D20m on the family BG gather.
Computed per candidate: N_b energy_ratio (diagnostic only; not detectability/SNR/
threshold), plus SVD diagnostics (relative singular values, numerical rank,
cutoff gap) on the target-event window. Diagnostic semantics only; G4 thresholds
remain null; no ranking or labels are derived here.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MT = ROOT / 'artifacts/research_checks/2026-09-26_eval_batch2d_mt_r1/records.json'
CO = ROOT / 'artifacts/research_checks/2026-09-26_eval_batch2d_co_r1/records.json'
OUT = ROOT / 'artifacts/research_checks/2026-09-26_eval_batch2d_co_analysis/mt_vs_co_comparison.json'

TARGET_EVENT = 'EV-C3m-D10m-W4m-T0.5m-E20-S0.02'
BGNEG_EVENTS = ['EV-C3m-BG-NEG-D5m', 'EV-C3m-BG-NEG-D10m', 'EV-C3m-BG-NEG-D20m']


def load(path):
    recs = json.loads(path.read_text(encoding='utf-8'))['records']
    return {(r['event_id'], r['candidate_id']): r for r in recs if r['availability'] == 'ran'}


mt, co = load(MT), load(CO)
candidates = sorted({cid for (_, cid) in mt} & {cid for (_, cid) in co})


def nb_energy(rec):
    return rec['metrics']['N_b']['energy_ratio']


rows = []
for cid in candidates:
    mt_t = mt.get((TARGET_EVENT, cid))
    co_t = co.get((TARGET_EVENT, cid))
    row = {
        'candidate_id': cid,
        'target_event': TARGET_EVENT,
        'nb_energy_ratio': {
            'mt_variable_offset': None if mt_t is None else nb_energy(mt_t),
            'co_common_offset': None if co_t is None else nb_energy(co_t),
        },
        'svd_diagnostics_co': None if co_t is None else co_t['svd_diagnostics'],
        'svd_diagnostics_mt': None if mt_t is None else mt_t['svd_diagnostics'],
    }
    rows.append(row)

bgneg = {}
for cid in candidates:
    mt_vals = [nb_energy(mt[(ev, cid)]) for ev in BGNEG_EVENTS if (ev, cid) in mt]
    co_vals = [nb_energy(co[(ev, cid)]) for ev in BGNEG_EVENTS if (ev, cid) in co]
    bgneg[cid] = {
        'events': BGNEG_EVENTS,
        'mt_median': sorted(mt_vals)[len(mt_vals) // 2] if mt_vals else None,
        'co_median': sorted(co_vals)[len(co_vals) // 2] if co_vals else None,
        'n_mt': len(mt_vals), 'n_co': len(co_vals),
    }

out = {
    'comparison_semantics': 'same frozen event windows applied unshifted; MT gather = 33-trace variable-offset single-shot (Tx fixed y15.35, offset 2.70->1.30->5.30 m); CO gather = 11-trace common-offset profile (offset 1.30 m, Rx y 12.65..22.65 m); N_b energy_ratio is a diagnostic residual ratio, not detectability/SNR/threshold',
    'target_event': TARGET_EVENT,
    'bgneg_events': BGNEG_EVENTS,
    'per_candidate': rows,
    'bgneg_median': bgneg,
}
Path(OUT).parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')

print(f'{"candidate":<16} {"N_b MT":>12} {"N_b CO":>12} {"CO/MT":>9}')
for r in rows:
    m, c = r['nb_energy_ratio']['mt_variable_offset'], r['nb_energy_ratio']['co_common_offset']
    ratio = '' if (m in (None, 0.0) or c is None) else f'{c / m:9.3f}'
    print(f"{r['candidate_id']:<16} {m!s:>12} {c!s:>12} {ratio:>9}")
print()
print('BG-NEG median energy_ratio (c3 D5/D10/D20 on BG):')
for cid in candidates:
    b = bgneg[cid]
    print(f"  {cid:<16} MT={b['mt_median']} (n={b['n_mt']})  CO={b['co_median']} (n={b['n_co']})")
print()
co_svd = next(r['svd_diagnostics_co'] for r in rows if r['candidate_id'] == 'B0_G1_BG')
mt_svd = next(r['svd_diagnostics_mt'] for r in rows if r['candidate_id'] == 'B0_G1_BG')
print('target-window SVD relative singular values (first 6):')
print('  CO:', None if co_svd is None else [round(v, 6) for v in co_svd['relative_singular_values'][:6]])
print('  MT:', None if mt_svd is None else [round(v, 6) for v in mt_svd['relative_singular_values'][:6]])
print('wrote', OUT)
