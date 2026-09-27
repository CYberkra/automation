"""Compare independent fsum calculations to the implementation's saved output."""
from pathlib import Path
import hashlib
import json
import math

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'START_HERE.md').exists())
OUT = ROOT / 'artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic'
root = json.loads((OUT / 'sfcw_pair_root_check.json').read_text(encoding='utf-8'))
actual = json.loads((OUT / 'results.json').read_text(encoding='utf-8'))
assert root['contract_sha256'] == actual['provenance']['contract_sha256_raw_bytes']
checks = []
def compare(a, b, label):
    if a is None or b is None:
        assert a is None and b is None, label
    else:
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-24), (label, a, b)
    checks.append(label)

index = {(r['geometry'], r['taper'], r['window_id'], r['operator']): r for r in actual['metrics']}
assert len(index) == len(root['rows']) == 48
for r in root['rows']:
    key = (r['geometry'], r['taper'], r['window'], r['operator'])
    a = index[key]
    for expected, observed in [('n_samples','sample_count'), ('first_time_s','first_sample_s'), ('last_time_s','last_sample_s'), ('background_norm','background_norm'), ('delta_norm','delta_norm')]:
        compare(r[expected], a[observed], str(key) + expected)
    for expected, observed in [('background_residual_ratio','background_residual_norm_ratio'), ('paired_change_error','paired_change_error_ratio'), ('delta_retention','delta_retained_norm_ratio')]:
        compare(r[expected], a[observed]['value'], str(key) + expected)
    compare(1/r['background_over_delta'], a['delta_over_background_norm']['value'], str(key)+'delta_over_background')

index = {(r['geometry'], r['window_id']): r for r in actual['tail_sensitivity']}
assert len(index) == len(root['taper_sensitivity']) == 12
for r in root['taper_sensitivity']:
    a = index[r['geometry'], r['window']]
    for expected, observed in [('background_taper_change','background_tail_change_over_window_no_taper_norm'), ('delta_taper_change','delta_tail_change_over_window_no_taper_norm')]:
        compare(r[expected], a[observed]['value'], str((r['geometry'],r['window']))+expected)
result = {'status':'passed_arithmetic_acceptance_only', 'checks':len(checks), 'metric_rows':48, 'raw_taper_rows':12,
          'relative_tolerance':1e-11, 'absolute_tolerance':1e-24,
          'tolerance_scope':'implementation agreement only; not physical error or reference certification',
          'results_sha256':hashlib.sha256((OUT/'results.json').read_bytes()).hexdigest(),
          'contract_sha256':root['contract_sha256'], 'checks_completed':checks,
          'solver_executed':False,'training_executed':False,'g4_relieved':False}
(OUT/'root_acceptance.json').write_text(json.dumps(result, indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='checks_completed'},indent=2))
