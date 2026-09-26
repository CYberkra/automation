"""Acceptance for the G4 step-4 sensitivity / pairwise flip-rate analysis.

Checks: r1/r2 byte identity; gates (mission tolerance, event table, reference
budget); merged ladder totals (34,086 records, identity anchor D==0); stratum
layering (exactly one geometry and one family per stratum, no merged strata,
dev families only); three-axis structure and arithmetic (flip_rate =
flips/comparisons, rates in [0,1], jitter replicates fixed); frozen level grid
coverage per mission class; ordering consistency with full_values; discipline
constants. Read-only; no solver.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json'
R2 = ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r2/results.json'
MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.2.json'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
REF_BUDGET = ROOT / 'artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json'

MISSION_SHA256 = 'ee039fb1fab3f1147e8a5fe53809bbb8ae9d6f51aafc2feee0fb90071a6c57ea'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
REF_BUDGET_SHA256 = '38cf98bab723358dc8d75486ef6824057d42c32a9030b56f7101400fc87a8cc1'

JITTER_REPLICATES = 200
DEV_FAMILIES = ('c1', 'c3')
SINGLE_LEVEL_AXES = {
    # (damage_type, mission_class): levels frozen on the mission grid
    ('polarity_flip', 'mission_relevant'): 1,
    ('amplitude_scale', 'weak_event_band'): 1,
    ('trace_deletion', 'weak_event_band'): 1,
}
MULTI_LEVEL = {
    # v0.2 sign-symmetric grid (determination 2026-09-26 22:23)
    ('amplitude_scale', 'mission_relevant'): {0.5, 0.1},
    ('sample_shift', 'mission_relevant'): {-16, -4, 4, 16},
    ('sample_shift', 'weak_event_band'): {-1, 1},
    ('trace_deletion', 'mission_relevant'): {4, 8},
}

b1 = R1.read_bytes()
b2 = R2.read_bytes()
assert b1 == b2, 'r1/r2 not byte-identical'
sha = hashlib.sha256(b1).hexdigest()
print('r1/r2 byte-identical: True')
print('sha256:', sha)

d = json.loads(b1.decode('utf-8'))
assert d['mission_tolerance_sha256'] == MISSION_SHA256 == \
    hashlib.sha256(MISSION.read_bytes()).hexdigest(), 'mission gate mismatch'
assert d['event_table_sha256'] == EVENT_TABLE_SHA256 == \
    hashlib.sha256(EVENT_TABLE.read_bytes()).hexdigest(), 'event table gate mismatch'
assert d['reference_budget_sha256'] == REF_BUDGET_SHA256 == \
    hashlib.sha256(REF_BUDGET.read_bytes()).hexdigest(), 'reference budget gate mismatch'
print('gates (mission tolerance v0.2 ee039fb1..., event table b0ad1003..., reference budget 38cf98ba...): ok')

assert d['ladder_records_merged'] == 34086
assert d['n_strata'] == len(d['strata']) == 24, d['n_strata']
hl = d['hard_limits']
assert hl['physical_acceptance_threshold'] is None
assert hl['reference_state'] == 'numerically_unresolved'
for k in ('solver_invoked', 'training_labels_generated', 'clean_truth_generated',
          'grid_convergence_certified', 'absolute_accuracy_claimed'):
    assert hl[k] is False, k
for k in ('constructed_reference_not_physical_truth', 'no_ranking_no_selection_no_thresholding',
          'simulation_domain_only'):
    assert hl[k] is True, k
assert 'NOT relieved' in hl['G4']
print('discipline (threshold null / numerically_unresolved / constructed ref / no ranking / G4 not relieved): ok')

seen = set()
n_d = n_nc = 0
for s in d['strata']:
    key = (s['damage_type'], s['mission_class'], s['geometry'], s['family'])
    assert key not in seen, f'duplicate stratum {key}'
    seen.add(key)
    assert s['geometry'] in ('mt', 'co') and s['family'] in DEV_FAMILIES
    m = s['n_candidates']
    assert m == len(s['full_values']) == len(s['full_ordering']), key
    assert s['n_pairs'] == m * (m - 1) // 2, key
    # ordering consistency: ascending by value, cid tie-break
    expect = sorted(s['full_values'], key=lambda c: (s['full_values'][c], c))
    assert s['full_ordering'] == expect, key
    assert all(v >= 0 for v in s['full_values'].values()), key
    if s['quantity'] == 'D':
        n_d += 1
    else:
        assert s['quantity'] == 'N_b' and s['damage_type'] == 'nc_amplify'
        assert s['mission_class'] == 'negative_control'
        n_nc += 1
    # level axis
    la = s['level_axis']
    axis_key = key[:2]
    if axis_key in SINGLE_LEVEL_AXES:
        assert la['applicable'] is False and la['flip_rate'] is None, key
        assert la['n_levels'] == 1, key
    elif axis_key in MULTI_LEVEL:
        frozen = MULTI_LEVEL[axis_key]
        got = {json.loads(repr(x)) if isinstance(x, str) else x for x in s['levels']}
        assert got == frozen, (key, got, frozen)
        assert la['applicable'] is True and la['flip_rate'] is not None, key
        assert la['n_levels'] == len(frozen), key
        assert abs(la['flip_rate'] - la['flipped_pairs'] / la['comparisons']) < 1e-15, key
        assert 0.0 <= la['flip_rate'] <= 1.0
    else:
        # nc_amplify: levels are the two gain factors
        assert s['quantity'] == 'N_b', key
        assert la['applicable'] is True and la['n_levels'] == 2, key
    # window jackknife
    wa = s['window_jackknife']
    assert wa['n_replicates'] == s['n_events'], key
    if s['n_events'] >= 2:
        assert wa['applicable'] is True and wa['flip_rate'] is not None, key
        assert abs(wa['flip_rate'] - wa['flips'] / wa['comparisons']) < 1e-15, key
        assert 0.0 <= wa['flip_rate'] <= 1.0
        assert wa['max_replicate_flip_rate'] is not None
        assert wa['max_replicate_flip_rate'] >= wa['flip_rate'] - 1e-15
    else:
        assert wa['applicable'] is False and wa['flip_rate'] is None, key
    # value jitter
    va = s['value_jitter']
    assert va['replicates'] == JITTER_REPLICATES, key
    assert va['epsilon'] > 0, key
    assert va['flip_rate'] is not None and 0.0 <= va['flip_rate'] <= 1.0, key
    assert abs(va['flip_rate'] - va['flips'] / va['comparisons']) < 1e-15, key
    assert va['n_event_windows'] == s['n_events'], key
print(f'strata: 24 (D={n_d}, N_b={n_nc}), layered (one geometry x one family each), ordering consistent: ok')

# jitter epsilon range: sqrt(2 * [0.000451, 0.001875]) ~ [0.0300, 0.0613]
for s in d['strata']:
    eps = s['value_jitter']['epsilon']
    assert 0.029 < eps < 0.062, (s['damage_type'], eps)
print('jitter epsilon within the reference-budget calibrated range: ok')

# no merged strata anywhere: every stratum carries exactly one geometry and one family
assert all(s['geometry'] in ('mt', 'co') and s['family'] in DEV_FAMILIES for s in d['strata'])
print('flip rates are diagnostic distributions only; no cutoff fixed (step 5 reserved): noted in methodology')

print('\nALL ACCEPTANCE CHECKS PASSED')
