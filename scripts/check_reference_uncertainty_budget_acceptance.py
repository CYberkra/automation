"""Acceptance for the G4 step-1 reference uncertainty budget.

Checks: r1/r2 byte identity; frozen event-table gate; per-window coverage
(43 windows, three metric classes each, amplitude class G1/G6-flagged);
preregistered judgment rule verbatim against the G4 plan; reference_state
and hard-limit discipline; and cross-checks of carried numbers against the
four archived source analyses. Read-only; no solver.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json'
R2 = ROOT / 'artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r2/results.json'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
PLAN = ROOT / 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md'
E6 = ROOT / 'artifacts/research_checks/2026-09-25_deep_zfine2_analysis_r1/results.json'
E7 = ROOT / 'artifacts/research_checks/2026-09-26_dep3d_gold_analysis_r1/results.json'
FINE2 = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1/results.json'
EB = ROOT / 'artifacts/research_checks/2026-09-25_error_budget/restoration_results.json'

RULE_ZH = '算子间指标差异 ≤ 对应窗参考不确定性 ⇒ 该差异记 undetermined，不得进入排序'

b1 = R1.read_bytes()
b2 = R2.read_bytes()
assert b1 == b2, 'r1/r2 results.json not byte-identical'
sha = hashlib.sha256(b1).hexdigest()
print('r1/r2 byte-identical: True')
print('sha256:', sha)

budget = json.loads(b1.decode('utf-8'))
et = json.loads(EVENT_TABLE.read_text(encoding='utf-8'))
e6 = json.loads(E6.read_text(encoding='utf-8'))
e7 = json.loads(E7.read_text(encoding='utf-8'))
f2 = json.loads(FINE2.read_text(encoding='utf-8'))
eb = json.loads(EB.read_text(encoding='utf-8'))

# 1. frozen event table gate
import subprocess
actual_et_sha = hashlib.sha256(EVENT_TABLE.read_bytes()).hexdigest()
assert actual_et_sha == EVENT_TABLE_SHA256 == budget['event_table_sha256'], 'event table gate mismatch'
print('event table gate: ok', EVENT_TABLE_SHA256[:16] + '...')

# 2. window coverage mirrors the frozen event table
wins = budget['per_event_window']
assert len(wins) == et['n_entries'] == budget['n_event_windows'] == 43
for w, e in zip(wins, et['entries']):
    assert w['event_id'] == e['event_id'], (w['event_id'], e['event_id'])
    assert w['family'] == e['family'] and w['role'] == e['role']
    for cls in ('arrival_time_class', 'shape_class', 'amplitude_class'):
        assert cls in w, (w['event_id'], cls)
    assert w['amplitude_class']['g1_g6_limited'] is True
    assert 'Warren 2016' in w['amplitude_class']['g1_g6_note']
print('window coverage: 43/43, three metric classes each, amplitude G1/G6 flagged: ok')

# 3. preregistered judgment rule verbatim vs the plan
assert budget['judgment_rule_preregistered']['zh'] == RULE_ZH
plan_text = PLAN.read_text(encoding='utf-8')
assert '算子间指标差异 ≤ 对应窗参考不确定性' in plan_text and '不得进入排序' in plan_text
print('judgment rule verbatim vs G4 plan §2: ok')

# 4. discipline: reference_state unchanged, threshold null, hard limits
assert budget['reference_state'] == 'numerically_unresolved'
assert budget['reference_state_unchanged_by_this_budget'] is True
assert budget['physical_acceptance_threshold'] is None
hl = budget['hard_limits']
for k in ('solver_invoked', 'training_labels_generated', 'clean_truth_generated',
          'grid_convergence_certified', 'absolute_accuracy_claimed'):
    assert hl[k] is False, k
print('discipline (numerically_unresolved / threshold null / hard limits): ok')

# 5. carried numbers cross-check against archived sources
fam = budget['per_family_direct_evidence']
p2 = f2['p2_direction_consistency']
for fam_id in ('c1', 'c3', 'c5', 'c8'):
    src = p2[fam_id.upper()]
    assert fam[fam_id]['envelope_peak_shift_ns']['fine2_minus_base_ns'] == \
        src['envelope_peak_direction']['fine2_minus_base_ns']
    assert fam[fam_id]['difference_spectrum_shape_correlation_main_200ns'] == \
        src['difference_spectrum_shape_correlation']['main_200ns']
gc = budget['scenario_level_uncertainty']['grid_chain_E6']
assert gc['chain_full_band_relative_L2']['FINE_to_ZFINE'] == \
    e6['zfine_vs_uniform_fine']['relative_L2_difference']
assert gc['chain_full_band_relative_L2']['ZFINE_to_ZFINE2'] == \
    e6['zfine2_vs_zfine']['relative_L2_difference']
cd = budget['scenario_level_uncertainty']['cross_dimension_E7']['pairs']
assert cd['B2D5CM_vs_DEP3D_5CM']['normalized_spectrum_shape_correlation'] == \
    e7['cross_dimension_same_spacing']['B2D5CM_vs_DEP3D_5CM']['normalized_spectrum_shape_correlation']
assert budget['scenario_level_uncertainty']['constructed_array_identifiability']['budget_max'] == \
    eb['statistics']['budget_max']
# per-window values are per-family fine2 values
for w in wins:
    f = fam[w['family']]
    assert w['arrival_time_class']['direct_evidence']['grid_tier_shift_ns'] == \
        f['envelope_peak_shift_ns']['fine2_minus_base_ns']
    assert w['shape_class']['direct_evidence']['per_frequency_sign_agreement_real'] == \
        f['per_frequency_sign_agreement_real']
print('carried numbers cross-check vs E6/E7/fine2/error-budget: ok')

# 6. C8 resolution attention note present
assert 'c8_resolution_attention' in budget['family_level_notes']
print('C8 resolution-attention note: ok')

print('\nALL ACCEPTANCE CHECKS PASSED')
