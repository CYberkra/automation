"""Acceptance for the G4 step-3 mission capability export.

Checks: r1/r2 byte identity; mission config + event table gates; merged
ladder totals (34,086 records, identity anchor D==0); mission/weak band
coverage against the frozen mapping; MT/CO and family layering (no merged
strata); discipline constants. Read-only; no solver.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-26_g4_mission_capability_v2_r1/results.json'
R2 = ROOT / 'artifacts/research_checks/2026-09-26_g4_mission_capability_v2_r2/results.json'
MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.1.json'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'

MISSION_SHA256 = '2d630c01cc95f3c13febbc19c4bf8f3b5a53b5fb82bbabcb2f68e7537ef78b28'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
RESULTS_SHA256 = 'c92929129718703fce7fa6f5b0922d48d22a6e075e9f130e37e51226103b21a0'

# level grid derived from the frozen mission tolerance config
# (ladder_level_mapping; v2 supersedes the v1 export whose hardcoded mapping
# also admitted sample_shift -4/-16/-1 outside the frozen grid)
MISSION_LEVELS = {
    ('amplitude_scale', 0.5), ('amplitude_scale', 0.1), ('polarity_flip', 1.0),
    ('sample_shift', 4), ('sample_shift', 16),
    ('trace_deletion', 4), ('trace_deletion', 8),
}
WEAK_LEVELS = {
    ('amplitude_scale', 0.9), ('sample_shift', 1), ('trace_deletion', 1),
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
print('gates (mission tolerance 2d630c01..., event table b0ad1003...): ok')

assert d['ladder_records_merged'] == 34086
assert d['identity_anchor']
hl = d['hard_limits']
assert hl['physical_acceptance_threshold'] is None
assert hl['reference_state'] == 'numerically_unresolved'
for k in ('solver_invoked', 'training_labels_generated', 'clean_truth_generated',
          'grid_convergence_certified', 'absolute_accuracy_claimed'):
    assert hl[k] is False, k
for k in ('constructed_reference_not_physical_truth', 'no_ranking_no_selection_no_thresholding',
          'simulation_domain_only'):
    assert hl[k] is True, k
print('discipline (threshold null / numerically_unresolved / constructed ref / no ranking): ok')

cap = d['capability']
assert d['n_capability_rows'] == len(cap) == 1579
assert sha == RESULTS_SHA256, 'unexpected capability v2 export hash'
seen = set()
for r in cap:
    key = (r['damage_type'], r['damage_level'])
    assert key in MISSION_LEVELS | WEAK_LEVELS, key
    expected = ('mission_relevant' if key in MISSION_LEVELS else 'weak_event_band')
    assert r['mission_class'] == expected, (key, r['mission_class'])
    assert r['geometry'] in ('mt', 'co') and r['family'] in ('c1', 'c3')
    assert set(r) >= {'D_median', 'D_p80', 'D_min', 'D_max', 'n'}
    strata = (r['candidate_id'], r['damage_type'], repr(r['damage_level']),
              r['mission_class'], r['geometry'], r['family'])
    assert strata not in seen, f'duplicate stratum {strata}'
    seen.add(strata)
print('capability rows: 1579, mission/weak classes per frozen mapping, strata unique: ok')
assert not [r for r in cap if r['damage_type'] == 'sample_shift' and r['damage_level'] < 0], \
    'negative sample_shift levels outside the frozen grid leaked back in'
print('frozen-grid guard (no sample_shift negative levels): ok')

nc = d['negative_control']
assert d['n_negative_control_rows'] == len(nc) == 324
for r in nc:
    assert r['nc_amplify_q'] in (2, 4)
    assert r['geometry'] in ('mt', 'co') and r['family'] in ('c1', 'c3')
print('negative-control rows: 324 (pure nc_amplify q in {2,4}, layered): ok')

# no merged strata anywhere: every row carries exactly one geometry and one family
assert all(r['geometry'] in ('mt', 'co') for r in cap + nc)
print('layering check (no MT/CO merge, no family merge): ok')

print('\nALL ACCEPTANCE CHECKS PASSED')
