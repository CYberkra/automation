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
R1 = ROOT / 'artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r1/results.json'
R2 = ROOT / 'artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r2/results.json'
MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.2.json'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'

MISSION_SHA256 = 'ee039fb1fab3f1147e8a5fe53809bbb8ae9d6f51aafc2feee0fb90071a6c57ea'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
RESULTS_SHA256 = '2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8'

# level grid derived from the frozen mission tolerance config v0.2
# (ladder_level_mapping; sign-symmetric explicit grid per determination
# 2026-09-26 22:23; v3 export gates on v0.2 and supersedes v1/v2)
MISSION_LEVELS = {
    ('amplitude_scale', 0.5), ('amplitude_scale', 0.1), ('polarity_flip', 1.0),
    ('sample_shift', -16), ('sample_shift', -4), ('sample_shift', 4), ('sample_shift', 16),
    ('trace_deletion', 4), ('trace_deletion', 8),
}
WEAK_LEVELS = {
    ('amplitude_scale', 0.9), ('sample_shift', -1), ('sample_shift', 1), ('trace_deletion', 1),
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
print('gates (mission tolerance v0.2 ee039fb1..., event table b0ad1003...): ok')

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
assert d['n_capability_rows'] == len(cap) == 2065
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
print('capability rows: 2065, mission/weak classes per frozen mapping v0.2, strata unique: ok')
# v0.2 grid is sign-symmetric: every shift magnitude must appear with both signs
shift_levels = {(r['damage_type'], r['damage_level']) for r in cap}
for mag in (4, 16):
    for sgn in (-1, 1):
        assert ('sample_shift', sgn * mag) in shift_levels, ('sample_shift', sgn * mag)
for sgn in (-1, 1):
    assert ('sample_shift', sgn * 1) in shift_levels, ('sample_shift', sgn * 1)
print('frozen-grid guard (v0.2 sign-symmetric shift grid present): ok')

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
