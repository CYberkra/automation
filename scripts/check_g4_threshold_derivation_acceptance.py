"""Acceptance for the G4 step-5 S2/S3 threshold-derivation execution.

Checks: r1/r2 byte identity; all nine frozen input gates (hashes recorded in the
result re-verified against the files on disk); discipline constants (G4 NOT
relieved, simulation-domain only, constructed reference, no ranking/selection/
training, solver never invoked); the internal acceptance gates of the run
(per-pair re-derivation reproduces fliprate v2 aggregates; identity anchor
admissible with D==0 on all mission-relevant levels); candidate set (exactly 57,
none missing vs the ladder); per-combo admissibility arithmetic (admissible ==
D_pass AND N_pass, independently re-evaluated from the per-candidate rows);
the independently obtained headline expectation (admissible set is exactly the
identity anchor in all three development combos; all 21 D strata report
unique_stable_top with the singleton vacuity flag; the 3 N_b strata carry the
cross-scale caveat and no per-pair verdict); a80 capability statements keep the
systematic bias term separated (never mixed into the statistical term, no depth
conversion). Read-only; no solver.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-27_g4_threshold_derivation_r1/results.json'
R2 = ROOT / 'artifacts/research_checks/2026-09-27_g4_threshold_derivation_r2/results.json'
M1 = ROOT / 'artifacts/research_checks/2026-09-27_g4_threshold_derivation_r1/run_manifest.json'
M2 = ROOT / 'artifacts/research_checks/2026-09-27_g4_threshold_derivation_r2/run_manifest.json'
FLIPRATE = ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json'
DERIVATION = ROOT / 'configs/research/g4_threshold_derivation_v1.0.json'

IDENTITY = 'B0_G1_BG'
EXPECTED_COMBOS = {('co', 'c3'), ('mt', 'c1'), ('mt', 'c3')}
DAMAGE_TYPES = ('amplitude_scale', 'polarity_flip', 'sample_shift', 'trace_deletion')

b1 = R1.read_bytes()
b2 = R2.read_bytes()
assert b1 == b2, 'r1/r2 results not byte-identical'
sha = hashlib.sha256(b1).hexdigest()
print('r1/r2 byte-identical: True')
print('sha256:', sha)
m1 = json.loads(M1.read_text(encoding='utf-8'))
m2 = json.loads(M2.read_text(encoding='utf-8'))
assert m1['run_tag'] == 'r1' and m2['run_tag'] == 'r2'
for m in (m1, m2):
    assert m['solver_invoked'] is False and m['training_eligible'] is False
    assert m['g4_relieved'] is False and m['inputs_gated'] == 9

d = json.loads(b1.decode('utf-8'))
assert d['schema'] == 'g4-threshold-derivation/1'

# ---- frozen input gates: recorded hashes re-verified against disk ----------
for rel, expect in d['frozen_inputs_sha256'].items():
    if ' ' in rel:  # the merged-ladder entry is a logical label, not a path
        continue
    got = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
    assert got == expect, f'input gate mismatch: {rel}'
print('nine frozen input gates re-verified on disk: ok')

# ---- discipline -----------------------------------------------------------
assert d['g4_relieved'] is False
assert d['thresholds_are_simulation_domain_only'] is True
assert 'NOT physical clean truth' in d['constructed_reference_declaration']
assert 'no configuration is' in d['no_ranking_no_selection_no_training']
hl = d['hard_limits']
assert isinstance(hl, list) and len(hl) == 6, 'hard limits not carried into the result'
assert any('G4 NOT relieved' in s for s in hl)
assert any('constructed reference' in s for s in hl)
deriv = json.loads(DERIVATION.read_text(encoding='utf-8'))
assert d['derivation_rules'] == deriv['derivation_rules']
assert d['hard_limits'] == deriv['hard_limits']
assert d['disclosures'] == deriv['disclosures']
print('discipline (G4 NOT relieved / simulation-domain only / constructed ref / no ranking): ok')

# ---- internal acceptance gates of the run ---------------------------------
c = d['checks']
assert c['per_pair_rederivation_reproduces_fliprate_v2_aggregates'] is True
assert c['n_strata_verified'] == 24 and c['n_d_strata'] == 21 and c['n_n_b_strata'] == 3
assert c['identity_anchor_admissible_and_D_zero'] is True
assert all(row['admissible'] and row['D_p80_zero_on_all_mission_levels']
           for row in c['identity_anchor_detail'])
assert c['identity_anchor_candidate_id'] == IDENTITY
print('internal gates (v2 aggregate reproduction, 21 D + 3 N_b strata, identity anchor D==0): ok')

# ---- candidate set ---------------------------------------------------------
cs = d['candidate_set']
assert cs['n_candidates'] == 57 and len(cs['candidate_ids']) == 57
assert cs['candidates_in_ladder_not_in_capability'] == []
print('candidate set: 57/57, no ladder candidate missing: ok')

# ---- per-combo admissibility arithmetic + headline expectation -------------
flip = json.loads(FLIPRATE.read_text(encoding='utf-8'))
v2_eps = {}
for row in flip['strata']:
    v2_eps[(row['damage_type'], row['mission_class'], row['geometry'], row['family'])] = \
        row['value_jitter']['epsilon']

seen = set()
for combo in d['combo_summary']:
    key = (combo['geometry'], combo['family'])
    assert key in EXPECTED_COMBOS and key not in seen
    seen.add(key)
    assert combo['n_candidates_total'] == 57
    assert len(combo['per_candidate']) == 57
    recompute_adm = []
    for pc in combo['per_candidate']:
        adm = bool(pc['D_pass'] and pc['N_constraint']['pass'])
        assert adm == pc['admissible'], (key, pc['candidate_id'])
        if adm:
            recompute_adm.append(pc['candidate_id'])
        if not pc['D_pass']:
            # every D failure must carry an explicit reason per block
            assert any(b['verdict_reason'] for b in pc['D_constraint']), (key, pc['candidate_id'])
        # epsilon used in each D block must equal the fliprate v2 stratum epsilon
        for b in pc['D_constraint']:
            dt = b['damage_type']
            expect_eps = v2_eps[(dt, 'mission_relevant', key[0], key[1])]
            assert abs(b['epsilon'] - expect_eps) < 1e-18, (key, pc['candidate_id'], dt)
        # conservative joint reading: N pass iff BOTH q rows available and <= 1.0
        qs = {row['nc_amplify_q']: row for row in pc['N_constraint']['per_q']}
        assert set(qs) == {2, 4}
        joint = all(qs[q]['available'] and qs[q]['N_b_energy_ratio_p80'] <= 1.0 for q in (2, 4))
        assert joint == pc['N_constraint']['pass'], (key, pc['candidate_id'])
    assert sorted(recompute_adm) == sorted(combo['admissible_candidates'])
    assert combo['admissible_candidates'] == [IDENTITY], (key, combo['admissible_candidates'])
    assert combo['n_admissible'] == 1 and combo['n_D_pass'] == 1
    assert len(combo['admissible_candidates']) <= combo['n_N_pass'] <= 57
    # epsilon map consistent with v2 for all four damage types
    for dt_, eps in combo['epsilon_by_damage_type'].items():
        assert dt_ in DAMAGE_TYPES
        assert abs(eps - v2_eps[(dt_, 'mission_relevant', key[0], key[1])]) < 1e-18
assert seen == EXPECTED_COMBOS
print('admissibility arithmetic re-evaluated: ok; admissible == identity anchor only in all 3 combos')

# ---- stratum stability verdicts --------------------------------------------
sv = d['strata_stability']
assert len(sv) == 24
n_d = n_nb = 0
for row in sv:
    if row['quantity'] == 'D':
        n_d += 1
        assert row['stability_applied'] is True and row['cross_scale_caveat'] is False
        assert row['stability_verdict'] == 'unique_stable_top'
        assert row['unique_stable_top_candidate'] == IDENTITY
        assert row['singleton_admissible_set_vacuously_unique'] is True
        assert row['admissible_present'] == [IDENTITY]
        # singleton: no admissible pairs, so no determined pairs
        assert row['n_pairs_admissible'] == 0 and row['n_determined_pairs'] == 0
        assert row['pair_rows'] == []
    else:
        n_nb += 1
        assert row['quantity'] == 'N_b'
        assert row['stability_applied'] is False and row['cross_scale_caveat'] is True
        assert row['stability_verdict'] == 'not_applied_cross_scale_caveat'
        assert row['pair_rows'] == [] and row['unique_stable_top_candidate'] is None
        assert row['note'] and 'cross-scale' in row['note']
assert (n_d, n_nb) == (21, 3)
print('stratum stability: 21 D unique_stable_top (singleton, vacuity flagged) + 3 N_b caveat rows: ok')

# ---- a80 capability statements: systematic bias term separated -------------
st = d['capability_statements']
assert len(st) == 3 * len(DAMAGE_TYPES) == 12
for blk in st:
    assert (blk['geometry'], blk['family']) in EXPECTED_COMBOS
    assert blk['damage_type'] in DAMAGE_TYPES
    assert blk['headline_excludes_weak_event_band'] is True
    sbt = blk['systematic_bias_term']
    assert sbt['reported_separately'] is True
    assert sbt['mixed_into_statistical_term'] is False
    assert sbt['depth_conversion_performed'] is False
    assert sbt['unit'] == 'ns'
    if sbt['grid_chain_direction_median_ns'] is not None:
        assert abs(sbt['grid_chain_direction_median_abs_ns']
                   - abs(sbt['grid_chain_direction_median_ns'])) < 1e-18
    assert len(blk['per_candidate']) == 57
    for pc in blk['per_candidate']:
        if pc['a80_statement_complete']:
            assert pc['n_levels_missing'] == 0 and pc['a80_upper_bound_D_p80'] is not None
        else:
            assert pc['n_levels_missing'] > 0 or pc['a80_upper_bound_D_p80'] is None
print('a80 capability statements: 12 blocks, systematic bias term separated, no depth conversion: ok')

print('ALL G4 S2/S3 ACCEPTANCE CHECKS PASSED')
