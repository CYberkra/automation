"""Freeze the G4 mission tolerance inputs v0.2 (deterministic, byte-reproducible).

v0.2 is DERIVED PROGRAMMATICALLY from configs/research/g4_mission_tolerance_v0.1.json
(the v0.1 file is loaded, deep-copied and patched key by key; no value is retyped
here, so no drift is possible). Re-running must produce a byte-identical file.

What v0.2 changes with respect to v0.1 (specification fix only):
  * the damage-ladder level grid becomes explicit and sign-symmetric:
    mission_relevant.sample_shift -> [-16, -4, 4, 16], weak_event_band.sample_shift
    -> [-1, 1], plus an explicit sentence in sign_convention;
  * provenance of the determination is recorded in the new key v02_rationale;
  * inputs gain the v0.1 file itself and the G4 step-4 flip-rate results.

Entry: user 2026-09-26 22:23 "这两个你好好调研研究确立一下呢？" — delegated investigation
and settlement of the two open G4 step-4 questions (level-grid sign symmetry; flip-rate
stability rule). This freeze produces NO threshold values and does NOT relieve G4.
"""
import copy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT / 'configs/research/g4_mission_tolerance_v0.1.json'
OUT = ROOT / 'configs/research/g4_mission_tolerance_v0.2.json'

# Frozen v0.1 input set, extended with the two determination evidence sources.
ADDITIONAL_INPUTS = {
    'previous_config': 'configs/research/g4_mission_tolerance_v0.1.json',
    'step4_fliprate': 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_r1/results.json',
}

SCHEMA_VERSION = 'g4_mission_tolerance_v0.2'
FROZEN_BY = 'kimi_autonomous_2026-09-26_22:23_user_delegated_determination'
USER_APPROVAL = (
    "user 2026-09-26 22:23 '这两个你好好调研研究确立一下呢？' — two open step-4 questions "
    "delegated for determination (shift-grid sign symmetry; flip-rate stability rule)"
)

MISSION_SHIFT_LEVELS = [-16, -4, 4, 16]
WEAK_SHIFT_LEVELS = [-1, 1]
SHIFT_SIGN_CONVENTION_SUFFIX = (
    "; sample_shift levels are explicit signed sample counts (negative = delay, "
    "positive = advance); the frozen grid is sign-symmetric by determination 2026-09-26 22:23"
)

V02_RATIONALE = {
    "change_scope": "specification fix only: the level grid is made explicit and sign-symmetric; the four mission inputs of v0.1 are unchanged",
    "supersedes": "configs/research/g4_mission_tolerance_v0.1.json (kept in git history; its positive-only grid was magnitude shorthand, not a sign exclusion — see commit d39065d message and decision_log stating 移位 ±4/±16)",
    "determination": "user 2026-09-26 22:23 '这两个你好好调研研究确立一下呢？' — delegated investigation and settlement of the two open G4 step-4 questions",
    "evidence_sign_symmetry": "empirical on damage-ladder r1 records: 484/486 candidate×geometry×family×|shift| comparisons have |median(D+) − median(D−)| ≤ 0.01; worst 0.012 < value-jitter epsilon ≈ 0.03; both signs already present in the ladder (no rerun needed)",
    "evidence_physics": "GPR landslide monitoring literature measures interface vertical displacement in both directions (bidirectional kinematics); at mission-tolerance magnitude the shift sign is a nuisance direction, not a distinct damage class",
    "evidence_intent": "v0.1 freeze record states mission band 移位 ±4/±16 and weak band ±1 (commit message + decision_log); the v0.1 JSON listed magnitudes only",
}

EXPECTED_EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'

# Top-level keys that must be deep-equal to v0.1 (no drift allowed).
UNCHANGED_KEYS = (
    'mission_inputs',
    'pod_mapping',
    'derivation_rules',
    'relation_to_gates',
    'status',
    'frozen_date',
    'plan',
)

# The complete set of paths allowed to differ between v0.1 and v0.2.
EXPECTED_DIFF_PATHS = {
    'schema_version',
    'frozen_by',
    'user_approval',
    'ladder_level_mapping.mission_relevant.sample_shift',
    'ladder_level_mapping.weak_event_band.sample_shift',
    'ladder_level_mapping.sign_convention',
    'v02_rationale (added)',
    'inputs.previous_config (added)',
    'inputs.step4_fliprate (added)',
    'inputs_sha256.previous_config (added)',
    'inputs_sha256.step4_fliprate (added)',
}


def sha256_file(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def load_v01():
    return json.loads(PREVIOUS.read_bytes().decode('utf-8'))


def diff_paths(old, new, prefix=''):
    """Paths where two nested structures differ; additions are tagged '(added)'."""
    found = []
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new)):
            child = f'{prefix}.{key}' if prefix else key
            if key not in old:
                found.append(f'{child} (added)')
            elif key not in new:
                found.append(f'{child} (removed)')
            else:
                found.extend(diff_paths(old[key], new[key], child))
    elif old != new:
        found.append(prefix or '<root>')
    return found


def derive(contract_v01):
    """Deep-copy v0.1 and apply the v0.2 determination patches."""
    contract = copy.deepcopy(contract_v01)
    contract['schema_version'] = SCHEMA_VERSION
    contract['frozen_by'] = FROZEN_BY
    contract['user_approval'] = USER_APPROVAL
    contract['v02_rationale'] = copy.deepcopy(V02_RATIONALE)
    ladder = contract['ladder_level_mapping']
    ladder['mission_relevant']['sample_shift'] = list(MISSION_SHIFT_LEVELS)
    ladder['weak_event_band']['sample_shift'] = list(WEAK_SHIFT_LEVELS)
    ladder['sign_convention'] = ladder['sign_convention'] + SHIFT_SIGN_CONVENTION_SUFFIX
    return contract


def assert_derivation(contract_v01, contract):
    # 1. per-key deep equality for every payload key v0.2 must not touch
    for key in UNCHANGED_KEYS:
        assert key in contract_v01, f'v0.1 is missing {key}'
        assert contract_v01[key] == contract[key], f'v0.2 must not alter {key}'
    # 2. ladder: only the two sample_shift lists and the sign_convention string may move
    old_ladder = contract_v01['ladder_level_mapping']
    new_ladder = contract['ladder_level_mapping']
    assert set(old_ladder) == set(new_ladder), 'ladder_level_mapping key set changed'
    for key in old_ladder:
        if key == 'sign_convention':
            assert new_ladder[key].startswith(old_ladder[key]), 'sign_convention must extend v0.1'
            assert new_ladder[key] == old_ladder[key] + SHIFT_SIGN_CONVENTION_SUFFIX
        elif key in ('mission_relevant', 'weak_event_band'):
            assert set(old_ladder[key]) == set(new_ladder[key]), f'{key} key set changed'
            for sub in old_ladder[key]:
                if sub == 'sample_shift':
                    assert new_ladder[key][sub] != old_ladder[key][sub], f'{key}.sample_shift not updated'
                else:
                    assert new_ladder[key][sub] == old_ladder[key][sub], f'{key}.{sub} must not change'
        else:
            assert new_ladder[key] == old_ladder[key], f'ladder_level_mapping.{key} must not change'
    # 3. new levels keep the v0.1 magnitudes, each sign made explicit (|grid| doubles)
    assert sorted(set(abs(v) for v in new_ladder['mission_relevant']['sample_shift'])) == sorted(old_ladder['mission_relevant']['sample_shift'])
    assert sorted(set(abs(v) for v in new_ladder['weak_event_band']['sample_shift'])) == sorted(old_ladder['weak_event_band']['sample_shift'])
    assert len(new_ladder['mission_relevant']['sample_shift']) == 2 * len(old_ladder['mission_relevant']['sample_shift'])
    assert len(new_ladder['weak_event_band']['sample_shift']) == 2 * len(old_ladder['weak_event_band']['sample_shift'])
    assert new_ladder['mission_relevant']['sample_shift'] == MISSION_SHIFT_LEVELS
    assert new_ladder['weak_event_band']['sample_shift'] == WEAK_SHIFT_LEVELS
    # 4. nothing beyond the declared diff set moved
    found = set(diff_paths(contract_v01, contract))
    assert found == EXPECTED_DIFF_PATHS, (
        f'unexpected v0.1 -> v0.2 differences: {sorted(found - EXPECTED_DIFF_PATHS)}; '
        f'missing expected ones: {sorted(EXPECTED_DIFF_PATHS - found)}'
    )


def main():
    contract_v01 = load_v01()
    inputs = dict(contract_v01['inputs'])
    for key, rel in ADDITIONAL_INPUTS.items():
        assert key not in inputs, f'v0.1 already declares input {key}'
        inputs[key] = rel

    actual = sha256_file(inputs['event_table'])
    assert actual == EXPECTED_EVENT_TABLE_SHA256, f'event table gate mismatch: {actual}'

    contract = derive(contract_v01)
    contract['inputs_sha256'] = {k: sha256_file(v) for k, v in inputs.items()}
    contract['inputs'] = inputs
    assert_derivation(contract_v01, contract)

    text = json.dumps(contract, indent=1, sort_keys=True, ensure_ascii=False)
    OUT.write_bytes(text.encode('utf-8'))
    print('wrote', OUT)
    print('sha256:', hashlib.sha256(OUT.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
