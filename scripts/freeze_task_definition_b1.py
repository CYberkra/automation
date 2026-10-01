"""Freeze the B1 task-definition contract v0.1 (deterministic, byte-reproducible).

Freezes ONLY the task definition of the contractual adaptive-ISP proposal
(docs/research/2026-10-01_contractual_adaptive_isp_proposal.md). This freeze
authorizes NO simulation and NO training; G4 stays unrelieved and the catalogue
keeps training_enabled=False / fdtd_execution_enabled=False.

Entry: user 2026-10-01 13:23 "确认，先试试效果" — sign-off of the B1 proposal
after review of the one-page review sheet (E:/automation_djh/B1_review_sheet.md).

Convention notes:
  * every payload value is typed here exactly once; per-key assertions verify
    structure, referenced-input hashes and the freeze semantics;
  * re-running must produce a byte-identical file (json indent=1, sort_keys);
  * EXPECTED_PROPOSAL_SHA256 gates on the exact proposal text frozen against.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/task_definition_contract_v0.1.json'

INPUTS = {
    'proposal': 'docs/research/2026-10-01_contractual_adaptive_isp_proposal.md',
    'operator_catalogue': 'configs/research/operator_catalogue_v0.2.json',
    'reward_tolerance_contract': 'configs/research/reward_tolerance_contract_v0.2.json',
    'reward_weights_contract': 'configs/research/reward_weights_contract_v0.3.json',
    'solver_execution_gate': 'configs/research/gprmax_v4_execution_gate.json',
}

# Gates on the exact proposal text frozen against (proposal drift detector).
EXPECTED_PROPOSAL_SHA256 = 'bc3db0dbadbd03655c954f49f681bdb136c46d71ffce29fbc1442851722eed04'

USER_APPROVAL = (
    'user 2026-10-01 13:23 "确认，先试试效果" — B1 proposal signed off after '
    'review of the one-page review sheet'
)

CONTRACT = {
    'schema': 'task_definition_contract/v0.1',
    'task_id': 'b1_contractual_adaptive_isp',
    'date': '2026-10-01',
    'frozen_by': 'kimi, freeze script scripts/freeze_task_definition_b1.py',
    'user_approval': USER_APPROVAL,
    'status': 'frozen_task_definition_only__no_simulation_no_training',
    'training_enabled': False,
    'fdtd_execution_enabled': False,
    'g4_relieved': False,
    'one_line_definition': (
        'per input B-scan, a lightweight policy network selects operator and '
        'discrete parameter tier stage by stage from operator_catalogue v0.2 '
        '(30 entries + STOP) to maximize the frozen contract reward'
    ),
    'mdp': {
        'state': 'current B-scan (dB scale, fixed size) + EC channels (stage counter t, operator one-hot history)',
        'action': '(operator_id, discrete_parameter_tier) or STOP; identity anchor B0_G1_BG always present',
        'transition': 'deterministic array operators (mean lambda / SVD k / local mean / RPCA lambda axis / gain class)',
        'termination': {'stop_action': True, 't_max': 4},
        'reward_timing': 'sparse terminal reward; per-stage reward zero',
        'gamma': 1.0,
    },
    'action_space': {
        'source': 'configs/research/operator_catalogue_v0.2.json',
        'enumerated_from_frozen_catalogue_only': True,
        'catalogue_entry_count': 30,
        'plus_stop': True,
        'identity_anchor': 'B0_G1_BG',
        'masking': 'illegality rules from operator contracts (truncation degeneracy refusal, mask limits, gain overflow) applied as -inf before softmax',
        'continuous_parameter_head': False,
        'discrete_tiers_only': True,
    },
    'reward': {
        'formula': 'R = R_contract - lambda_c * cost(a_1:T)',
        'r_contract_computed_by': [
            'configs/research/reward_tolerance_contract_v0.2.json',
            'configs/research/reward_weights_contract_v0.3.json',
        ],
        'weights_contract_note': (
            'v0.3 supersedes v0.2 with the background-class section verbatim '
            'unchanged (its own supersedes block asserts this); binding to v0.3 '
            'is equivalent to the proposal text "权重契约 v0.2" for the '
            'background class and additionally covers the gain class'
        ),
        'dB_identity_coefficients_fixed_at_one_not_fittable': True,
        'hard_constraint_priority': 'violations clamp R below the identity anchor ("worse than doing nothing")',
        'anchors_online_per_episode': ['identity reads -due', 'oracle reads +due'],
        'lambda_c': 'inference-time knob; switching accuracy/latency mode needs no retrain',
        'cost_initial_values': 'from existing batch telemetry',
    },
    'data_discipline': {
        'dev_families': ['C1', 'C3'],
        'test_families': ['C5', 'C8'],
        'test_zero_contact_until_one_shot_evaluation': True,
    },
    'training_method': {
        'algorithm': 'PPO (reward need not be differentiable)',
        'episode': 'sample one simulated gather from dev families -> rollout <= t_max stages -> contract reward settlement',
        'reproducibility': 'fixed seeds, deterministic operators, r1/r2 double pass; training config frozen as its own contract JSON before any start',
        'monitoring': ['train/val reward curves', 'per-operator selection frequency entropy (collapse alarm)', 'repeat runs report std'],
    },
    'goodhart_guardrails': [
        'reward coefficients stay dB-identity-constructed (=1); never back-fitted from policy training performance',
        'tau_A / tau_D / epsilon update only via mechanism-ladder recalibration (pending: S1/S3), never from training results',
        'identity/oracle anchors evaluated every episode as online calibration of the reward scale',
        'selection-frequency entropy below bound -> diagnose, not judge',
        'test families evaluated once (S5-style); no threshold re-tuning on test',
    ],
    'extension_spec_v0_1': {
        'five_step_new_family_entry': [
            'catalogue version bump with composition-schema review',
            'mechanism-ladder calibration producing discrete tiers',
            'scope_limits registration of uncovered regions',
            'action-space column extension synced to the frozen catalogue version; old weights keep working (new columns small-constant init)',
            'evaluation-rule coverage for the new family registered in the next tolerance/weights contract version',
        ],
        'explicit_non_reservations': [
            'no continuous parameter head interface',
            'no multi-objective reward structure (lambda_c*cost covers the accuracy/latency trade-off)',
        ],
    },
    'acceptance_criteria_draft': [
        'one-shot synthetic test evaluation: significantly above identity baseline',
        'not below the protocol first-run global winners (S2X/S2TZX->B9, C3mX->B7) on the same families',
        'no selection-frequency collapse alarm',
        'any failure -> undetermined, not negative; no test-family iteration',
    ],
    'prerequisites_before_any_training': [
        'A1 G4 relief (user sign-off)',
        'A2 S1/S3 tolerance recalibration',
        'B2/B3 sampling design freeze + batch simulation budget',
        'B4 label/reward generation chain acceptance (r1/r2 byte-identical)',
        'C1 PyTorch environment locked (wheel hashes archived, import verified)',
    ],
}


def sha256_file(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def assert_contract(contract, inputs_sha256):
    # 1. referenced inputs exist and hashes were computed
    assert set(inputs_sha256) == set(INPUTS)
    # 2. proposal drift gate
    if EXPECTED_PROPOSAL_SHA256 is not None:
        assert inputs_sha256['proposal'] == EXPECTED_PROPOSAL_SHA256, (
            f"proposal changed since freeze: {inputs_sha256['proposal']}"
        )
    # 3. catalogue structural facts the contract binds to
    cat = json.loads((ROOT / INPUTS['operator_catalogue']).read_bytes().decode('utf-8'))
    assert cat['candidate_count'] == 30 and len(cat['configurations']) == 30
    assert contract['action_space']['catalogue_entry_count'] == 30
    ids = [c['id'] for c in cat['configurations']]
    assert contract['action_space']['identity_anchor'] in ids, 'identity anchor missing from catalogue'
    assert cat['training_enabled'] is False and cat['fdtd_execution_enabled'] is False
    # 4. freeze semantics: this contract itself must not enable anything
    assert contract['training_enabled'] is False
    assert contract['fdtd_execution_enabled'] is False
    assert contract['g4_relieved'] is False
    # 5. weights contract supersedes relation (v0.3 embeds v0.2 background class unchanged)
    w = json.loads((ROOT / INPUTS['reward_weights_contract']).read_bytes().decode('utf-8'))
    assert w['status'] == 'frozen'
    assert w['supersedes']['path'] == 'reward_weights_contract_v0.2.json'
    assert 'background-class section unchanged' in w['supersedes']['relation']
    # 6. tolerance contract is frozen
    t = json.loads((ROOT / INPUTS['reward_tolerance_contract']).read_bytes().decode('utf-8'))
    assert t['status'] == 'frozen'
    # 7. reward refs match the loaded files
    assert contract['reward']['r_contract_computed_by'] == [
        INPUTS['reward_tolerance_contract'], INPUTS['reward_weights_contract'],
    ]
    # 8. guardrails / extension spec / acceptance present and complete
    assert len(contract['goodhart_guardrails']) == 5
    assert len(contract['extension_spec_v0_1']['five_step_new_family_entry']) == 5
    assert len(contract['extension_spec_v0_1']['explicit_non_reservations']) == 2
    assert contract['acceptance_criteria_draft'][:3] and len(contract['prerequisites_before_any_training']) == 5
    # 9. action-space discipline flags
    assert contract['action_space']['continuous_parameter_head'] is False
    assert contract['action_space']['discrete_tiers_only'] is True
    assert contract['action_space']['enumerated_from_frozen_catalogue_only'] is True


def main():
    inputs_sha256 = {k: sha256_file(v) for k, v in INPUTS.items()}
    contract = dict(CONTRACT)
    contract['inputs'] = dict(INPUTS)
    contract['inputs_sha256'] = inputs_sha256
    assert_contract(contract, inputs_sha256)

    text = json.dumps(contract, indent=1, sort_keys=True, ensure_ascii=False)
    OUT.write_bytes(text.encode('utf-8'))
    print('wrote', OUT)
    print('sha256:', hashlib.sha256(OUT.read_bytes()).hexdigest())
    print('proposal_sha256:', inputs_sha256['proposal'])


if __name__ == '__main__':
    main()
