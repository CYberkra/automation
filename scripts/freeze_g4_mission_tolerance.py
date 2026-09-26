"""Freeze the G4 mission tolerance inputs v0.1 (deterministic, byte-reproducible).

Generates configs/research/g4_mission_tolerance_v0.1.json from this script's
frozen content plus SHA-256 of registered inputs (G4 plan, frozen event table,
reference uncertainty budget, damage ladder acceptance evidence, Wang 2026
reading note). Re-running must produce a byte-identical file.

Entry: user 2026-09-26 21:38 "同意" — approves the four G4 step-3 mission
inputs at the proposed defaults (anchored to Wang et al. 2026, doi:
10.1109/MGRS.2026.3702334). This freeze locks the INPUTS and the derivation
mapping; it produces no threshold values and does not relieve G4.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/g4_mission_tolerance_v0.1.json'

INPUTS = {
    'g4_plan': 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md',
    'event_table': 'configs/research/batch2d_v1_event_table_v0.1.json',
    'reference_budget': 'artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json',
    'damage_ladder_acceptance_evidence': 'scripts/check_damage_ladder_acceptance.py',
    'wang2026_reading_note': 'docs/research/2026-09-26_wang2026_uav_gpr_landslide_reading.md',
    's4_baseline_group': 'configs/research/s4_baseline_group_v0.1.json',
}

CONTRACT = {
    "schema_version": "g4_mission_tolerance_v0.1",
    "frozen_date": "2026-09-26",
    "frozen_by": "kimi_autonomous_2026-09-26_21:38_user_agreed",
    "status": "frozen",
    "plan": "docs/research/2026-09-26_g4_calibration_plan_v0.1.md",
    "user_approval": "user 2026-09-26 21:38 '同意' — four mission inputs fixed at the proposed defaults",
    "mission_inputs": {
        "minimum_recognizable_damage": {
            "statement": "界面型电性对比异常（10 m 深、4 m 宽、epsilon_r=20 对围岩，锚点模型 B2D-C3m-D10m-W4m-T0.5m-E20-S0.02 及其深度变体 D5m/D20m）在冻结事件窗内的形状/极性保真",
            "anchor": "Wang et al. 2026 field-validated mission class: bedrock/sliding interface at 10-20 m depth, borehole depth error <1 m absolute / <10% relative (best 0.05 m / 0.27%)",
        },
        "false_call_cost": {
            "statement": "误报 = 鬼界面/水平振铃带被当成真实界面（Wang 2026 自列局限：离轴反射 ghost interfaces 误导解释与钻孔布设）；负控 N_b 对准鬼界/振铃带能量比",
            "negative_control_policy": "nc_zero / bg_absent / off_path 行独立评价，永不并入目标行",
        },
        "conservativeness": {
            "statement": "能力值 + 显式系统偏差项的保守表述（a80 式）：能力值取 80% 保守界限（D 的 p80 作为抹除上界），系统偏差项单列（如速度偏差 +/-10% 的到时-深度换算项），统计项与系统项不混合",
            "template_anchor": "Wang 2026 probabilistic DP: 80% band = statistical interval + systematic velocity-bias term (their eq. 26-30)",
        },
        "family_scope": {
            "statement": "容差分覆盖层厚度族给出（C1/C3/C5/C8 各自一族），不四族统一",
            "anchor": "Wang 2026 validates two geologically distinct sites separately, no pooling across sites",
            "dev_test_split": "能力值导出仅限开发族 {C1,C3}；测试族 {C5,C8} 维持 S4 冻结，阈值锁定后一次性确认评价",
        },
    },
    "ladder_level_mapping": {
        "note": "operationalization of minimum_recognizable_damage onto the frozen damage-ladder level grid (damage ladder v0.1); adjustable only by re-freeze (versioned), never by candidate performance",
        "mission_relevant": {
            "amplitude_scale": [0.5, 0.1],
            "polarity_flip": ["all"],
            "sample_shift": [4, 16],
            "trace_deletion": [4, 8],
        },
        "weak_event_band": {
            "amplitude_scale": [0.9],
            "sample_shift": [1],
            "trace_deletion": [1],
            "role": "protected weak events: feed tau calibration and undetermined-band checks, never enter the capability headline",
        },
        "rationale": "mission-relevant = perturbations at or above the Wang-2026 mission class (major interface response loss >=50%, polarity reversal, decimeter-class displacement, loss of lateral continuity across >=4 traces); weak band = the finest ladder steps, kept for sensitivity (G4 step 4) and R-rule tau calibration",
        "sign_convention": "amplitude_scale factor = residual amplitude after damage (0.5 = 50% energy loss); sample_shift in samples; trace_deletion in traces",
    },
    "pod_mapping": {
        "defect_size_a": "ladder level (amplitude factor / shift samples / deleted traces)",
        "signal_response_a_hat": "evaluation metric per class (arrival-time / shape / amplitude, amplitude class G1/G6-limited)",
        "detection_threshold_a_hat_th": "mission tolerance (damage upper limit) + minimum effective gain; derived ONLY from ladder data x these mission inputs",
        "a90_95_capability": "capability value + explicit systematic bias term (a80-style statement per conservativeness input)",
        "false_call_rate": "negative-control N_b constraint, reported alongside detection (never merged)",
    },
    "derivation_rules": [
        "capability numbers come from damage-ladder records (constructed reference, not physical truth) on dev families {C1,C3} only",
        "MT and CO geometries are layered separately and never merged or cross-paired",
        "no ranking, no selection, no pass/fail thresholding in this export; undetermined determination uses the reference-uncertainty ruler (reference_uncertainty_budget_v0_1) at G4 step 4",
        "numerical uncertainty (step 1 budget) and mission tolerance (this freeze) are accounted separately",
        "all capability claims carry simulation-domain validity only; no field extrapolation (G6/G7 open)",
    ],
    "relation_to_gates": "G4 NOT relieved by this freeze; physical_acceptance_threshold stays null; candidates stay undetermined; training_eligible=false",
}

EXPECTED_EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'


def sha256_file(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def main():
    actual = sha256_file(INPUTS['event_table'])
    assert actual == EXPECTED_EVENT_TABLE_SHA256, f'event table gate mismatch: {actual}'
    contract = dict(CONTRACT)
    contract['inputs_sha256'] = {k: sha256_file(v) for k, v in INPUTS.items()}
    contract['inputs'] = INPUTS
    text = json.dumps(contract, indent=1, sort_keys=True, ensure_ascii=False)
    OUT.write_bytes(text.encode('utf-8'))
    print('wrote', OUT)
    print('sha256:', hashlib.sha256(OUT.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
