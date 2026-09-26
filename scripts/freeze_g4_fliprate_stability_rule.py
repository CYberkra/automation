"""Freeze the G4 flip-rate stability rule v0.1 (deterministic, byte-reproducible).

Generates configs/research/g4_fliprate_stability_rule_v0.1.json from this script's
frozen content plus SHA-256 of registered inputs (G4 step-4 flip-rate results,
evaluation protocol v0.2, G4 calibration plan, mission tolerance v0.2). Re-running
must produce a byte-identical file.

The rule freezes the boundary of the preregistered reading "翻转频繁 ⇒ 保留集合 /
partial_only": a pairwise candidate ordering inside one stratum counts as determined
only under ZERO observed flips on all three axes of the frozen flip-rate analysis.
Numeric content of the statistical certificate (rule-of-three bound and the implied
separation in units of the reference-uncertainty magnitude) is COMPUTED here, not
hand-transcribed. No numerical threshold is derived from candidate performance,
this freeze makes no ranking verdict, and it does not relieve G4.
"""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/g4_fliprate_stability_rule_v0.1.json'

INPUTS = {
    'step4_fliprate': 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_r1/results.json',
    'evaluation_protocol': 'docs/research/2026-09-24_evaluation_and_labels_v0.2.md',
    'g4_plan': 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md',
    'mission_tolerance': 'configs/research/g4_mission_tolerance_v0.2.json',
}

# Certificate parameters (frozen inputs to the computation, not tuned).
CONFIDENCE_LEVEL = 0.95
K_REPLICATES = 200


def certificate_numbers():
    """Return the numbers written into statistical_certificate, computed here.

    Rule of three: observing zero flips in K replicates bounds the per-replicate
    flip probability at p <= 1 - (1 - confidence)^(1/K) with the stated confidence.
    Inverting the triangular-difference flip model p(d) = (2 - d)^2 / 8 (separation d
    in units of eps) gives the implied minimum separation d = 2 - sqrt(8 p).
    """
    p_upper = 1.0 - math.exp(math.log(1.0 - CONFIDENCE_LEVEL) / K_REPLICATES)
    # internal consistency: zero flips survive with exactly the stated confidence at p_upper
    assert abs((1.0 - p_upper) ** K_REPLICATES - (1.0 - CONFIDENCE_LEVEL)) < 1e-12
    d_over_eps = 2.0 - math.sqrt(8.0 * p_upper)
    assert 0.0 < d_over_eps < 2.0, 'implied separation outside the triangular support'
    return p_upper, d_over_eps


P_UPPER_95, D_OVER_EPS = certificate_numbers()

CONTRACT = {
    "schema_version": "g4_fliprate_stability_rule_v0.1",
    "frozen_date": "2026-09-26",
    "frozen_by": "kimi_autonomous_2026-09-26_22:23_user_delegated_determination",
    "status": "frozen",
    "plan": "docs/research/2026-09-26_g4_calibration_plan_v0.1.md",
    "user_approval": "user 2026-09-26 22:23 '这两个你好好调研研究确立一下呢？' — two open step-4 questions delegated for determination (shift-grid sign symmetry; flip-rate stability rule)",
    "unit": "pairwise candidate ordering within one stratum of the frozen flip-rate analysis (stratum = quantity × damage_type × mission_class × geometry × family); relations are never established cross-stratum, cross-geometry or cross-family",
    "determined_order_rule": {
        "statement_zh": "一对候选的排序记为 determined 当且仅当在冻结翻转率分析的全部三条轴上观察到零翻转：级别轴在任何两个级别间无排序符号变化；window_jackknife 在任何一个留一事件 replicate 对全量排序无翻转；value_jitter 在 K=200 个 replicate 的任何一个中无翻转。任一轴观察到翻转 ⇒ 该对排序 undetermined（与预注册'翻转频繁 ⇒ 保留集合/partial_only'衔接：本规则是全有或全无的保守读法，零翻转=该对差异在全部扰动下稳定）。tie（精确 float 相等）按翻转率口径排除并单独计数，不构成 determined 排序。",
        "axes": {
            "level_axis": "zero sign changes across all frozen ladder level pairs",
            "window_jackknife": "zero flips vs the full-stratum ordering across all leave-one-event-out replicates",
            "value_jitter": "zero flips across all K replicates against the full-stratum ordering"
        },
        "all_or_nothing_conservative_reading": "the rule implements the step-1 preregistered ruler (差异与参考不确定性同量级 ⇒ undetermined) in its maximal form: ANY observed instability under reference-uncertainty-magnitude perturbation demotes the pair to undetermined"
    },
    "statistical_certificate": {
        "jitter_model": "independent uniform jitter U(-eps, +eps) on each candidate value; for pair separation d the per-replicate flip probability is p(d) = max(0, (2*eps - d)^2 / (8*eps^2)) (triangular difference distribution)",
        "zero_flips_bound": "0 flips in K=200 replicates certifies p <= 1 - 0.05^(1/200) ≈ {rule_of_three:.4f} at 95% confidence (rule of three)".format(rule_of_three=P_UPPER_95),
        "equivalent_separation": "the certificate implies pair separation d >= {d_over_eps:.3f} * eps, i.e. the pair difference exceeds ~1.65x the reference-uncertainty magnitude — the quantitative form of the step-1 ruler".format(d_over_eps=D_OVER_EPS),
        "ruler_link": "quantifies the preregistered ruler of reference_uncertainty_budget_v0_1: differences within ~1.65x the reference uncertainty are demoted to undetermined; the factor is certified by the zero-flip observation itself, not tuned"
    },
    "stratum_level_rule": "at G4 step 5 a candidate may be reported as the unique stable top of a stratum only if its ordering against EVERY other candidate in that stratum is determined; otherwise the stratum keeps set/partial_only semantics per evaluation protocol v0.2 §7 (unique_admissible vs partial_only). This rule governs ordering stability ONLY; admissibility additionally needs the step-5 frozen threshold constraints.",
    "anti_overfit": "rule frozen BEFORE any application that produces verdicts; the zero-flip cutoff is not tuned to candidate outcomes; changes only by versioned re-freeze; Cawley-Talbot precedent (evaluation protocol v0.2 §9): the selection procedure itself must be included in any performance evaluation",
    "what_this_rule_does_not_do": [
        "produces no threshold values",
        "makes no candidate selection or ranking verdict",
        "generates no labels",
        "amplitude class remains G1/G6-limited",
        "applies only to the D and N_b diagnostic quantities (constructed reference, simulation-domain only)"
    ],
    "application": "executed once at G4 step 5 on a frozen flip-rate analysis (inputs hash chain); step-4 artifacts remain diagnostics and this freeze alone does not apply the rule"
}


def sha256_file(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def main():
    contract = dict(CONTRACT)
    contract['inputs_sha256'] = {k: sha256_file(v) for k, v in INPUTS.items()}
    contract['inputs'] = INPUTS
    text = json.dumps(contract, indent=1, sort_keys=True, ensure_ascii=False)
    OUT.write_bytes(text.encode('utf-8'))
    print('wrote', OUT)
    print('sha256:', hashlib.sha256(OUT.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
