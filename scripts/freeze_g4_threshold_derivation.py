"""Freeze configs/research/g4_threshold_derivation_v1.0.json (G4 step 5, stage S1).

User approval: 2026-09-27 "按你建议，开始做吧" — the five review points of
docs/research/2026-09-26_g4_step5_lock_program_draft.md resolved at the
draft's recommended values (D_th = 1x epsilon_stratum, N_th = 1.0,
four-prerequisite relief timing, a80 systematic term separately accounted,
all 57 candidates in scope).

This freeze LOCKS THE PROGRAM ONLY. It derives no threshold values; threshold
numbers appear for the first time at S3 (development-group execution).
Deterministic: rerunning must be byte-identical. CPU read-only, no solver.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "configs" / "research" / "g4_threshold_derivation_v1.0.json"
DATE = "2026-09-27"

EXPECTED = {
    "event_table": ("configs/research/batch2d_v1_event_table_v0.1.json",
                    "b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c"),
    "mission_tolerance_v02": ("configs/research/g4_mission_tolerance_v0.2.json",
                              "ee039fb1fab3f1147e8a5fe53809bbb8ae9d6f51aafc2feee0fb90071a6c57ea"),
    "stability_rule": ("configs/research/g4_fliprate_stability_rule_v0.1.json",
                       "8ee61d0b7777f188d84fd754e3302ade9fe678ed0eae1b7f3f0ade4c451b65f6"),
    "s4_baseline_split": ("configs/research/s4_baseline_group_v0.1.json",
                          "8676f37d6aee603b7f2481779ffe361d2878735316863da0bed0db82e16b5593"),
    "reference_budget_r1": ("artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json",
                            "38cf98bab723358dc8d75486ef6824057d42c32a9030b56f7101400fc87a8cc1"),
    "capability_v3_r1": ("artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r1/results.json",
                         "2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8"),
    "fliprate_v2_r1": ("artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json",
                       "f96d317a3e046dfcc9efb01da55d49f78ab09c6fdb562c37e4ce93c6480e558f"),
    "fliprate_v1_r1_superseded": ("artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_r1/results.json",
                                  "b30f787763aa353fd21ada0e51111b802adfa9f87bd4c1b5e004261e44a8eec9"),
}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ladder_merged_sha() -> str:
    """Recompute the damage-ladder merged stripped hash (same algorithm as
    scripts/check_damage_ladder_acceptance.py: chunk order, strip per-record
    resource, sort_keys canonical JSON)."""
    recs = []
    for k in range(8):
        d = ROOT / f"artifacts/research_checks/2026-09-26_damage_ladder_r1_c{k:02d}"
        rs = json.loads((d / "records.json").read_text(encoding="utf-8"))["records"]
        for r in rs:
            r.pop("resource", None)
        recs.extend(rs)
    assert len(recs) == 34086, len(recs)
    return hashlib.sha256(json.dumps(recs, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def main():
    inputs_sha256 = {k: sha(ROOT / p) for k, (p, _) in EXPECTED.items()}
    for k, (_, want) in EXPECTED.items():
        assert inputs_sha256[k] == want, f"{k}: {inputs_sha256[k]} != {want}"
    ladder_sha = ladder_merged_sha()
    assert ladder_sha.startswith("6bc13f63"), ladder_sha

    doc = {
        "schema_version": "g4_threshold_derivation_v1.0",
        "status": "frozen_program_only_no_threshold_values_derived",
        "frozen_date": DATE,
        "frozen_by": "kimi (S1 of the step-5 lock program)",
        "user_approval": ("user 2026-09-27 '按你建议，开始做吧' — approved execution of the step-5 lock "
                          "program with the five review points resolved at the draft's recommended values"),
        "review_points_resolution": {
            "1_D_th_margin": "D_th = 1x epsilon_stratum (NOT 2x); margin may only change before S3, never after seeing candidate performance",
            "2_N_th": "N_th = 1.0 (negative-control no-amplification; identity anchor constant, not a tuned parameter)",
            "3_relief_timing": "four prerequisites (S1 frozen program + S3 dev execution + S5 one-shot test evaluation + user sign-off); test evaluation BEFORE relief",
            "4_a80_systematic_term": "systematic bias term accounted SEPARATELY from the statistical a80 statement (Wang 2026 template; statistical and systematic terms never mixed)",
            "5_candidate_set": "all 57 configurations in scope (27 catalogue + 27 frozen L + 3 frozen R); unavailable configurations recorded as-is, never scored zero",
        },
        "locked_objects": {
            "metrics": "D (in-window nrmse vs constructed reference) and N_b (negative-control residual ratio), evaluation-runner v0.1 definitions; A/H/rho stay unavailable-as-recorded under reference_state=mixed; no new metrics",
            "threshold_derivation_program": "this file",
            "mask_program": "frozen event windows + damage rectangle (central 1/2 samples x central min(8, n_traces) traces, applied as-is, no shifting) — damage ladder v0.1 convention",
            "candidate_set": "all 57 configurations present in capability export v3; unavailable as-is",
            "features": "in-window waveforms + N_b diagnostic scalars only; NO training features defined (training eligibility is an independent gate)",
            "split": "S4 frozen: development {B2D-C1m, B2D-C3m} / test {B2D-C5m, B2D-C8m}; MT/CO stratified throughout, never merged, never cross-paired",
        },
        "inputs_sha256": inputs_sha256,
        "damage_ladder_merged_stripped_sha256": ladder_sha,
        "derivation_rules": {
            "D_constraint": {
                "formula": "candidate X passes the mission-tolerance constraint in stratum (geometry x family) iff for EVERY damage_type in {amplitude_scale, polarity_flip, sample_shift, trace_deletion} and EVERY frozen mission_relevant level of that damage type, capability_v3 row (X, damage_type, level, family, geometry, mission_relevant).D_p80 <= epsilon[geometry, family, damage_type, mission_relevant]",
                "epsilon_source": "fliprate v2 stratum value_jitter.epsilon (= sqrt(2 x median deviation_from_unity) of the stratum's event windows, calibrated from the step-1 budget shape class; 0.0300-0.0341) — the same ruler as the whole chain",
                "levels_source": "configs/research/g4_mission_tolerance_v0.2.json ladder_level_mapping.mission_relevant (sign-symmetric explicit grid)",
                "no_relaxation": "if no candidate passes in a stratum, the verdict is 'no admissible candidate in this stratum'; D_th is NEVER relaxed (protocol v0.2 step 3; POD preregistration: threshold frozen first)",
            },
            "N_constraint": {
                "formula": "candidate X passes the false-alarm constraint in (geometry x family) iff for EACH q in {2, 4}, negative_control row (X, family, geometry, nc_amplify_q=q).N_b_energy_ratio_p80 <= 1.0",
                "merged_q_scope_operationalization": "the draft's 'q in {2,4} merged scope' is frozen as the CONSERVATIVE JOINT reading: both q rows must pass (equivalent to max over q of the p80 <= 1.0); this is pinned here, before any candidate performance is consulted",
                "semantics": "N_th = 1.0 = identity anchor (no amplification on the ghost-interface/ringing negative-control band); amplifiers worse than doing nothing are inadmissible",
                "parallel_discipline": "N_b reported alongside detection, never merged into target rows",
            },
            "capability_statement": {
                "form": "per stratum x candidate: D_p80 capability value (a80-style: upper-bound statement) + EXPLICIT systematic bias term (time-of-arrival to depth conversion under velocity deviation, simulation-domain grid-chain direction magnitude), statistical and systematic terms separately accounted",
                "weak_event_band_protection": "weak-event-band levels NEVER enter the capability headline (protection clause)",
                "reported_regardless": "capability values reported whether or not constraints pass",
            },
            "stability_rule_application": {
                "rule": "configs/research/g4_fliprate_stability_rule_v0.1.json applied ONCE, on fliprate v2",
                "unit": "pairwise candidate ordering within one stratum (quantity x damage_type x mission_class x geometry x family); never cross-stratum",
                "per_pair_rederivation": ("fliprate v2 stores aggregate flip counts only; per-pair flips are re-derived at S3 from the "
                                          "frozen v2 methodology (deterministic jitter seed protocol Random(repr(stratum_key)|replicate|candidate_id), "
                                          "level-axis over the frozen signed grid, leave-one-event jackknife). ACCEPTANCE GATE: the per-pair "
                                          "re-derivation MUST reproduce every v2 aggregate (flips/comparisons/ties per stratum per axis) exactly; "
                                          "any mismatch aborts S3"),
                "determined_pair": "zero flips on ALL THREE axes (ties excluded and counted separately, per the frozen rule)",
                "stratum_verdict": "unique stable top only if determined-better against EVERY other admissible candidate in the stratum; otherwise set/partial_only per protocol v0.2 section 7; structure reported only, no deployment selection",
            },
        },
        "disclosures": [
            ("stability-rule input chain contains the SUPERSEDED fliprate v1 (b30f7877...): the rule was frozen "
             "2026-09-26 before v2 existed, so the anti-overfit discipline (rule frozen before any verdict) is intact; "
             "the rule is applied to v2 aggregates at S3 and this disclosure satisfies the 2026-09-27 audit observation"),
            ("cross-scale caveat: value_jitter epsilon is calibrated on the shape-class NRMSE scale and is applied in v2 "
             "to BOTH D (linear amplitude scale) and N_b (energy ratio, quadratic scale); determined relations in the 3 N_b "
             "strata therefore carry cross_scale_caveat=true and are reported separately from D-stratum relations"),
        ],
        "hard_limits": [
            "this freeze derives NO threshold values; threshold numbers first appear at S3",
            "G4 NOT relieved by this freeze; relief requires the four prerequisites including user sign-off",
            "no physical/training labels; training_eligible stays false; G4 relief does not unlock training",
            "constructed reference (identity of damaged input), simulation domain only; no field extrapolation",
            "amplitude class remains G1/G6-limited; reference_state stays numerically_unresolved",
            "test families {C5,C8} untouched until S4/S5; S5 executes the frozen protocol exactly once, no feedback tuning",
        ],
        "version_discipline": "any post-lock change requires a new version (v1.1) + independent verification; the test family is never re-evaluated on a new version",
    }

    OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True), encoding="utf-8", newline="\n")
    print(json.dumps({"frozen": str(OUT.relative_to(ROOT)),
                      "sha256": sha(OUT),
                      "ladder_merged_stripped_sha256": ladder_sha}, indent=1))


if __name__ == "__main__":
    main()
