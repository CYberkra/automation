"""Adversarial metric/label checks, using constructed arrays and fixture-only limits."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from research_evaluation_contract import waveform_metrics, energy_metrics, interval, configuration_labels
from research_operator_contract import apply_configuration, catalogue


def run():
    checks, examples = [], {}

    def record(name, condition, **values):
        if not condition:
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **values})

    def metric(y, s, mask=None, **kwargs):
        mask = np.ones(s.shape, bool) if mask is None else mask
        return waveform_metrics(y, s, mask, reference_kind=kwargs.pop("reference_kind", "complete_clean"),
                                state=kwargs.pop("state", "isolated"),
                                direction_error_bound=kwargs.pop("direction_error_bound", 0.), **kwargs)

    def row(key, d=None, r=None, available=True):
        return {"id": key, "available": available,
                "metrics": {"D": None if d is None else interval(d), "R": None if r is None else interval(r)}}

    s = np.zeros((16, 8)); s[10, :2] = [1., -1.]
    identity = metric(s, s)["metrics"]
    record("E01_identity", identity["nrmse"] == 0 and identity["amplitude_factor"] == 1)
    zero = metric(np.zeros_like(s), s)["metrics"]
    record("E02_zero_is_not_perfect", zero["nrmse"] == 1 and zero["amplitude_error"] == 1
           and zero["signed_cosine"] is None)
    flip = metric(-s, s)["metrics"]
    record("E03_polarity_loss", flip["nrmse"] == 2 and flip["signed_cosine"] == -1)
    shifted = metric(np.roll(s, 1, axis=0), s)["metrics"]
    record("E04_shift_not_aligned_away", shifted["nrmse"] > 1 and shifted["signed_cosine"] == 0)
    small = metric(.01 * s, s)["metrics"]
    record("E05_scale_invariant_shape_not_fidelity", np.isclose(small["signed_cosine"], 1)
           and np.isclose(small["shape_residual"], 0) and np.isclose(small["nrmse"], .99), **small)
    n = np.zeros_like(s); n[2, :2] = [1., -1.]
    full = metric(s + n, s, state="mixed", scope="complete")["metrics"]
    half = metric(.5 * (s + n), s, state="mixed", scope="complete")["metrics"]
    record("E06_mse_can_reward_uniform_shrinkage", half["nrmse"] < full["nrmse"]
           and half["amplitude_error"] > full["amplitude_error"],
           raw_nrmse=full["nrmse"], half_nrmse=half["nrmse"], half_amplitude_factor=half["amplitude_factor"])
    strong_weak = s + .01 * n
    weak_mask = np.zeros(s.shape, bool); weak_mask[2] = True
    overall = metric(s, strong_weak)["metrics"]["nrmse"]
    weak_error = metric(s, strong_weak, weak_mask)["metrics"]["nrmse"]
    record("E07_weak_event_not_hidden_by_global_error", overall < .011 and weak_error == 1,
           global_nrmse=overall, weak_event_nrmse=weak_error)
    absent = metric(np.zeros_like(s), np.zeros_like(s), state="absent")
    record("E08_zero_reference_masked", not absent["available"] and absent["metrics"] is None)
    zero_energy = energy_metrics(s, np.zeros_like(s), np.ones(s.shape, bool))
    record("E09_zero_baseline_uses_absolute_energy", zero_energy["energy_ratio"] is None
           and zero_energy["output_mean_square"] > 0)
    state_results = {state: metric(s, s, reference_kind="isolated_event", state=state)
                     for state in ("mixed", "ambiguous", "numerically_unresolved")}
    record("E10_uncertain_event_not_given_clean_label", all(not v["available"] for v in state_results.values()))
    mixed_complete = metric(s, s, state="mixed", scope="complete")
    record("E11_known_complete_reference_can_evaluate_mixture", mixed_complete["absolute_preservation_eligible"])
    common = np.zeros_like(s); common[5, :] = .2
    f0, f1 = common, common + s
    z0 = apply_configuration(f0, "B3_G1_BG")["output"]
    z1 = apply_configuration(f1, "B3_G1_BG")["output"]
    contrast = metric(z1 - z0, f1 - f0, reference_kind="paired_contrast", state="mixed", scope="contrast")
    common_mask = np.zeros(s.shape, bool); common_mask[5] = True
    common_error = metric(z1, common, common_mask)["metrics"]["nrmse"]
    record("E12_contrast_cannot_certify_common_layer", contrast["metrics"]["nrmse"] == 0
           and not contrast["absolute_preservation_eligible"] and common_error == 1,
           contrast_nrmse=0., deleted_common_layer_nrmse=common_error)
    gain = np.power(4., np.arange(16) / 15)[:, None]
    negative = np.zeros_like(s); negative[-1] = .01
    growth = energy_metrics(gain * negative, negative, np.ones(s.shape, bool))["energy_ratio"]
    record("E13_gain_negative_control", growth == 16, energy_growth=growth)
    attenuated = s / gain
    recovered = metric(gain * attenuated, s)["metrics"]
    overgain = metric(gain * s, s)["metrics"]
    record("E14_gain_target_must_be_declared", recovered["nrmse"] < 1e-12 and overgain["nrmse"] > 1,
           known_attenuation_recovery_error=recovered["nrmse"], unattenuated_overgain_error=overgain["nrmse"])

    # All numerical limits below are fixture assertions, not GPR acceptance values.
    limits, objectives = {"D": .1}, ["D", "R"]
    bad_zero = configuration_labels([row("identity", 0, 1), row("zero", 1, 0)], limits=limits, objectives=objectives)
    record("L01_constraints_reject_delete_all", bad_zero["unique_label"] == "identity")
    ties = configuration_labels([row("a", 0, .2), row("b", 0, .2)], limits=limits, objectives=objectives)
    record("L02_ties_not_forced_to_id", ties["status"] == "multiple_admissible" and ties["unique_label"] is None)
    tradeoff = configuration_labels([row("a", .01, .5), row("b", .05, .1)], limits=limits, objectives=objectives)
    record("L03_tradeoff_not_called_equivalence", tradeoff["status"] == "multiple_admissible"
           and not tradeoff["robust_preference_pairs"] and tradeoff["incomparable_is_not_equivalent"])
    uncertain = row("uncertain", .09, .1); uncertain["metrics"]["D"] = interval(.09, .07, .12)
    partial = configuration_labels([row("known", .01, .5), uncertain], limits=limits, objectives=objectives)
    record("L04_interval_crossing_constraint_abstains", partial["status"] == "partial_only" and partial["unique_label"] is None)
    no_gate = configuration_labels([row("a", 0, 0)], limits={"D": None}, objectives=objectives)
    record("L05_unfrozen_threshold_cannot_certify_winner", no_gate["status"] == "partial_only")
    missing = configuration_labels([row("a", 0, .1), row("missing", None, .01)], limits=limits, objectives=objectives)
    record("L06_missing_reference_not_zero_loss", missing["status"] == "partial_only" and missing["unique_label"] is None)
    none = configuration_labels([row("bad", 1, 0)], limits=limits, objectives=objectives)
    record("L07_no_safe_action_is_not_automatic_identity_label", none["status"] == "no_admissible_candidate" and none["unique_label"] is None)
    failure = row("numeric_failure", available=False)
    failed = configuration_labels([row("valid", .01, .1), failure], limits=limits, objectives=objectives)
    record("L08_failure_kept_separate", failed["candidate_records"][1]["feasibility"] == "unavailable")
    dominated = configuration_labels([row("good", .01, .1), row("worse", .03, .5)], limits=limits, objectives=objectives)
    record("L09_pareto_preference", dominated["unique_label"] == "good"
           and dominated["robust_preference_pairs"] == [{"better": "good", "worse": "worse"}])
    overlapping = [row("a", .02, .2), row("b", .03, .3)]
    overlapping[0]["metrics"]["D"] = interval(.02, .01, .04)
    robust = configuration_labels(overlapping, limits=limits, objectives=objectives)
    record("L10_nominal_order_not_robust_order", robust["status"] == "multiple_admissible")

    # End-to-end seven-background example with an isolated zero-mean lateral target.
    clutter = np.zeros_like(s); clutter[2, :] = 1.
    target_mask = np.zeros(s.shape, bool); target_mask[9:12] = True
    clutter_mask = np.zeros(s.shape, bool); clutter_mask[1:4] = True
    rows = []
    for c in catalogue():
        if c["end_gain"] != 1:
            continue
        z = apply_configuration(s + clutter, c["id"])["output"]
        d = metric(z, s, target_mask)["metrics"]["nrmse"]
        r = energy_metrics(z, s + clutter, clutter_mask)["energy_ratio"]
        # Exact controls snap only floating point identity residue at 1e-12.
        rows.append(row(c["id"], 0. if d < 1e-12 else d, 0. if r < 1e-12 else r))
    example = configuration_labels(rows, limits=limits, objectives=objectives)
    record("L11_operator_effects_generate_set_label", set(example["nondominated_supported_set"])
           == {"B3_G1_BG", "B4_G1_BG"} and example["unique_label"] is None)
    raw_record = row("raw", full["nrmse"], 1.)
    half_record = row("half", half["nrmse"], .25)
    raw_record["metrics"]["A"] = interval(full["amplitude_error"])
    half_record["metrics"]["A"] = interval(half["amplitude_error"])
    amplitude_guard = configuration_labels([raw_record, half_record], limits={"D": 2., "A": .1}, objectives=objectives)
    record("L12_amplitude_guard_stops_mse_shrinkage_reward", amplitude_guard["unique_label"] == "raw")
    examples.update({"background_known_components": example, "uncertain_constraint": partial,
                     "no_frozen_gate": no_gate, "tradeoff": tradeoff,
                     "contrast_only": contrast, "unknown_event_states": state_results,
                     "amplitude_guard": amplitude_guard})
    sources = [Path(__file__), Path(__file__).with_name("research_evaluation_contract.py"),
               Path(__file__).with_name("research_operator_contract.py")]
    return {"schema": "evaluation-label-checks/0.2", "evidence_level": "constructed_arrays_only",
            "field_data_used": False, "fdtd_executed": False, "training_executed": False,
            "numpy_version": np.__version__, "fixture_only_nrmse_limit": .1,
            "physical_acceptance_thresholds_frozen": False,
            "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            "checks": checks, "examples": examples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"passed": len(result["checks"]), "output": str(args.output),
                      "candidate_label_example": result["examples"]["background_known_components"]["nondominated_supported_set"]}))
