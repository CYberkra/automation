"""Meaningful operator/reference counterexamples on constructed arrays only.

Run: python scripts/check_operator_contract.py --output <new-json-path>
No raw field files, solver, training, or output-dependent evaluation windows.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from research_operator_contract import (
    VERSION, ConfigUnavailable, apply_configuration, catalogue, local_mean_control,
)


def run():
    checks = []

    def record(name, condition, **values):
        if not condition:
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **values})

    def rejected(fn, exception=ValueError, reason=None):
        try:
            fn()
        except exception as exc:
            return reason is None or reason in str(exc)
        return False

    rng = np.random.default_rng(2026092403)
    x = rng.normal(size=(16, 8))
    original = x.copy()
    configs = catalogue()
    record("O01_catalogue", len(configs) == 27 and len({c["id"] for c in configs}) == 27
           and sum(c["order"] == "GB" for c in configs) == 6,
           entries=27, background_then_gain=21, gain_then_background=6)
    results = [apply_configuration(x, c["id"]) for c in configs]
    record("O02_all_configurations_preserve_input_axes_finiteness",
           np.array_equal(x, original) and all(r["output"].shape == x.shape
           and np.isfinite(r["output"]).all() and len(r["steps"]) == 2 for r in results))
    record("O03_identity_exact", np.array_equal(results[0]["output"], x))
    layer = np.outer(np.array([0., 1., -0.5, 0.]), np.ones(8))
    partial = apply_configuration(layer, "B1_G1_BG")["output"]
    record("O04_partial_mean_damages_real_layer", np.allclose(partial, .75 * layer),
           protected_layer_amplitude_ratio=.75, signed_nrmse=.25)
    gain = apply_configuration(x, "B0_G4_BG")
    record("O05_gain_shared_axis_and_endpoints", gain["gain_curve"][0] == 1
           and gain["gain_curve"][-1] == 4
           and np.allclose(gain["output"], x * gain["gain_curve"][:, None]))
    max_commutation_error = max(float(np.max(np.abs(
        apply_configuration(gain["output"], f"B{b}_G1_BG")["output"] -
        apply_configuration(x, f"B{b}_G4_BG")["output"]))) for b in (1, 2, 3))
    record("O06_all_partial_means_commute_with_shared_gain", max_commutation_error < 1e-12,
           maximum_absolute_difference=max_commutation_error)
    gb = apply_configuration(x, "B4_G4_GB")
    record("O07_reverse_order_has_no_actual_pre_gain_background_result", gb["pre_gain"] is None
           and gb["steps"][0]["operator"] == "G" and np.allclose(
               gb["audit_before_gain"] * gb["gain_curve"][:, None], gb["output"]))
    tied = np.diag([2., 2., 1., .2])
    record("O08_cutoff_tie_rejected_but_complete_cluster_allowed",
           rejected(lambda: apply_configuration(tied, "B4_G1_BG"), ConfigUnavailable, "gap_unresolved")
           and np.allclose(apply_configuration(tied, "B5_G1_BG")["output"], np.diag([0., 0., 1., .2])))
    near = np.diag([2., 2. * (1 - 1e-10), 1., .2])
    record("O09_near_tie_guard_and_sensitivity", rejected(lambda: apply_configuration(near, "B4_G1_BG"),
           ConfigUnavailable, "gap_unresolved") and np.isfinite(
               apply_configuration(near, "B4_G1_BG", gap_rtol=1e-12)["output"]).all(),
           default_gap_rtol=1e-8, deliberately_relaxed_guard=1e-12,
           note="guard changes availability; no physical tolerance established")
    redundant = apply_configuration(layer, "B5_G1_BG")
    record("O10_numerical_nullspace_does_not_create_arbitrary_labels",
           redundant["steps"][0]["diagnostics"]["status"] == "rank_redundant"
           and np.allclose(redundant["output"], apply_configuration(layer, "B4_G1_BG")["output"]))
    record("O11_zero_input_is_valid_but_not_informative", all(np.array_equal(
           apply_configuration(np.zeros_like(x), c["id"])["output"], np.zeros_like(x)) for c in configs))
    mask = np.ones(x.shape, dtype=bool); mask[2, 3] = False
    kept = apply_configuration(x, "B0_G1_BG", mask=mask)
    record("O12_mask_preserved_or_explicitly_rejected", np.array_equal(kept["mask"], mask)
           and np.array_equal(kept["output"], x) and rejected(
               lambda: apply_configuration(x, "B1_G1_BG", mask=mask), ConfigUnavailable)
           and rejected(lambda: apply_configuration(x, "B0_G2_BG", mask=mask), ConfigUnavailable))
    invalid = x.copy(); invalid[0, 0] = np.nan
    record("O13_bad_input_rejected", all(rejected(fn) for fn in [
        lambda: apply_configuration(invalid, "B0_G1_BG"),
        lambda: apply_configuration(x.astype(complex), "B0_G1_BG"),
        lambda: apply_configuration(x.astype(int), "B0_G1_BG"),
        lambda: apply_configuration(x, "B0_G1_BG", domain="envelope"),
        lambda: apply_configuration(x.T[None, :], "B0_G1_BG"),
        lambda: apply_configuration(x, "unknown"),
        lambda: apply_configuration(x, "B0_G1_BG", mask=np.ones(x.shape)),
    ]))
    record("O14_rank_dimension_not_silently_clamped", rejected(
           lambda: apply_configuration(np.ones((3, 8)), "B6_G1_BG"), ConfigUnavailable))
    record("O15_gain_overflow_not_clipped", rejected(
           lambda: apply_configuration(np.full((4, 4), 1e308), "B0_G4_BG"),
           ConfigUnavailable, "numeric_failure"))
    permutation = rng.permutation(x.shape[1])
    max_equivariance_error = max(float(np.max(np.abs(
        apply_configuration(2.5 * x[:, permutation], c["id"])["output"] -
        2.5 * apply_configuration(x, c["id"])["output"][:, permutation]))) for c in configs)
    record("O16_scale_and_trace_permutation_equivariance", max_equivariance_error < 1e-11,
           maximum_absolute_difference=max_equivariance_error)
    constant = np.ones((4, 11))
    record("O17_local_mean_edges_normalized", np.array_equal(local_mean_control(constant, 5),
           np.zeros_like(constant)), note="also deletes a laterally constant protected layer")
    impulse = np.zeros((4, 21)); impulse[2, 10] = 1
    local = local_mean_control(impulse, 5)
    record("O18_local_mean_creates_flanking_residuals", np.allclose(local[2, 8:13], [-.2, -.2, .8, -.2, -.2]),
           centre=.8, four_neighbour_residuals=-.2, note="not new physical reflectors")

    # States are declared from construction and independent reference evidence.
    # They are NOT inferred by inspecting the processed output.
    signal = np.zeros((16, 8)); signal[10, :] = .1
    clutter = np.zeros_like(signal); clutter[2, :] = 1
    roi = np.zeros(signal.shape, dtype=bool); roi[9:12, :] = True
    observed = signal + clutter
    filtered = apply_configuration(observed, "B3_G1_BG")["output"]
    error = float(np.linalg.norm((filtered - signal)[roi]) / np.linalg.norm(signal[roi]))
    record("R01_isolated_allows_defined_event_error", error == 1,
           state="isolated", waveform_nrmse=error, roi_defined_before_processing=True,
           evidence="known additive arrays, temporally disjoint support")
    mixed = signal + 2 * signal
    record("R02_mixed_no_pure_event_label", np.count_nonzero(mixed) > 0,
           state="mixed", waveform_nrmse=None, eligible=False,
           reason="common-support contributions; retain total-output diagnostic only")
    # Two semantic worlds have the same observation and demand different outputs.
    world_a = {"signal": signal, "clutter": clutter}
    world_b = {"signal": clutter, "clutter": signal}
    record("R03_ambiguous_identical_observation_different_targets",
           np.array_equal(world_a["signal"] + world_a["clutter"], world_b["signal"] + world_b["clutter"])
           and not np.array_equal(world_a["signal"], world_b["signal"]), state="ambiguous",
           waveform_nrmse=None, eligible=False, reason="algebraic non-identifiability, not a geological simulation")
    negative = np.zeros_like(signal); negative[-1, :] = .01
    enlarged = apply_configuration(negative, "B0_G4_BG")["output"]
    before_energy, after_energy = float(np.sum(negative ** 2)), float(np.sum(enlarged ** 2))
    record("R04_absent_uses_energy_not_zero_reference_nrmse", np.isclose(after_energy / before_energy, 16),
           state="absent", waveform_nrmse=None, residual_energy_before=before_energy,
           residual_energy_after=after_energy, energy_ratio=after_energy / before_energy,
           reason="predeclared target-free control; fixed last-sample window")
    reference_amplitude, discrepancy = .01, .02
    record("R05_unresolved_is_not_absence", discrepancy > reference_amplitude,
           state="numerically_unresolved", waveform_nrmse=None, eligible=False,
           reference_amplitude=reference_amplitude, constructed_discrepancy=discrepancy,
           reason="constructed credibility failure; no solver convergence experiment")
    root = Path(__file__).resolve().parent
    return {"schema": "operator-contract-checks/0.1", "operator_version": VERSION,
            "seed": 2026092403, "numpy_version": np.__version__,
            "evidence_level": "deterministic_array_contract_and_counterexamples_only",
            "field_data_used": False, "physical_simulation": False, "training_run": False,
            "file_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            [Path(__file__), root / "research_operator_contract.py"]},
            "checks": checks}


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
                      "evidence_level": result["evidence_level"]}))
