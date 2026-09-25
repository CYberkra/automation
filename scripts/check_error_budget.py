"""Error budget certificate checks on constructed arrays; no field data, FDTD, or training.

Run: python scripts/check_error_budget.py --output <new-json-path>
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from research_operator_contract import apply_configuration, catalogue, ConfigUnavailable
from research_evaluation_contract import waveform_metrics
from research_error_budget import (
    configuration_output_budget, difference_budget, svd_perturbation_bound, _half_ulp,
)


def _twoprod(a, b):
    """Dekker TwoProd: exact error of the float64 product a*b (no FMA needed)."""
    p = a * b
    splitter = 134217729.0  # 2**27 + 1
    ah = (splitter * a) - (splitter * a - a)
    al = a - ah
    bh = (splitter * b) - (splitter * b - b)
    bl = b - bh
    err = ((ah * bh - p) + ah * bl + al * bh) + al * bl
    return p, err


def run():
    checks = []

    def check(name, ok, **details):
        if not bool(ok):
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **details})

    rng = np.random.default_rng(2026092501)
    x = rng.normal(size=(48, 24))

    # C01: identity chain on an exact constructed input certifies a zero budget.
    identity = configuration_output_budget(x, "B0_G1_BG")
    check("identity_zero_budget", identity["available"] and identity["budget_full"] == 0.0
          and "input:exact_construction" in identity["certificates"])

    # C02: mean certificate bounds the actual error against an fsum reference.
    mean_cfg = configuration_output_budget(x, "B2_G1_BG")
    y = mean_cfg["output"]
    m_ref = np.array([math.fsum(row.tolist()) / x.shape[1] for row in x])[:, None]
    y_ref = x - 0.5 * m_ref
    measured = float(np.linalg.norm(y - y_ref))
    slack = float(np.linalg.norm(_half_ulp(y_ref)))
    check("mean_certificate_conservative", measured <= mean_cfg["budget_full"] + slack
          and mean_cfg["budget_full"] > 0, measured=measured, budget=mean_cfg["budget_full"])

    # C03: the certificate scales with the input units.
    scaled = configuration_output_budget(1e6 * x, "B2_G1_BG")
    check("mean_budget_scales_with_units",
          np.isclose(scaled["budget_full"], 1e6 * mean_cfg["budget_full"], rtol=1e-12))

    # C04: gain certificate bounds the exact product error (Dekker TwoProd).
    gain_cfg = configuration_output_budget(x, "B0_G2_BG")
    curve = np.power(2.0, np.arange(x.shape[0]) / (x.shape[0] - 1))[:, None]
    exact_err = np.empty_like(x)
    for index in np.ndindex(x.shape):
        p, err = _twoprod(float(curve[index[0], 0]), float(x[index]))
        exact_err[index] = abs(gain_cfg["output"][index] - (p + err))
    check("gain_certificate_conservative",
          float(np.linalg.norm(exact_err)) <= gain_cfg["budget_full"],
          measured=float(np.linalg.norm(exact_err)), budget=gain_cfg["budget_full"])

    # C05/C09: SVD certificate present and conservative on a well-separated array.
    u_r, _ = np.linalg.qr(rng.normal(size=(48, 3)))
    v_r, _ = np.linalg.qr(rng.normal(size=(24, 3)))
    xs_well = u_r @ np.diag([3.0, 1.0, 0.3]) @ v_r.T + 1e-3 * rng.normal(size=(48, 24))
    svd_cfg = configuration_output_budget(xs_well, "B4_G1_BG")
    check("svd_budget_available_large_gap", svd_cfg["available"] and svd_cfg["budget_full"] > 0)
    xs = xs_well / float(np.max(np.abs(xs_well)))
    uu, ss, vt = np.linalg.svd(xs, full_matrices=False)
    removed_ref = np.empty_like(xs)
    for index in np.ndindex(xs.shape):
        removed_ref[index] = math.fsum([float(uu[index[0], 0] * ss[0]) * float(vt[0, index[1]])])
    y_ref = (xs - removed_ref) * float(np.max(np.abs(xs_well)))
    measured = float(np.linalg.norm(svd_cfg["output"] - y_ref))
    check("svd_certificate_conservative", measured <= svd_cfg["budget_full"],
          measured=measured, budget=svd_cfg["budget_full"])

    # C06: perturbation formula b_i = r_i + 2*sqrt(2)*sigma_i*sinTheta_i, sinTheta = sqrt(2) r/gap.
    bound, missing = svd_perturbation_bound([1e-12, 2e-12], [3.0, 1.0, 0.3], 2)
    expected_0 = 1e-12 + 2 * math.sqrt(2) * 3.0 * (math.sqrt(2) * 1e-12 / (3.0 - 1.0))
    check("svd_perturbation_formula", missing is None
          and np.isclose(bound["components"][0]["output_perturbation_bound"], expected_0)
          and np.isclose(bound["components"][1]["gap"], 0.7))

    # C07: vacuous bound (gap too small for the residual) is missing, never forced.
    bound, missing = svd_perturbation_bound([1e-3], [1.0, 0.999], 1)
    check("svd_vacuous_bound_missing", bound is None and missing == "svd_subspace_bound_vacuous")

    # C08: operator-level unresolved gap propagates as a missing budget with reason.
    tied = u_r @ np.diag([1.0, 1.0 - 5e-9, 0.5]) @ v_r.T
    tied_cfg = configuration_output_budget(tied, "B4_G1_BG")
    check("svd_gap_unavailable_missing", not tied_cfg["available"]
          and tied_cfg["reason"].startswith("operator_unavailable:svd_cutoff_gap_unresolved"))

    # C10: two-step chain (mean then gain x2) is conservative against an exact reference.
    chain = configuration_output_budget(x, "B1_G2_BG")
    y = chain["output"]
    m_ref = np.array([math.fsum(row.tolist()) / x.shape[1] for row in x])[:, None]
    stage = x - 0.25 * m_ref
    exact_err = np.empty_like(x)
    for index in np.ndindex(x.shape):
        p, err = _twoprod(float(curve[index[0], 0]), float(stage[index]))
        exact_err[index] = abs(y[index] - (p + err))
    slack = float(np.linalg.norm(0.25 * _half_ulp(m_ref) * curve)
                  + np.linalg.norm(_half_ulp(stage) * curve))  # fsum-path slacks
    check("chain_composition_conservative",
          float(np.linalg.norm(exact_err)) <= chain["budget_full"] + slack + 1e-300,
          measured=float(np.linalg.norm(exact_err)), budget=chain["budget_full"])

    # C11/C12: difference budget is the triangle sum and reports amplification.
    diff = difference_budget(1e-3, 2e-3, 1.0, 1.0, 0.01)
    check("difference_triangle", diff["budget"] == 3e-3
          and diff["certificate"] == "triangle_inequality")
    check("difference_amplification", np.isclose(diff["cancellation_amplification"], 200.0))

    # C13: restored cosine on a constructed scene matches the direct value.
    traces = np.arange(24)
    flat = np.exp(-((np.arange(48)[:, None] - 24.0) / 3.0) ** 2) * np.ones((1, 24))
    clutter = 2.0 * np.exp(-((np.arange(48)[:, None] - 8.0) / 3.0) ** 2) * np.ones((1, 24))
    scene_x, scene_s = flat + clutter, flat
    restored = configuration_output_budget(scene_x, "B1_G1_BG")
    full_mask = np.ones(scene_x.shape, dtype=bool)
    with_budget = waveform_metrics(restored["output"], scene_s, full_mask,
                                   reference_kind="complete_clean", state="mixed",
                                   scope="complete",
                                   direction_error_bound=restored["budget_full"])
    yn, sn = restored["output"] / 1.0, scene_s / 1.0
    direct = float((yn.ravel() @ sn.ravel())
                   / np.sqrt((yn.ravel() @ yn.ravel()) * (sn.ravel() @ sn.ravel())))
    check("restored_cosine_matches_direct", with_budget["metrics"]["signed_cosine"] is not None
          and np.isclose(with_budget["metrics"]["signed_cosine"], direct),
          cosine=with_budget["metrics"]["signed_cosine"], budget=restored["budget_full"])

    # C14: outputs at or below the certified budget stay missing with reasons.
    zero = configuration_output_budget(scene_x, "B3_G1_BG")
    zero_mask = waveform_metrics(zero["output"], scene_s, full_mask,
                                 reference_kind="complete_clean", state="mixed",
                                 scope="complete", direction_error_bound=zero["budget_full"])
    check("cancelled_output_stays_missing", zero_mask["metrics"]["signed_cosine"] is None
          and zero_mask["metric_reasons"]["signed_cosine"]
          in ("zero_output", "output_direction_numerically_unresolved"),
          reason=zero_mask["metric_reasons"]["signed_cosine"],
          budget=zero["budget_full"],
          output_norm=float(np.linalg.norm(zero["output"])))

    # C15/C16: FDTD-sourced and unknown provenance never receive a budget.
    fdtd = configuration_output_budget(scene_x, "B1_G1_BG", input_provenance="fdtd_unresolved")
    check("fdtd_provenance_missing", not fdtd["available"]
          and fdtd["reason"] == "fdtd_numerically_unresolved")
    unknown = configuration_output_budget(scene_x, "B1_G1_BG", input_provenance="mystery")
    check("unknown_provenance_missing", not unknown["available"]
          and unknown["reason"] == "unknown_input_provenance")

    # C17: unit rescaling x1e6 changes no availability decision or cosine value.
    big = configuration_output_budget(1e6 * scene_x, "B1_G1_BG")
    with_budget_big = waveform_metrics(big["output"], 1e6 * scene_s, full_mask,
                                       reference_kind="complete_clean", state="mixed",
                                       scope="complete",
                                       direction_error_bound=big["budget_full"])
    check("unit_scale_invariant_decision",
          (with_budget_big["metrics"]["signed_cosine"] is None)
          == (with_budget["metrics"]["signed_cosine"] is None)
          and np.isclose(with_budget_big["metrics"]["signed_cosine"],
                         with_budget["metrics"]["signed_cosine"]))

    # C18: window budget never exceeds the full budget.
    window = np.zeros(scene_x.shape, dtype=bool)
    window[16:32, :] = True
    win = configuration_output_budget(scene_x, "B1_G1_BG", mask=window)
    check("window_budget_within_full", 0 <= win["budget_window"] <= win["budget_full"])

    # C19: invalid declared input bounds are rejected.
    try:
        configuration_output_budget(x, "B0_G1_BG", input_error_bound=float("nan"))
        rejected = False
    except ValueError:
        rejected = True
    check("invalid_input_bound_rejected", rejected)

    # C20: nonzero budgets always carry a nonempty certificate chain.
    all_certified = all(
        configuration_output_budget(x, c["id"])["certificates"]
        for c in catalogue() if c["end_gain"] == 1
        and configuration_output_budget(x, c["id"])["available"])
    check("certificates_nonempty", all_certified)

    root = Path(__file__).parent
    return {"schema": "error-budget-checks/1", "checks": checks,
            "evidence_level": "constructed_array_certificates_only",
            "numpy": np.__version__, "physical_thresholds_frozen": False,
            "fdtd_executed": False, "training_executed": False,
            "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                              [Path(__file__), root / "research_error_budget.py",
                               root / "research_operator_contract.py",
                               root / "research_evaluation_contract.py"]}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"passed": len(result["checks"]), "output": str(args.output)}))
