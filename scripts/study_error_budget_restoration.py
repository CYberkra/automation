"""Restore direction diagnostics on the constructed damage-pilot arrays.

Budgets come from computed certificates (exact construction + fsum
recomputation + Wedin/Davis-Kahan for SVD). No field data, FDTD, or training.
Writes one JSON with per-row budgets, certificate types, and restored or
still-missing signed cosines. Run: python scripts/study_error_budget_restoration.py
--output <new-json-path>
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import study_damage_pilot
from research_error_budget import configuration_output_budget
from research_evaluation_contract import waveform_metrics

MISSING_REASONS = {"zero_output", "output_direction_numerically_unresolved"}


def run():
    checks, rows = [], []
    pilot, arrays = study_damage_pilot.run()

    def check(name, ok, **details):
        if not bool(ok):
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **details})

    def evaluate(y, s, mask, budget):
        return waveform_metrics(y, s, mask, reference_kind="complete_clean", state="mixed",
                                scope="complete", direction_error_bound=budget)

    scenes = sorted({r["scene"] for r in pilot["candidates"]})
    for scene in scenes:
        x, s = arrays[f"{scene}_input"], arrays[f"{scene}_reference"]
        event_masks = [arrays[f"{scene}_mask{i}"] for i in range(100)
                       if f"{scene}_mask{i}" in arrays]
        for row in [r for r in pilot["candidates"] if r["scene"] == scene]:
            cid = row["configuration"]
            if not row["available"]:
                rows.append({"scene": scene, "configuration": cid, "available": False,
                             "reason": row["reason"]})
                continue
            budget = configuration_output_budget(x, cid)
            assert budget["available"]
            y = budget["output"]
            full_mask = np.ones(x.shape, dtype=bool)
            whole = evaluate(y, s, full_mask, budget["budget_full"])
            entry = {"scene": scene, "configuration": cid, "available": True,
                     "budget_full": budget["budget_full"],
                     "certificate": " + ".join(sorted(set(
                         c if isinstance(c, str) else "svd_components" for c in
                         budget["certificates"]))),
                     "whole": {**whole["metrics"],
                               "cosine_missing_reason": whole["metric_reasons"]["signed_cosine"],
                               "output_l2": float(np.linalg.norm(y)),
                               "budget": budget["budget_full"]},
                     "events": []}
            for mask in event_masks:
                windowed = configuration_output_budget(x, cid, mask=mask)
                metrics = evaluate(y, s, mask, windowed["budget_window"])
                entry["events"].append(
                    {**metrics["metrics"],
                     "cosine_missing_reason": metrics["metric_reasons"]["signed_cosine"],
                     "output_l2": float(np.linalg.norm(y[mask])),
                     "budget": windowed["budget_window"]})
            rows.append(entry)

    evaluated = [r for r in rows if r["available"]]
    check("all_rows_budgeted", len(evaluated) == 42 and all(
        np.isfinite(r["budget_full"]) and r["budget_full"] >= 0 for r in evaluated))
    allowed = {"input:exact_construction", "identity_exact_copy",
               "mean_exact_recomputation_fsum", "gain_by_exact_one",
               "gain_counted_half_ulp", "svd_arith_fsum+wedin_davis_kahan_sin_theta",
               "svd_components"}
    check("certificate_types_allowed", all(
        set(r["certificate"].split(" + ")) <= allowed for r in evaluated))

    restored = missing = 0
    for r in evaluated:
        for section in [r["whole"], *r["events"]]:
            cosine, reason = section["signed_cosine"], section["cosine_missing_reason"]
            if cosine is None:
                missing += 1
                if reason not in MISSING_REASONS:
                    raise AssertionError(f"unexpected_missing_reason:{reason}")
            else:
                restored += 1
                if not -1 - 1e-12 <= cosine <= 1 + 1e-12:
                    raise AssertionError("cosine_out_of_range")
    check("cosine_values_valid", True, restored=restored, still_missing=missing)

    # Whole-mask restoration reproduces the direct cosine computed without a budget.
    r = evaluated[0]
    x, s = arrays[f"{r['scene']}_input"], arrays[f"{r['scene']}_reference"]
    y = configuration_output_budget(x, r["configuration"])["output"]
    direct = float((y.ravel() @ s.ravel())
                   / np.sqrt((y.ravel() @ y.ravel()) * (s.ravel() @ s.ravel())))
    check("restored_cosine_matches_direct", np.isclose(r["whole"]["signed_cosine"], direct))

    # Unit rescaling changes no decision on a sampled scene.
    big = configuration_output_budget(1e6 * x, r["configuration"])
    big_metrics = evaluate(big["output"], 1e6 * s, np.ones(x.shape, bool), big["budget_full"])
    check("unit_scale_invariant_decisions",
          (big_metrics["metrics"]["signed_cosine"] is None)
          == (r["whole"]["signed_cosine"] is None)
          and np.isclose(big_metrics["metrics"]["signed_cosine"], r["whole"]["signed_cosine"]))

    budgets = [r["budget_full"] for r in evaluated]
    result = {"schema": "error-budget-restoration/1",
              "evidence_level": "constructed_array_certificates_only",
              "field_data_used": False, "fdtd_executed": False, "training_executed": False,
              "budget_source": "exact_construction+fsum_recomputation+wedin_davis_kahan",
              "numpy": np.__version__, "checks": checks, "rows": evaluated,
              "statistics": {
                  "rows": len(evaluated),
                  "cosine_evaluations": restored + missing,
                  "cosines_restored": restored, "cosines_still_missing": missing,
                  "budget_min": min(budgets), "budget_max": max(budgets)},
              "limitations": ["budgets certify numerical identifiability only, not physical accuracy",
                              "window budgets carry conservative full-array earlier-step terms",
                              "constructed arrays only; no FDTD or field-data statement"]}
    root = Path(__file__).parent
    result["source_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                               [Path(__file__), root / "research_error_budget.py",
                                root / "study_damage_pilot.py",
                                root / "research_evaluation_contract.py",
                                root / "research_operator_contract.py"]}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"checks_passed": len(result["checks"]),
                      "statistics": result["statistics"], "output": str(args.output)}))
