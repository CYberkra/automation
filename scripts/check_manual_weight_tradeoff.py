"""Illustrate hand-weighted scoring on a known-component array, not a GPR benchmark."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from research_operator_contract import apply_configuration, catalogue


def run():
    # Both events are coherent along the trace axis but have disjoint sample support.
    # No operator here has the semantic/window information supplied to this evaluator.
    signal = np.zeros((16, 8)); signal[10, :] = .1
    clutter = np.zeros_like(signal); clutter[2, :] = 1.
    x = signal + clutter
    rows = []
    for config in catalogue():
        if config["end_gain"] != 1:
            continue  # Establish background evaluation before assigning gain utility.
        result = apply_configuration(x, config["id"])["output"]
        distortion = float(np.sum((result[10] - signal[10]) ** 2) / np.sum(signal[10] ** 2))
        residual = float(np.sum(result[2] ** 2) / np.sum(clutter[2] ** 2))
        rows.append({"configuration": config["id"], "protected_squared_nrmse": distortion,
                     "clutter_residual_energy_ratio": residual})
    sweeps = []
    for weight in [0., .2, .5, .8, 1.]:
        scores = [weight * r["protected_squared_nrmse"] + (1 - weight) * r["clutter_residual_energy_ratio"] for r in rows]
        best = min(scores)
        indices = [i for i, score in enumerate(scores) if abs(score - best) < 1e-12]
        sweeps.append({"structure_weight": weight, "clutter_weight": 1 - weight,
                       "minimum_score": best,
                       "equivalent_minimisers": [rows[i]["configuration"] for i in indices],
                       "winner_squared_nrmse": rows[indices[0]]["protected_squared_nrmse"],
                       "winner_clutter_residual_energy_ratio": rows[indices[0]]["clutter_residual_energy_ratio"]})
    assert "B3_G1_BG" in sweeps[1]["equivalent_minimisers"]
    assert sweeps[2]["equivalent_minimisers"] == ["B2_G1_BG"]
    assert sweeps[3]["equivalent_minimisers"] == ["B1_G1_BG"]
    assert sweeps[4]["equivalent_minimisers"] == ["B0_G1_BG"]
    # Illustrative limit only, not a proposed geological acceptance standard.
    protected_nrmse_limit = .1
    feasible = [r["configuration"] for r in rows if r["protected_squared_nrmse"] <= protected_nrmse_limit ** 2 + 1e-12]
    assert feasible == ["B0_G1_BG"]
    return {"schema": "manual-weight-tradeoff/0.1", "date": "2026-09-24",
            "evidence_level": "constructed_known_component_array_only", "field_data_used": False,
            "fdtd_executed": False, "training_executed": False, "numpy_version": np.__version__,
            "source_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                              [Path(__file__), Path(__file__).with_name("research_operator_contract.py")]},
            "score": "w*protected_squared_nrmse+(1-w)*clutter_residual_energy_ratio",
            "evaluated_background_configurations": rows, "weight_sweeps": sweeps,
            "illustrative_gate": {"protected_nrmse_limit": protected_nrmse_limit,
                                  "feasible_configurations": feasible, "physical_standard": False},
            "checks_passed": 5,
            "conclusion": "Weight changes can reward deletion; a constraint can reveal an inadequate action set. No selector can exceed the same set's valid oracle."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"checks_passed": result["checks_passed"], "output": str(args.output),
                      "weight_sweeps": result["weight_sweeps"]}, ensure_ascii=False))
