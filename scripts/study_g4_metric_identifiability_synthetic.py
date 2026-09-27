"""Pure algebraic check of reference-metric identifiability; no field data used."""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from research_evaluation_contract import energy_metrics, waveform_metrics


OUT = ROOT / "artifacts/research_checks/2026-09-27_g4_s6_metric_identifiability"
MASK = np.ones((2, 1), dtype=bool)


def main():
    if OUT.exists():
        raise SystemExit(f"refusing_to_overwrite:{OUT}")
    s = np.array([[1.0], [1.0]])
    c = np.array([[1.0], [-1.0]])
    x = s + c
    ps = (s @ s.T) / (s.T @ s).item()
    pc = (c @ c.T) / (c.T @ c).item()
    zero = np.zeros((2, 2), dtype=float)
    operators = {
        "identity": np.eye(2),
        "P_S": ps,
        "P_C": pc,
        "zero": zero,
        "2P_S": 2.0 * ps,
        "A_target_plus_2_background": ((s + 2.0 * c) @ s.T) / (s.T @ s).item(),
    }

    rows = []
    for name, operator in operators.items():
        y = operator @ x
        negative_output = operator @ c
        d_x = waveform_metrics(
            y, x, MASK, reference_kind="complete_clean", state="mixed", scope="complete"
        )
        d_s = waveform_metrics(
            y, s, MASK, reference_kind="complete_clean", state="mixed", scope="complete"
        )
        n_b = energy_metrics(negative_output, c, MASK)
        rows.append({
            "operator": name,
            "matrix": operator.tolist(),
            "output_for_X": y.ravel().tolist(),
            "output_for_C_negative_control": negative_output.ravel().tolist(),
            "D_X_nrmse": d_x["metrics"]["nrmse"],
            "D_S_nrmse": d_s["metrics"]["nrmse"],
            "N_b_energy_ratio": n_b["energy_ratio"],
            "waveform_metric_reference_contract": "complete_clean / mixed / complete",
        })

    by_name = {row["operator"]: row for row in rows}
    close = np.isclose
    assert close(by_name["P_S"]["D_X_nrmse"], by_name["P_C"]["D_X_nrmse"])
    assert close(by_name["P_S"]["D_S_nrmse"], 0.0)
    assert close(by_name["P_C"]["D_S_nrmse"], np.sqrt(2.0))
    assert close(by_name["P_S"]["N_b_energy_ratio"], 0.0)
    assert close(by_name["P_C"]["N_b_energy_ratio"], 1.0)
    assert close(by_name["P_S"]["D_X_nrmse"], by_name["A_target_plus_2_background"]["D_X_nrmse"])
    assert close(by_name["P_S"]["N_b_energy_ratio"], by_name["A_target_plus_2_background"]["N_b_energy_ratio"])
    assert close(by_name["P_S"]["D_S_nrmse"], 0.0)
    assert close(by_name["A_target_plus_2_background"]["D_S_nrmse"], 2.0)

    contract_path = ROOT / "scripts/research_evaluation_contract.py"
    result = {
        "schema_version": "g4_s6_metric_identifiability_synthetic_v1",
        "purpose": "mathematical_identifiability_demonstration_only",
        "data_scope": "newly_constructed_arrays_only; no simulation, measured, or test-family data",
        "inputs": {"S_target": s.ravel().tolist(), "C_background_negative_control": c.ravel().tolist(),
                   "X_mixture": x.ravel().tolist(), "mask_shape": list(MASK.shape)},
        "operator_scope": "each fixed matrix is applied identically to X and C; operators are oracle algebraic controls, not deployable algorithms or electromagnetic decompositions",
        "metrics": {"D_X": "waveform_metrics(output, X, full_mask, complete_clean/mixed/complete).nrmse",
                    "D_S": "waveform_metrics(output, S, full_mask, complete_clean/mixed/complete).nrmse",
                    "N_b": "energy_metrics(operator(C), C, full_mask).energy_ratio"},
        "operators": rows,
        "checks": {
            "P_S_vs_P_C_same_D_X_different_D_S_and_N_b": True,
            "P_S_vs_A_same_D_X_and_N_b_different_D_S": True,
            "values_are_diagnostics_only_no_threshold_application": True,
        },
        "provenance": {
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "metric_contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
            "numpy_version": np.__version__,
        },
    }
    OUT.mkdir(parents=True)
    (OUT / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT / "results.json"), "operators": len(rows), "checks": result["checks"]}, sort_keys=True))


if __name__ == "__main__":
    main()
