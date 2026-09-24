"""Reproduce algebraic counterexamples for the two-operator research scope.

Inputs: deterministic, dimensionless numeric arrays constructed in this file.
Output: JSON of checks and measured values, never a GPR performance benchmark.
No field data, FDTD solver, model training, or physical depth interpretation.
Run: python scripts/check_processing_algebra.py --output <new-json-path>
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def nrmse(output, reference):
    norm = np.linalg.norm(reference)
    if norm == 0:
        raise ValueError("NRMSE is undefined for a zero reference")
    return float(np.linalg.norm(output - reference) / norm)


def remove_first_singular_component(x):
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    return x - s[0] * np.outer(u[:, 0], vt[0])


def run():
    rng = np.random.default_rng(20260924)
    sample = np.arange(64, dtype=np.float64)
    pulse = np.exp(-((sample - 28) / 4) ** 2) * np.cos((sample - 28) * 0.8)
    layer = np.outer(pulse, np.ones(16))
    x = layer + 0.1 * rng.normal(size=layer.shape)
    checks = []

    def record(name, condition, **values):
        if not condition:
            raise AssertionError(name)
        checks.append({"id": name, "passed": True, **values})

    record("C01_identity", nrmse(layer.copy(), layer) == 0,
           nrmse=nrmse(layer.copy(), layer))
    removed = layer - layer.mean(axis=1, keepdims=True)
    record("C02_mean_removes_constant_lateral_layer", np.linalg.norm(removed) < 1e-12,
           protected_layer_nrmse=nrmse(removed, layer))
    record("C03_zero_and_polarity_error", nrmse(np.zeros_like(layer), layer) == 1
           and nrmse(-layer, layer) == 2,
           zero_nrmse=nrmse(np.zeros_like(layer), layer),
           polarity_reversal_nrmse=nrmse(-layer, layer))
    shifted = np.zeros_like(layer)
    shifted[1:] = layer[:-1]
    shift_error = nrmse(shifted, layer)
    record("C04_one_sample_shift", shift_error > 0, nrmse=shift_error)

    noise = rng.normal(size=layer.shape) * 0.03
    snr = lambda s, n: float(10 * np.log10(np.sum(s * s) / np.sum(n * n)))
    before, after = snr(layer, noise), snr(100 * layer, 100 * noise)
    record("C05_constant_gain_preserves_component_snr", abs(after - before) < 1e-12,
           before_db=before, after_db=after, gain=100)

    mean_subtract = lambda a: a - a.mean(axis=1, keepdims=True)
    gain = np.linspace(1, 8, x.shape[0])[:, None]
    commutation_error = float(np.max(np.abs(gain * mean_subtract(x) - mean_subtract(gain * x))))
    record("C06_fixed_shared_gain_mean_commute", commutation_error < 1e-12,
           maximum_absolute_difference=commutation_error)

    diagonal = np.diag([2.0, 1.0])
    g = np.diag([1.0, 3.0])
    bg = g @ remove_first_singular_component(diagonal)
    gb = remove_first_singular_component(g @ diagonal)
    record("C07_refitted_svd_order_changes_result",
           np.allclose(bg, np.diag([0.0, 3.0])) and np.allclose(gb, np.diag([2.0, 0.0])),
           background_then_gain=bg.tolist(), gain_then_background=gb.tolist(),
           frobenius_difference=float(np.linalg.norm(bg - gb)))

    a, b = np.diag([2.0, 0.0]), np.diag([0.0, 1.0])
    nonlinear_error = float(np.linalg.norm(remove_first_singular_component(a + b)
                            - remove_first_singular_component(a) - remove_first_singular_component(b)))
    record("C08_refitted_svd_is_not_additive", abs(nonlinear_error - 1) < 1e-12,
           frobenius_difference=nonlinear_error)

    attenuation = np.exp(-np.linspace(0, 3, layer.shape[0]))[:, None]
    inverse = 1 / attenuation
    clean_error = nrmse(inverse * attenuation * layer, layer)
    amplified_noise_ratio = float(np.linalg.norm(inverse * noise) / np.linalg.norm(noise))
    record("C09_known_attenuator_inverse_amplifies_added_noise",
           clean_error < 1e-12 and amplified_noise_ratio > 1,
           noiseless_recovery_nrmse=clean_error,
           added_noise_l2_amplification=amplified_noise_ratio,
           maximum_inverse_gain=float(inverse.max()))

    rejected = False
    try:
        nrmse(layer, np.zeros_like(layer))
    except ValueError:
        rejected = True
    record("C10_zero_reference_nrmse_undefined", rejected)
    return {
        "schema": "processing-algebra-checks/1",
        "evidence_level": "algebra_and_deterministic_numeric_examples_only",
        "physical_simulation": False, "field_data_used": False, "training_run": False,
        "seed": 20260924, "numpy_version": np.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "roundoff_tolerance": 1e-12,
        "tolerance_role": "floating_point_identity_checks_not_acceptance_thresholds",
        "checks": checks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"passed": len(result["checks"]), "output": str(args.output),
                      "evidence_level": result["evidence_level"]}))
