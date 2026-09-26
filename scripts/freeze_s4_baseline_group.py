"""Freeze the S4 baseline group contract v0.1 (deterministic, byte-reproducible).

Generates configs/research/s4_baseline_group_v0.1.json from this script's frozen
content plus SHA-256 of registered inputs (proposal doc, event table, operator
contract, evaluation contract). Re-running must produce a byte-identical file.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/s4_baseline_group_v0.1.json'

INPUTS = {
    'proposal': 'docs/research/2026-09-26_s4_baseline_group_proposal.md',
    'event_table': 'configs/research/batch2d_v1_event_table_v0.1.json',
    'operator_contract': 'scripts/research_operator_contract.py',
    'evaluation_contract': 'scripts/research_evaluation_contract.py',
}

CONTRACT = {
    "schema_version": "s4_baseline_group_v0.1",
    "frozen_date": "2026-09-26",
    "frozen_by": "kimi_autonomous_2026-09-26_19:10_continue",
    "status": "frozen",
    "proposal": "docs/research/2026-09-26_s4_baseline_group_proposal.md",
    "selection_status": {
        "F_config": "undetermined_G4",
        "R_threshold_tau": "undetermined_G4",
        "L_lambda_width_gain": "undetermined_G4",
        "O": "unavailable_no_computable_reference_no_approximation",
        "note": "G4 thresholds all null + reference_state numerically_unresolved: no quality-based selection is permitted; this freeze locks spaces, protocol and group split only",
    },
    "identity": {
        "candidate_id": "B0_G1_BG",
        "records_mt": "artifacts/research_checks/2026-09-26_eval_batch2d_mt_r1",
        "records_co": "artifacts/research_checks/2026-09-26_eval_batch2d_co_r1",
    },
    "F_protocol": {
        "candidate_set": "27-entry operator catalogue (scripts/research_operator_contract.py)",
        "gates_order": ["protective_constraints", "protected_event_fidelity", "frozen_utility_partial_order"],
        "protected_events": ["E12-S0.005", "D20m", "T0.25m", "W1m", "W2m"],
        "execution": "only_after_G4_resolution; dev group only; single confirmatory test-group evaluation",
    },
    "L_grid": {
        "operator": "local_mean_control (scripts/research_operator_contract.py :115)",
        "strength_lambda": [0.25, 0.5, 1.0],
        "window_width_odd": [5, 11, 21],
        "end_gain_q": [1, 2, 4],
        "legality_rule": "width <= n_traces else unavailable_for_this_window (no auto-shortening)",
        "order": "local_mean_then_fixed_shared_gain (BG-like; no GB variant)",
        "parameter_status": "development starting grid, not calibrated optimum",
    },
    "R_rule": {
        "input": "input window singular spectrum",
        "gaps": "d_k=(sigma_k-sigma_{k+1})/sigma_1 for k in {1,2,3}",
        "selection": "k*=argmax d_k if max(d_k) >= tau else identity",
        "threshold_candidates": [0.05, 0.1, 0.2],
        "threshold_status": "uncalibrated development candidates; tau choice undetermined (G4)",
        "numerical_rejection": "components tied below numerical rank floor are not resolved (same口径 as gap_rtol=1e-8)",
        "gain": "q=1 ablation first; full-config q choice undetermined with tau",
    },
    "groups": {
        "split_unit": "group_id (mother-model family); all variants of a family in the same group",
        "dev": ["B2D-C1m", "B2D-C3m"],
        "test": ["B2D-C5m", "B2D-C8m"],
        "rationale": "dev covers the most protected-event classes with the minimum two families for the stability-difference judgment; test retains deep/absorption-pressure families (both contain D20m, C8 contains T0.25m)",
        "inert_until": "G4 resolution; any change requires re-freeze with trace",
        "no_generalization_claim": "4 families total do not support any generalization claim",
    },
    "evaluation": {
        "same_functions_as": "scripts/research_evaluation_contract.py",
        "windows": "same frozen event table, same unshifted mapping, same time-axis convention",
        "metrics": ["N_b (energy_metrics)", "gain risk diagnostics", "svd spectrum diagnostics", "resource"],
        "not_computed": ["D/A/H/rho (reference_state numerically_unresolved)", "R_c (no legitimate pure-clutter window)"],
        "gather_semantics": "per-record version carried; MT and CO layered separately, never merged or cross-paired",
        "arrays": "persist .npy to git-ignored dir with sha256; never inline large arrays into JSON",
        "feasibility": "configuration_labels with all-null limits -> partial_only; no ranking, no unique label",
    },
    "declarations": {
        "reference_state": "numerically_unresolved",
        "physical_thresholds_null": True,
        "tolerances_null": True,
        "training_eligible": False,
        "training_labels_generated": False,
        "physical_label_eligible": False,
        "clean_truth_generated": False,
        "isolated_event_granted": False,
        "cross_family_differencing": False,
        "cross_tier_differencing": False,
        "negative_controls_merged_into_targets": False,
        "grid_tier": "BASE",
        "grid_convergence_certified": False,
        "fdtd_solver_invoked_by_baseline_runner": False,
        "gates_open": {"G1": True, "G4": True, "G5": True},
        "g1_note": "A0 3D contract still awaits user approval",
    },
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    inputs = {k: sha256_file(ROOT / v) for k, v in INPUTS.items()}
    doc = dict(CONTRACT)
    doc['inputs'] = {INPUTS[k]: v for k, v in inputs.items()}
    text = json.dumps(doc, indent=1, ensure_ascii=False, sort_keys=True) + '\n'
    OUT.write_bytes(text.encode('utf-8'))  # binary write: no CRLF translation
    print(json.dumps({'output': str(OUT.relative_to(ROOT)),
                      'sha256': hashlib.sha256(text.encode('utf-8')).hexdigest(),
                      'inputs': inputs}, indent=1))


if __name__ == '__main__':
    main()
