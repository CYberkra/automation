"""Freeze the a0_3d_v1 contract package (C3 anchor A0 3D validation pair, 5 cm isotropic).

Derives the two .in files programmatically from the measured dep3d_gold_v1 5 cm
inputs (SHA-pinned), changing only the title and (TGT) the target box z range
9.75-10.25 -> 19.75-20.25 (D20m class -> D10m anchor, the studied variable).
Writes configs/research/a0_3d_v1/{*.in, cases.json, groups.json,
static_check.json, budget.json} and stages the gate with
approved_to_simulate=false (user ratification flips it in a separate step).
Deterministic: rerunning must be byte-identical. Never invokes a solver.
"""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = 299792458.0
BATCH = "a0_3d_v1"
DATE = "2026-09-27"
OUT = ROOT / "configs" / "research" / "a0_3d_v1"
DEP = ROOT / "configs" / "research" / "dep3d_gold_v1"

DEP_BG_SHA = "8705969e425886cd01ee12dab8b9cab5d21a1bdf1cdab7c2b7fa2d4f96ee2343"
DEP_TGT_SHA = "00d370d8590366269517a5b2090ba18a5b12ea91ddad7fe692d51326421ab848"

RUN_BG = "B2D-C3m-BG-3D5CM"
RUN_TGT = "B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM"

# dep3d measured supervision values (read at freeze time for provenance)
DEP_MEAS = {}
for r in ("DEP3D_5CM_BG", "DEP3D_5CM_TGT"):
    s = json.loads((ROOT / f"artifacts/research_checks/2026-09-26_{r}/supervision.json").read_text(encoding="utf-8"))
    assert s["reason"] == "completed"
    DEP_MEAS[r] = {"wall_s": s["wall_s"], "peak_job_commit_bytes": s["peak_job_commit_bytes"]}


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(case_id: str) -> int:
    h = hashlib.sha256(f"{BATCH}|B2D-C3m|{case_id}|3d5cm".encode()).digest()
    return int.from_bytes(h[:8], "big")


def build_inputs():
    bg_src = (DEP / "DEP3D_5CM_BG.in").read_bytes()
    tgt_src = (DEP / "DEP3D_5CM_TGT.in").read_bytes()
    assert sha(bg_src) == DEP_BG_SHA, "dep3d BG source changed"
    assert sha(tgt_src) == DEP_TGT_SHA, "dep3d TGT source changed"

    bg = bg_src.decode("utf-8")
    old_title_bg = "#title: DEP3D_5CM_BG DEP family 20m target, true 3D isotropic 5cm, background"
    assert bg.count(old_title_bg) == 1
    bg = bg.replace(old_title_bg,
                    "#title: B2D-C3m-BG-3D5CM C3 family A0 pair, true 3D isotropic 5cm, background")

    tgt = tgt_src.decode("utf-8")
    old_title_tgt = "#title: DEP3D_5CM_TGT DEP family 20m target, true 3D isotropic 5cm, with target"
    old_box = "#box: 3 3 9.75 5 7 10.25 wet"
    assert tgt.count(old_title_tgt) == 1 and tgt.count(old_box) == 1
    tgt = tgt.replace(old_title_tgt,
                      "#title: B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM C3 family A0 pair, true 3D isotropic 5cm, target D10m")
    tgt = tgt.replace(old_box, "#box: 3 3 19.75 5 7 20.25 wet")
    return {RUN_BG: bg, RUN_TGT: tgt}


def static_check():
    d = 0.05
    nx, ny, nz = 160, 200, 1000
    cells = nx * ny * nz
    assert cells == 32_000_000 and cells <= 64_000_000
    dt = 1.0 / (C * math.sqrt(3 * (1 / d) ** 2))
    tw = 1200e-9
    time_steps = int(math.ceil(tw / dt)) + 1
    # cross-validate against the measured dep3d 5 cm runs (same domain/grid/window)
    assert abs(dt - 9.629166007732353e-11) < 1e-24
    assert time_steps == 12464

    def divisible(x):
        return abs(x / d - round(x / d)) < 1e-9

    coords = {
        "domain_m": [8.0, 10.0, 50.0], "surface_z_m": 30.0, "cover_bottom_z_m": 27.0,
        "antenna_z_m": 45.0, "target_z_m": [19.75, 20.25], "target_y_m": [3.0, 7.0],
        "antenna_y_m": [4.35, 5.65], "target_x_m": [3.0, 5.0], "tx_rx_x_m": 4.0,
        "pml_thickness_m": 1.0,
    }
    flat = [coords["domain_m"][0], coords["domain_m"][1], coords["domain_m"][2],
            coords["surface_z_m"], coords["cover_bottom_z_m"], coords["antenna_z_m"],
            *coords["target_z_m"], *coords["target_y_m"], *coords["antenna_y_m"],
            *coords["target_x_m"], coords["tx_rx_x_m"], coords["pml_thickness_m"]]
    assert all(divisible(x) for x in flat)
    margins = {
        "bottom_pml_to_target_z_m": 19.75 - 1.0,
        "side_x_pml_to_target_m": 3.0 - 1.0,
        "side_y_pml_to_target_m": 3.0 - 1.0,
        "top_pml_to_antenna_z_m": 49.0 - 45.0,
    }
    return {
        "schema": "a0_3d_v1_static_check/1",
        "checked_date": DATE,
        "derivation": "inputs derived from dep3d_gold_v1 5 cm pair (SHA-pinned); only title and target z range differ",
        "source_inputs_sha256": {"DEP3D_5CM_BG.in": DEP_BG_SHA, "DEP3D_5CM_TGT.in": DEP_TGT_SHA},
        "grid": {"dx_m": d, "dy_m": d, "dz_m": d, "nx": nx, "ny": ny, "nz": nz, "cells": cells,
                 "cell_cap_3d": 64_000_000, "dt_formula_s": dt, "time_steps": time_steps,
                 "dt_crosscheck": "equals dep3d 5 cm measured dt 9.629166007732353e-11 s and Iterations=12464"},
        "time_window_s": tw,
        "pml": {"cells": [20, 20, 20, 20, 20, 20], "formulation": "HORIPML",
                "physical_thickness_m": 1.0, "pml_cfs": "identical to dep3d 5 cm pair"},
        "coordinates": coords,
        "all_coordinates_on_grid": True,
        "margins_m": margins,
        "diagnostic_window_ns": [240.0, 400.0],
        "diagnostic_window_basis": "same +/-80 ns rule around the D10m self-computed two-way time 320.2 ns (air 30 m/c + cover 6 m/(c/4) + rock 14 m/(c/3)); reference only, not a physical timing calibration",
        "materials": {"rock": [9.0, 0.001], "cover": [16.0, 0.01], "wet_target": [20.0, 0.02]},
    }


def budget():
    walls = [DEP_MEAS[r]["wall_s"] for r in ("DEP3D_5CM_BG", "DEP3D_5CM_TGT")]
    peaks = [DEP_MEAS[r]["peak_job_commit_bytes"] for r in ("DEP3D_5CM_BG", "DEP3D_5CM_TGT")]
    return {
        "schema": "a0_3d_v1_budget/1",
        "date": DATE,
        "method": "dep3d_gold_v1 5 cm measured supervision values only; no extrapolation to other grids/domains/machines",
        "dep3d_measured": DEP_MEAS,
        "grid": {"dx_dy_dz_m": [0.05, 0.05, 0.05], "cells": 32_000_000,
                 "dt_s": 9.629166007732353e-11, "time_steps": 12464, "time_window_ns": 1200},
        "single_case": {"measured_wall_s_envelope": [min(walls), max(walls)],
                        "measured_job_peak_commit_bytes": [min(peaks), max(peaks)]},
        "batch_totals": {"n_cases": 2, "expected_wall_s_envelope": [2 * min(walls), 2 * max(walls)]},
        "hard_stops": {"wall_minutes_per_case": 40, "job_commit_GiB": 20, "output_GiB": 1,
                       "retries": 0, "max_fdtd_runs": 2},
        "preflight": {"minimum_available_RAM_GiB": 24, "minimum_available_VRAM_GiB": 8},
        "execution": {"threads": 8, "backend": "CUDA", "precision": "double", "device_id": 0,
                      "cpu_fallback": False},
        "known_limitations": [
            "dep3d wall clock includes unmeasured vctip residency; no basis to lower the 3D wall envelope",
            "budgets are magnitude references, not ETA promises; supervisor records govern",
        ],
    }


def cases_and_groups(input_sha):
    cases, groups = [], []
    specs = [
        (RUN_BG, "B2D-C3m-BG", "bg"),
        (RUN_TGT, "B2D-C3m-D10m-W4m-T0.5m-E20-S0.02", "target"),
    ]
    for run_id, mother, role in specs:
        seed = seed_for(mother)
        cases.append({
            "case_id": mother, "run_id": run_id, "file": f"{run_id}.in",
            "input_sha256": input_sha[run_id], "batch_id": BATCH, "grid_tier": "3D5CM",
            "family": "C3", "cover_thickness_m": 3.0, "role": role,
            "mother_model_id": mother, "group_id": "B2D-C3m", "seed": seed,
            "variant_tag": "3d5cm",
            "materials": {"rock": [9.0, 0.001], "cover": [16.0, 0.01], "wet_target": [20.0, 0.02]},
            "target_box_m": None if role == "bg" else [3.0, 19.75, 5.0, 20.25],
            "target_box_note_3d": None if role == "bg" else
                "3D box x 3-5 m, y 3-7 m (= 2D y 14-18 m relocated by y_3D = y_2D - 11), z 19.75-20.25 m; x-limited 2 m finite target",
            "diagnostic_window_ns": [240.0, 400.0],
        })
        groups.append({"case_id": mother, "group_id": "B2D-C3m", "run_id": run_id,
                       "mother_model_hash": input_sha[run_id], "seed": seed,
                       "variant_tag": "3d5cm"})
    cases_doc = {
        "batch_id": BATCH, "grid_tier": "3D5CM",
        "spec": "docs/research/2026-09-26_a0_3d_validation_contract_proposal.md (user ratified 2026-09-27) + docs/research/2026-09-26_3d_validation_subset_decision.md D3",
        "n_cases": 2, "n_exceptions": 0,
        "derivation": "dep3d_gold_v1 5 cm pair, target z 9.75-10.25 -> 19.75-20.25 (D20m class -> D10m anchor); BG otherwise byte-derived from DEP3D_5CM_BG.in",
        "cases": cases,
    }
    groups_doc = {
        "batch_id": BATCH, "grid_tier": "3D5CM",
        "convention": "spec 4.2: any split is by group_id only; the 2D A0 mother and its 3D validation pair share group B2D-C3m",
        "groups": groups,
    }
    return cases_doc, groups_doc


def write_json(path: Path, doc):
    path.write_text(json.dumps(doc, ensure_ascii=False, indent=1) + "\n", encoding="utf-8", newline="\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = build_inputs()
    input_sha = {}
    for run_id, text in inputs.items():
        b = text.encode("utf-8")
        (OUT / f"{run_id}.in").write_bytes(b)
        input_sha[run_id] = sha(b)

    sc = static_check()
    write_json(OUT / "static_check.json", sc)
    write_json(OUT / "budget.json", budget())
    cases_doc, groups_doc = cases_and_groups(input_sha)
    write_json(OUT / "cases.json", cases_doc)
    write_json(OUT / "groups.json", groups_doc)

    # stage the gate: contracts written, approved_to_simulate stays False
    gate_path = ROOT / "configs" / "research" / "gprmax_v4_execution_gate.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    assert gate["approved_to_simulate"] is False and gate["approved_run_ids"] == []
    launcher_sha = sha((ROOT / "scripts" / "run_approved_a0_3d.py").read_bytes())
    supervisor_sha = sha((ROOT / "scripts" / "bounded_windows_process.py").read_bytes())
    rt = "artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json"
    cu = "artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json"
    rt_sha = sha((ROOT / rt).read_bytes())
    cu_sha = sha((ROOT / cu).read_bytes())
    assert rt_sha == "14abf893d12a6db11e264f324ee83ec66194b9983efbcbc46cab5b79d74c9f38"
    assert cu_sha == "65c7c240db5146e1bf37f10a171f572f8ac599933d229a2272c44af0b780e9df"
    contracts = []
    for run_id in (RUN_BG, RUN_TGT):
        contracts.append({
            "packet_id": "A0-3D-V1",
            "run_id": run_id,
            "input_path": f"configs/research/a0_3d_v1/{run_id}.in",
            "input_sha256": input_sha[run_id],
            "runtime_identity_path": rt, "runtime_identity_sha256": rt_sha,
            "cuda_identity_path": cu, "cuda_identity_sha256": cu_sha,
            "launcher_sha256": launcher_sha,
            "supervisor_sha256": supervisor_sha,
            "attempt_record": f"artifacts/research_checks/{DATE}_{run_id}_attempt.json",
            "run_directory": f"artifacts/simulations/{DATE}_{run_id}",
            "continuation": "Serial in frozen order (BG then TGT); any failure aborts the batch, later attempts untouched. No retries.",
        })
    gate["updated_date"] = DATE
    gate["batch_id"] = f"{BATCH}_pending_agreement"
    gate["approved_to_simulate"] = False
    gate["approved_run_ids"] = []
    gate["approved_run_ids_pending_agreement"] = [RUN_BG, RUN_TGT]
    gate["approved_execution_contracts"] = contracts
    gate["approved_compute_budget"] = {
        "wall_minutes": 40, "job_commit_GiB": 20, "minimum_available_RAM_GiB": 24,
        "minimum_available_VRAM_GiB": 8, "output_GiB": 1, "threads": 8,
        "backend": "CUDA", "precision": "double", "retries": 0, "max_fdtd_runs": 2,
        "device_id": 0,
        "VRAM_limit_semantics": "Preflight availability threshold; not a hard device-memory quota",
    }
    gate["note"] = ("a0_3d_v1 frozen 2026-09-27 (2 runs: C3 anchor A0 3D validation pair, 5 cm isotropic, "
                    "derived from dep3d_gold_v1 5 cm inputs with only title and target z 9.75-10.25 -> 19.75-20.25 changed). "
                    "approved_to_simulate stays false until the user agrees to execute; G1 stays open until analysis; "
                    "no physical/training labels.")
    gate_path.write_text(json.dumps(gate, ensure_ascii=False, indent=1) + "\n", encoding="utf-8", newline="\n")

    print(json.dumps({"batch": BATCH, "runs": [RUN_BG, RUN_TGT],
                      "input_sha256": input_sha,
                      "static_check_sha256": sha((OUT / "static_check.json").read_bytes()),
                      "budget_sha256": sha((OUT / "budget.json").read_bytes()),
                      "cases_sha256": sha((OUT / "cases.json").read_bytes()),
                      "groups_sha256": sha((OUT / "groups.json").read_bytes()),
                      "launcher_sha256": launcher_sha,
                      "gate_staged": "approved_to_simulate=false"}, indent=1))


if __name__ == "__main__":
    main()
