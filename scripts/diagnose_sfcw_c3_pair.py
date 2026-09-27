"""Fixed-window paired diagnostics for the archived C3 SFCW development proxy."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs/research/sfcw_c3_pair_diagnostic_v1.json"
OUTPUT_DIR = ROOT / "artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic"
EXPECTED_MANIFEST_SHA = "00ecdf7c62fabd62b45af7c43aaed125161f9456ea95c20d077c933bb3c6125a"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_contract(path: Path) -> tuple[dict, str]:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    return json.loads(raw.decode("utf-8")), digest


def ratio(numerator: float, denominator: float, denominator_label: str) -> dict:
    if denominator == 0.0:
        return {"value": None, "reason": f"{denominator_label}_is_exactly_zero"}
    return {"value": float(numerator / denominator), "reason": None}


def apply_operator(values: np.ndarray, operator: str) -> np.ndarray:
    if operator == "identity":
        return values.copy()
    if operator == "trace_mean_removal":
        return values - np.mean(values, axis=1, keepdims=True)
    raise ValueError(f"unrecognized frozen operator: {operator}")


def synthetic_check() -> dict:
    identical_target = np.full((8, 4), 2.0 + 3.0j, dtype=np.complex128)
    removed = apply_operator(identical_target, "trace_mean_removal")
    assert np.max(np.abs(removed)) == 0.0
    zero = ratio(0.0, 0.0, "synthetic_denominator")
    assert zero == {"value": None, "reason": "synthetic_denominator_is_exactly_zero"}
    complex_input = np.array([[1.0 + 2.0j, 3.0 - 4.0j]], dtype=np.complex128)
    complex_output = apply_operator(complex_input, "trace_mean_removal")
    assert np.iscomplexobj(complex_output)
    assert np.allclose(complex_output, [[-1.0 + 3.0j, 1.0 - 3.0j]], rtol=0, atol=1e-15)
    e0 = np.array([[1.0 + 2j, 2.0 - 1j], [4.0 - 2j, 3.0 + 5j]])
    e1 = np.array([[2.0 + 4j, -1.0 + 3j], [1.0 - 2j, 8.0 + 4j]])
    delta = e1 - e0
    linear_error = apply_operator(e1, "trace_mean_removal") - apply_operator(e0, "trace_mean_removal") - apply_operator(delta, "trace_mean_removal")
    assert np.max(np.abs(linear_error)) < 1e-14
    return {
        "status": "passed",
        "identical_cross_trace_target_removed_exactly": True,
        "zero_denominator_result": zero,
        "complex_values_preserved": True,
        "linearity_max_abs_error": float(np.max(np.abs(linear_error))),
    }


def load_gather(z: np.lib.npyio.NpzFile, rows: list[dict], family: str, role: str, taper: str):
    selected = [r for r in rows if r["family"] == family and r["role"] == role and r["taper"] == taper]
    if family == "MT":
        assert len(selected) == 1
        row = selected[0]
        key = row["official_reconstruction"]["output_keys"]
        data = np.asarray(z[key["complex_envelope"]])
        if data.ndim == 1:
            data = data[:, None]
        positions = np.asarray(row["receiver_positions_xyz_m"], dtype=float)
        names = list(row["receiver_names"])
        paths = [row["input_path"]] * len(names)
        axis = np.asarray(z[key["envelope_time_s"]], dtype=float)
    else:
        assert len(selected) == 11
        # CO t01..t11 paths are not lexically ordered by geometry; sort on recorded y.
        selected.sort(key=lambda r: float(r["receiver_positions_xyz_m"][0][1]))
        data_columns, positions_list, names, paths = [], [], [], []
        axis = None
        for row in selected:
            key = row["official_reconstruction"]["output_keys"]
            values = np.asarray(z[key["complex_envelope"]])
            if values.ndim == 1:
                values = values[:, None]
            assert values.shape[1] == 1
            data_columns.append(values[:, 0])
            positions_list.append(row["receiver_positions_xyz_m"][0])
            names.append(row["receiver_names"][0])
            paths.append(row["input_path"])
            current_axis = np.asarray(z[key["envelope_time_s"]], dtype=float)
            if axis is None:
                axis = current_axis
            else:
                np.testing.assert_array_equal(axis, current_axis)
        data = np.column_stack(data_columns)
        positions = np.asarray(positions_list, dtype=float)
    assert data.shape == (2004, 33 if family == "MT" else 11)
    assert positions.shape == (data.shape[1], 3)
    assert len(names) == data.shape[1]
    assert np.all(np.diff(positions[:, 1]) > 0), f"{family} traces not ordered by increasing y"
    assert np.iscomplexobj(data) and np.isfinite(data).all()
    source_positions = np.asarray([r["source_position_xyz_m"] for r in selected], dtype=float)
    if family == "MT":
        source_positions = np.repeat(source_positions, data.shape[1], axis=0)
    assert source_positions.shape == positions.shape
    return data, axis, positions, source_positions, names, paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    check = synthetic_check()
    if args.self_test_only:
        print(json.dumps(check, indent=2))
        return
    if OUTPUT_DIR.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {OUTPUT_DIR}")

    contract, contract_sha = parse_contract(CONTRACT_PATH)
    assert contract["status"] == "fixed_before_diagnostic_output"
    assert contract["representation"] == "complex_envelope"
    assert contract["operators"] == ["identity", "trace_mean_removal"]
    assert contract["operator_application"].startswith("whole gather before windowing")
    windows = contract["windows"]
    assert len(windows) == 6
    for prev, cur in zip(windows, windows[1:]):
        assert round(prev["stop_s"] * 1e12) == round(cur["start_s"] * 1e12), "fixed boundaries differ by more than 1 ps"

    manifest_path = ROOT / contract["input_manifest"]
    manifest_sha = sha256(manifest_path)
    assert manifest_sha == EXPECTED_MANIFEST_SHA == contract["input_manifest_sha256"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    archive_path = manifest_path.parent / manifest["array_archive"]["path"]
    archive_sha = sha256(archive_path)
    assert archive_sha == manifest["array_archive"]["sha256"]
    records = manifest["records"]
    assert len(records) == 48 and manifest["inputs_n"] == 24
    input_hashes = {}
    for row in records:
        input_path = ROOT / row["input_path"]
        assert "C3m" in input_path.name and "C5m" not in str(input_path) and "C8m" not in str(input_path)
        actual = input_hashes.setdefault(str(input_path), sha256(input_path))
        assert actual == row["input_sha256"]
    assert len(input_hashes) == 24

    metric_rows, sensitivity_rows, per_trace = [], [], []
    with np.load(archive_path, allow_pickle=False) as z:
        for family in ("CO", "MT"):
            taper_data = {}
            for taper in ("no_taper", "tail_200ns"):
                bg, time_axis, bg_pos, bg_srcpos, bg_names, bg_paths = load_gather(z, records, family, "BG", taper)
                tgt, tgt_time, tgt_pos, tgt_srcpos, tgt_names, tgt_paths = load_gather(z, records, family, "TGT", taper)
                np.testing.assert_array_equal(time_axis, tgt_time)
                np.testing.assert_array_equal(bg_pos, tgt_pos)
                np.testing.assert_array_equal(bg_srcpos, tgt_srcpos)
                assert bg_names == tgt_names
                assert bg.shape == tgt.shape
                assert np.allclose(bg_pos[:, 2], 45.0, rtol=0, atol=1e-9)
                assert np.allclose(bg_srcpos[:, 2], 45.0, rtol=0, atol=1e-9)
                if family == "CO":
                    assert np.allclose(bg_pos[:, 1] - bg_srcpos[:, 1], 1.30, rtol=0, atol=1e-9)
                else:
                    assert np.allclose(bg_srcpos[:, 1], 15.35, rtol=0, atol=1e-9)
                    assert bg_names == [f"mt{k:02d}" for k in range(1, 34)]
                taper_data[taper] = (bg, tgt, time_axis, bg_pos, bg_srcpos, bg_names, bg_paths, tgt_paths)

            bg0, tgt0, time_axis, positions, source_positions, names, bg_paths0, tgt_paths0 = taper_data["no_taper"]
            bg1, tgt1, time_axis1, positions1, source_positions1, names1, bg_paths1, tgt_paths1 = taper_data["tail_200ns"]
            np.testing.assert_array_equal(time_axis, time_axis1)
            np.testing.assert_array_equal(positions, positions1)
            np.testing.assert_array_equal(source_positions, source_positions1)
            assert names == names1
            delta0, delta1 = tgt0 - bg0, tgt1 - bg1
            window_masks = {}
            sample_coverage = np.zeros(time_axis.shape, dtype=np.int8)
            for window in windows:
                mask = (time_axis >= window["start_s"]) & (time_axis < window["stop_s"])
                window_masks[window["id"]] = mask
                sample_coverage += mask.astype(np.int8)
            assert np.all(sample_coverage == 1), "fixed windows must assign every saved sample exactly once"
            for window in windows:
                mask = window_masks[window["id"]]
                tsel = time_axis[mask]
                for taper, (bg, tgt, _, _, _, _, bg_paths, tgt_paths) in taper_data.items():
                    delta = tgt - bg
                    bg_norm = float(np.linalg.norm(bg[mask]))
                    delta_norm = float(np.linalg.norm(delta[mask]))
                    raw_ratio = ratio(delta_norm, bg_norm, "background_norm")
                    for operator in contract["operators"]:
                        tb = apply_operator(bg, operator)
                        tt = apply_operator(tgt, operator)
                        td = apply_operator(delta, operator)
                        changed = tt - tb
                        input_norm = float(np.linalg.norm(delta[mask]))
                        changed_norm = float(np.linalg.norm(changed[mask]))
                        paired_error = float(np.linalg.norm((changed - delta)[mask]))
                        linear_error = float(np.linalg.norm((changed - td)[mask]))
                        metric_rows.append({
                            "geometry": family,
                            "taper": taper,
                            "operator": operator,
                            "window_id": window["id"],
                            "start_s": float(window["start_s"]),
                            "stop_s_exclusive": float(window["stop_s"]),
                            "sample_count": int(mask.sum()),
                            "first_sample_s": float(tsel[0]),
                            "last_sample_s": float(tsel[-1]),
                            "trace_count": int(bg.shape[1]),
                            "background_norm": bg_norm,
                            "delta_norm": input_norm,
                            "delta_over_background_norm": raw_ratio,
                            "background_residual_norm_ratio": ratio(float(np.linalg.norm(tb[mask])), bg_norm, "background_norm"),
                            "paired_change_error_ratio": ratio(paired_error, input_norm, "delta_norm"),
                            "delta_retained_norm_ratio": ratio(changed_norm, input_norm, "delta_norm"),
                            "operator_linearity_error_norm": linear_error,
                            "operator_linearity_error_over_input_delta_norm": ratio(linear_error, input_norm, "delta_norm"),
                            "background_input_paths": bg_paths,
                            "target_input_paths": tgt_paths,
                        })

                # Keep input background and delta tail sensitivity separate. Same-window
                # no-taper denominators are primary; full-period denominators are also explicit.
                bg_change = float(np.linalg.norm((bg1 - bg0)[mask]))
                delta_change = float(np.linalg.norm((delta1 - delta0)[mask]))
                sensitivity_rows.append({
                    "geometry": family,
                    "window_id": window["id"],
                    "sample_count": int(mask.sum()),
                    "background_tail_change_over_window_no_taper_norm": ratio(
                        bg_change, float(np.linalg.norm(bg0[mask])), "same_window_no_taper_background_norm"),
                    "delta_tail_change_over_window_no_taper_norm": ratio(
                        delta_change, float(np.linalg.norm(delta0[mask])), "same_window_no_taper_delta_norm"),
                    "background_tail_change_over_full_period_no_taper_norm": ratio(
                        bg_change, float(np.linalg.norm(bg0)), "full_period_no_taper_background_norm"),
                    "delta_tail_change_over_full_period_no_taper_norm": ratio(
                        delta_change, float(np.linalg.norm(delta0)), "full_period_no_taper_delta_norm"),
                })

                for trace_i, (name, pos) in enumerate(zip(names, positions)):
                    no = delta0[mask, trace_i]
                    tail = delta1[mask, trace_i]
                    per_trace.append({
                        "geometry": family,
                        "window_id": window["id"],
                        "trace_index_zero_based": trace_i,
                        "receiver_name": name,
                        "receiver_position_xyz_m": pos.tolist(),
                        "no_taper_delta_norm": float(np.linalg.norm(no)),
                        "tail_200ns_delta_norm": float(np.linalg.norm(tail)),
                        "tail_change_over_no_taper_trace_delta_norm": ratio(
                            float(np.linalg.norm(tail - no)), float(np.linalg.norm(no)), "no_taper_trace_delta_norm"),
                        "background_input_paths": [bg_paths0[trace_i], bg_paths1[trace_i]],
                        "target_input_paths": [tgt_paths0[trace_i], tgt_paths1[trace_i]],
                    })

    OUTPUT_DIR.mkdir(parents=True)
    script_path = Path(__file__).resolve()
    result = {
        "schema": "sfcw-c3-pair-diagnostic-results/1",
        "status": "completed_development_diagnostic_only",
        "scope": contract["scope"],
        "prohibited_interpretations": contract["prohibited_interpretations"],
        "units": contract["units"],
        "operator_definition": "identity: T(X)=X; trace_mean_removal: T(X)[t,j]=X[t,j]-mean_j(X[t,j]) using complex values across the full gather at each time",
        "operator_application": contract["operator_application"],
        "pair_definition": "E0=matched BG complex_envelope; E1=D10 TGT complex_envelope; delta=E1-E0 before operator application",
        "window_rule": contract["window_rule"],
        "windows": windows,
        "provenance": {
            "contract_path": str(CONTRACT_PATH.relative_to(ROOT)).replace("\\", "/"),
            "contract_sha256_raw_bytes": contract_sha,
            "contract_encoding": "UTF-8; hash covers unchanged raw bytes",
            "input_manifest_path": str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
            "input_manifest_sha256": manifest_sha,
            "input_npz_path": str(archive_path.relative_to(ROOT)).replace("\\", "/"),
            "input_npz_sha256": archive_sha,
            "source_h5_files_n": len(input_hashes),
            "source_h5_sha256": {Path(k).relative_to(ROOT).as_posix(): v for k, v in input_hashes.items()},
            "diagnostic_script_sha256": sha256(script_path),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
        },
        "synthetic_check": check,
        "metrics": metric_rows,
        "tail_sensitivity": sensitivity_rows,
        "per_trace_delta": per_trace,
    }
    out = OUTPUT_DIR / "results.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "metric_rows": len(metric_rows), "tail_sensitivity_rows": len(sensitivity_rows), "per_trace_rows": len(per_trace), "result_sha256": sha256(out)}, indent=2))


if __name__ == "__main__":
    main()
