"""Decompose the archived 200 ns receiver taper for the existing C3 CO BG gather."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import (
    apply_tail_taper,
    direct_frequency_response,
    load_receiver,
    load_source,
    reconstruct_time_response,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs/research/sfcw_co_tail_decomposition_v1.json"
OUTPUT_DIR = ROOT / "artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ratio(num: float, den: float, name: str) -> dict:
    return {"value": None, "reason": f"{name}_is_exactly_zero"} if den == 0.0 else {"value": float(num / den), "reason": None}


def center(values: np.ndarray) -> np.ndarray:
    return values - np.mean(values, axis=1, keepdims=True)


def synthetic_check() -> dict:
    x = np.array([[1 + 2j, 3 - 1j, 5 + 4j], [2 - 3j, 4 + 1j, 1 + 5j]], dtype=np.complex128)
    w = np.array([1.0, 0.25, 0.0])
    w2 = np.array([1.0, 0.25])[:, None]
    removed = (1.0 - w2) * x
    after = w2 * x
    assert np.array_equal(x, after + removed)
    assert np.allclose(center(x), center(after) + center(removed), rtol=0, atol=1e-15)
    z = ratio(0.0, 0.0, "synthetic_denominator")
    assert z == {"value": None, "reason": "synthetic_denominator_is_exactly_zero"}
    assert np.iscomplexobj(center(x))
    return {
        "status": "passed",
        "complex_centering_linear_decomposition": True,
        "zero_denominator": z,
    }


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

    contract_raw = CONTRACT_PATH.read_bytes()
    contract_sha = hashlib.sha256(contract_raw).hexdigest()
    contract = json.loads(contract_raw.decode("utf-8"))
    assert contract["status"] == "fixed_before_output"
    parent_contract_path = ROOT / contract["parent_window_contract"]
    assert sha256(parent_contract_path) == contract["parent_window_contract_sha256"]
    manifest_path = ROOT / contract["manifest"]
    manifest_sha = sha256(manifest_path)
    assert manifest_sha == contract["manifest_sha256"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    archive_path = manifest_path.parent / manifest["array_archive"]["path"]
    archive_sha = sha256(archive_path)
    assert archive_sha == manifest["array_archive"]["sha256"]

    rows = [r for r in manifest["records"] if r["family"] == "CO" and r["role"] == "BG"]
    no_rows = [r for r in rows if r["taper"] == "no_taper"]
    taper_rows = [r for r in rows if r["taper"] == "tail_200ns"]
    assert len(no_rows) == len(taper_rows) == 11
    # Sort from the recorded physical receiver y coordinate, not by t01..t11 path strings.
    no_rows.sort(key=lambda r: float(r["receiver_positions_xyz_m"][0][1]))
    taper_rows.sort(key=lambda r: float(r["receiver_positions_xyz_m"][0][1]))

    sources, receivers, positions, source_positions, h5_hashes = [], [], [], [], {}
    with np.load(archive_path, allow_pickle=False) as arrays:
        x_columns, before_columns, after_columns = [], [], []
        for no, tw in zip(no_rows, taper_rows):
            assert no["input_path"] == tw["input_path"] and no["input_sha256"] == tw["input_sha256"]
            path = ROOT / no["input_path"]
            assert path.is_file() and "C3m-BG-CO11" in path.name and "C5m" not in str(path) and "C8m" not in str(path)
            h5_hashes[no["input_path"]] = sha256(path)
            assert h5_hashes[no["input_path"]] == no["input_sha256"]
            source = load_source(path)
            receiver = load_receiver(path, receiver_path="name:measurement", component="Ex")
            with h5py.File(path, "r") as h5:
                rx_attrs = dict(h5["rxs/rx1"].attrs)
                src_attrs = dict(h5["srcs/src1"].attrs)
            pos = np.asarray(rx_attrs["Position"], dtype=float)
            src_pos = np.asarray(src_attrs["Position"], dtype=float)
            np.testing.assert_array_equal(pos, np.asarray(no["receiver_positions_xyz_m"][0], dtype=float))
            np.testing.assert_array_equal(src_pos, np.asarray(no["source_position_xyz_m"], dtype=float))
            x_columns.append(receiver.samples)
            sources.append(source)
            receivers.append(receiver)
            positions.append(pos)
            source_positions.append(src_pos)
            k_no = no["official_reconstruction"]["output_keys"]
            k_tw = tw["official_reconstruction"]["output_keys"]
            before_columns.append(np.asarray(arrays[k_no["complex_envelope"]]))
            after_columns.append(np.asarray(arrays[k_tw["complex_envelope"]]))
            np.testing.assert_array_equal(arrays[k_no["frequency_hz"]], arrays[k_tw["frequency_hz"]])

        # Source history, lattice and source scaling must be identical for this common-offset gather.
        source0, receiver0 = sources[0], receivers[0]
        for source, receiver in zip(sources[1:], receivers[1:]):
            np.testing.assert_array_equal(source.samples, source0.samples)
            assert (source.dt, source.time_offset, source.spatial_scale, source.quantity, source.units, source.source_type) == (
                source0.dt, source0.time_offset, source0.spatial_scale, source0.quantity, source0.units, source0.source_type)
            assert (receiver.dt, receiver.time_offset, len(receiver.samples)) == (
                receiver0.dt, receiver0.time_offset, len(receiver0.samples))
        X = np.column_stack(x_columns)
        E_no = np.column_stack(before_columns)
        E_tail = np.column_stack(after_columns)
        assert X.shape == (20352, 11) and E_no.shape == E_tail.shape == (2004, 11)
        assert np.all(np.diff(np.asarray(positions)[:, 1]) > 0)
        assert np.allclose(np.asarray(positions)[:, 1], 12.65 + np.arange(11), rtol=0, atol=1e-12)
        assert np.allclose(np.asarray(positions)[:, 2], 45.0, rtol=0, atol=1e-12)
        assert np.allclose(np.asarray(source_positions)[:, 1], np.asarray(positions)[:, 1] - 1.30, rtol=0, atol=1e-12)

        dt = receiver0.dt
        t_raw = receiver0.times
        raw_centered = center(X)
        raw_residual_norm_by_time = np.linalg.norm(raw_centered, axis=1)
        taper_fraction = float(taper_rows[0]["tail_taper_fraction"])
        assert all(float(r["tail_taper_fraction"]) == taper_fraction for r in taper_rows)
        expected_fraction = (round(200e-9 / dt) - 0.25) / len(X)
        assert taper_fraction == expected_fraction
        W = apply_tail_taper(np.ones(len(X), dtype=np.float64), taper_fraction)
        assert W.shape == (len(X),) and np.all((W >= 0) & (W <= 1))
        U = X - W[:, None] * X
        X_after = W[:, None] * X
        raw_decomposition_closure = float(np.max(np.abs(X - X_after - U)))
        raw_scale = float(np.max(np.abs(X)))
        assert raw_decomposition_closure <= np.finfo(np.float64).eps * max(1.0, raw_scale) * 4

        # Apply the same official source-normalised frequency response and inverse transform to U.
        frequencies = np.asarray(manifest["frequency_hz"], dtype=np.float64)
        rx_u = replace(receiver0, samples=U)
        fu = direct_frequency_response(source0, rx_u, frequencies, source_floor_db=-100.0)
        assert bool(np.all(fu.source_valid)) and np.isfinite(fu.response).all()
        tu = reconstruct_time_response(fu, window="rectangular", zero_pad_factor=4, time_shift=0.0)
        np.testing.assert_array_equal(tu.time, arrays[k_no["envelope_time_s"]])
        F_U = tu.complex_envelope
        np.testing.assert_allclose(E_no - E_tail, F_U, rtol=1e-10, atol=1e-15)

        R_U = center(F_U)
        R_before, R_after = center(E_no), center(E_tail)
        R_difference = R_before - R_after
        full_raw_sq = float(np.sum(np.abs(raw_centered) ** 2))
        raw_assign = np.zeros(len(t_raw), dtype=np.int8)
        raw_windows = []
        for window in contract["raw_windows"]:
            mask = (t_raw >= window["start_s"]) & (t_raw < window["stop_s"])
            raw_assign += mask.astype(np.int8)
            values = raw_centered[mask]
            sq = float(np.sum(np.abs(values) ** 2))
            raw_windows.append({
                "window_id": window["id"], "start_s": window["start_s"], "stop_s_exclusive": window["stop_s"],
                "sample_count": int(mask.sum()),
                "first_sample_s": float(t_raw[mask][0]) if np.any(mask) else None,
                "last_sample_s": float(t_raw[mask][-1]) if np.any(mask) else None,
                "centered_l2": float(np.linalg.norm(values)),
                "centered_squared_norm": sq,
                "share_of_full_raw_centered_squared_norm": ratio(sq, full_raw_sq, "full_raw_centered_squared_norm"),
            })
        assert np.all(raw_assign <= 1)
        unassigned = raw_assign == 0
        unassigned_sq = float(np.sum(np.abs(raw_centered[unassigned]) ** 2))

        support = W < 1.0
        assert np.any(support)
        support_sq = float(np.sum(np.abs(raw_centered[support]) ** 2))
        support_u_sq = float(np.sum(np.abs(U[support]) ** 2))
        raw_global = {
            "sample_count": int(len(X)),
            "first_sample_s": float(t_raw[0]),
            "last_sample_s": float(t_raw[-1]),
            "dt_s": float(dt),
            "centered_l2": float(np.linalg.norm(raw_centered)),
            "centered_squared_norm": full_raw_sq,
            "per_time_centered_l2_array": "diagnostics.npz:raw_centered_residual_l2_by_time",
            "actual_taper_nonzero_support": {
                "definition": "official W<1; U=(1-W)X is not identically zero by construction, independent of observed amplitude",
                "sample_count": int(support.sum()),
                "first_index_zero_based": int(np.flatnonzero(support)[0]),
                "last_index_zero_based": int(np.flatnonzero(support)[-1]),
                "first_sample_s": float(t_raw[support][0]),
                "last_sample_s": float(t_raw[support][-1]),
                "centered_squared_norm_share_of_full_raw": ratio(support_sq, full_raw_sq, "full_raw_centered_squared_norm"),
                "removed_U_squared_norm": support_u_sq,
            },
            "fixed_windows_unassigned_samples": {
                "count": int(unassigned.sum()),
                "first_sample_s": float(t_raw[unassigned][0]) if np.any(unassigned) else None,
                "last_sample_s": float(t_raw[unassigned][-1]) if np.any(unassigned) else None,
                "centered_squared_norm_share_of_full_raw": ratio(unassigned_sq, full_raw_sq, "full_raw_centered_squared_norm"),
            },
        }

        per_trace = []
        for j, (row, pos, src_pos) in enumerate(zip(no_rows, positions, source_positions)):
            rj = raw_centered[:, j]
            per_trace.append({
                "input_path": row["input_path"], "input_sha256": row["input_sha256"],
                "receiver_name": row["receiver_names"][0], "receiver_y_m": float(pos[1]),
                "source_y_m": float(src_pos[1]),
                "distance_to_lateral_pml_inner_y1_m": float(pos[1] - 1.0),
                "distance_to_lateral_pml_inner_y31_m": float(31.0 - pos[1]),
                "raw_mean_v_per_m": float(np.mean(X[:, j])),
                "raw_last_sample_v_per_m": float(X[-1, j]),
                "last_sample_time_s": float(t_raw[-1]),
                "centered_l2_full_record": float(np.linalg.norm(rj)),
                "centered_l2_actual_taper_support": float(np.linalg.norm(rj[support])),
            })

        reconstruction_windows = []
        recon_assign = np.zeros(len(tu.time), dtype=np.int8)
        for window in contract["reconstruction_windows"]:
            mask = (tu.time >= window["start_s"]) & (tu.time < window["stop_s"])
            recon_assign += mask.astype(np.int8)
            before, after, removed = R_before[mask], R_after[mask], R_U[mask]
            nb, na, nu = (float(np.linalg.norm(v)) for v in (before, after, removed))
            closure = before - after - removed
            cross_term = float(2.0 * np.real(np.vdot(after, removed)))
            reconstruction_windows.append({
                "window_id": window["id"], "start_s": window["start_s"], "stop_s_exclusive": window["stop_s"],
                "sample_count": int(mask.sum()),
                "first_sample_s": float(tu.time[mask][0]) if np.any(mask) else None,
                "last_sample_s": float(tu.time[mask][-1]) if np.any(mask) else None,
                "centered_before_norm": nb,
                "centered_after_taper_norm": na,
                "centered_removed_tail_norm": nu,
                "relative_linear_closure": ratio(float(np.linalg.norm(closure)), nb, "centered_before_norm"),
                "after_over_before_norm": ratio(na, nb, "centered_before_norm"),
                "squared_norm_cross_term": cross_term,
                "squared_norm_identity_residual": float(np.sum(np.abs(before) ** 2) - np.sum(np.abs(after) ** 2) - np.sum(np.abs(removed) ** 2) - cross_term),
            })
        assert np.all(recon_assign == 1)

        freq_keys = no_rows[0]["official_reconstruction"]["output_keys"]
        saved_no_resp = np.column_stack([
            arrays[r["official_reconstruction"]["output_keys"]["response"]] for r in no_rows
        ])
        saved_tail_resp = np.column_stack([
            arrays[r["official_reconstruction"]["output_keys"]["response"]] for r in taper_rows
        ])
        response_closure = fu.response - (saved_no_resp - saved_tail_resp)
        envelope_closure = F_U - (E_no - E_tail)
        processing_closure = {
            "frequency_response_removed_vs_archived_difference_max_abs": float(np.max(np.abs(response_closure))),
            "frequency_response_removed_relative_l2": ratio(float(np.linalg.norm(response_closure)), float(np.linalg.norm(saved_no_resp - saved_tail_resp)), "archived_frequency_difference_norm"),
            "complex_envelope_removed_vs_archived_difference_max_abs": float(np.max(np.abs(envelope_closure))),
            "complex_envelope_removed_relative_l2": ratio(float(np.linalg.norm(envelope_closure)), float(np.linalg.norm(E_no - E_tail)), "archived_envelope_difference_norm"),
            "centered_decomposition_max_abs": float(np.max(np.abs(R_before - R_after - R_U))),
        }
        arrays_out = {
            "raw_time_s": t_raw,
            "raw_centered_residual_l2_by_time": raw_residual_norm_by_time,
            "raw_centered_residual": raw_centered,
            "official_tail_weight_W": W,
            "removed_tail_U": U,
            "frequency_hz": frequencies,
            "removed_tail_frequency_response": fu.response,
            "removed_tail_complex_envelope": F_U,
            "archived_no_taper_complex_envelope": E_no,
            "archived_tail_200ns_complex_envelope": E_tail,
            "reconstruction_time_s": tu.time,
        }

    OUTPUT_DIR.mkdir(parents=True)
    npz_path = OUTPUT_DIR / "diagnostics.npz"
    np.savez_compressed(npz_path, **arrays_out)
    processing_path = Path(__import__("gprMax.toolboxes.SFCW.processing", fromlist=["__file__"]).__file__).resolve()
    script_path = Path(__file__).resolve()
    result = {
        "schema": "sfcw-co-tail-decomposition-results/1",
        "status": "completed_development_diagnostic_only",
        "scope": contract["scope"],
        "contract_path": str(CONTRACT_PATH.relative_to(ROOT)).replace("\\", "/"),
        "contract_sha256": contract_sha,
        "parent_window_contract_sha256": sha256(parent_contract_path),
        "input_manifest_path": str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
        "input_manifest_sha256": manifest_sha,
        "input_npz_path": str(archive_path.relative_to(ROOT)).replace("\\", "/"),
        "input_npz_sha256": archive_sha,
        "input_h5_sha256": h5_hashes,
        "official_processing_source_path": str(processing_path),
        "official_processing_source_sha256": sha256(processing_path),
        "script_path": str(script_path.relative_to(ROOT)).replace("\\", "/"),
        "script_sha256": sha256(script_path),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scope_note": "C3 CO BG only; raw E-field variation and deterministic taper contribution; not PML/finite-domain cause attribution, clean reference or physical error budget",
        "source_metadata": {
            "source_type": source0.source_type, "quantity": source0.quantity, "units": source0.units,
            "spatial_scale": float(source0.spatial_scale), "dt_s": float(source0.dt),
            "time_offset_s": float(source0.time_offset), "sample_count": int(len(source0.samples)),
            "history_sha256_by_input": {r["input_path"]: hashlib.sha256(source.samples.tobytes()).hexdigest() for r, source in zip(no_rows, sources)},
        },
        "receiver_metadata": {
            "component": "Ex", "units": receiver0.units, "dt_s": float(dt),
            "time_offset_s": float(receiver0.time_offset), "sample_count": int(len(X)),
            "trace_count": int(X.shape[1]), "trace_sort": "ascending numeric receiver y from manifest position metadata",
            "source_lattice_and_history_identical_across_traces": True,
        },
        "taper": {
            "source_unchanged": True,
            "official_function": "gprMax.toolboxes.SFCW.processing.apply_tail_taper",
            "duration_nominal_ns": 200.0,
            "fraction": taper_fraction,
            "weights_array": "diagnostics.npz:official_tail_weight_W",
            "removed_signal_definition": "U = X - W*X; official receiver-only raised-cosine W, no amplitude-selected support",
            "raw_linear_decomposition_max_abs_roundoff": raw_decomposition_closure,
        },
        "fixed_raw_windows": raw_windows,
        "raw_global": raw_global,
        "per_trace": per_trace,
        "reconstruction": {
            "frequency_start_hz": float(frequencies[0]), "frequency_stop_hz": float(frequencies[-1]),
            "frequency_count": int(len(frequencies)), "frequency_step_hz": float(frequencies[1] - frequencies[0]),
            "method": "official direct_frequency_response + reconstruct_time_response",
            "source_floor_db": -100.0, "time_window": "full saved receiver history",
            "reconstruct_window": "rectangular", "zero_pad_factor": 4, "time_shift_s": 0.0,
            "norm_units": "unweighted discrete complex-envelope L2; not calibrated energy",
            "decomposition": "C(E_no)=C(E_tail)+C(F(U)); squared norms include cross term and are not independent shares",
            "processing_closure": processing_closure,
            "windows": reconstruction_windows,
        },
        "synthetic_check": check,
        "diagnostics_npz": {"path": npz_path.name, "sha256": sha256(npz_path), "keys": list(arrays_out)},
        "limits": [
            "tail decomposition proves only the deterministic linear processing contribution under identical source/history normalization",
            "no physical attribution to PML, finite domain, or truncation cause is made",
            "air-only prior tail/PML findings are not transferred as a correction or acceptance threshold to layered CO data",
            "raw E-field and source-normalized transfer response are diagnostics, not clean field truth",
        ],
    }
    result_path = OUTPUT_DIR / "results.json"
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(OUTPUT_DIR), "npz_sha256": result["diagnostics_npz"]["sha256"], "results_sha256": sha256(result_path), "raw_samples": len(X), "traces": X.shape[1], "raw_windows": len(raw_windows), "reconstruction_windows": len(reconstruction_windows)}, indent=2))


if __name__ == "__main__":
    main()
