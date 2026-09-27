"""Read-only frozen-window F0/delta diagnostics for the existing C3 CO pair."""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "artifacts/research_checks/2026-09-27_c3_co_pair_window_contamination"
EVENT_TABLE = ROOT / "configs/research/batch2d_v1_event_table_v0.1.json"
EVENT_TABLE_SHA256 = "b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c"
RECORDS = ROOT / "artifacts/research_checks/2026-09-26_eval_batch2d_co_r1/records.json"
MANIFEST = ROOT / "artifacts/research_checks/2026-09-26_eval_batch2d_co_r1/run_manifest.json"
N_SAMPLES = 20352
N_TRACES = 11
RX_Y = np.array([12.65 + k for k in range(N_TRACES)], dtype=float)
OFFSET_M = 1.30
EVENT_IDS = (
    "EV-C3m-COV",
    "EV-C3m-D10m-W4m-T0.5m-E20-S0.02",
    "EV-C3m-BG-NEG-D5m",
    "EV-C3m-BG-NEG-D10m",
    "EV-C3m-BG-NEG-D20m",
)
IDENTITY = "B0_G1_BG"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_rows():
    rows = json.loads(RECORDS.read_text(encoding="utf-8"))["records"]
    return [r for r in rows if r.get("candidate_id") == IDENTITY
            and r.get("family", "").upper() == "C3" and r.get("event_id") == "EV-C3m-COV"
            and r.get("case") in {"BG", "TGT"}]


def extract_sources(rows):
    sources = {}
    for row in rows:
        case = row["case"]
        for run_id, pair in row["provenance"]["case_h5_files"].items():
            record = {"path": pair[0].replace("\\", "/"), "sha256": pair[1]}
            if run_id in sources:
                assert sources[run_id] == record, f"conflicting source provenance: {run_id}"
            sources[run_id] = record
            expected_mother = "B2D-C3m-BG-CO11-" if case == "BG" else "B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-CO11-"
            assert run_id.startswith(expected_mother), f"wrong case provenance: {run_id}"
    assert len(sources) == 2 * N_TRACES, f"expected 22 paired trace files, found {len(sources)}"
    return sources


def load_gather(case, sources, sample_interval):
    prefix = "B2D-C3m-BG-CO11-" if case == "BG" else "B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-CO11-"
    columns, metadata = [], []
    for k, expected_y in enumerate(RX_Y, start=1):
        run_id = f"{prefix}t{k:02d}"
        source = sources[run_id]
        path = ROOT / source["path"]
        actual_sha = sha256(path)
        assert actual_sha == source["sha256"], f"H5 SHA mismatch: {run_id}"
        with h5py.File(path, "r") as h5:
            rx = h5["rxs/rx1"]
            ex = rx["Ex"]
            x = ex[:].astype(np.float64)
            rx_pos = np.asarray(rx.attrs["Position"], dtype=float)
            tx_pos = np.asarray(h5["srcs/src1"].attrs["Position"], dtype=float)
            dt = float(ex.attrs["SampleInterval"])
            offset = float(ex.attrs["TimeSampleOffset"])
        assert x.shape == (N_SAMPLES,), f"unexpected sample count: {run_id} {x.shape}"
        assert np.isfinite(x).all(), f"non-finite data: {run_id}"
        assert abs(rx_pos[1] - expected_y) <= 1e-9, f"Rx y mismatch: {run_id}"
        assert abs(tx_pos[1] - (expected_y - OFFSET_M)) <= 1e-9, f"Tx/Rx offset mismatch: {run_id}"
        assert dt == sample_interval and offset == 0.0, f"time axis mismatch: {run_id}"
        columns.append(x)
        metadata.append({"run_id": run_id, "file": source["path"], "sha256": actual_sha,
                         "rx_position_m": rx_pos.tolist(), "tx_position_m": tx_pos.tolist(),
                         "sample_interval_s": dt, "time_sample_offset_s": offset,
                         "n_samples": int(x.size)})
    return np.stack(columns, axis=1), metadata


def window_metrics(f0, delta, lo, hi):
    f0w = f0[lo:hi + 1, :]
    dw = delta[lo:hi + 1, :]

    def vector_stats(a):
        ss = float(np.sum(a * a, dtype=np.float64))
        n = int(a.size)
        return {"l2_norm": float(np.sqrt(ss)), "sum_squares": ss,
                "mean_square": ss / n, "max_abs": float(np.max(np.abs(a))),
                "nonzero_sample_count": int(np.count_nonzero(a)),
                "exact_zero": bool(np.all(a == 0.0)), "n_values": n}

    traces = []
    for k in range(N_TRACES):
        fs = vector_stats(f0w[:, k])
        ds = vector_stats(dw[:, k])
        traces.append({"trace": f"t{k + 1:02d}", "f0": fs, "delta": ds,
                       "amplitude_norm_ratio_f0_over_delta": None if ds["l2_norm"] == 0.0 else fs["l2_norm"] / ds["l2_norm"],
                       "ratio_reason": "zero_delta_no_epsilon_added" if ds["l2_norm"] == 0.0 else None,
                       "delta_nonzero": bool(np.any(dw[:, k] != 0.0))})
    whole_f0 = vector_stats(f0w)
    whole_delta = vector_stats(dw)
    ratio = None if whole_delta["l2_norm"] == 0.0 else whole_f0["l2_norm"] / whole_delta["l2_norm"]
    return {"window": {"index_lo": lo, "index_hi_inclusive": hi,
                        "n_samples": int(hi - lo + 1), "n_traces": N_TRACES,
                        "mask": "all_samples_by_all_11_traces; no shift/alignment"},
            "whole_window_merged_across_samples_and_traces": {
                "f0": whole_f0, "delta": whole_delta,
                "amplitude_norm_ratio_f0_over_delta": ratio,
                "ratio_reason": "zero_delta_no_epsilon_added" if ratio is None else None,
                "delta_nonzero": bool(np.any(dw != 0.0)),
                "delta_nonzero_sample_count": int(np.count_nonzero(dw)),
                "delta_max_abs": float(np.max(np.abs(dw))),
                "f0_max_abs": float(np.max(np.abs(f0w)))},
            "per_trace": traces}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    if out.exists():
        raise SystemExit(f"refusing_to_overwrite:{out}")

    assert sha256(EVENT_TABLE) == EVENT_TABLE_SHA256, "event table SHA mismatch"
    table = json.loads(EVENT_TABLE.read_text(encoding="utf-8"))
    assert table["status"] == "frozen"
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    co_rows = source_rows()
    assert {r["case"] for r in co_rows} == {"BG", "TGT"}
    sources = extract_sources(co_rows)
    for case in ("BG", "TGT"):
        for run_id, item in manifest["co_cases"][case]["h5"].items():
            assert sources[run_id] == {"path": item["path"].replace("\\", "/"), "sha256": item["sha256"]}, f"manifest/provenance mismatch: {run_id}"
    entries = {e["event_id"]: e for e in table["entries"] if e["event_id"] in EVENT_IDS}
    assert set(entries) == set(EVENT_IDS)
    dt_s = float(table["time_axis"]["dt_s"])

    f0, f0_meta = load_gather("BG", sources, dt_s)
    f1, f1_meta = load_gather("TGT", sources, dt_s)
    assert f0.shape == f1.shape == (N_SAMPLES, N_TRACES)
    assert all(np.array_equal(a["rx_position_m"], b["rx_position_m"])
               and np.array_equal(a["tx_position_m"], b["tx_position_m"])
               for a, b in zip(f0_meta, f1_meta)), "paired Tx/Rx geometry differs"
    delta = f1 - f0
    first_nonzero = np.flatnonzero(np.any(delta != 0.0, axis=1))

    role_note = {
        "EV-C3m-COV": "cover_interface role; no cover-toggle pair; diagnostic is target-pair change inside frozen cover window",
        "EV-C3m-D10m-W4m-T0.5m-E20-S0.02": "target role; target-pair difference window",
        "EV-C3m-BG-NEG-D5m": "negative_control role retained; target-pair difference sampled in this BG-NEG window, not relabeled target",
        "EV-C3m-BG-NEG-D10m": "negative_control role retained; exact same frozen indices as D10 target window, not independent evidence",
        "EV-C3m-BG-NEG-D20m": "negative_control role retained; target-pair difference sampled in this BG-NEG window, not relabeled target",
    }
    window_rows = []
    for event_id in EVENT_IDS:
        e = entries[event_id]
        result = window_metrics(f0, delta, int(e["index_lo"]), int(e["index_hi"]))
        window_rows.append({"event_id": event_id, "role": e["role"],
                            "reference_kind_metadata": e["reference_type"],
                            "reference_state_metadata": e["reference_state"],
                            "window_hash": e["window_hash"], "mask_hash": e["mask_hash"],
                            "role_interpretation": role_note[event_id], **result})

    paired_inputs = {
        "F0": {"case": "BG", "source_files": f0_meta},
        "F1": {"case": "TGT", "source_files": f1_meta},
    }
    result = {
        "schema_version": "c3_co_pair_window_contamination_v1",
        "scope": "existing C3 CO BG/TGT paired 11-trace gathers only; CPU read-only; no solver/evaluator invocation",
        "definition": "delta = F1 - F0; amplitude contamination ratio = ||M F0||_2 / ||M delta||_2; signal sum_squares/mean_square are discrete-field-sample proxies, not joules",
        "source_files": {"event_table": str(EVENT_TABLE.relative_to(ROOT)).replace("\\", "/"),
                         "event_table_sha256": sha256(EVENT_TABLE),
                         "records": str(RECORDS.relative_to(ROOT)).replace("\\", "/"),
                         "records_sha256": sha256(RECORDS),
                         "run_manifest": str(MANIFEST.relative_to(ROOT)).replace("\\", "/"),
                         "run_manifest_sha256": sha256(MANIFEST),
                         "evaluation_script": "scripts/run_eval_batch2d_co.py",
                         "evaluation_script_sha256": sha256(ROOT / "scripts/run_eval_batch2d_co.py"),
                         "analysis_script": "scripts/analyze_batch2d_co.py",
                         "analysis_script_sha256": sha256(ROOT / "scripts/analyze_batch2d_co.py"),
                         "diagnostic_script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
                         "diagnostic_script_sha256": sha256(__file__)},
        "axis_and_pair_checks": {"dt_s": dt_s, "time_sample_offset_s": 0.0,
                                  "sample_count_per_trace": N_SAMPLES, "trace_count": N_TRACES,
                                  "rx_y_m": RX_Y.tolist(), "tx_rx_offset_m": OFFSET_M,
                                  "paired_positions_equal": True,
                                  "source_provenance_sha_matched": True,
                                  "global_delta_nonzero": bool(np.any(delta != 0.0)),
                                  "first_nonzero_sample_index_zero_based": None if first_nonzero.size == 0 else int(first_nonzero[0]),
                                  "leading_all_trace_zero_delta_samples": N_SAMPLES if first_nonzero.size == 0 else int(first_nonzero[0]),
                                  "anchor_t05_check_reused_from_existing_evaluation_manifest": manifest["co_cases"]["BG"]["mother_check"]["anchor_max_abs_difference"] == 0.0 and manifest["co_cases"]["TGT"]["mother_check"]["anchor_max_abs_difference"] == 0.0},
        "paired_inputs": paired_inputs,
        "windows": window_rows,
        "interpretation_limits": ["No threshold or pass/fail judgment is applied.",
                                   "No clean/isolated event reference eligibility is granted.",
                                   "BG-NEG rows retain their negative_control role; the pair delta is diagnostic only.",
                                   "COV metadata is not a cover-toggle difference; this computes the target-pair change within that window.",
                                   "D10 target and BG-NEG-D10 use identical frozen indices and are role-distinct, not independent evidence."]
    }
    out.mkdir(parents=True)
    (out / "results.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out / "results.json"), "windows": len(window_rows),
                      "global_delta_nonzero": result["axis_and_pair_checks"]["global_delta_nonzero"],
                      "first_nonzero_sample": result["axis_and_pair_checks"]["first_nonzero_sample_index_zero_based"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
