"""Read-only audit of supplied materials; derived metadata goes under artifacts/.

Run: python scripts/audit_materials.py
Requires numpy, pandas, h5py, pypdf. No simulator or ML training is invoked.
"""
from __future__ import annotations

import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "探地雷达背景资料"
OUT = ROOT / "artifacts" / "initial_audit"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def audit_csv(path):
    with path.open(encoding="utf-8-sig") as stream:
        header = [next(stream).strip() for _ in range(4)]
    values = [float(re.search(r"=\s*([\d.]+)", line)[1]) for line in header]
    ns, window, nt, spacing = values
    ns, nt = int(ns), int(nt)
    count = nonfinite = metadata_changes = 0
    lo, hi = np.full(5, np.inf), np.full(5, -np.inf)
    traces = []
    for chunk in pd.read_csv(path, skiprows=4, header=None, chunksize=ns * 256, dtype=float):
        a = chunk.to_numpy()
        if a.shape[1] != 5:
            raise ValueError(f"Expected five columns: {path}, got {a.shape}")
        count += len(a)
        nonfinite += int((~np.isfinite(a)).sum())
        lo, hi = np.minimum(lo, np.nanmin(a, axis=0)), np.maximum(hi, np.nanmax(a, axis=0))
        complete = len(a) // ns
        if complete:
            cube = a[:complete * ns].reshape(complete, ns, 5)
            coords = cube[:, :, [0, 1, 2, 4]]
            metadata_changes += int(np.any(coords != coords[:, :1, :], axis=(1, 2)).sum())
            traces.append(cube[:, 0, [0, 1, 2, 4]])
    xyz = np.concatenate(traces)
    dlat = np.deg2rad(np.diff(xyz[:, 1]))
    dlon = np.deg2rad(np.diff(xyz[:, 0]))
    lat = np.deg2rad(xyz[:, 1])
    hav = np.sin(dlat / 2) ** 2 + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(dlon / 2) ** 2
    distance = 6371008.8 * 2 * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))
    return {
        "header": header, "samples_per_trace": ns, "traces_declared": nt,
        "time_window_ns": window, "trace_interval_header_m": spacing,
        "payload_rows": count, "expected_rows": ns * nt, "row_count_matches": count == ns * nt,
        "nonfinite_values": nonfinite, "traces_with_varying_position_metadata": metadata_changes,
        "column_min": lo.tolist(), "column_max": hi.tolist(),
        "position_columns": ["longitude", "latitude", "surface_elevation_m_documented", "flight_height_m_documented"],
        "coordinate_distance_note": "Spherical horizontal distance diagnostic; CRS, datum, alignment not validated",
        "coordinate_path_length_m": float(distance.sum()),
        "coordinate_step_median_m": float(np.median(distance)),
        "header_path_length_m": (nt - 1) * spacing,
        "domain": "exported_real_ascan; upstream transform/calibration unverified",
    }


def pdf_text(blob):
    reader = PdfReader(io.BytesIO(blob))
    return "\n\n".join(f"=== PAGE {i + 1} ===\n{page.extract_text()}" for i, page in enumerate(reader.pages)), len(reader.pages)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    texts = OUT / "document_text"
    texts.mkdir(exist_ok=True)
    records = []
    for path in sorted(SOURCE.rglob("*")):
        if not path.is_file():
            continue
        record = {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)}
        print("AUDIT", record["path"], flush=True)
        suffix = path.suffix.lower()
        if suffix == ".csv":
            record["csv"] = audit_csv(path)
        elif suffix == ".h5":
            with h5py.File(path, "r") as f:
                record["attributes"] = {k: np.asarray(v).tolist() for k, v in f.attrs.items()}
                data = f["data"]
                counts = {}
                for plane in range(data.shape[0]):
                    for start in range(0, data.shape[1], 256):
                        ids, nums = np.unique(data[plane, start:start + 256, :], return_counts=True)
                        for identifier, number in zip(ids, nums):
                            counts[str(int(identifier))] = counts.get(str(int(identifier)), 0) + int(number)
                record["geometry"] = {"shape": list(data.shape), "dtype": str(data.dtype), "material_id_counts": counts,
                    "classification": "geometry material grid; not receiver waveform output"}
        elif suffix == ".pdf":
            body, pages = pdf_text(path.read_bytes())
            target = texts / (path.stem + ".txt")
            target.write_text(body, encoding="utf-8", errors="backslashreplace")
            record["pdf"] = {"pages": pages, "extracted_text": target.relative_to(ROOT).as_posix()}
        elif suffix == ".zip":
            members = []
            with zipfile.ZipFile(path) as z:
                for item in z.infolist():
                    if item.is_dir():
                        continue
                    name = item.filename
                    if not item.flag_bits & 0x800:
                        try:
                            name = name.encode("cp437").decode("gbk")
                        except UnicodeError:
                            pass
                    content = z.read(item)
                    member = {"name": name, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
                    if name.lower().endswith(".pdf"):
                        body, pages = pdf_text(content)
                        # Only write derived text under a known directory, never extract ZIP paths.
                        target = texts / (Path(name).stem + ".txt")
                        target.write_text(body, encoding="utf-8", errors="backslashreplace")
                        member.update(pages=pages, extracted_text=target.relative_to(ROOT).as_posix())
                    members.append(member)
            record["zip_members"] = members
        records.append(record)
    payload = {"audit_utc": datetime.now(timezone.utc).isoformat(), "source": SOURCE.relative_to(ROOT).as_posix(),
               "script_sha256": sha256(Path(__file__)), "files": records}
    (OUT / "materials_inventory.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    rows = []
    for record in records:
        if "csv" in record:
            c = record["csv"]
            rows.append({"file": record["path"], "samples": c["samples_per_trace"], "traces": c["traces_declared"],
                         "window_ns": c["time_window_ns"], "rows_ok": c["row_count_matches"],
                         "nonfinite": c["nonfinite_values"], "metadata_changes": c["traces_with_varying_position_metadata"],
                         "height_min_m": c["column_min"][4], "height_max_m": c["column_max"][4],
                         "sha256": record["sha256"]})
    pd.DataFrame(rows).to_csv(OUT / "csv_summary.csv", index=False, encoding="utf-8-sig")
    print("DONE", len(records), "files", len(rows), "CSV", flush=True)


if __name__ == "__main__":
    main()
