"""Read-only numeric audit of the public Schmidt et al. Zenodo S2P archive.

This is a research source audit, not permittivity inversion or private GPR import.
Only the archive's declared Hz / S / RI / 50-ohm format is accepted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np


def audit(path):
    files = []
    with zipfile.ZipFile(path) as archive:
        for name in sorted(archive.namelist()):
            if not name.lower().endswith(".s2p"):
                continue
            raw = archive.read(name)
            lines = raw.decode("utf-8").splitlines()
            headers = [line.split("!")[0].split() for line in lines if line.strip().startswith("#")]
            if len(headers) != 1:
                raise ValueError(f"Expected one format header: {name}")
            header = headers[0]
            if len(header) != 6 or [x.lower() for x in header[:5]] != ["#", "hz", "s", "ri", "r"] or float(header[5]) != 50:
                raise ValueError(f"Unsupported source format: {name}: {header}")
            tokens = []
            for line in lines:
                line = line.split("!")[0].strip()
                if not line or line.startswith("#"):
                    continue
                tokens.extend(float(x) for x in line.split())
            if not tokens or len(tokens) % 9:
                raise ValueError(f"Invalid two-port row count: {name}")
            array = np.array(tokens).reshape(-1, 9)
            if not np.isfinite(array).all() or not (np.diff(array[:, 0]) > 0).all():
                raise ValueError(f"Nonfinite data or unordered frequencies: {name}")
            files.append({
                "file": name, "sha256": hashlib.sha256(raw).hexdigest(),
                "header": " ".join(header), "samples": len(array),
                "min_frequency_hz": float(array[0, 0]),
                "max_frequency_hz": float(array[-1, 0]),
                "count_in_provisional_20_170_mhz": int(((array[:, 0] >= 20e6) & (array[:, 0] <= 170e6)).sum()),
                "all_values_finite": True, "frequency_strictly_increasing": True,
            })
    if not files:
        raise ValueError("Archive contains no S2P data")
    return {
        "schema": "public-sparameter-source-audit/1", "dataset_doi": "10.5281/zenodo.15473270",
        "record_version": "1.0", "zip_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "archive_extracted": False, "private_field_data_used": False, "material_inversion_run": False,
        "record_description_start_hz": 300e6, "paper_acquisition_start_hz": 300e3,
        "paper_analysis_start_hz": 1e6, "files": files,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"files_checked": len(result["files"]), "output": str(args.output)}))
