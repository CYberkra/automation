"""Run existing array checks without field data, a solver, or source downloads."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import uuid

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CHECKS = (
    "check_processing_algebra.py",
    "check_operator_contract.py",
    "check_manual_weight_tradeoff.py",
    "check_evaluation_labels.py",
    "study_damage_pilot.py",
    "check_direction_diagnostics.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="New directory; never overwrite an earlier run")
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = args.output_dir or ROOT / "artifacts" / "local_checks" / f"{stamp}-{uuid.uuid4().hex[:8]}"
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for name in CHECKS:
        script = ROOT / "scripts" / name
        is_pilot = name == "study_damage_pilot.py"
        destination = output / "damage_pilot" / "results.json" if is_pilot else output / f"{script.stem}.json"
        output_args = ["--output-dir", str(destination.parent)] if is_pilot else ["--output", str(destination)]
        # Ignore PYTHONOPTIMIZE/PYTHONPATH: one historical check uses assertions.
        result = subprocess.run(
            [sys.executable, "-E", str(script), *output_args],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
            timeout=120,
        )
        if result.returncode:
            print(result.stdout, end="")
            print(result.stderr, end="", file=sys.stderr)
            raise SystemExit(result.returncode)
        evidence = json.loads(destination.read_text(encoding="utf-8"))
        if "checks" in evidence:
            count = len(evidence["checks"])
            if not count or not all(item["passed"] is True for item in evidence["checks"]):
                raise RuntimeError(f"Incomplete or failed checks: {name}")
        else:
            count = evidence["checks_passed"]
            if not isinstance(count, int) or count <= 0:
                raise RuntimeError(f"Invalid check count: {name}")
        rows.append({"script": f"scripts/{name}", "checks_passed": count,
                     "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
                     "result_file": destination.relative_to(output).as_posix(),
                     "result_sha256": hashlib.sha256(destination.read_bytes()).hexdigest()})
        print(f"PASS {name}: {count}")
    report = {"schema": "workspace-verification/1", "utc": stamp,
              "python": platform.python_version(), "platform": platform.platform(),
              "numpy": np.__version__, "total_checks_passed": sum(r["checks_passed"] for r in rows),
              "evidence_level": "constructed_array_checks_only",
              "field_data_used": False, "fdtd_executed": False, "training_executed": False,
              "runs": rows}
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"PASS total: {report['total_checks_passed']}; results: {output}")


if __name__ == "__main__":
    main()
