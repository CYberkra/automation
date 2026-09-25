"""Attest that the three-grid FDTD-template direction diagnostics stay missing.

Reads the archived joint-roles results (42 rows x layer/target = 84 cosines),
confirms every cosine is null, and emits a new file refining the machine
reason to fdtd_numerically_unresolved per the error budget contract. The
archive itself is never modified. No solver, training, or field data.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "artifacts" / "research_checks" / "2026-09-25_review_fixes" / "joint_roles_results.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    entries = []
    for row in archived["rows"]:
        if not row["available"]:
            continue
        for role in ("layer", "target"):
            section = row[role]
            if section["signed_cosine"] is not None:
                raise AssertionError(f"unexpected restored cosine: {row['config_id']} {role}")
            entries.append({"config_id": row["config_id"], "grid": row["grid"],
                            "profile": row["profile"], "role": role,
                            "signed_cosine": None,
                            "previous_reason": section["metric_reasons"]["signed_cosine"],
                            "reason": "fdtd_numerically_unresolved"})
    if len(entries) != 84:
        raise AssertionError(f"expected 84 missing cosines, got {len(entries)}")
    if any(e["previous_reason"] != "direction_error_bound_not_declared" for e in entries):
        raise AssertionError("unexpected previous reason")
    result = {"schema": "three-grid-direction-attestation/1",
              "source_archive": str(ARCHIVE.relative_to(ROOT)),
              "source_archive_sha256": hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
              "missing_cosines": len(entries), "refined_reason": "fdtd_numerically_unresolved",
              "budget_contract": "configs/research/error_budget_contract_v0.1.json",
              "note": "FDTD templates remain numerically unresolved; no budget is declared "
                      "and no cosine is restored. This attestation modifies no archive.",
              "entries": entries}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"missing_cosines": len(entries), "output": str(args.output)}))


if __name__ == "__main__":
    main()
