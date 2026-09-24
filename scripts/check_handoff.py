"""Read-only handoff integrity checks; no solver, field data, or network."""
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    tracked = [p for p in subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0") if p]
    pattern = re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})")
    missing_tracked_files = [p for p in tracked if not (ROOT/p).exists()]
    secret_files = [p for p in tracked if (ROOT/p).is_file() and pattern.search((ROOT/p).read_bytes())]
    gate = json.loads((ROOT/"configs/research/gprmax_v4_execution_gate.json").read_text(encoding="utf-8"))
    registry = json.loads((ROOT/"docs/research/research_artifact_registry.json").read_text(encoding="utf-8"))
    missing_artifacts = []
    for run in registry["runs"]:
        if run.get("tracked", True):
            p = run.get("directory", run.get("record"))
            if not (ROOT/p).exists():
                missing_artifacts.append(p)
    entries = ["START_HERE.md", "docs/research/2026-09-24_gprmax_v4_review.md",
               "docs/research/2026-09-24_damage_pilot_findings.md", "docs/research/decision_log.md"]
    links, missing_links = 0, []
    for name in entries:
        path = ROOT/name
        for link in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if link.startswith(("https://", "http://", "#")):
                continue
            links += 1
            if not (path.parent/link.split("#")[0]).exists():
                missing_links.append({"file": name, "link": link})
    result_dir = ROOT/"artifacts/research_checks/2026-09-24_damage_pilot_r2"
    evidence = json.loads((result_dir/"results.json").read_text(encoding="utf-8"))
    bad_hashes = [name for name, expected in evidence["artifact_sha256"].items()
                  if hashlib.sha256((result_dir/name).read_bytes()).hexdigest() != expected]
    for name, expected in evidence["source_sha256"].items():
        source = ROOT/("docs/research" if name.endswith(".md") else "scripts")/name
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            bad_hashes.append(str(source.relative_to(ROOT)))
    # A recorded user approval is valid; this checker must not silently demand
    # that every future run stay unapproved. Validate the exact bound contract.
    safe_gate = gate["approved_to_simulate"] is False and gate["approved_run_ids"] == []
    gate_errors = []
    if gate["approved_to_simulate"]:
        try:
            contract = gate["approved_execution_contract"]
            budget = gate["approved_compute_budget"]
            safe_gate = (gate["approved_run_ids"] == ["M00_x_3d"]
                         and budget["max_fdtd_runs"] == 1 and budget["retries"] == 0
                         and all(budget[k] > 0 for k in ("wall_minutes", "job_commit_GiB", "output_GiB")))
            for path_key, hash_key in (("input_path", "input_sha256"),
                                       ("runtime_identity_path", "runtime_identity_sha256")):
                if hashlib.sha256((ROOT/contract[path_key]).read_bytes()).hexdigest() != contract[hash_key]:
                    gate_errors.append(path_key)
            for filename, key in (("run_approved_airborne.py", "launcher_sha256"),
                                  ("bounded_windows_process.py", "supervisor_sha256")):
                if hashlib.sha256((ROOT/"scripts"/filename).read_bytes()).hexdigest() != contract[key]:
                    gate_errors.append(filename)
            safe_gate = safe_gate and not gate_errors
            if gate.get('execution_policy',{}).get('required_solver_backend')=='cuda':
                safe_gate = safe_gate and budget['backend']=='CUDA' and gate['execution_policy']['allow_cpu_solver_fallback'] is False
                if hashlib.sha256((ROOT/contract['cuda_identity_path']).read_bytes()).hexdigest()!=contract['cuda_identity_sha256']:
                    gate_errors.append('cuda_identity_path')
                safe_gate = safe_gate and not gate_errors
        except (KeyError, OSError, TypeError) as exc:
            safe_gate = False
            gate_errors.append(type(exc).__name__)
    report = {"tracked_files": len(tracked), "entry_links_checked": links, "missing_links": missing_links,
              "missing_tracked_artifacts": missing_artifacts, "current_evidence_hash_mismatches": bad_hashes,
              "credential_pattern_files": secret_files, "missing_tracked_files": missing_tracked_files,
              "simulation_authorization_consistent": safe_gate, "authorization_errors": gate_errors,
              "checker_invoked_solver": False}
    report["passed"] = not (missing_links or missing_artifacts or bad_hashes or secret_files or missing_tracked_files) and safe_gate
    print(json.dumps(report, ensure_ascii=False))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
