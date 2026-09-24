"""Read-only source fingerprint and optional official-doc caching; never imports gprMax."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen


READS = {
    "gprMax/_version.py": "version declaration",
    "setup.py": "Python bounds and dependency metadata; not executed",
    "gprMax/__init__.py": "public run/Scene/study exports",
    "gprMax/gprMax.py": "run signature/defaults and CLI options",
    "gprMax/waveforms.py": "Waveform coefficients, Ricker and impulse evaluation",
    "gprMax/sources.py": "HertzianDipole waveform sampling and electric update only",
    "gprMax/fields_outputs.py": "source excitation and receiver timing/identity writers only",
    "gprMax/grid/fdtd_grid.py": "calculate_dt and reduced-domain declarations only",
    "README.rst": "installation, repository layout and Python requirements",
    "docs/source/migration_v3_v4.rst": "migration sections 1-6 and checklist",
    "docs/source/input_hash_cmds.rst": "domain/domain_mode, spatial/time steps, waveform, dipole, receiver commands",
    "docs/source/output.rst": "root metadata, receiver identity, fields and source timing sections",
    "docs/source/sources_ports.rst": "source families and fields/ports/time axes",
    "docs/source/studies.rst": "families, compatibility, restart and outputs",
    "docs/source/accelerators.rst": "optional runtimes, solver precision and feature restrictions",
    "docs/source/features.rst": "2D, source/study, material and subgrid sections",
    "gprMax/toolboxes/SFCW/README.rst": "full toolbox README, including finite-record validity",
    "examples/gpr/basic/cylinder_Ascan_2D.in": "full small official syntax example; not executed",
    "examples/gpr/basic/cylinder_Bscan_2D_study.csv": "acquisition header and example rows; not executed",
}
WEB = [
    ("index", "", "navigation and displayed 4.0.0 label"),
    ("migration", "migration_v3_v4.html", "migration guide"),
    ("installation", "inc_README.html", "installation and environment requirements"),
    ("input", "input_hash_cmds.html", "selected domain, source, receiver and time commands"),
    ("output", "output.html", "selected identity/timing/source metadata sections"),
    ("sources_ports", "sources_ports.html", "source families, physical quantities and timing"),
    ("studies", "studies.html", "study families and reuse restrictions"),
    ("sfcw", "inc_SFCW.html", "SFCW synthesis and finite-record limitations"),
    ("accelerators", "accelerators.html", "precision and backend restrictions"),
    ("modelling", "gprmodelling.html", "grid, coordinates and boundary guidance"),
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--docs-cache", type=Path, help="Explicitly enables downloading the listed official HTML pages")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root = args.source_root.resolve()
    if not (root / "gprMax/_version.py").is_file():
        raise ValueError("source_root_must_contain_gprMax_version_file")
    files = []
    for name, scope in READS.items():
        data = (root/name).read_bytes()
        files.append({"path": name, "bytes": len(data), "sha256": digest(data), "reading_scope": scope})
    match = re.search(r'__version__\s*=\s*[\"\']([^\"\']+)', (root/"gprMax/_version.py").read_text(encoding="utf-8"))
    # Fingerprint selected evidence, not an assertion about the entire distribution.
    tree = "\n".join(f"{x['path']}\t{x['sha256']}" for x in sorted(files, key=lambda x:x["path"]))
    web = []
    if args.docs_cache:
        args.docs_cache.mkdir(parents=True, exist_ok=False)

        def fetch(entry):
            name, suffix, scope = entry
            url = "https://docs.gprmax.com/en/latest/"+suffix
            row = {"id": name, "url": url, "reading_scope": scope}
            try:
                with urlopen(Request(url, headers={"User-Agent": "gpr-research-source-archive/1"}), timeout=25) as response:
                    body = response.read()
                    row.update({"resolved_url": response.geturl(), "content_type": response.headers.get("Content-Type")})
                path = args.docs_cache/f"{name}.html"
                path.write_bytes(body)
                row.update({"status": "cached", "bytes": len(body), "sha256": digest(body), "cache_filename": path.name})
            except Exception as exc:
                row.update({"status": "cache_unavailable", "error_type": type(exc).__name__, "message": str(exc)})
            return row

        with ThreadPoolExecutor(max_workers=4) as executor:
            web = list(executor.map(fetch, WEB))
    else:
        web = [{"id": n, "url": "https://docs.gprmax.com/en/latest/"+s, "reading_scope": r,
                "status": "not_downloaded_by_this_command"} for n,s,r in WEB]
    report = {"schema": "gprmax-v4-source-audit/1", "audited_utc": datetime.now(timezone.utc).isoformat(),
              "source_root_on_audit_machine": str(root), "declared_version": match.group(1) if match else None,
              "source_git_metadata_present": (root/".git").exists(), "upstream_commit": None,
              "official_archive_byte_equivalence_verified": False,
              "selected_evidence_fingerprint_sha256": digest(tree.encode()), "selected_files": files,
              "official_docs": web, "latest_is_mutable": True,
              "source_imported": False, "solver_invoked": False, "official_tests_executed": False,
              "installed_v4_runtime_verified": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"version": report["declared_version"], "selected_files": len(files),
                      "docs_cached": sum(x["status"]=="cached" for x in web), "solver_invoked": False}))


if __name__ == "__main__":
    main()
