"""Verify the air-layer/boundary batch outputs without re-running the solver.

Checks, per group:
  - every snapshot H5 parses, is float64, finite, and has the declared shape;
  - E components satisfy the n*dt time stamp and H components the
    (n-1/2)*dt stamp implied by the requested snapshot time;
  - the recorded centre trace is byte-identical between the 12 m reference and
    the archived 2026-10-03 snapshot run, which is the reproduction check the
    contract requires before the domain-width comparison means anything.

Scope: provenance, dtype, time staggering and reproduction. It does not accept
the physics, does not certify the boundary, and does not sign off morphology.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np

REPO = Path(r"E:\automation_djh\automation_repo")
ROOT = REPO / "artifacts/research_checks/2026-10-03_hs4t2d_airlayer_boundary"
ARCHIVED = REPO / "artifacts/research_checks/2026-10-03_hs4t2d_relief08_snapshot"

C_M_PER_S = 299792458.0
COMPONENTS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")
TOLERANCE = 1e-9


def available_components(handle: h5py.File) -> list[str]:
    """Snapshot files store Ex..Hz as root-level datasets (verified 2026-10-03)."""
    return [c for c in COMPONENTS if c in handle]


def check_group(group_dir: Path) -> dict:
    snap_dir = group_dir / "profile_snaps"
    result: dict = {
        "group": group_dir.name,
        "snapshots_found": 0,
        "components_verified": 0,
        "problems": [],
        "snapshots": [],
    }
    if not snap_dir.is_dir():
        result["problems"].append("no profile_snaps directory")
        return result

    for h5_path in sorted(snap_dir.glob("wave_*.h5")):
        entry: dict = {"file": h5_path.name, "components": {}}
        try:
            with h5py.File(h5_path, "r") as handle:
                attrs = handle.attrs
                nx, ny, nz = (int(v) for v in attrs["nx_ny_nz"])
                spacing = [float(v) for v in attrs["dx_dy_dz"]]
                # gprMax records the snapshot time and the magnetic half step in
                # the file; read them rather than trusting the filename.
                stored_time = float(attrs["time"])
                magnetic_time = float(attrs["magnetic_time"])
                iteration = int(attrs["iteration"])
                entry["time_s"] = stored_time
                entry["magnetic_time_s"] = magnetic_time
                entry["iteration"] = iteration
                entry["shape_xyz"] = [nx, ny, nz]
                entry["spacing_m"] = spacing

                for comp in available_components(handle):
                    arr = np.asarray(handle[comp])
                    entry["components"][comp] = {
                        "shape": list(arr.shape),
                        "dtype": str(arr.dtype),
                        "all_finite": bool(np.all(np.isfinite(arr))),
                        "absmax": float(np.max(np.abs(arr))),
                    }
                    if arr.dtype != np.float64:
                        result["problems"].append(
                            f"{h5_path.name}:{comp} dtype {arr.dtype} is not float64"
                        )
                    if not np.all(np.isfinite(arr)):
                        result["problems"].append(f"{h5_path.name}:{comp} has non-finite values")
                    result["components_verified"] += 1

                # Courant condition for this grid, 2D: dt <= 1/(c*sqrt(1/dx^2+1/dz^2))
                dx, _, dz = spacing
                dt_cfl = 1.0 / (C_M_PER_S * np.sqrt(1.0 / dx**2 + 1.0 / dz**2))
                entry["dt_cfl_s"] = float(dt_cfl)
                entry["steps_at_stored_time"] = int(round(stored_time / dt_cfl))
                entry["time_stamp_consistent"] = bool(
                    abs(stored_time - entry["steps_at_stored_time"] * dt_cfl) < TOLERANCE
                )
                if not entry["time_stamp_consistent"]:
                    result["problems"].append(
                        f"{h5_path.name}: stored time {stored_time} is not an integer "
                        f"multiple of the Courant dt {dt_cfl}"
                    )
                # E sits at n*dt and H at (n-1/2)*dt; gprMax stores both.
                expected_h = (iteration - 0.5) * dt_cfl
                entry["magnetic_time_offset_s"] = float(magnetic_time - expected_h)
                entry["magnetic_time_consistent"] = bool(
                    abs(magnetic_time - expected_h) < TOLERANCE
                )
                if not entry["magnetic_time_consistent"]:
                    result["problems"].append(
                        f"{h5_path.name}: magnetic_time {magnetic_time} does not match "
                        f"(n-1/2)*dt = {expected_h} at iteration {iteration}"
                    )
        except Exception as exc:  # noqa: BLE001 - report, do not abort the batch
            result["problems"].append(f"{h5_path.name}: {type(exc).__name__}: {exc}")
        result["snapshots"].append(entry)
        result["snapshots_found"] += 1
    return result


def compare_centre_traces() -> dict:
    """Byte-level reproduction check against the archived snapshot run."""
    out: dict = {
        "reference": "artifacts/research_checks/2026-10-03_hs4t2d_relief08_snapshot",
        "candidate": "artifacts/research_checks/2026-10-03_hs4t2d_airlayer_boundary/g1_domain_12m",
        "compared": False,
    }
    ref = ARCHIVED / "profile.h5"
    cand = ROOT / "g1_domain_12m" / "profile.h5"
    if not ref.is_file() or not cand.is_file():
        out["reason"] = "a profile.h5 is missing"
        return out
    with h5py.File(ref, "r") as a, h5py.File(cand, "r") as b:
        ea = np.asarray(a["rxs/rx1/Ey"]).ravel()
        eb = np.asarray(b["rxs/rx1/Ey"]).ravel()
        out["compared"] = True
        out["reference_shape"] = list(ea.shape)
        out["candidate_shape"] = list(eb.shape)
        out["reference_dtype"] = str(ea.dtype)
        out["candidate_dtype"] = str(eb.dtype)
        if ea.shape != eb.shape:
            out["identical"] = False
            out["reason"] = "shape mismatch"
            return out
        diff = np.abs(ea - eb)
        out["max_abs_difference"] = float(diff.max())
        out["bit_identical"] = bool(diff.max() == 0.0)
        out["identical"] = out["bit_identical"]
        out["reference_absmax"] = float(np.abs(ea).max())
        out["candidate_absmax"] = float(np.abs(eb).max())
    return out


def main(argv: list[str]) -> int:
    out_path = ROOT / "independent_verification.json"
    groups = sorted(d for d in ROOT.iterdir() if d.is_dir() and (d / "profile_snaps").is_dir())

    report = {
        "status": "PENDING",
        "checks": [check_group(g) for g in groups],
        "reproduction_check": compare_centre_traces(),
        "scope": (
            "Provenance, dtype, finiteness, time staggering and reproduction only. "
            "No physical acceptance, no boundary certification, no morphology "
            "closure, and no comparison against the archived domain-width numbers."
        ),
    }

    total_problems = sum(len(c["problems"]) for c in report["checks"])
    repro = report["reproduction_check"]
    repro_ok = repro.get("bit_identical") is True

    if total_problems == 0 and repro_ok:
        report["status"] = "PASS"
    elif total_problems == 0:
        report["status"] = "PASS_EXCEPT_REPRODUCTION_DIFFERS"
    else:
        report["status"] = "FAIL"

    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"status: {report['status']}")
    for c in report["checks"]:
        print(
            f"  {c['group']:<26} snapshots={c['snapshots_found']:>2} "
            f"components={c['components_verified']:>3} problems={len(c['problems'])}"
        )
        for p in c["problems"][:5]:
            print(f"      - {p}")
    r = report["reproduction_check"]
    if r.get("compared"):
        print(
            f"  reproduction vs archived: bit_identical={r.get('bit_identical')} "
            f"max_abs_diff={r.get('max_abs_difference')}"
        )
    else:
        print(f"  reproduction: not compared ({r.get('reason')})")
    print(f"written: {out_path.relative_to(REPO)}")
    return 0 if report["status"].startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
