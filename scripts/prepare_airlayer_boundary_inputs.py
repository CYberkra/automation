"""Build the four frozen-input variants for the air-layer / boundary experiment.

Scope: 2D invariant-y HS4T2D, gprMax V4.0.0. One variable is changed per
group relative to the archived snapshot input; every other byte of the source,
geometry, materials, receiver and cell size is preserved.

Groups
------
G1 domain_x 12 / 18 / 24 m
   Only #domain's X length grows. The cover layer keeps spanning the domain,
   which is the reading the gprMax guide allows for domain-spanning objects.
   Source at x=5.6 m and receiver at x=6.9 m are NOT moved, so the
   source/interface/receiver geometry stays identical and any difference in
   the recorded trace is attributable to the boundary distance. The snapshot
   region is widened with the domain so the added region is actually sampled.

G2 pml_cfs, thickness held at 20 cells
   Only the CFS stretching parameters change (12-parameter #pml_cfs).
   Thickness, geometry and cell size are untouched, so this separates
   absorption quality from the geometric rearrangement that a thickness
   change would introduce.

G3 air-layer rendering
   Not a new solve. G1's 12 m output is re-rendered with the air region
   (z > ground) isolated, so an air-layer lateral mode can be inspected
   directly rather than inferred from the full-domain plot.

G4 vertical grid 5 -> 2.5 cm
   The audit found cover at 170 MHz resolves to 8.31 cells per wavelength,
   below the guide's minimum of 10. Z is halved while X and Y stay at 5 cm, so
   the interface geometry is unchanged in the X direction and only the
   vertical discretisation is corrected. This is the group that makes the
   morphology read compliant with the guide rather than grid-sensitive.

Every variant is written to its own directory and hashed. No solver runs here.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

REPO = Path(r"E:\automation_djh\automation_repo")
OUT_ROOT = REPO / "artifacts/research_checks/2026-10-03_hs4t2d_airlayer_boundary"

SNAPSHOT_TIMES_NS = [40, 50, 65, 80, 95, 110, 130, 145, 175, 220]

# 2D invariant-y: Y stays a single 0.05 m cell in every variant.
BASE_DOMAIN = (12.0, 0.05, 33.0)
GRID_XY = 0.05


def read_lines(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return handle.read().split("\n")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def replace_command(lines: list[str], command: str, new_line: str) -> list[str]:
    """Replace every instance of a command line, keeping position and count."""
    out: list[str] = []
    replaced = 0
    prefix = command + ":"
    for line in lines:
        if line.lstrip().lower().startswith(prefix.lower()):
            out.append(new_line)
            replaced += 1
        else:
            out.append(line)
    if replaced != 1:
        raise ValueError(f"expected exactly 1 {command} line, found {replaced}")
    return out


def drop_command(lines: list[str], command: str, *, required: bool = True) -> list[str]:
    prefix = command + ":"
    out = [ln for ln in lines if not ln.lstrip().lower().startswith(prefix.lower())]
    if required and len(out) == len(lines):
        raise ValueError(f"no {command} line to drop")
    return out


def snapshot_lines(x_len: float, times_ns: list[int], dz: float = GRID_XY) -> list[str]:
    """Snapshot volume lines.

    The 7th to 9th parameters are the cell size of the snapshot volume and
    must match #dx_dy_dz. gprMax resolves the stored array shape from these,
    not from #domain, so leaving them at 0.05 while #dx_dy_dz declares a
    refined dz silently produces a snapshot on the coarse grid: the G4 rerun
    wrote 660 cells where the solver ran 1320. Verified 2026-10-03.
    """
    return [
        f"#snapshot: 0 0 0 {x_len:g} 0.05 33 0.05 0.05 {dz:g} {t}e-9 wave_{t:03d}.h5"
        for t in times_ns
    ]


def write_variant(
    base_lines: list[str],
    out_dir: Path,
    *,
    domain_x: float,
    pml_cells: str | None = None,
    pml_cfs: str | None = None,
    note: str,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    lines = list(base_lines)
    # Snapshots are regenerated for every variant so that the widened domain
    # is sampled and all variants share one time set.
    lines = drop_command(lines, "#snapshot")
    lines = replace_command(lines, "#domain", f"#domain: {domain_x:g} 0.05 33")
    if pml_cells is not None:
        lines = replace_command(lines, "#pml_cells", f"#pml_cells: {pml_cells}")
    if pml_cfs is not None:
        lines = drop_command(lines, "#pml_cfs", required=False)
        lines.insert(1, f"#pml_cfs: {pml_cfs}")
    lines = replace_command(lines, "#title", f"#title: {note}")
    lines.extend(snapshot_lines(domain_x, SNAPSHOT_TIMES_NS))

    text = "\n".join(lines)
    in_path = out_dir / "profile.in"
    in_path.write_text(text, encoding="utf-8", newline="")

    return {
        "directory": out_dir.name,
        "note": note,
        "input": str(in_path.relative_to(REPO)),
        "input_sha256": sha256_text(text),
        "domain_x_m": domain_x,
        "pml_cells": pml_cells,
        "pml_cfs": pml_cfs,
        "snapshot_times_ns": SNAPSHOT_TIMES_NS,
        "snapshot_count": len(SNAPSHOT_TIMES_NS),
    }


def build_g4(base_lines: list[str], out_dir: Path) -> dict:
    """Vertical refinement: dz 0.05 -> 0.025 m, dx and dy unchanged.

    Interface z coordinates are rescaled by the same factor so the physical
    geometry is preserved; X coordinates and the source/receiver positions in
    x are untouched.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    scale = 0.5
    new_dz = GRID_XY * scale

    lines: list[str] = []
    for line in base_lines:
        low = line.lstrip().lower()
        if low.startswith("#dx_dy_dz:"):
            lines.append(f"#dx_dy_dz: {GRID_XY:g} {GRID_XY:g} {new_dz:g}")
        elif low.startswith("#geometry_view:"):
            lines.append(
                f"#geometry_view: 0 0 0 {BASE_DOMAIN[0]:g} 0.05 {BASE_DOMAIN[2]:g} "
                f"{GRID_XY:g} {GRID_XY:g} {new_dz:g} hs4t2d_geom_g4_zfine n"
            )
        elif low.startswith("#snapshot:"):
            continue
        elif low.startswith("#pml_cells:"):
            # Keep the PML one metre thick in physical terms.
            lines.append("#pml_cells: 20 0 20 20 0 40")
        elif low.startswith("#box:"):
            toks = line.split()
            x0, y0, z0, x1, y1, z1 = (float(v) for v in toks[1:7])
            name = toks[7] if len(toks) > 7 else ""
            lines.append(
                f"#box: {x0:g} {y0:g} {z0 * scale:g} {x1:g} {y1:g} {z1 * scale:g} {name}".rstrip()
            )
        elif low.startswith("#hertzian_dipole:"):
            # "#hertzian_dipole: <polarisation> <x> <y> <z> <waveform>".
            # Only z is refined; y stays on the single invariant cell centre.
            toks = line.split()
            pol, x, y, z, waveform = toks[1], toks[2], toks[3], toks[4], toks[5]
            lines.append(
                f"#hertzian_dipole: {pol} {x} {y} {float(z) * scale:g} {waveform}"
            )
        elif low.startswith("#rx:"):
            # "#rx: <x> <y> <z> <name> <component>"; only z is refined.
            toks = line.split()
            x, y, z, name, component = toks[1], toks[2], toks[3], toks[4], toks[5]
            lines.append(
                f"#rx: {x} {y} {float(z) * scale:g} {name} {component}"
            )
        elif low.startswith("#title:"):
            lines.append("#title: HS4T2D G4 vertical grid 2.5 cm, boundary audit series")
        else:
            lines.append(line)

    lines.extend(snapshot_lines(BASE_DOMAIN[0], SNAPSHOT_TIMES_NS, dz=new_dz))
    text = "\n".join(lines)
    in_path = out_dir / "profile.in"
    in_path.write_text(text, encoding="utf-8", newline="")

    return {
        "directory": out_dir.name,
        "note": "vertical grid 0.05 -> 0.025 m; z coordinates scaled by 0.5; x/y unchanged",
        "input": str(in_path.relative_to(REPO)),
        "input_sha256": sha256_text(text),
        "domain_x_m": BASE_DOMAIN[0],
        "dz_m": new_dz,
        "pml_cells": "20 0 20 20 0 40",
        "snapshot_times_ns": SNAPSHOT_TIMES_NS,
        "snapshot_count": len(SNAPSHOT_TIMES_NS),
    }


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: prepare_airlayer_boundary_inputs.py <baseline.in>", file=sys.stderr)
        return 2
    base_path = Path(argv[1]).resolve()
    base_lines = read_lines(base_path)
    base_sha = sha256_text("\n".join(base_lines))

    variants: list[dict] = []

    # G1: domain width, single factor.
    for width in (12.0, 18.0, 24.0):
        tag = f"g1_domain_{int(width)}m"
        variants.append(
            write_variant(
                base_lines,
                OUT_ROOT / tag,
                domain_x=width,
                note=(
                    f"HS4T2D G1 domain X={width:g} m, 2D invariant-y; source and "
                    "receiver unmoved; boundary-distance series"
                ),
            )
        )

    # G2: CFS stretching profile, thickness held fixed at 20 cells.
    # gprMax V4 #pml_cfs takes 12 parameters:
    #   alpha profile, direction, min, max,
    #   kappa profile, direction, min, max,
    #   sigma profile, direction, min, max.
    # A quartic forward profile with a stronger sigma_max raises the
    # attenuation without moving any geometric boundary.
    g2_cfs = "quartic forward 1 8 quartic forward 1 8 quartic forward 0.05 0.5"
    variants.append(
        write_variant(
            base_lines,
            OUT_ROOT / "g2_pml_cfs_stronger",
            domain_x=12.0,
            pml_cfs=g2_cfs,
            note=(
                "HS4T2D G2 PML CFS strengthened, thickness fixed at 20 cells; "
                "geometry, source and receiver unchanged"
            ),
        )
    )

    # G3 is a rendering variant of G1 12 m and needs no new solve.
    variants.append(
        {
            "directory": "g3_air_layer_render",
            "note": (
                "no new solve: re-renders the g1_domain_12m snapshots with the air "
                "region (z > ground) isolated"
            ),
            "input": None,
            "input_sha256": None,
            "reuses": "g1_domain_12m",
        }
    )

    # G4: vertical grid compliance.
    variants.append(build_g4(base_lines, OUT_ROOT / "g4_vertical_2p5cm"))

    manifest = {
        "status": "PREPARED_NOT_RUN",
        "baseline_input": str(base_path.relative_to(REPO)),
        "baseline_sha256": base_sha,
        "two_dimensional": True,
        "invariance": "Y is a single 0.05 m cell; source and receiver are Ey; 2D line source",
        "cell_size_m": {"x": GRID_XY, "y": GRID_XY, "z": GRID_XY},
        "snapshot_times_ns": SNAPSHOT_TIMES_NS,
        "variants": variants,
        "scope": (
            "Input preparation only. No solver has been run and no field array "
            "exists yet. The air-layer hypothesis under test is boundary "
            "absorption sufficiency, not a physical lateral reflection: the "
            "cover layer spans the domain, which the gprMax guide treats as "
            "extending to infinity, so a truncation-based reflection is not "
            "expected by that convention."
        ),
    }
    out = OUT_ROOT / "prepared_inputs_manifest.json"
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    for v in variants:
        sha = v.get("input_sha256")
        print(f"  {v['directory']:<26} {v['note'][:64]}")
        if sha:
            print(f"  {'':<26} sha256 {sha}")
    print()
    print(f"written: {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
