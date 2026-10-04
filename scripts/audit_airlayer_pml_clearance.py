"""Pre-solver audit: distances from geometry, source and receivers to the PML.

The gprMax V4 modelling guide states two conditions for trustworthy absorbing
boundaries (docs.gprmax.com/en/latest/gprmodelling.html, read 2026-10-03):

  1. "all sources and targets are kept at least 15 cells away from them [the
     ABCs]";
  2. "free space (i.e. air) should be always included above a source for at
     least 15-20 cells".

It also states that objects spanning the computational domain "are assumed to
extend to infinity", so boundary terminations of such objects are not
themselves reflection sources - only ABC imperfection is.

This script computes those distances from the frozen input so the
domain-width experiment can be interpreted correctly: if any object, the
source or a receiver sits closer than 15 cells to a PML inner face, the
"spanning object extends to infinity" reading no longer holds for that
configuration and a domain-width difference cannot be attributed to the ABC
alone.

Reads only. No solver is invoked and no field array is produced.
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

REPO = Path(r"E:\automation_djh\automation_repo")

# Official guide thresholds, cells.
MIN_PML_SEPARATION_CELLS = 15
MIN_FREE_SPACE_CELLS = 15
FREE_SPACE_PREFERRED_CELLS = 20

# Shortest wavelength in the slowest (highest permittivity) declared medium,
# used for the cells-per-wavelength check at the top of the SFCW band.
SFCW_MAX_FREQ_HZ = 170e6
SFCW_MIN_FREQ_HZ = 20e6
GUIDE_MIN_CELLS_PER_WAVELENGTH = 10
GUIDE_PREFERRED_CELLS_PER_WAVELENGTH = 20
C_M_PER_S = 299792458.0


def parse_commands(in_path: Path) -> dict:
    cmds: dict[str, list[str]] = {}
    with in_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if not line.startswith("#"):
                raise ValueError(f"unexpected non-command line: {line!r}")
            parts = line[1:].split(None, 1)
            name = parts[0]
            # gprMax commands are written as "#name: value"; strip the colon.
            if name.endswith(":"):
                name = name[:-1]
            name = name.lower()
            value = parts[1] if len(parts) > 1 else ""
            cmds.setdefault(name, []).append(value)
    return cmds


def parse_boxes(cmds: dict) -> list[tuple[float, float, float, float, float, float]]:
    boxes = []
    for value in cmds.get("box", []):
        tokens = value.split()
        if len(tokens) < 6:
            raise ValueError(f"#box needs 6 coordinates, got {value!r}")
        # A #box line may carry a trailing material name after the 6 coordinates.
        # Order is p1=(x0,y0,z0) then p2=(x1,y1,z1), per
        # user_objects/cmds_geometry/box.py::Box. Return named corners so the
        # x/z indices cannot be mixed up again downstream.
        x0, _y0, z0, x1, _y1, z1 = (float(t) for t in tokens[:6])
        boxes.append((x0, z0, x1, z1))
    return boxes


def cells_per_wavelength(eps_r: float, freq_hz: float, cell: float) -> float:
    """Cells per wavelength for a non-magnetic medium.

    lambda = c / (f * sqrt(eps_r)), so cells-per-wavelength = lambda / cell.
    The guide's rule of thumb is >=10 cells per wavelength at the shortest
    wavelength, 20 preferred.
    """
    wavelength_m = C_M_PER_S / (freq_hz * math.sqrt(eps_r))
    return wavelength_m / cell


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: audit_airlayer_pml_clearance.py <input.in>", file=sys.stderr)
        return 2
    in_path = Path(argv[1]).resolve()
    out_path = in_path.parent / "pml_clearance_audit.json"

    cmds = parse_commands(in_path)
    domain = [float(v) for v in cmds["domain"][0].split()]
    x_len, _, z_len = domain[0], domain[1], domain[2]
    dxyz = [float(v) for v in cmds["dx_dy_dz"][0].split()]
    dx, _, dz = dxyz

    pml_raw = [int(v) for v in cmds["pml_cells"][0].split()]
    if len(pml_raw) == 1:
        pml_raw = pml_raw * 6
    if len(pml_raw) != 6:
        raise ValueError(f"#pml_cells needs 1 or 6 values, got {cmds['#pml_cells']!r}")
    # gprMax V4 order per user_objects/cmds_singleuse.py::PMLThickness and
    # grid/fdtd_grid.py::_validate_pml_thickness: x0, y0, z0, xmax, ymax, zmax.
    pml_x0, pml_y0, pml_z0, pml_xmax, pml_ymax, pml_zmax = pml_raw

    formulation = cmds.get("pml_formulation", ["(unset)"])[0]
    cfs_present = "pml_cfs" in cmds
    cfs_lines = cmds.get("pml_cfs", [])

    boxes = parse_boxes(cmds)

    # Inner faces of the PML, in metres.
    x_pml_inner_min = pml_x0 * dx
    x_pml_inner_max = x_len - pml_xmax * dx
    z_pml_inner_min = pml_z0 * dz
    z_pml_inner_max = z_len - pml_zmax * dz

    box_x_min = min(b[0] for b in boxes)
    box_x_max = max(b[2] for b in boxes)
    box_z_min = min(b[1] for b in boxes)
    box_z_max = max(b[3] for b in boxes)

    src = cmds["hertzian_dipole"][0].split()
    src_x, _, src_z = float(src[1]), float(src[2]), float(src[3])
    receivers = []
    for value in cmds.get("rx", []):
        f = value.split()
        receivers.append(
            {
                "name": f[3] if len(f) > 3 else None,
                "x_m": float(f[0]),
                "z_m": float(f[2]),
                "component": f[4] if len(f) > 4 else None,
                "clearance_x_min_m": float(f[0]) - x_pml_inner_min,
                "clearance_x_max_m": x_pml_inner_max - float(f[0]),
            }
        )

    # Air layer above the source: the source is at z=27 m, ground at z=12 m.
    air = [b for b in boxes if b[1] >= box_z_max - 1e-9]
    ground_z = min((b[1] for b in boxes if b[1] > 0), default=None)
    air_above_source_m = (src_z - ground_z) if ground_z is not None else None

    # Cells per wavelength in every declared medium across the SFCW band.
    media = []
    for value in cmds.get("material", []):
        toks = value.split()
        # "#material: <er> <se> <murt> <muz> <name>", per
        # user_objects/cmds_multiuse.py::Material (er, se, murt, muz).
        er = float(toks[0])
        se = float(toks[1]) if len(toks) > 1 else 0.0
        media.append(
            {
                "name": toks[4] if len(toks) > 4 else f"material_{len(media)}",
                "eps_r": er,
                "conductivity_S_per_m": se,
                "cells_per_wavelength_at_20MHz": cells_per_wavelength(er, SFCW_MIN_FREQ_HZ, dx),
                "cells_per_wavelength_at_170MHz": cells_per_wavelength(er, SFCW_MAX_FREQ_HZ, dx),
            }
        )

    nearest_rx_to_xmax = min(r["clearance_x_max_m"] for r in receivers)

    checks = []

    def add(name: str, value: dict, ok: bool | None, note: str) -> None:
        checks.append({"check": name, **value, "pass": ok, "note": note})

    add(
        "source_to_pml_inner_x_min",
        {"metres": src_x - x_pml_inner_min, "cells": (src_x - x_pml_inner_min) / dx},
        (src_x - x_pml_inner_min) / dx >= MIN_PML_SEPARATION_CELLS,
        f"guide requires >= {MIN_PML_SEPARATION_CELLS} cells from source to ABC",
    )
    add(
        "source_to_pml_inner_x_max",
        {"metres": x_pml_inner_max - src_x, "cells": (x_pml_inner_max - src_x) / dx},
        (x_pml_inner_max - src_x) / dx >= MIN_PML_SEPARATION_CELLS,
        f"guide requires >= {MIN_PML_SEPARATION_CELLS} cells from source to ABC",
    )
    add(
        "nearest_receiver_to_pml_inner_x_max",
        {"metres": nearest_rx_to_xmax, "cells": nearest_rx_to_xmax / dx},
        nearest_rx_to_xmax / dx >= MIN_PML_SEPARATION_CELLS,
        "closest receiver to the right-hand PML; governs the usable trace window",
    )
    add(
        "air_above_source",
        {
            "metres": air_above_source_m,
            "cells": (air_above_source_m / dz) if air_above_source_m is not None else None,
        },
        (air_above_source_m / dz) >= FREE_SPACE_PREFERRED_CELLS
        if air_above_source_m is not None
        else None,
        f"guide requires >= {MIN_FREE_SPACE_CELLS}-{FREE_SPACE_PREFERRED_CELLS} cells of air above source",
    )
    add(
        "geometry_spans_domain_x",
        {
            "box_x_min_m": box_x_min,
            "box_x_max_m": box_x_max,
            "domain_x_m": x_len,
            "touches_x_min": abs(box_x_min) < 1e-9,
            "touches_x_max": abs(box_x_max - x_len) < 1e-9,
        },
        abs(box_x_min) < 1e-9 and abs(box_x_max - x_len) < 1e-9,
        "guide: objects spanning the domain are assumed to extend to infinity, "
        "so their truncation is not itself a reflection source; this reading "
        "depends on the clearance checks above also passing",
    )
    add(
        "pml_cfs_explicit",
        {"present": cfs_present, "lines": len(cfs_lines), "formulation": formulation},
        None,
        "CFS parameters are left at the formulation default; recorded because the "
        "boundary-absorption experiment varies these while holding thickness fixed",
    )
    for m in media:
        add(
            f"cells_per_wavelength_{m['name']}_at_170MHz",
            {
                "eps_r": m["eps_r"],
                "cells_per_wavelength": m["cells_per_wavelength_at_170MHz"],
                "guide_minimum": GUIDE_MIN_CELLS_PER_WAVELENGTH,
                "guide_preferred": GUIDE_PREFERRED_CELLS_PER_WAVELENGTH,
            },
            m["cells_per_wavelength_at_170MHz"] >= GUIDE_MIN_CELLS_PER_WAVELENGTH,
            "below the guide minimum the band is under-resolved and arrival times "
            "carry grid sensitivity; this is independent of any boundary effect",
        )

    report = {
        "status": "READ_ONLY_AUDIT",
        "input": str(in_path),
        "domain_m": {"x": x_len, "y": domain[1], "z": z_len},
        "cell_size_m": {"x": dx, "z": dz},
        "pml_cells": {
            "x0": pml_x0,
            "xmax": pml_xmax,
            "y0": pml_y0,
            "ymax": pml_ymax,
            "z0": pml_z0,
            "zmax": pml_zmax,
            "formulation": formulation,
        },
        "pml_inner_faces_m": {
            "x_min": x_pml_inner_min,
            "x_max": x_pml_inner_max,
            "z_min": z_pml_inner_min,
            "z_max": z_pml_inner_max,
        },
        "geometry_extent_m": {
            "x_min": box_x_min,
            "x_max": box_x_max,
            "z_min": box_z_min,
            "z_max": box_z_max,
            "box_count": len(boxes),
            "air_layer_box_count": len(air),
        },
        "source_m": {"x": src_x, "z": src_z, "ground_z": ground_z},
        "receivers": receivers,
        "media": media,
        "checks": checks,
        "scope": (
            "Distance and resolution audit against the gprMax V4 modelling guide. "
            "No solver run, no field array, no physical acceptance. Passing the "
            "clearance checks does not certify the boundary; it only establishes "
            "whether the 'spanning geometry extends to infinity' reading is "
            "available for interpreting a domain-width comparison."
        ),
    }

    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"input: {in_path}")
    print(f"PML inner faces: x=[{x_pml_inner_min:.3f}, {x_pml_inner_max:.3f}] m")
    print(f"geometry extent: x=[{box_x_min:.3f}, {box_x_max:.3f}] m, {len(boxes)} boxes")
    print()
    width = max(len(c["check"]) for c in checks)
    for c in checks:
        flag = {True: "PASS", False: "FAIL", None: "INFO"}[c["pass"]]
        detail = ""
        if "metres" in c and c["metres"] is not None:
            detail = f"{c['metres']:.3f} m"
            if c.get("cells") is not None:
                detail += f" / {c['cells']:.1f} cells"
        elif "cells_per_wavelength" in c:
            detail = f"{c['cells_per_wavelength']:.2f} cells/wavelength"
        elif "present" in c:
            detail = f"present={c['present']} lines={c['lines']}"
        print(f"  {flag:4}  {c['check']:<{width}}  {detail}")
    print()
    print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
