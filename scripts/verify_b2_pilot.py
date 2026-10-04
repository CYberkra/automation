# -*- coding: utf-8 -*-
"""B2 pilot batch post-run verification: NC-BG bitwise pairing + h5 integrity."""
import json
from pathlib import Path

import h5py
import numpy as np

REPO = Path(r"E:\automation_djh\automation_repo")
SIMS = REPO / "artifacts/simulations"

# (nc_mother, nc_date, bg_mother, bg_date)
PAIRS = [
    ("B2D-C1mX-NC", "2026-10-01", "B2D-C1mX-BG", "2026-10-01"),
    ("B2D-C3mX-NC", "2026-10-01", "B2D-C3mX-BG", "2026-09-28"),
    ("B2D-C3mS2X-NC", "2026-10-01", "B2D-C3mS2X-BG", "2026-09-28"),
    ("B2D-C3mS2TZX-NC", "2026-10-01", "B2D-C3mS2TZX-BG", "2026-09-28"),
]

cases = json.loads((REPO / "configs/research/batch2d_b2_pilot_co/cases.json").read_text(encoding="utf-8"))["cases"]
missing, bad = [], []
for c in cases:
    rid = c["run_id"]
    h5 = SIMS / f"2026-10-01_{rid}" / f"{rid}.h5"
    if not h5.exists():
        missing.append(rid)
        continue
    try:
        with h5py.File(h5, "r") as h:
            assert abs(float(h.attrs["dt"]) - 5.896635841874211e-11) < 1e-20
            n_rx = len(list(h["rxs"].items()))
            assert n_rx == 1
    except Exception as e:
        bad.append((rid, str(e)))
print(f"h5 integrity: {len(cases)} cases, missing={len(missing)}, bad={len(bad)}")
if missing: print("MISSING:", missing[:5])
if bad: print("BAD:", bad[:5])

print("\nNC-BG 配对差分（逐位比较 Ex）:")
for nc_m, nc_d, bg_m, bg_d in PAIRS:
    worst = 0.0
    n_bit_identical = 0
    for k in range(33):
        rid_nc = f"{nc_m}-CO33-t{k+1:02d}"
        rid_bg = f"{bg_m}-CO33-t{k+1:02d}"
        with h5py.File(SIMS / f"{nc_d}_{rid_nc}" / f"{rid_nc}.h5", "r") as h:
            ex_nc = list(h["rxs"].items())[0][1]["Ex"][:]
        with h5py.File(SIMS / f"{bg_d}_{rid_bg}" / f"{rid_bg}.h5", "r") as h:
            ex_bg = list(h["rxs"].items())[0][1]["Ex"][:]
        if ex_nc.shape == ex_bg.shape and np.array_equal(ex_nc, ex_bg):
            n_bit_identical += 1
        else:
            d = float(np.abs(ex_nc - ex_bg).max()) if ex_nc.shape == ex_bg.shape else float("nan")
            worst = max(worst, d)
    status = "逐位恒零 ✓" if n_bit_identical == 33 else f"异常! 仅{n_bit_identical}/33, max|diff|={worst}"
    print(f"  {nc_m} ({nc_d}) vs {bg_m} ({bg_d}): {n_bit_identical}/33 {status}")
