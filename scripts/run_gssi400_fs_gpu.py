"""GSSI 400 MHz official antenna model in free space — pipeline validation run.

Mirrors the official example
E:/gprMax-v.4.0.0/gprMax-v.4.0.0/examples/gpr/antennas/antenna_like_GSSI_400_fs.py
verbatim in geometry/source (2 mm, domain 0.340x0.340x0.318 m, 15 ns,
Stadler et al. 2022 updated model). Only changes: explicit output path,
GPU backend selection, gpu_precision kept at the official default 'single'.

Purpose: validate that this local gprMax v4 install + official antenna
framework + GPU backend produce a sane standard forward sample. It does NOT
represent the project SFCW device (20-170 MHz); the GSSI 400 band is far
above the project band and results must not be mixed into project samples.
"""

from pathlib import Path

import gprMax
from gprMax.toolboxes.GPRAntennaModels.GSSI import antenna_like_GSSI_400

OUT_DIR = Path("E:/automation_djh/automation_repo/artifacts/simulations/2026-09-27_GSSI400FS-VALID")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE = OUT_DIR / "GSSI400FS-VALID.h5"

# Discretisation (official example value)
dl = 0.002

# Domain (official example values)
x = 0.340
y = 0.340
z = 0.318

scene = gprMax.Scene()

title = gprMax.Title(name="GSSI400FS-VALID")
domain = gprMax.Domain(p1=(x, y, z))
dxdydz = gprMax.Discretisation(p1=(dl, dl, dl))
time_window = gprMax.TimeWindow(time=15e-9)

scene.add(title)
scene.add(domain)
scene.add(dxdydz)
scene.add(time_window)

# Import antenna model and add to model (official example position)
ant_pos = (0.170, 0.170, 0.100)
gssi_objects = antenna_like_GSSI_400(ant_pos[0], ant_pos[1], ant_pos[2], resolution=dl)
for obj in gssi_objects:
    scene.add(obj)

gprMax.run(scenes=[scene], geometry_only=False, outputfile=OUT_FILE, gpu=[0])
print("GSSI400FS-VALID done ->", OUT_FILE)
