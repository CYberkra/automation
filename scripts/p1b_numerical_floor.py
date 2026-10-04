"""P1b numerical-floor diagnostic: duplicate solve of the P1 rough scene.

Purpose: the rough-minus-halfspace contrast shows a -44 dB residual at the
surface-reflection time (~106 ns), which is non-physical (strict causality:
no influence from the bedrock relief can return before ~118 ns even at c).
Hypothesis: rough and halfspace grids contain different smoothed-material ID
sets, so the GPU update sequence is not bit-identical between roles and the
matched-background subtraction has a numerical floor. This run solves the
IDENTICAL rough scene a second time; diff(rough_a, rough_b) measures the pure
run-to-run numerical floor of the solver itself (expected ~0 if deterministic).

Declared: same scene construction as P1 rough (same script hash inputs),
different output directory; 95 s class run; diagnostic only.
"""
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import gprMax
from hs_capsule_identity import sha256
from p1_antenna_gssi400x4_3d import build_scene, manifest, audit_h5

OUT = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p1b'
P1 = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p1'


def main():
    if OUT.exists():
        raise SystemExit(f'{OUT} exists; frozen, no rerun')
    OUT.mkdir(parents=True)
    scene, boxes, adapt_log, ant_xyz = build_scene('rough')
    man = manifest('rough', boxes, adapt_log, ant_xyz)
    man['purpose'] = 'duplicate of P1 rough for run-to-run numerical floor'
    (OUT / 'input_manifest.json').write_text(json.dumps(man, ensure_ascii=False,
                                                       indent=1, default=str) + '\n',
                                             encoding='utf-8')
    gprMax.run(scenes=[scene], geometry_only=False, outputfile=str(OUT / 'p1b.h5'),
               gpu=[0], gpu_precision='double', subgrid=False)
    info = audit_h5(OUT / 'p1b.h5')
    with h5py.File(P1 / 'rough/p1.h5') as h:
        ya = h['rxs/rx1/Ey'][:]
        dt = float(h.attrs['dt'])
    with h5py.File(OUT / 'p1b.h5') as h:
        yb = h['rxs/rx1/Ey'][:]
    d = np.abs(ya - yb)
    info['rerun_diff'] = {
        'max_abs': float(d.max()), 'argmax_ns': float(np.argmax(d) * dt * 1e9),
        'nonzero_samples': int(np.count_nonzero(d)),
        'relative_to_signal_peak': float(d.max() / np.abs(ya).max())}
    (OUT / 'audit.json').write_text(json.dumps(info, ensure_ascii=False, indent=1) + '\n',
                                    encoding='utf-8')
    print(json.dumps(info['rerun_diff']))


if __name__ == '__main__':
    main()
