"""P2 batch: 13-station 3D B-scan with the x4-scaled GSSI-400-like antenna.

User (2026-10-05): B-scan with the real (gprMax built-in, x4-scaled) antenna,
rough + halfspace roles, same 3D domain/grid/PML/materials as the frozen P1
probe (scripts/p1_antenna_gssi400x4_3d.py). Stations x_abs = 15.25 + 0.5 k m,
k = 0..12, matching the archived 2D dense layout. GPU double, no subgrid.

Freeze discipline: batch root created once (exist_ok=False); per-trace input
manifest written before its solve; audit FAIL stops the batch (no retry).
"""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import gprMax
from antenna_gssi400_x4 import scaled_antenna, coarse_adapt
from hs_capsule_identity import sha256
from p1_antenna_gssi400x4_3d import (DL, DOMAIN, PML, TIME_WINDOW, ROCK, COVER,
                                     DEBYE, GROUND_ZP, SKID_ZP, X_SHIFT,
                                     geometry, manifest, audit_h5)

OUT = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p2'
N_STATIONS = 13
ROLES = ('rough', 'halfspace')


def build_scene(role, k):
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name=f'gssi400x4_p2_{role}_k{k:02d}'))
    scene.add(gprMax.Domain(p1=DOMAIN))
    scene.add(gprMax.Discretisation(p1=(DL, DL, DL)))
    scene.add(gprMax.TimeWindow(time=TIME_WINDOW))
    scene.add(gprMax.PMLFormulation(formulation='HORIPML'))
    scene.add(gprMax.PMLThickness(thickness=PML))
    scene.add(gprMax.Material(**ROCK))
    scene.add(gprMax.Material(**COVER))
    scene.add(gprMax.AddDebyeDispersion(**DEBYE))
    boxes = geometry(role)
    for b in boxes:
        scene.add(gprMax.Box(**b))
    ant_xyz = (15.25 + 0.5 * k - X_SHIFT, 0.5 * DOMAIN[1], SKID_ZP)
    ant_objs, adapt_log = coarse_adapt(scaled_antenna(*ant_xyz), DL)
    for o in ant_objs:
        scene.add(o)
    if k == 6:  # centre station geometry view only (vtkhdf volume cost)
        scene.add(gprMax.GeometryView(p1=(0, 0, 0), p2=DOMAIN, dl=(DL, DL, DL),
                                      output_type='n', filename=f'p2_{role}_geom'))
    return scene, boxes, adapt_log, ant_xyz


def main():
    if OUT.exists():
        raise SystemExit(f'{OUT} exists; frozen batch, no rerun')
    OUT.mkdir(parents=True)
    index = []
    for role in ROLES:
        for k in range(N_STATIONS):
            tag = f'{role}_k{k:02d}'
            tdir = OUT / tag
            tdir.mkdir()
            scene, boxes, adapt_log, ant_xyz = build_scene(role, k)
            man = manifest(role, boxes, adapt_log, ant_xyz)
            man['station_k'] = k
            man['station_x_abs'] = 15.25 + 0.5 * k
            man['batch_script_sha256'] = sha256(Path(__file__))
            (tdir / 'input_manifest.json').write_text(
                json.dumps(man, ensure_ascii=False, indent=1, default=str) + '\n',
                encoding='utf-8')
            gprMax.run(scenes=[scene], geometry_only=False,
                       outputfile=str(tdir / 'p2.h5'), gpu=[0],
                       gpu_precision='double', subgrid=False)
            info = audit_h5(tdir / 'p2.h5')
            info['station_k'] = k
            info['role'] = role
            (tdir / 'audit.json').write_text(json.dumps(info, ensure_ascii=False,
                                                        indent=1) + '\n', encoding='utf-8')
            index.append({'tag': tag, 'status': info['status'],
                          'h5_sha256': info['h5_sha256'],
                          'rx_Position': info['rx']['Position'],
                          'src_Position': info['src']['Position']})
            print(tag, info['status'], 'peak', f"{info['rx']['peak_abs']:.3e}")
            if info['status'] != 'PASS':
                (OUT / 'batch_index.json').write_text(
                    json.dumps(index, indent=1) + '\n', encoding='utf-8')
                raise SystemExit(f'{tag} audit FAIL; batch stopped, no retry')
    (OUT / 'batch_index.json').write_text(json.dumps(index, indent=1) + '\n',
                                          encoding='utf-8')
    print('P2 batch complete:', len(index), 'traces')


if __name__ == '__main__':
    main()
