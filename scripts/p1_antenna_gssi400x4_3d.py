"""P1 probe: x4-scaled GSSI-400-like antenna over the 3D HS4 model, centre station.

User (2026-10-05): real antenna (gprMax built-in GSSI-400, scaled x4 to the
~100 MHz class) replacing the Hertzian dipole, 15 m altitude, rough/halfspace
pair at the centre station to measure cost and behaviour before the 13-station
B-scan batch (P2). GPU double precision, no subgrid (V4 subgrids are CPU-only).

Declared coordinate shifts (grid alignment with the 2D archive):
    x' = x_abs - 11.6 m, z' = z_abs - 6.0 m; y unchanged.
    Domain 12.0 x 3.2 x 24.0 m at 0.04 m -> 300 x 80 x 600 cells.
Declared approximations:
    - antenna thin details quantized to the 4 cm grid (coarse_adapt log);
    - relief extruded invariantly along y from the archived 2D transect
      (interface_relief of the centre_rough profile), z quantized to 0.04 m;
    - halfspace role = uniform cover subsurface (no bedrock), matching the
      archived 2D halfspace semantics;
    - PML 50 cells in x, 25 in z (matching the 2 m / 1 m PML of the 2D runs),
      10 cells in y (y is the invariant direction; declared choice);
    - antenna skid bottom at z_abs = 27 m (15 m altitude), centre station
      k=6 -> x_abs = 18.25 m; source/receiver snap to the 4 cm grid, actual
      GridPosition is read back from the H5 and recorded.
Freeze discipline: each role directory is created once (exist_ok=False); the
input manifest (script/module hashes, constants, per-object hashes) is written
before solving. One fresh attempt; no automatic retry.
"""
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import gprMax
from antenna_gssi400_x4 import scaled_antenna, coarse_adapt, CASE_SIZE
from hs_capsule_identity import sha256
from plot_hs4_permittivity_bscans import interface_relief

OUT = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p1'
DL = 0.04
X_SHIFT, Z_SHIFT = 11.6, 6.0
DOMAIN = (12.0, 3.2, 24.0)          # shifted frame
GROUND_ZP = 12.0 - Z_SHIFT          # 6.0
SKID_ZP = 27.0 - Z_SHIFT            # 21.0  (15 m altitude)
PML = (50, 10, 25, 50, 10, 25)      # x0 y0 z0 x1 y1 z1 cells
TIME_WINDOW = 600e-9
STATION_K = 6
STATION_X_ABS = 15.25 + 0.5 * STATION_K   # 18.25
ROCK = dict(er=9, se=0.001, mr=1, sm=0, id='rock')
COVER = dict(er=18.017, se=0.003, mr=1, sm=0, id='cover')
DEBYE = dict(poles=1, er_delta=(7.878,), tau=(6.4567e-9,), material_ids=['cover'])


def relief_z_abs(x_abs):
    """Piecewise-constant interface bottom z (abs frame) from the 2D archive."""
    xs, zs = interface_relief()
    i = np.searchsorted(xs, x_abs, side='right') - 1
    i = min(max(i, 0), len(zs) - 1)
    return float(zs[i])


def snap(v):
    return round(round(v / DL) * DL, 10)


def geometry(role):
    """Subsurface boxes in the shifted frame; yields dict records for audit."""
    nx = int(round(DOMAIN[0] / DL))
    boxes = []
    for i in range(nx):
        x0a = X_SHIFT + i * DL
        xc = x0a + 0.5 * DL
        zif = relief_z_abs(xc) if role == 'rough' else 12.0
        zif_p = snap(zif - Z_SHIFT)
        if role == 'rough' and zif_p > 0:
            boxes.append(dict(p1=(x0a - X_SHIFT, 0.0, 0.0), p2=(x0a - X_SHIFT + DL, DOMAIN[1], zif_p),
                              material_id='rock'))
        if role == 'halfspace':
            zif_p = 0.0
        boxes.append(dict(p1=(x0a - X_SHIFT, 0.0, zif_p), p2=(x0a - X_SHIFT + DL, DOMAIN[1], GROUND_ZP),
                          material_id='cover'))
    return boxes


def build_scene(role):
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name=f'gssi400x4_p1_{role}_k{STATION_K}'))
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
    ant_xyz = (STATION_X_ABS - X_SHIFT, 0.5 * DOMAIN[1], SKID_ZP)
    ant_objs, adapt_log = coarse_adapt(scaled_antenna(*ant_xyz), DL)
    for o in ant_objs:
        scene.add(o)
    scene.add(gprMax.GeometryView(p1=(0, 0, 0), p2=DOMAIN, dl=(DL, DL, DL),
                                  output_type='n', filename=f'p1_{role}_geom'))
    return scene, boxes, adapt_log, ant_xyz


def manifest(role, boxes, adapt_log, ant_xyz):
    def h(s):
        return hashlib.sha256(s.encode('utf-8')).hexdigest()
    obj_hashes = [h(json.dumps(b, sort_keys=True)) for b in boxes]
    return {'role': role, 'station_k': STATION_K, 'station_x_abs': STATION_X_ABS,
            'shifts': {'x': X_SHIFT, 'z': Z_SHIFT}, 'domain_m': DOMAIN, 'dl_m': DL,
            'pml_cells': PML, 'time_window_s': TIME_WINDOW,
            'ground_z_abs': 12.0, 'skid_z_abs': 27.0,
            'antenna_anchor_shifted': ant_xyz, 'antenna_case_m': CASE_SIZE,
            'materials': {'rock': ROCK, 'cover': COVER, 'debye': DEBYE},
            'n_subsurface_boxes': len(boxes),
            'subsurface_boxes_sha256_merkle': h(''.join(obj_hashes)),
            'coarse_adapt_log': adapt_log,
            'script_sha256': sha256(Path(__file__)),
            'antenna_module_sha256': sha256(ROOT / 'scripts/antenna_gssi400_x4.py'),
            'relief_source': 'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre/centre_rough/profile.in'}


def audit_h5(h5path):
    info = {'h5_sha256': sha256(h5path)}
    with h5py.File(h5path) as h:
        info['version'] = str(h.attrs['gprMax'])
        info['nx_ny_nz'] = h.attrs['nx_ny_nz'].tolist()
        info['dx_dy_dz'] = h.attrs['dx_dy_dz'].tolist()
        info['dt'] = float(h.attrs['dt'])
        info['iterations'] = int(h.attrs['Iterations'])
        y = h['rxs/rx1/Ey'][:]
        info['rx'] = {'dtype': str(y.dtype), 'len': int(y.shape[0]),
                      'finite': bool(np.isfinite(y).all()),
                      'peak_abs': float(np.max(np.abs(y))),
                      'GridPosition': h['rxs/rx1'].attrs['GridPosition'].tolist(),
                      'Position': h['rxs/rx1'].attrs['Position'].tolist()}
        info['src'] = {'GridPosition': h['srcs/src1'].attrs['GridPosition'].tolist(),
                       'Position': h['srcs/src1'].attrs['Position'].tolist()}
        ok = (info['version'] == '4.0.0' and y.dtype == np.float64 and info['rx']['finite']
              and info['rx']['len'] == info['iterations']
              and info['nx_ny_nz'] == [300, 80, 600]
              and np.allclose(info['dx_dy_dz'], [DL, DL, DL])
              and info['rx']['peak_abs'] > 0)
        info['status'] = 'PASS' if ok else 'FAIL'
    return info


def main():
    if OUT.exists():
        raise SystemExit(f'{OUT} exists; frozen batch, no rerun')
    OUT.mkdir(parents=True)
    for role in ('rough', 'halfspace'):
        rdir = OUT / role
        rdir.mkdir()
        scene, boxes, adapt_log, ant_xyz = build_scene(role)
        (rdir / 'input_manifest.json').write_text(
            json.dumps(manifest(role, boxes, adapt_log, ant_xyz), ensure_ascii=False,
                       indent=1, default=str) + '\n', encoding='utf-8')
        gprMax.run(scenes=[scene], geometry_only=False,
                   outputfile=str(rdir / 'p1.h5'), gpu=[0], gpu_precision='double',
                   subgrid=False)
        info = audit_h5(rdir / 'p1.h5')
        (rdir / 'audit.json').write_text(json.dumps(info, ensure_ascii=False, indent=1) + '\n',
                                         encoding='utf-8')
        print(role, json.dumps(info, ensure_ascii=False)[:600])
        if info['status'] != 'PASS':
            raise SystemExit(f'{role} audit FAIL; batch left for inspection, no retry')


if __name__ == '__main__':
    main()
