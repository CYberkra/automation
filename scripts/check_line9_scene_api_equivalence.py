"""Equivalence checks between the hash-command (.in) chain and the Scene API chain.

Two approved checks (S3 verification scope, no long solve):
1. twin: a small 2D TM model is built twice -- once as an .in file, once through the
   Scene API -- with identical domain, grid, PML, waveform, dipole, rx and material;
   both run 2000 iterations at FP64 and the Ez traces must agree bit-exactly.
2. geometry: the real x190 H0 model is built through the API with geometry_only=True
   (GeometryObjectsRead against the real H5 + JSON database), verifying the API
   geometry import path on the production geometry without a solve.

Runs in the verified v401 env. Writes equivalence.json into --out.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

TWIN_IN = '''#title: api twin small model; equivalence check only
#domain_mode: TM
#domain: 2 2 inf
#dx_dy_dz: 0.025 0.025 0.025
#time_window: {time_window}
#omp_threads: 1
#pml_cells: 10 10 0 10 10 0
#pml_formulation: HORIPML
#dispersive_averaging: y
#material: 11 0.001 1 0 cover
#box: 0 0 0 2 1 0.025 cover
#waveform: impulse 1 1 impulse
#hertzian_dipole: z 0.5 1.5 0.0125 impulse
#rx: 1.5 1.5 0.0125 twin_rx Ez
'''
ITERATIONS = 2000
DL = 0.025
DT = DL / (299792458. * np.sqrt(2))


def build_twin_scene():
    import gprMax
    inf = float('inf')
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name='api twin small model; equivalence check only'))
    scene.add(gprMax.DomainMode(mode='TM'))
    scene.add(gprMax.Domain(p1=(2.0, 2.0, inf)))
    scene.add(gprMax.Discretisation(p1=(DL, DL, DL)))
    scene.add(gprMax.TimeWindow(iterations=ITERATIONS))
    scene.add(gprMax.OMPThreads(1))
    scene.add(gprMax.PMLFormulation(formulation='HORIPML'))
    scene.add(gprMax.PMLThickness(thickness=(10, 10, 0, 10, 10, 0)))
    scene.add(gprMax.DispersiveAveraging(enabled=True))
    scene.add(gprMax.Material(er=11, se=0.001, mr=1, sm=0, id='cover'))
    scene.add(gprMax.Box(p1=(0, 0, 0), p2=(2.0, 1.0, DL), material_id='cover'))
    scene.add(gprMax.Waveform(wave_type='impulse', amp=1, freq=1, id='impulse'))
    scene.add(gprMax.HertzianDipole(p1=(0.5, 1.5, 0.0125), polarisation='z', waveform_id='impulse'))
    scene.add(gprMax.Rx(p1=(1.5, 1.5, 0.0125), id='twin_rx', outputs=['Ez']))
    return scene


def read_ez(path):
    with h5py.File(path) as h:
        x = h['rxs/rx1/Ez'][:]
        dt = float(h.attrs['dt'])
        it = int(h.attrs['Iterations'])
    return x, dt, it


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--geofile', type=Path, default=None,
                   help='real H0.h5 for the geometry_only check (database JSON beside it)')
    a = p.parse_args()
    a.out.mkdir(parents=True)
    assert not any(a.out.iterdir()), 'out must be empty'
    import gprMax

    # --- twin check: identical small model through both chains ---
    tw = float((ITERATIONS - 1.5) * DT)  # .in chain does ceil(tw/dt)+1; this gives ITERATIONS
    in_path = a.out / 'twin.in'
    in_path.write_text(TWIN_IN.format(time_window=repr(tw)), encoding='utf-8')
    gprMax.run(inputfile=str(in_path), cpu_precision='double', hide_progress_bars=True)
    scene = build_twin_scene()
    gprMax.run(scenes=[scene], outputfile=a.out / 'twin_api', cpu_precision='double',
               hide_progress_bars=True)
    x_in, dt_in, it_in = read_ez(a.out / 'twin.h5')
    x_api, dt_api, it_api = read_ez(a.out / 'twin_api.h5')
    assert it_in == it_api == ITERATIONS and dt_in == dt_api
    max_abs = float(np.max(np.abs(x_in - x_api)))
    result = dict(twin=dict(iterations=ITERATIONS, bitwise_equal=bool(np.array_equal(x_in, x_api)),
                            max_abs_diff=max_abs, rms_in=float(np.sqrt(np.mean(x_in ** 2)))))

    # --- geometry check: real model through API, geometry_only ---
    if a.geofile is not None:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from build_line9_scene_api import build_scene
        scene = build_scene('x190_H0', Path(a.geofile).resolve(), iterations=10)
        gprMax.run(scenes=[scene], outputfile=a.out / 'api_geom_H0', geometry_only=True,
                   hide_progress_bars=True)
        result['geometry_only_real_model'] = 'BUILD_OK'

    (a.out / 'equivalence.json').write_text(
        json.dumps(dict(status='PASS' if result['twin']['bitwise_equal'] else 'DIFFERS',
                        **result), indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
