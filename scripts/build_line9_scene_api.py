"""Build the model_current x190 AGL5 model through the official gprMax v4 Python Scene API.

Parameter mapping from model_current/cases/x190_<role>/profile.in (hash-command chain):

    #title: ...                     -> gprMax.Title(name=...)
    #domain_mode: TM                -> gprMax.DomainMode(mode='TM')
    #domain: 210 42.5 inf           -> gprMax.Domain(p1=(210, 42.5, inf))
    #dx_dy_dz: 0.025 0.025 0.025    -> gprMax.Discretisation(p1=(0.025,)*3)
    #time_window: 2.4e-06           -> gprMax.TimeWindow(time=2.4e-6)  (or iterations= for short equivalence runs)
    #omp_threads: 2                 -> gprMax.OMPThreads(2)
    #pml_cells: 80 80 0 80 80 0     -> gprMax.PMLThickness(thickness=(80, 80, 0, 80, 80, 0))
    #pml_formulation: HORIPML       -> gprMax.PMLFormulation(formulation='HORIPML')
    #dispersive_averaging: y        -> gprMax.DispersiveAveraging(enabled=True)
    #waveform: impulse 1 1 impulse  -> gprMax.Waveform(wave_type='impulse', amp=1, freq=1, id='impulse')
    #hertzian_dipole: z x y 0.0125 impulse -> gprMax.HertzianDipole(p1=(x, y, 0.0125), polarisation='z', waveform_id='impulse')
    #rx: x y 0.0125 <id> Ez         -> gprMax.Rx(p1=(x, y, 0.0125), id=<id>, outputs=['Ez'])
    #geometry_objects_read: 0 0 0 <geo.h5> line9_materials_v0_2 y
                                    -> gprMax.GeometryObjectsRead(p1=(0, 0, 0), geofile=<abs path>,
                                         material_database='line9_materials_v0_2', averaging='y')

API facts verified against the installed 4.0.1 package on ROG (2026-10-11):
- Waveform kwargs: wave_type/amp/freq/id; docstring states an impulse is one nonzero
  sample on the source update lattice and its positive frequency is unused.
- GeometryObjectsRead resolves the JSON material database beside the HDF5 file and
  checks the file's MaterialDatabase attribute against the requested database name.
- Official API example pattern: gprMax._examples.features.impedance_surface.virtual_waveguide_2d
  (DomainMode/Domain(inf)/Discretisation/TimeWindow/PMLThickness/OMPThreads/run).
"""
import argparse
import json
from pathlib import Path

TX_M = (169.35, 35.95, 0.0125)
RX_M = (170.65, 36.025, 0.0125)
DOMAIN_M = (210.0, 42.5, float('inf'))
DL = 0.025
TIME_WINDOW_S = 2.4e-6
PML_CELLS = (80, 80, 0, 80, 80, 0)
DATABASE = 'line9_materials_v0_2'


def build_scene(case_id, geofile, iterations=None):
    import gprMax
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name=f'{case_id}; Scene API build of model_current x190 AGL5'))
    scene.add(gprMax.DomainMode(mode='TM'))
    scene.add(gprMax.Domain(p1=DOMAIN_M))
    scene.add(gprMax.Discretisation(p1=(DL, DL, DL)))
    if iterations is None:
        scene.add(gprMax.TimeWindow(time=TIME_WINDOW_S))
    else:
        scene.add(gprMax.TimeWindow(iterations=int(iterations)))
    scene.add(gprMax.OMPThreads(2))
    scene.add(gprMax.PMLFormulation(formulation='HORIPML'))
    scene.add(gprMax.PMLThickness(thickness=PML_CELLS))
    scene.add(gprMax.DispersiveAveraging(enabled=True))
    scene.add(gprMax.Waveform(wave_type='impulse', amp=1, freq=1, id='impulse'))
    scene.add(gprMax.HertzianDipole(p1=TX_M, polarisation='z', waveform_id='impulse'))
    scene.add(gprMax.Rx(p1=RX_M, id=f'{case_id}_rx1', outputs=['Ez']))
    scene.add(gprMax.GeometryObjectsRead(p1=(0, 0, 0), geofile=str(geofile),
                                         material_database=DATABASE, averaging='y'))
    return scene


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['run'])
    p.add_argument('--case', required=True, help='case id, e.g. x190_H0')
    p.add_argument('--geofile', type=Path, required=True, help='absolute path to H0.h5/H1.h5 (database JSON beside it)')
    p.add_argument('--outputfile', type=Path, required=True, help='output base path without .h5')
    p.add_argument('--iterations', type=int, default=None, help='override time window for short equivalence runs')
    p.add_argument('--geometry-only', action='store_true')
    p.add_argument('--gpu', type=int, default=None)
    a = p.parse_args()
    import gprMax
    scene = build_scene(a.case, a.geofile.resolve(), a.iterations)
    kwargs = dict(scenes=[scene], outputfile=a.outputfile.resolve(),
                  hide_progress_bars=True, geometry_only=a.geometry_only)
    if a.gpu is not None:
        kwargs.update(gpu=a.gpu, gpu_precision='double')
    gprMax.run(**kwargs)
    print(json.dumps(dict(status='SCENE_API_RUN_RETURNED', case=a.case,
                          geometry_only=a.geometry_only, iterations=a.iterations)))


if __name__ == '__main__':
    main()
