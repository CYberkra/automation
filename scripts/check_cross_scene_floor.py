"""Cross-scene subtraction floor: minimal reproducers (2026-10-05).

Companion to docs/research/2026-10-05_cross_scene_subtraction_floor.md.
Three configurations of the same 1.6 x 1.6 x 3.2 m @4 cm scene (Hertzian
dipole y at (0.8,0.8,2.2), rx at (0.92,0.8,2.2), impulse, HORIPML, 100 ns):

- CASE 'same'      : identical scenes, one with an extra UNUSED material
                     definition -> expected bit-identical (solver deterministic).
- CASE 'pml_touch' : bottom layer z 0.4-0.6 m (touches bottom PML), rock vs dry
                     -> first rx difference at step ~19, far earlier than the
                     physical two-way time (~step 139).
- CASE 'interior'  : central block (0.6-1.0 m in x/y, z 0.8-1.0 m) clear of
                     every PML slab, rock vs dry -> first difference ~step 72,
                     still earlier than the physical two-way (~step 104),
                     consistent with the grid-speed precursor channel.

Diagnostic only; not a physical acceptance run.
"""
import argparse
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import gprMax
from hs_capsule_identity import sha256

CASES = ('same', 'pml_touch', 'interior')
DT_EXPECT = 7.703332806185881e-11


def build(case, variant):
    s = gprMax.Scene()
    s.add(gprMax.Domain(p1=(1.6, 1.6, 3.2)))
    s.add(gprMax.Discretisation(p1=(0.04, 0.04, 0.04)))
    s.add(gprMax.TimeWindow(time=100e-9))
    s.add(gprMax.PMLFormulation(formulation='HORIPML'))
    s.add(gprMax.Waveform(wave_type='impulse', amp=1.0, freq=100e6, id='imp'))
    s.add(gprMax.HertzianDipole(polarisation='y', p1=(0.8, 0.8, 2.2), waveform_id='imp'))
    s.add(gprMax.Rx(p1=(0.92, 0.8, 2.2), id='rx1', outputs=['Ey']))
    s.add(gprMax.Material(er=9, se=0.001, mr=1, sm=0, id='rock'))
    s.add(gprMax.Material(er=4.0, se=0.0005, mr=1, sm=0, id='dry'))
    if case == 'same':
        if variant == 'b':
            s.add(gprMax.Material(er=5.0, se=0.0, mr=1, sm=0, id='dummy_unused'))
        s.add(gprMax.Box(p1=(0, 0, 0.4), p2=(1.6, 1.6, 0.6), material_id='rock'))
    elif case == 'pml_touch':
        s.add(gprMax.Box(p1=(0, 0, 0.4), p2=(1.6, 1.6, 0.6),
                         material_id='rock' if variant == 'a' else 'dry'))
    else:
        s.add(gprMax.Box(p1=(0.6, 0.6, 0.8), p2=(1.0, 1.0, 1.0),
                         material_id='rock' if variant == 'a' else 'dry'))
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--gpu', action='store_true', help='use GPU double (default CPU double)')
    args = ap.parse_args()
    out = args.out
    if out.exists():
        raise SystemExit(f'{out} exists; new output directory required')
    out.mkdir(parents=True)
    run_kw = (dict(gpu=[0], gpu_precision='double', subgrid=False) if args.gpu
              else dict(cpu_precision='double'))
    report = {'case_order': [], 'script_sha256': sha256(Path(__file__))}
    for case in CASES:
        ys = {}
        for variant in ('a', 'b'):
            fn = out / f'{case}_{variant}.h5'
            gprMax.run(scenes=[build(case, variant)], geometry_only=False,
                       outputfile=str(fn), **run_kw)
            with h5py.File(fn) as h:
                ys[variant] = h['rxs/rx1/Ey'][:]
                dt = float(h.attrs['dt'])
        d = np.abs(ys['a'] - ys['b'])
        nz = np.nonzero(d > 0)[0]
        report['case_order'].append(case)
        report[case] = {
            'dt': dt, 'h5_a_sha256': sha256(out / f'{case}_a.h5'),
            'h5_b_sha256': sha256(out / f'{case}_b.h5'),
            'first_nonzero_diff_step': int(nz[0]) if len(nz) else None,
            'first_nonzero_diff_ns': float(nz[0] * dt * 1e9) if len(nz) else None,
            'max_abs_diff': float(d.max()),
            'nonzero_samples': int(len(nz))}
        print(case, json.dumps(report[case]))
    (out / 'floor_repro_report.json').write_text(
        json.dumps(report, indent=1) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
