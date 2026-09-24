"""Generate review-only V4 TMz inputs and estimates using stdlib; never import a solver.

Usage: python scripts/prepare_v4_calibration_packet.py --output <new directory>
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    runs = []
    for case in ('M00', 'M01'):
        for variant in ('base', 'fine', 'boundary', 'long'):
            dx = .0125 if variant == 'fine' else .025
            shift = 2 if variant == 'boundary' else 0
            lx, ly = 16+2*shift, 10+2*shift
            ground = 7+shift
            sx, sy, rx = 7.8+shift, 7.5+shift, 8.2+shift
            ns = 480 if variant == 'long' else 240
            pml = round(.3/dx)
            name = f'{case}_{variant}'
            text = (f'#title: REVIEW ONLY - V4 calibration {name}\n'
                    '#domain_mode: TM\n'
                    f'#domain: {lx} {ly} inf\n'
                    f'#dx_dy_dz: {dx} {dx} {dx}\n'
                    f'#time_window: {ns}e-9\n'
                    '#omp_threads: 8\n'
                    f'#pml_cells: {pml} {pml} 0 {pml} {pml} 0\n'
                    '#pml_formulation: HORIPML\n'
                    '#waveform: ricker 1 80e6 pulse80\n'
                    f'#hertzian_dipole: z {sx:g} {sy:g} inf pulse80\n'
                    f'#rx: {rx:g} {sy:g} inf measurement Ez\n')
            if case == 'M01':
                text += ('#material: 9 0 1 0 ideal_halfspace\n'
                         f'#box: 0 0 0 {lx} {ground} inf ideal_halfspace\n')
            # Geometry-only consistency checks, not the gprMax parser.
            for coordinate in (lx, ly, ground, sx, sy, rx):
                assert math.isclose(coordinate/dx, round(coordinate/dx), abs_tol=1e-9)
            assert .3 < sx < lx-.3 and .3 < rx < lx-.3 and .3 < sy < ly-.3
            assert math.isclose(rx-sx, .4, abs_tol=1e-12)
            assert sy-ground == .5 and math.isclose(pml*dx, .3, abs_tol=1e-12)
            path = args.output/f'{name}.in'
            path.write_text(text, encoding='utf-8', newline='\n')
            cells = round(lx/dx)*round(ly/dx)
            dt = dx/(299792458*math.sqrt(2))
            steps = math.ceil(ns*1e-9/dt)+1
            # A transparent lower-bound proxy, NOT allocated/RSS memory.
            six_field_bytes = 6*(round(lx/dx)+1)*(round(ly/dx)+1)*2*8
            runs.append(dict(id=name, input=path.name,
                             input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             case=case, variant=variant, dx_m=dx, domain_xy_m=[lx,ly],
                             ground_y_m=ground if case=='M01' else None,
                             source_xy_m=[sx,sy], receiver_xy_m=[rx,sy],
                             pml_cells_active=pml, requested_time_ns=ns,
                             active_xy_cells=cells, cfl_dt_estimate_s=dt,
                             steps_upper_estimate=steps,
                             xy_cell_step_proxy=cells*steps,
                             six_double_field_arrays_proxy_bytes=six_field_bytes))
    manifest = dict(schema='v4-calibration-packet/1', packet_id='V4-CAL-01',
                    status='superseded_draft_pending_sfcw_redesign', approved_to_simulate=False,
                    solver_imported=False, solver_executed=False,
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    run_count=len(runs), runs=runs,
                    proposed_backend='CPU, double, 8 OpenMP threads, sequential, no MPI',
                    proposed_limits=dict(per_run_wall_minutes=30, total_solver_wall_minutes=120,
                                         process_tree_rss_GiB=8, output_budget_GiB=1,
                                         retries=0, automatic_expansion=False),
                    limits_enforced=False,
                    prerequisites=['explicit user agreement on packet and run IDs',
                                   'separate Python 3.11-3.13 environment with verified V4 build',
                                   'runtime/source identity and full dependency record',
                                   'execution supervisor enforcing resource caps'],
                    estimates_note='CFL estimates and six-field storage proxies, not measured dt, runtime or total RAM',
                    syntax_basis='local V4 docs/source/input_hash_cmds.rst and examples_simple_2D.rst',
                    source_timing='built-in Ricker chi=sqrt(2)/80e6; current samples at (n+0.5)dt; no extra start delay')
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'packet':str(args.output),'inputs':len(runs),'solver_executed':False,
                      'total_xy_cell_steps':sum(r['xy_cell_step_proxy'] for r in runs)}))


if __name__ == '__main__':
    main()
