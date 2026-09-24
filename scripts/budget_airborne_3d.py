"""Static V4 dense-grid budget, stdlib only. No gprMax import/allocation/solve."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--source-root', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise SystemExit('Refusing overwrite')
    cases = [('calibration', (24,24,24), 4, 400, 9),
             ('calibration_expanded', (32,32,28), 6, 400, 9),
             ('depth20_candidate', (40,40,50), 25, 1000, 16)]
    rows = []
    for name, dimensions, ground, window_ns, er in cases:
        for dx in (.1, .05, .025):
            n = [round(d/dx) for d in dimensions]
            cells = math.prod(n)
            nodes = math.prod(i+1 for i in n)
            # fdtd_grid.py initialise_geometry_arrays + initialise_field_arrays:
            # solid uint32; rigidE/H (12+6) int8; ID 6 uint32; 6 float64 fields.
            core_bytes = 22*cells + 72*nodes
            dt = dx/(299792458*math.sqrt(3))
            steps = math.ceil(window_ns*1e-9/dt)+1
            ytx = dimensions[1]/2-.7
            yrx = ytx+1.3
            assert math.isclose(yrx-ytx,1.3,abs_tol=1e-12)
            for v in (*dimensions,ground,ground+15,dimensions[0]/2,ytx,yrx):
                assert math.isclose(v/dx,round(v/dx),abs_tol=1e-8)
            rows.append(dict(id=f'{name}_dx{dx}', role=name, domain_xyz_m=dimensions,
                ground_z_m=ground, tx_xyz_m=[dimensions[0]/2,ytx,ground+15],
                rx_xyz_m=[dimensions[0]/2,yrx,ground+15], grid_m=dx,
                pml_proposed_thickness_m=1, pml_cells=round(1/dx),
                assumed_max_er=er, requested_time_ns=window_ns,
                cells=cells, nodes=nodes, core_array_bytes=core_bytes,
                core_array_GiB=core_bytes/2**30, cfl_dt_estimate_ps=dt*1e12,
                steps_estimate=steps, cell_steps_per_solve=cells*steps,
                cells_per_wavelength_at_170MHz=299792458/(170e6*math.sqrt(er)*dx),
                two_polarization_sequential_work=cells*steps*2,
                peak_memory_note='one solve at a time; two polarizations double work, not this peak'))
    result=dict(schema='airborne-3d-static-budget/1', solver_executed=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        v4_grid_source_sha256=hashlib.sha256((a.source_root/'gprMax/grid/fdtd_grid.py').read_bytes()).hexdigest(),
        core_formula_bytes='22*nx*ny*nz + 72*(nx+1)*(ny+1)*(nz+1) for double precision',
        exclusions=['PML auxiliary arrays','temporary/build arrays','dispersive poles',
                    'Python/runtime overhead','outputs','GPU/host duplication'],
        runtime_seconds=None, throughput_measured=False,
        resource_fact_from_prior_readonly_inventory_physical_RAM_bytes=67871469568,
        rows=rows, geometry_is_proposal=True,
        coordinate_convention='X flight; Y cross-track; Z up; planar ground',
        baseline_midpoint_note='Y centre minus 0.05m keeps both endpoints on all compared grids',
        checks={'baseline_length_and_grid_alignment':True,'case_grid_combinations':len(rows)})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
    for r in rows:
        print(f"{r['id']}: {r['cells']/1e6:.3f}M cells; core {r['core_array_GiB']:.2f} GiB; {r['steps_estimate']} steps; {r['cell_steps_per_solve']/1e12:.3f}e12 cell-steps; {r['cells_per_wavelength_at_170MHz']:.2f} cells/wavelength")


if __name__ == '__main__':
    main()
