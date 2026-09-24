"""Static PML budget and four review-only 3D inputs. Does not import/run gprMax."""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def pml_values(nx, ny, nz, p):
    # pml.py initialise_field_arrays, one CFS term, two slabs per axis.
    x = (p+1)*ny*(nz+1)+(p+1)*(ny+1)*nz+p*(ny+1)*nz+p*ny*(nz+1)
    y = nx*(p+1)*(nz+1)+(nx+1)*(p+1)*nz+(nx+1)*p*nz+nx*p*(nz+1)
    z = nx*(ny+1)*(p+1)+(nx+1)*ny*(p+1)+(nx+1)*ny*p+nx*(ny+1)*p
    return 2*(x+y+z)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    a = parser.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    previous = ROOT/'artifacts/research_checks/2026-09-24_airborne_3d_budget.json'
    rows = json.loads(previous.read_text())['rows']
    for r in rows:
        n = [round(d/r['grid_m']) for d in r['domain_xyz_m']]
        p = r['pml_cells']
        phi = 8*pml_values(*n,p)
        # 4 electric and 4 magnetic vectors per slab, conservative E length p+1.
        coeff = 6*8*(4*(p+1)+4*p)
        known = r['core_array_bytes']+phi+coeff
        r.update(pml_phi_bytes=phi, pml_coefficient_upper_bytes=coeff,
                 enumerated_arrays_GiB=known/2**30,
                 planning_allowance_GiB=(1.5*known+2*2**30)/2**30,
                 allowance_note='50% plus 2GiB contingency; chosen allowance, not proven peak bound')
    # Independent cubic reduction; no dense allocation.
    for n,p in ((4,1),(480,20)):
        assert pml_values(n,n,n,p)==6*(4*p+2)*n*(n+1)
    files=[]
    for case in ('M00','M01'):
        for axis in ('x','y'):
            name=f'{case}_{axis}_3d'
            text=(f'#title: REVIEW ONLY airborne 3D {name}\n'
                  '#domain_mode: 3D\n#domain: 24 24 24\n'
                  '#dx_dy_dz: 0.05 0.05 0.05\n#time_window: 400e-9\n'
                  '#omp_threads: 8\n#pml_cells: 20 20 20 20 20 20\n'
                  '#pml_formulation: HORIPML\n#waveform: impulse 1 1 impulse\n'
                  f'#hertzian_dipole: {axis} 12 11.3 19 impulse\n'
                  f'#rx: 12 12.6 19 measurement E{axis}\n')
            if case=='M01':
                text+='#material: 9 0 1 0 ideal_halfspace\n#box: 0 0 0 24 24 4 ideal_halfspace\n'
            path=a.output/f'{name}.in'
            path.write_text(text,encoding='utf-8',newline='\n')
            files.append(dict(id=name,path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                              status='review_only_not_parser_validated'))
    sources=['gprMax/grid/fdtd_grid.py','gprMax/pml.py','build_config.py','setup.py','pyproject.toml']
    result=dict(schema='airborne-first-run-preparation/1',packet_id='V4-AIR3D-DRAFT-01',
        approved_to_simulate=False,solver_executed=False,runtime_verified=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parent_budget_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
        official_source_sha256={s:hashlib.sha256((a.source_root/s).read_bytes()).hexdigest() for s in sources},
        pml_assumptions='six full-face external slabs, one default CFS term, double, no internal PML',
        exclusions=['unknown transient build peak','allocator overhead beyond chosen allowance',
                    'dispersion','GPU duplicate buffers','full-volume output'],
        rows=rows,inputs=files,
        proposed_initial_run_ids=['M00_x_3d'],
        later_unapproved_run_ids=['M01_x_3d','M00_y_3d','M01_y_3d'],
        proposed_initial_limits={'threads':8,'backend':'CPU','precision':'double',
            'wall_minutes':30,'process_tree_rss_GiB':24,'minimum_available_RAM_GiB':32,
            'output_GiB':1,'retries':0},
        limits_enforced=False,
        stop_rules=['wall/RSS/output cap','nonzero exit','missing required source/receiver metadata',
                    'nonfinite or all-zero response','wrong version/precision/geometry'],
        continue_rule='no automatic progression to later IDs; inspect first result and agree scope',
        environment_readonly_findings={'python_candidate':'bundled Python 3.12.14',
            'VS_BuildTools':'2022','MSVC_directory_version':'14.44.35207',
            'gprMax_V4_installed_in_isolated_env':False,
            'compiler_invocation_tested':False},
        unresolved=['isolated V4 build and import identity','runtime resource supervisor',
                    'actual Yee source/receiver positions','grid convergence beyond this base'],
        validation={'cubic_PML_formula_two_sizes':True})
    (a.output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps([{'id':r['id'],'arrays_GiB':round(r['enumerated_arrays_GiB'],2),
                      'planning_GiB':round(r['planning_allowance_GiB'],2)} for r in rows],indent=2))


if __name__=='__main__':
    main()
