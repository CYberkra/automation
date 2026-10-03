"""Freeze a 1/60m centre pair to test contraction after the 2.5cm pilot."""
import argparse
import json
from pathlib import Path

from prepare_hs4t2d_joint_grid import ROOT, BASE, create_input
from hs_capsule_identity import sha256


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--centre', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise ValueError('new study required')
    pilot_path = args.centre / 'execution_contract.json'
    pilot = json.loads(pilot_path.read_text('utf-8'))
    check = json.loads((args.centre/'completed_verification.json').read_text('utf-8'))
    if check['status'] != 'PASS' or check['contract_sha256'] != sha256(pilot_path) or pilot['stage'] != 'centre':
        raise ValueError('completed verified2.5cm centre required')
    args.out.mkdir(parents=True)
    groups = []
    for role in ('rough', 'halfspace'):
        source = BASE / ('x40wide_'+role) / 'profile.in'
        text = create_input(source, 61, 1)
        text = text.replace('#dx_dy_dz: 0.025 0.05 0.025', '#dx_dy_dz: 0.016666666666666666 0.05 0.016666666666666666')
        text = text.replace('#pml_cells: 80 0 40 80 0 40', '#pml_cells: 120 0 60 120 0 60')
        text = text.replace('#geometry_view: 0 0 0 36 0.05 33 0.025 0.05 0.025 hs4t2d_geom n',
                            '#geometry_view: 0 0 0 36 0.05 33 0.016666666666666666 0.05 0.016666666666666666 hs4t2d_geom n')
        name = 'centre_'+role
        folder = args.out/name
        folder.mkdir()
        path = folder/'profile.in'
        path.write_text(text, encoding='utf-8')
        groups.append({'id':name, 'input':str(path.resolve()), 'sha256':sha256(path),
                       'source_input':str(source.relative_to(ROOT)), 'source_sha256':sha256(source),
                       'baseline_group':'x40wide_'+role, 'halfspace':role=='halfspace',
                       'traces':1, 'old_station_indices':[61], 'station_x_m':[17.6],
                       'spacing_m':[1/60,.05,1/60], 'pml_cells':[120,0,60,120,0,60],
                       'minimum_free_VRAM_GiB':4.2})
    code = ['scripts/prepare_hs4t2d_joint_grid_third.py', 'scripts/prepare_hs4t2d_joint_grid.py',
            'scripts/run_hs4t2d_joint_grid_third.py', 'scripts/run_hs4t2d_joint_grid_third.cmd',
            'scripts/check_hs4t2d_joint_grid_third.py', 'scripts/gprmax_cached_cuda_entry.py']
    contract = dict(pilot)
    contract.update({'resolution_level':'third_1over60m', 'groups':groups, 'max_runs':2, 'max_wall_s':1800,
        'purpose':'Centre contraction diagnostic across5cm/2.5cm/1over60m at identical15m physical model and boundary widths.',
        'min_available_RAM_GiB':2.0,
        'approval_basis':'User 做吧 approving joint-grid checking; follow-up contraction diagnostic after the completed2.5cm scan. Frozen before any third-grid attempt.',
        'prerequisite':{'directory':str(args.centre.resolve()), 'contract_sha256':sha256(pilot_path),
                        'completed_verification_sha256':sha256(args.centre/'completed_verification.json')},
        'numerical_factor':'X/Z 1/60m centre only, grid2160x1x1980, dx/dz5cm divided3. dy=.05 unchanged; sides2m/Z1m PML held; actual2D CFL dt changes.',
        'expected_geometry':'Third-resolution voxel map must equal original5cm map repeated3x along X/Z. No smoothing or physical geometry/source shift.',
        'no_duplicate_attempts':'New grid only; coarse and2.5cm centre outputs reused. No1.25cm run (estimated above currently free device memory at36m width).',
        'scope':'Adjacent-grid contraction at centre only. No full13-station third-grid convergence or absolute numerical error bound.',
        'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in code}})
    (args.out/'execution_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    print('Frozen third-resolution centre pair')


if __name__ == '__main__':
    main()
