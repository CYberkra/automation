"""Freeze a centre pilot or the remaining stations of joint X/Z refinement."""
import argparse
import json
from pathlib import Path

from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/research_checks/2026-10-03_hs4t2d_boundary_controls'


def create_input(source, first_index, traces):
    result = []
    offset = .05 * (first_index - 1)
    for line in source.read_text('utf-8').splitlines():
        words = line.split()
        if line.startswith(('#src_steps:', '#rx_steps:')):
            continue
        if line.startswith('#dx_dy_dz:'):
            line = '#dx_dy_dz: 0.025 0.05 0.025'
        elif line.startswith('#pml_cells:'):
            line = '#pml_cells: 80 0 40 80 0 40'
        elif line.startswith(('#hertzian_dipole:', '#rx:')):
            start = 2 if words[0] == '#hertzian_dipole:' else 1
            words[start] = f'{float(words[start]) + offset:.12g}'
            line = ' '.join(words)
        elif line.startswith('#geometry_view:'):
            line = '#geometry_view: 0 0 0 36 0.05 33 0.025 0.05 0.025 hs4t2d_geom n'
        result.append(line)
    if traces > 1:
        result.extend(['#src_steps: 0.5 0 0', '#rx_steps: 0.5 0 0'])
    return '\n'.join(result) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['centre', 'remaining'], required=True)
    parser.add_argument('--centre', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new study required; consumed attempts cannot be reused')
    prerequisite = None
    if args.stage == 'remaining':
        if args.centre is None:
            raise ValueError('completed centre pilot required')
        pilot = json.loads((args.centre / 'completed_verification.json').read_text('utf-8'))
        if pilot['status'] != 'PASS' or pilot['contract_sha256'] != sha256(args.centre / 'execution_contract.json'):
            raise ValueError('pilot not independently verified')
        prerequisite = {'directory': str(args.centre.resolve()),
                        'contract_sha256': sha256(args.centre / 'execution_contract.json'),
                        'completed_verification_sha256': sha256(args.centre / 'completed_verification.json')}
    old = json.loads((BASE / 'execution_contract.json').read_text('utf-8'))
    specs = [('centre', [61])] if args.stage == 'centre' else [
        ('left', list(range(1, 52, 10))), ('right', list(range(71, 122, 10)))]
    args.out.mkdir(parents=True)
    groups = []
    for segment, indices in specs:
        for role in ('rough', 'halfspace'):
            source = BASE / ('x40wide_' + role) / 'profile.in'
            name = segment + '_' + role
            folder = args.out / name
            folder.mkdir()
            path = folder / 'profile.in'
            path.write_text(create_input(source, indices[0], len(indices)), encoding='utf-8')
            groups.append({'id': name, 'input': str(path.resolve()), 'sha256': sha256(path),
                           'source_input': str(source.relative_to(ROOT)), 'source_sha256': sha256(source),
                           'baseline_group': 'x40wide_' + role, 'halfspace': role == 'halfspace',
                           'traces': len(indices), 'old_station_indices': indices,
                           'spacing_m': [.025, .05, .025], 'pml_cells': [80, 0, 40, 80, 0, 40],
                           'station_x_m': [14.6 + .05*(i-1) for i in indices],
                           'minimum_free_VRAM_GiB': 3.5})
    code = ['scripts/prepare_hs4t2d_joint_grid.py', 'scripts/run_hs4t2d_joint_grid.py',
            'scripts/run_hs4t2d_joint_grid.cmd', 'scripts/check_hs4t2d_joint_grid.py',
            'scripts/gprmax_cached_cuda_entry.py']
    contract = {'status': 'FROZEN_APPROVED', 'stage': args.stage,
                'approval_basis': 'User: 做吧, approving joint X/Z refinement at actual15m height and subsequent spatial investigation. Frozen before solver attempt on2026-10-04.',
                'purpose': '5cm ->2.5cm X/Z sensitivity with exactly preserved staircase/material/source locations and physical PML widths. No full-convergence or physical event-identity claim.',
                'backend': 'CUDA', 'gpu_device': 0, 'precision': 'double', 'python': old['python'],
                'max_runs': sum(g['traces'] for g in groups),
                'max_wall_s': 900 if args.stage == 'centre' else 5400,
                'min_available_RAM_GiB': 1.5, 'groups': groups,
                'source_identities': old['source_identities'], 'prerequisite': prerequisite,
                'baseline_directory': str(BASE.resolve()), 'baseline_contract_sha256': sha256(BASE / 'execution_contract.json'),
                'invariants': '36x.05x33m; original .8m full-profile relief; Z12 ground;15m height;1.3m offset; invariant Y .05m line-source scale; impulse I1;600ns; HORIPML/default singleCFS. Physical PML sides2m, Zfaces1m.',
                'numerical_factor': 'X/Z half cell size, dt follows actual2D CFL, PML cell counts doubled to hold widths; grid sampling of fixed source changes with dt and must be normalized from saved complex source spectrum.',
                'expected_geometry': 'Fine voxel material map must equal coarse map repeated2x in X and Z. Original 5cm-aligned staircase is preserved, no smoother surrogate surface.',
                'no_duplicate_attempts': 'Centre61 executed only in pilot; remaining1..51/71..121 assembled with pilot. Existing coarse13 reused without solver replay.',
                'processing': 'actual saved source-normalized complex501 frequencies20-170MHz; Hann/8x/200ns tail; rough-halfspace complex subtraction before reconstruction; fixed85-120 and160-220ns windows; noSVD/AGC/per-trace scaling.',
                'scope': 'Development diagnostic, not3D/field/clean truth or training qualification; two grids do not certify an absolute numerical error bound.',
                'code_identities': {str(ROOT / name): sha256(ROOT / name) for name in code}}
    (args.out / 'execution_contract.json').write_text(json.dumps(contract, indent=2) + '\n', encoding='utf-8')
    print('Frozen', args.stage, contract['max_runs'], 'traces')


if __name__ == '__main__':
    main()
