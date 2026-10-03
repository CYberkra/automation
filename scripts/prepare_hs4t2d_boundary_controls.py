"""Freeze geometry-preserving side/top/bottom PML controls before FDTD."""
import argparse
import json
from pathlib import Path

from hs_capsule_identity import sha256
from prepare_hs4t2d_cause_controls import ROOT, DESIGN, create_input


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new study required')
    base = json.loads((DESIGN / 'execution_contract_cuda_cached_att4.json').read_text('utf-8'))
    specs = [('x40wide', [40, 0, 20, 40, 0, 20], 13),
             ('x60wide', [60, 0, 20, 60, 0, 20], 1),
             ('bottom40wide', [20, 0, 40, 20, 0, 20], 1),
             ('top40wide', [20, 0, 20, 20, 0, 40], 1)]
    args.out.mkdir(parents=True)
    groups = []
    for factor, pml, traces in specs:
        for role in ('rough', 'halfspace'):
            source = DESIGN / 'relief2p0' / ('hs4t2d_t61.in' if traces == 1 else 'hs4t2d_t01.in')
            text = create_input(source, 36, 15, .05, role == 'halfspace', traces)
            before = '#pml_cells: 20 0 20 20 0 20'
            if text.count(before) != 1:
                raise ValueError('unexpected baseline boundary')
            text = text.replace(before, '#pml_cells: ' + ' '.join(map(str, pml)))
            name = factor + '_' + role
            folder = args.out / name
            folder.mkdir()
            target = folder / 'profile.in'
            target.write_text(text, encoding='utf-8')
            groups.append({'id': name, 'input': str(target.resolve()), 'sha256': sha256(target),
                           'source_input': str(source.relative_to(ROOT)), 'source_sha256': sha256(source),
                           'width_m': 36, 'x_translation_m': 12, 'height_m': 15, 'dz_m': .05,
                           'halfspace': role == 'halfspace', 'traces': traces, 'pml_cells': pml,
                           'old_station_indices': [61] if traces == 1 else list(range(1, 122, 10)),
                           'minimum_free_VRAM_GiB': 3.5})
    contract = {'status': 'FROZEN_APPROVED',
                'approval_basis': 'User active goal: continue, read official manual and identify cause; existing autonomous GPU authorization. This bounded contract is frozen before any attempt.',
                'purpose': 'Change only native PML thickness, with identical material geometry/source/receiver/grid. Compare 13 stations and centre; no physical PASS or training labels.',
                'backend': 'CUDA', 'gpu_device': 0, 'precision': 'double', 'python': base['python'],
                'max_runs': 32, 'max_wall_s': 1800, 'min_available_RAM_GiB': 1,
                'groups': groups, 'source_identities': base['source_identities'],
                'baseline': 'artifacts/research_checks/2026-10-03_hs4t2d_cause_controls/wide36_{rough,halfspace}',
                'invariants': '36 x .05 x 33 m; .05 m grid; ground Z12; 15 m standoff; offset1.3; impulse I1; same endpoint-padded material map; HORIPML with default single CFS term; 600ns raw time.',
                'factors': 'x40: both X faces20->40 cells, 13 pairs; x60: X20->60 centre pair; bottom40: Z0 only20->40 centre pair; top40: Zmax only20->40 centre pair.',
                'clearance': 'X60 PML inner faces3/33m; even all transect stations TX14.6..20.6 RX15.9..21.9 would remain outside PML by >=11.1m; thick top PML inner Z31, sourceZ27 clearance4m.',
                'processing': 'exact saved source-normalized complex501 spectrum20-170MHz; Hann/8x/200ns tail; complex rough-halfspace subtraction; additionally raw full-record and fixed windows85-120/160-220/300-450ns.',
                'claims': 'Boundary sensitivity only, finite comparisons; thickness invariance is not a full boundary-error bound; 5cm grid is still coarse.',
                'code_identities': {str(ROOT / name): sha256(ROOT / name) for name in
                    ['scripts/prepare_hs4t2d_boundary_controls.py', 'scripts/prepare_hs4t2d_cause_controls.py',
                     'scripts/run_hs4t2d_cause_controls.py', 'scripts/run_hs4t2d_cause_controls.cmd',
                     'scripts/gprmax_cached_cuda_entry.py']}}
    (args.out / 'execution_contract.json').write_text(json.dumps(contract, indent=2) + '\n', encoding='utf-8')
    print('Frozen', len(groups), 'groups /', sum(g['traces'] for g in groups), 'traces')


if __name__ == '__main__':
    main()
