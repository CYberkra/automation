"""Freeze the PML sensitivity segment: S2 BG anchor trace t17, PML thickness sweep.

Basis: user 2026-09-28 — "可以，顺便调研确立好我们的色散参数" (the "可以" answers the
2026-09-28 proposal to schedule a PML sensitivity experiment after the CO B-scan
display work found the interface echo only 3-8 dB above a trace-invariant floor).

Question: is the broadband floor (median envelope ~-65 dB vs direct, trace-invariant)
dominated by PML boundary reflections, or is it physical/geometric?

Design (all runs reuse the mother anchor geometry Tx y=15.35 / Rx y=16.65, z=45,
i.e. exactly the archived CO33-t17 pair; the archived t17 run is the pml40 baseline
and is NOT rerun):

- B2D-C3mS2-BG-PML20-T17: mother text, #pml_cells 0 40 40 0 40 40 -> 0 20 20 0 20 20
- B2D-C3mS2-BG-PML80-T17: mother text, #pml_cells -> 0 80 80 0 80 80 (2.0 m)
- B2D-FS-PML40-T17: free-space control; #material:/#box: lines removed, pml40 kept.
  Floor level here isolates PML+numerical floor with no subsurface at all.

Readout (post-run, separate analysis): median reconstructed envelope in 130-170 ns
and 240-280 ns relative to trace max, plus the interface-window peak. Ordering
pml20 > pml40 > pml80 of the floor confirms PML origin quantitatively; FS floor
near the BG floor confirms non-physical (non-subsurface) origin.

Mirrors scripts/freeze_batch2d_slope_t2_co.py structure; --verify-only rebuilds in
memory and asserts byte equality with disk without touching the gate.
Read-only over archived contracts; no solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SLOPE_BASE = ROOT / 'configs/research/batch2d_slope_t2'
OUT = ROOT / 'configs/research/pml_sensitivity_s2bg'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
TIME_WINDOW_S = 1200e-9
BATCH_ID = 'pml_sensitivity_s2bg_t17'
RUN_DATE = '2026-09-28'
MOTHER = 'B2D-C3mS2-BG'
PML_LINE = '#pml_cells: 0 40 40 0 40 40'
BASELINE_RUN = '2026-09-27_B2D-C3mS2-BG-CO33-t17'  # archived pml40 anchor, not rerun

VARIANTS = [
    # (run_id, pml replacement or None, drop_materials, note)
    # run_id PML20b: the first PML20 attempt (2026-09-28) was consumed by a
    # launch-infrastructure failure (non-ASCII em-dash in the generated #title line
    # crashed gprMax's GBK input reader before any solving); per the no-rerun rule a
    # corrected input runs under a new run_id. The failed attempt record is kept.
    ('B2D-C3mS2-BG-PML20b-T17', '#pml_cells: 0 20 20 0 20 20', False,
     '0.5 m PML (20 cells), HORIPML + CFS unchanged'),
    ('B2D-C3mS2-BG-PML80-T17', '#pml_cells: 0 80 80 0 80 80', False,
     '2.0 m PML (80 cells), HORIPML + CFS unchanged'),
    ('B2D-FS-PML40-T17', None, True,
     'free-space control: materials/boxes removed, pml40 kept'),
]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(run_id: str) -> int:
    txt = f'{BATCH_ID}|{run_id}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def build_input(mother_text: str, run_id: str, pml_line, drop_materials: bool, note: str) -> str:
    lines = mother_text.splitlines()
    title_idx = [i for i, l in enumerate(lines) if l.startswith('#title:')]
    assert len(title_idx) == 1
    assert '#hertzian_dipole: x inf 15.35 45 impulse' in lines
    assert '#rx: inf 16.65 45 measurement Ex' in lines
    new_lines = []
    for l in lines:
        if l.startswith('#title:'):
            new_lines.append(f'#title: {run_id} {BATCH_ID} PML sensitivity (mother {MOTHER}) - {note}')
        elif pml_line is not None and l == PML_LINE:
            new_lines.append(pml_line)
        elif drop_materials and (l.startswith('#material:') or l.startswith('#box:')):
            continue
        else:
            new_lines.append(l)
    text = '\n'.join(new_lines) + '\n'
    if pml_line is not None:
        assert pml_line in text and PML_LINE not in text
    if drop_materials:
        assert '#material:' not in text and '#box:' not in text and PML_LINE in text
    text.encode('ascii')  # gprMax reads inputs with the system codec (GBK here); ASCII-only
    return text


def collect():
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    mother_text = (SLOPE_BASE / (MOTHER + '.in')).read_text(encoding='utf-8')
    assert mother_text.count(PML_LINE) == 1

    static = {
        'schema': 'pml_sensitivity_s2bg_static_check/1',
        'checked_date': RUN_DATE,
        'question': 'origin of the trace-invariant broadband floor seen in the CO B-scan '
                    '(median envelope ~-65 dB vs direct): PML reflection or physical/geometric',
        'mother': {'id': MOTHER, 'base': 'configs/research/batch2d_slope_t2',
                   'sha256': sha256_bytes(mother_text.encode('utf-8'))},
        'anchor': {'tx_y_m': 15.35, 'rx_y_m': 16.65, 'z_m': 45.0,
                   'equals': 'archived CO33-t17 pair; archived run is the pml40 baseline',
                   'baseline_run_directory': f'artifacts/simulations/{BASELINE_RUN}'},
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'pml_variants': [{'run_id': r, 'note': n} for r, _, _, n in VARIANTS],
        'pml_cfs_unchanged': '#pml_cfs: constant forward 0 0 constant forward 1 1 '
                             'quartic forward 0 0.21235349838321013',
        'formulation': 'HORIPML (unchanged from mother)',
        'interior_check': 'pml80 leaves interior y in (2,30) m, z in (2,48) m; antennas '
                          '(y 15.35/16.65, z 45), surface (z 30) and slope interface '
                          '(z 25.0-26.65) all remain inside the interior',
        'readout_windows_ns': {'mid': [130, 170], 'deep': [240, 280],
                               'interface_approx': [199, 223]},
    }
    assert cells == 2560000

    inputs, cases_out = {}, []
    for run_id, pml_line, drop_materials, note in VARIANTS:
        text = build_input(mother_text, run_id, pml_line, drop_materials, note)
        inputs[run_id] = text
        cases_out.append({
            'run_id': run_id, 'file': run_id + '.in',
            'input_sha256': sha256_bytes(text.encode('utf-8')),
            'mother_model_id': MOTHER, 'case_id': run_id, 'group_id': BATCH_ID,
            'segment': 'pml_sensitivity', 'variant_tag': note,
            'seed': seed_for(run_id),
            'anchor': {'tx_y_m': 15.35, 'rx_y_m': 16.65, 'z_m': 45.0},
        })

    budget = {
        'schema': 'pml_sensitivity_s2bg_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'measured_co_wall_s': 19.3,
                        'note': 'same grid/window as the archived CO t17 anchor '
                                '(measured 17.8-23.6 s); pml80 adds cells only inside '
                                'the PML region of an unchanged domain'},
        'batch_totals': {'n_cases': len(VARIANTS), 'measured_estimate_s': len(VARIANTS) * 25,
                         'conservative_budget_min': 10},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 1, 'retries': 0, 'max_fdtd_runs': len(VARIANTS)},
        'known_limitations': ['estimates are for cap-setting only, not ETA'],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': 'pml_sensitivity', 'grid_tier': 'BASE',
                 'spec': 'decision_log 2026-09-28 (CO B-scan floor observability) + '
                         'user 2026-09-28 "可以" scheduling the PML sensitivity experiment',
                 'n_cases': len(VARIANTS), 'n_exceptions': 0,
                 'derivation': 'mother text byte-identical except #title, and #pml_cells '
                               '(PML20/PML80) or dropped #material/#box lines (FS control); '
                               'pml40 baseline is the archived CO33-t17 run, not rerun',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'pml_sensitivity', 'grid_tier': 'BASE',
                  'convention': 'single diagnostic group; not part of any dev/validation split',
                  'groups': [{'case_id': BATCH_ID, 'group_id': BATCH_ID,
                              'run_id_pattern': 'B2D-*-T17', 'n_traces': len(VARIANTS)}]}
    return inputs, cases_doc, groups_doc, budget, static


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    inputs, cases_doc, groups_doc, budget, static = collect()

    if args.verify_only:
        for run_id, text in inputs.items():
            assert (OUT / (run_id + '.in')).read_text(encoding='utf-8') == text, run_id
        assert json.loads((OUT / 'cases.json').read_text(encoding='utf-8')) == cases_doc
        assert json.loads((OUT / 'groups.json').read_text(encoding='utf-8')) == groups_doc
        assert json.loads((OUT / 'budget.json').read_text(encoding='utf-8')) == budget
        assert json.loads((OUT / 'static_check.json').read_text(encoding='utf-8')) == static
        print(f'verify-only: all {len(inputs)} PML-sensitivity inputs + contracts byte-identical')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_pml_sensitivity_s2bg.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'PML-SENS-S2BG-T17',
        'run_id': rid,
        'input_path': f'configs/research/pml_sensitivity_s2bg/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (PML20, PML80, FS); any failure aborts the '
                        'batch, later attempts untouched. No retries.',
    } for rid in inputs]

    (OUT / 'cases.json').write_text(json.dumps(cases_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'groups.json').write_text(json.dumps(groups_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'budget.json').write_text(json.dumps(budget, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'static_check.json').write_text(json.dumps(static, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')

    old_gate = json.loads(GATE.read_text(encoding='utf-8'))
    new_gate = {
        'schema': old_gate['schema'],
        'updated_date': RUN_DATE,
        'authority': old_gate['authority'],
        'user_instruction': old_gate['user_instruction'],
        'target_solver': old_gate['target_solver'],
        'source_root_on_current_machine': old_gate['source_root_on_current_machine'],
        'batch_id': BATCH_ID,
        'approved_to_simulate': False,
        'approved_run_ids': [],
        'approved_run_ids_pending_agreement': [c['run_id'] for c in contracts],
        'approved_compute_budget': {
            'wall_minutes': 20, 'job_commit_GiB': 4, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 1, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': len(contracts),
            'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-28: "可以" — schedules the PML sensitivity experiment proposed after '
            'the CO B-scan display found the interface echo only a few dB above a '
            'trace-invariant broadband floor',
            'decision_log 2026-09-28: floor origin (PML vs physical) explicitly marked '
            'uncertified; this segment is the certifying experiment',
            'cost basis: measured 17.8-23.6 s/case on identical grid; 3 cases conservative '
            '10 min total; per-case wall cap 20 min',
        ],
        'permitted_preparation': old_gate['permitted_preparation'],
        'requires_agreement_before_execution': old_gate['requires_agreement_before_execution'],
        'approved_execution_contracts': contracts,
        'execution_policy': old_gate['execution_policy'],
        'execution_outcome': {'status': 'pending_user_agreement'},
        'last_completed_execution_contract': {
            'batch_id': old_gate['batch_id'],
            'execution_outcome': old_gate.get('execution_outcome'),
            'note': old_gate.get('note'),
        },
        'note': 'pml_sensitivity_s2bg_t17 frozen 2026-09-28 (S2 BG anchor t17: PML20, PML80, '
                'free-space control; pml40 baseline archived as CO33-t17, not rerun). '
                'approved_to_simulate stays false until the recorded basis is applied; '
                'G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
