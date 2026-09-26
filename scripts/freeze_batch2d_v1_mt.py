"""Freeze the batch2d_v1 multitrace (MT) segment: 22 cases x 33 receivers (proposal 2026-09-26 §9).

Generates .in inputs (mother model + 33 explicit #rx lines, mt01..mt33, Ex only), cases/groups
contracts with seeds per spec §4.1, static_check.json (spec §5.3 step 2), budget.json (spec §3.1
formulas), the runner clone, and rewrites the global execution gate with approved_to_simulate=false
(step 4 of the freeze flow). Execution still requires the user's explicit agreement (step 5-6).

Read-only over archived contracts; no solver is launched here.
"""

import hashlib
import json
import math
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'configs/research/batch2d_v1'
OUT = ROOT / 'configs/research/batch2d_v1_mt'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
N_TRACES = 33
DY_RX = 0.25
Y_FIRST = 12.65
RX_Z = 45.0
ANCHOR_INDEX = 16  # 0-based; mt17 == legacy single-trace Rx (y16.65)
TIME_WINDOW_S = 1200e-9
VARIANT = 'multitrace_33rx_0p25m'
BATCH_ID = 'batch2d_v1'
SEGMENT = 'mt'

# Proposal §3.2 execution order: family BG precedes its targets; C3 single-factor block last.
MT_ORDER = [
    'B2D-C1m-BG', 'B2D-C1m-D10m-W4m-T0.5m-E20-S0.02',
    'B2D-C3m-BG', 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02',
    'B2D-C5m-BG', 'B2D-C5m-D10m-W4m-T0.5m-E20-S0.02',
    'B2D-C8m-BG', 'B2D-C8m-D10m-W4m-T0.5m-E20-S0.02',
    'B2D-C1m-D10m-NC', 'B2D-C3m-D10m-NC', 'B2D-C5m-D10m-NC', 'B2D-C8m-D10m-NC',
    'B2D-C3m-D10m-W1m-T0.5m-E20-S0.02', 'B2D-C3m-D10m-W2m-T0.5m-E20-S0.02',
    'B2D-C3m-D10m-W4m-T0.25m-E20-S0.02', 'B2D-C3m-D10m-W4m-T1m-E20-S0.02',
    'B2D-C3m-D10m-W4m-T0.5m-E12-S0.005', 'B2D-C3m-D10m-W4m-T0.5m-E28-S0.05',
    'B2D-C3m-D10m-W4m-T0.5m-E20-S0.005', 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.05',
    'B2D-C3m-D20m-W4m-T0.5m-E20-S0.02', 'B2D-C3m-D5m-W4m-T0.5m-E20-S0.02',
]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(group_id: str, case_id: str, variant: str) -> int:
    txt = f'{BATCH_ID}|{group_id}|{case_id}|{variant}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def build_mt_input(mother_text: str, run_id: str) -> str:
    lines = mother_text.splitlines()
    rx_idx = [i for i, l in enumerate(lines) if l.startswith('#rx:')]
    assert len(rx_idx) == 1, run_id
    rx_i = rx_idx[0]
    title_idx = next(i for i, l in enumerate(lines) if l.startswith('#title:'))
    new_lines = list(lines)
    new_lines[title_idx] = f'#title: {run_id} batch2d_v1 MT 33rx variable-offset receiver array'
    rx_block = [f'#rx: inf {Y_FIRST + k * DY_RX:.2f} {RX_Z} mt{k + 1:02d} Ex' for k in range(N_TRACES)]
    new_lines[rx_i:rx_i + 1] = rx_block
    return '\n'.join(new_lines) + '\n'


def main():
    assert len(MT_ORDER) == 22 and len(set(MT_ORDER)) == 22
    OUT.mkdir(parents=True, exist_ok=True)
    base_cases = json.loads((BASE / 'cases.json').read_text(encoding='utf-8'))
    base_groups = json.loads((BASE / 'groups.json').read_text(encoding='utf-8'))
    by_id = {c['run_id']: c for c in base_cases['cases']}
    missing = [r for r in MT_ORDER if r not in by_id]
    assert not missing, missing

    # ---------- static geometry self-checks (spec §5.3 step 2) ----------
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    static = {
        'schema': 'batch2d_v1_mt_static_check/1',
        'checked_date': '2026-09-26',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT,
                 'dt_formula_over_contract': dt_formula / DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'time_steps_contract': 20352,
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'receiver_array': {
            'n_traces': N_TRACES, 'spacing_m': DY_RX, 'first_y_m': Y_FIRST,
            'last_y_m': Y_FIRST + (N_TRACES - 1) * DY_RX, 'z_m': RX_Z,
            'anchor_index': ANCHOR_INDEX, 'anchor_y_m': Y_FIRST + ANCHOR_INDEX * DY_RX,
            'trace_ids': [f'mt{k + 1:02d}' for k in range(N_TRACES)],
        },
        'trace_checks': [], 'materials_unchanged_vs_mother': [], 'source_unchanged_vs_mother': [],
    }
    assert cells == 2560000
    assert abs(main_fields_gib - 0.3437627702951431) < 1e-9
    assert static['receiver_array']['last_y_m'] == 20.65
    assert static['receiver_array']['anchor_y_m'] == 16.65

    cases_out, groups_out, contracts = [], [], []
    runner_src = ROOT / 'scripts/run_approved_batch2d.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_mt.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1'", f"BATCH = '{BATCH_ID}_mt'"),
                          encoding='utf-8', newline='\n')
    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    for mother_id in MT_ORDER:
        m = by_id[mother_id]
        run_id = mother_id + '-MT33'
        mother_text = (BASE / m['file']).read_text(encoding='utf-8')
        mt_text = build_mt_input(mother_text, run_id)
        in_path = OUT / (run_id + '.in')
        in_path.write_text(mt_text, encoding='utf-8', newline='\n')
        in_sha = sha256_bytes(in_path.read_bytes())

        # per-trace static checks
        y_cells, ids = [], []
        for k in range(N_TRACES):
            y = Y_FIRST + k * DY_RX
            yc, zc = y / DY, RX_Z / DZ
            assert abs(yc - round(yc)) < 1e-9 and abs(zc - round(zc)) < 1e-9
            assert 1.0 < y < 31.0, f'rx inside PML: {y}'
            y_cells.append(int(round(yc)))
            ids.append(f'mt{k + 1:02d}')
        assert y_cells == [506 + 10 * k for k in range(N_TRACES)]
        assert y_cells[ANCHOR_INDEX] == 666  # identical cell to the archived single-trace Rx
        assert len(set(ids)) == N_TRACES
        rx_lines = [l for l in mt_text.splitlines() if l.startswith('#rx:')]
        assert len(rx_lines) == N_TRACES and all(l.endswith(' Ex') for l in rx_lines)
        # everything except title and rx block byte-identical to the mother input
        strip = lambda s: [l for l in s.splitlines() if not l.startswith('#rx:') and not l.startswith('#title:')]
        assert strip(mother_text) == strip(mt_text), run_id
        static['trace_checks'].append({'run_id': run_id, 'n_rx_lines': len(rx_lines),
                                       'y_cells_first_last': [y_cells[0], y_cells[-1]],
                                       'all_ex_only': True, 'ids_unique': True,
                                       'outside_pml': True, 'anchor_cell_666': True})
        static['materials_unchanged_vs_mother'].append(run_id)
        static['source_unchanged_vs_mother'].append(run_id)

        entry = dict(m)
        entry.update(run_id=run_id, file=run_id + '.in', input_sha256=in_sha,
                     segment=SEGMENT, variant_tag=VARIANT,
                     seed=seed_for(m['group_id'], m['case_id'], VARIANT),
                     multitrace={
                         'n_traces': N_TRACES, 'trace_spacing_m': DY_RX,
                         'trace_span_m': (N_TRACES - 1) * DY_RX,
                         'trace_axis_y_range_m': [Y_FIRST, Y_FIRST + (N_TRACES - 1) * DY_RX],
                         'trace_receiver_z_m': RX_Z, 'anchor_trace_index': ANCHOR_INDEX,
                         'anchor_trace_y_m': 16.65, 'trace_ids': ids,
                         'trace_order': 'construction',
                         'acquisition_geometry': 'single-shot variable-offset receiver array along y (offset 2.70-5.30 m)',
                         'rx_output_components': ['Ex'],
                     })
        cases_out.append(entry)
        groups_out.append({'case_id': m['case_id'], 'group_id': m['group_id'],
                           'run_id': run_id,
                           'mother_model_hash': next(g['mother_model_hash'] for g in base_groups['groups']
                                                     if g['group_id'] == m['group_id']),
                           'seed': entry['seed'], 'variant_tag': VARIANT})
        contracts.append({
            'packet_id': 'BATCH2D-V1-MT',
            'input_path': f'configs/research/batch2d_v1_mt/{run_id}.in',
            'input_sha256': in_sha,
            'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
            'runtime_identity_sha256': runtime_sha,
            'launcher_sha256': launcher_sha,
            'supervisor_sha256': supervisor_sha,
            'attempt_record': f'artifacts/research_checks/2026-09-26_{run_id}_attempt.json',
            'run_directory': f'artifacts/simulations/2026-09-26_{run_id}',
            'continuation': 'Serial in frozen order; family BG precedes its targets; any failure aborts batch, later attempts untouched. No retries.',
            'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
            'cuda_identity_sha256': cuda_sha,
            'run_id': run_id,
        })

    # ---------- budget.json (spec §3.1 formulas + proposal §4 increments) ----------
    yfine_ref = {'wall_s': 1088.781, 'cells': 40960000, 'dt_s': 1.474e-11}
    wall_scaled = yfine_ref['wall_s'] * (cells / yfine_ref['cells']) * (yfine_ref['dt_s'] / DT_CONTRACT)
    budget = {
        'schema': 'batch2d_v1_mt_budget/1',
        'date': '2026-09-26',
        'formulas': 'spec §3.1: cells=(32/dy)*(50/dz); dt=1/(c*sqrt(dy^-2+dz^-2)); '
                    'main_fields_ID_GiB=72B*2*(ny+1)(nz+1)/2^30; wall=ref*(cells/cells_ref)*(dt_ref/dt)',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'dt_formula_s': dt_formula, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'wall_scaled_from_YFINE_s': wall_scaled,
                        'measured_base_batch_mean_s': 48.3,
                        'measured_base_batch_range_s_excl_first': [29.859, 59.218],
                        'planning_envelope_min': 2.5},
        'multitrace_increment': {
            'output_per_case_bytes': N_TRACES * 20352 * 8,
            'output_per_case_MiB': N_TRACES * 20352 * 8 / 2 ** 20,
            'device_rx_array_bytes': 6 * 20352 * N_TRACES * 8,
            'device_rx_array_MiB': 6 * 20352 * N_TRACES * 8 / 2 ** 20,
            'job_peak_increment_GiB': 0.034,
            'wall_increment': 'expected negligible (unmeasured); not used to relax budget',
        },
        'batch_totals': {
            'n_cases': 22,
            'planning_envelope_s': 22 * 150,
            'measured_mean_estimate_s': 22 * 48.3,
            'conservative_budget_min': 65,
        },
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 22},
        'known_limitations': [
            'working scaling is not a guaranteed runtime',
            'main-field array excludes PML, update coefficients, CUDA context and host temporaries',
            'Job committed memory is not GPU VRAM',
            'estimates are for ranking and cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-26_multitrace_batch_proposal.md (user-approved 2026-09-26) + '
                         'docs/research/2026-09-26_batch_2d_spec_v1.md',
                 'n_cases': 22, 'n_exceptions': 0,
                 'derivation': 'mother_model_id identical to batch2d_v1 base cases; receiver-array derivative only',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                  'convention': 'spec §4.2: any split is by group_id only, never by case_id, window or trace',
                  'groups': groups_out}
    (OUT / 'cases.json').write_text(json.dumps(cases_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'groups.json').write_text(json.dumps(groups_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'budget.json').write_text(json.dumps(budget, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'static_check.json').write_text(json.dumps(static, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')

    # ---------- gate rewrite (spec §5.3 step 4): approved_to_simulate=false ----------
    old_gate = json.loads(GATE.read_text(encoding='utf-8'))
    new_gate = {
        'schema': old_gate['schema'],
        'updated_date': '2026-09-26',
        'authority': old_gate['authority'],
        'user_instruction': old_gate['user_instruction'],
        'target_solver': old_gate['target_solver'],
        'source_root_on_current_machine': old_gate['source_root_on_current_machine'],
        'batch_id': f'{BATCH_ID}_mt',
        'approved_to_simulate': False,
        # check_handoff safety rule requires approved_run_ids == [] while the flag is false
        # (spec §5.2 cites scripts/check_handoff.py:53); the frozen order is kept in
        # approved_run_ids_pending_agreement and copied into approved_run_ids on agreement.
        'approved_run_ids': [],
        'approved_run_ids_pending_agreement': [c['run_id'] for c in contracts],
        'approved_compute_budget': {
            'wall_minutes': 20, 'job_commit_GiB': 4, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': 22, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'multitrace batch proposal v0.1 (docs/research/2026-09-26_multitrace_batch_proposal.md), user-approved 2026-09-26: lifts G2 data-shape block so the 27 operator candidates and baselines L/R can run on [sample, trace] inputs',
            'eval runner design v0.1 (docs/research/2026-09-26_eval_runner_and_selector_design_v0.1.md) §A.3 gate 2: 27/27 candidates rejected on single-trace archived data; G2 is the total downstream blocker',
            'cost basis: batch2d_v1 coarse measured 29.9-59.2 s/case excl. first, job peak ~2.72 GiB; 33-rx increment ~0.034 GiB and negligible wall increment (proposal §4.1); batch cap 65 min',
            'input package verified 2026-09-26: 22 .in derived from frozen base mothers with byte-identical geometry/materials/source; 33 rx lines per case, Ex only, ids mt01-mt33, all y integer cells (506+10k), anchor mt17 == archived single-trace Rx cell 666; static_check.json + budget.json on disk',
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
        'note': 'batch2d_v1_mt frozen 2026-09-26 (22 cases x 33 rx, BASE, single-shot variable-offset array, proposal user-approved). approved_to_simulate stays false until the user agrees to execute; G2 lifted in data-shape sense only; G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: 22 .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: batch2d_v1_mt, approved_to_simulate=false, 22 contracts')
    print('launcher_sha256', launcher_sha)
    print('anchor mt17 cell', y_cells[ANCHOR_INDEX], '; dt_formula/contract', dt_formula / DT_CONTRACT)


if __name__ == '__main__':
    main()
