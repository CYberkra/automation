"""Freeze the batch2d_slope_t2 multitrace (MT) segment: 6 slope cases x 33 receivers.

Basis: the base batch2d_slope_t2 batch (single-trace) completed 2026-09-27 and is single-offset
by design; the user asked to see B-scans of the finished batch ("跑完的那批先给我看看bscan",
2026-09-27), which a single trace cannot provide. The MT derivative (same 33-receiver
variable-offset array convention as batch2d_v1_mt, frozen mothers = the 6 slope inputs) is the
minimal way to produce B-scans for the slope tier-T2 variants. Acquisition geometry is unchanged
(scheme A): Tx (15.35, 45), rx array y = 12.65..20.65 step 0.25 at z = 45, Ex only, ids
mt01..mt33, anchor mt17 == the single-trace receiver cell (y=16.65).

Mirrors scripts/freeze_batch2d_v1_mt.py: generates .in inputs (mother + 33 explicit #rx lines),
cases/groups contracts, static_check.json, budget.json, the runner clone, and rewrites the global
execution gate with approved_to_simulate=false. --verify-only rebuilds in memory and asserts byte
equality with disk without touching the gate.

Read-only over archived contracts; no solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'configs/research/batch2d_slope_t2'
OUT = ROOT / 'configs/research/batch2d_slope_t2_mt'
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
ANCHOR_INDEX = 16  # 0-based; mt17 == single-trace Rx (y16.65)
TIME_WINDOW_S = 1200e-9
VARIANT = 'multitrace_33rx_0p25m'
BATCH_ID = 'batch2d_slope_t2'
SEGMENT = 'mt'
RUN_DATE = '2026-09-27'

MT_ORDER = [
    'B2D-C3mS2-BG', 'B2D-C3mS2-NC', 'B2D-C3mS2-D10m-W4m-T0.5m-E20-S0.02',
    'B2D-C3mS2TZ-BG', 'B2D-C3mS2TZ-NC', 'B2D-C3mS2TZ-D10m-W4m-T0.5m-E20-S0.02',
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
    new_lines[title_idx] = f'#title: {run_id} {BATCH_ID} MT 33rx variable-offset receiver array'
    rx_block = [f'#rx: inf {Y_FIRST + k * DY_RX:.2f} {RX_Z} mt{k + 1:02d} Ex' for k in range(N_TRACES)]
    new_lines[rx_i:rx_i + 1] = rx_block
    return '\n'.join(new_lines) + '\n'


def collect():
    """Build all in-memory artifacts. Returns (inputs, cases_doc, groups_doc, budget, static)."""
    base_cases = json.loads((BASE / 'cases.json').read_text(encoding='utf-8'))
    by_id = {c['run_id']: c for c in base_cases['cases']}
    missing = [r for r in MT_ORDER if r not in by_id]
    assert not missing, missing

    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    static = {
        'schema': 'batch2d_slope_t2_mt_static_check/1',
        'checked_date': RUN_DATE,
        'base_batch': 'batch2d_slope_t2 (mothers are the frozen slope inputs; only #rx/#title differ)',
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

    inputs, cases_out, groups_out = {}, [], []
    for mother_id in MT_ORDER:
        m = by_id[mother_id]
        run_id = mother_id + '-MT33'
        mother_text = (BASE / m['file']).read_text(encoding='utf-8')
        mt_text = build_mt_input(mother_text, run_id)
        inputs[run_id] = mt_text

        y_cells, ids = [], []
        for k in range(N_TRACES):
            y = Y_FIRST + k * DY_RX
            yc, zc = y / DY, RX_Z / DZ
            assert abs(yc - round(yc)) < 1e-9 and abs(zc - round(zc)) < 1e-9
            assert 1.0 < y < 31.0, f'rx inside PML: {y}'
            y_cells.append(int(round(yc)))
            ids.append(f'mt{k + 1:02d}')
        assert y_cells == [506 + 10 * k for k in range(N_TRACES)]
        assert y_cells[ANCHOR_INDEX] == 666  # identical cell to the single-trace Rx
        assert len(set(ids)) == N_TRACES
        rx_lines = [l for l in mt_text.splitlines() if l.startswith('#rx:')]
        assert len(rx_lines) == N_TRACES and all(l.endswith(' Ex') for l in rx_lines)
        strip = lambda s: [l for l in s.splitlines() if not l.startswith('#rx:') and not l.startswith('#title:')]
        assert strip(mother_text) == strip(mt_text), run_id
        static['trace_checks'].append({'run_id': run_id, 'n_rx_lines': len(rx_lines),
                                       'y_cells_first_last': [y_cells[0], y_cells[-1]],
                                       'all_ex_only': True, 'ids_unique': True,
                                       'outside_pml': True, 'anchor_cell_666': True})
        static['materials_unchanged_vs_mother'].append(run_id)
        static['source_unchanged_vs_mother'].append(run_id)

        entry = dict(m)
        entry.update(run_id=run_id, file=run_id + '.in',
                     input_sha256=sha256_bytes(mt_text.encode('utf-8')),
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
                           'mother_model_hash': next(g['mother_model_hash'] for g in
                                                     json.loads((BASE / 'groups.json').read_text(encoding='utf-8'))['groups']
                                                     if g['run_id'] == m['run_id']),
                           'seed': entry['seed'], 'variant_tag': VARIANT})

    yfine_ref = {'wall_s': 1088.781, 'cells': 40960000, 'dt_s': 1.474e-11}
    wall_scaled = yfine_ref['wall_s'] * (cells / yfine_ref['cells']) * (yfine_ref['dt_s'] / DT_CONTRACT)
    budget = {
        'schema': 'batch2d_slope_t2_mt_budget/1',
        'date': RUN_DATE,
        'formulas': 'spec 3.1: cells=(32/dy)*(50/dz); dt=1/(c*sqrt(dy^-2+dz^-2)); '
                    'main_fields_ID_GiB=72B*2*(ny+1)(nz+1)/2^30; wall=ref*(cells/cells_ref)*(dt_ref/dt)',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'dt_formula_s': dt_formula, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'wall_scaled_from_YFINE_s': wall_scaled,
                        'measured_slope_base_batch_range_s': [17.812, 23.61],
                        'measured_slope_base_note': 'post vctip-warden walls; pre-warden first case was '
                                                    '911.875 s (vctip.exe lingering in the job object)'},
        'multitrace_increment': {
            'output_per_case_bytes': N_TRACES * 20352 * 8,
            'output_per_case_MiB': N_TRACES * 20352 * 8 / 2 ** 20,
            'device_rx_array_bytes': 6 * 20352 * N_TRACES * 8,
            'device_rx_array_MiB': 6 * 20352 * N_TRACES * 8 / 2 ** 20,
            'job_peak_increment_GiB': 0.034,
            'wall_increment': 'expected negligible (measured on batch2d_v1_mt); not used to relax budget',
        },
        'batch_totals': {'n_cases': 6, 'planning_envelope_s': 6 * 150,
                         'measured_mean_estimate_s': 6 * 60,
                         'conservative_budget_min': 20},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 6},
        'known_limitations': [
            'working scaling is not a guaranteed runtime',
            'main-field array excludes PML, update coefficients, CUDA context and host temporaries',
            'Job committed memory is not GPU VRAM',
            'estimates are for ranking and cap-setting only, not ETA',
            'vctip.exe job-lingering can add ~15 min per case if the out-of-job singleton lapses',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-27_slope_family_design_draft.md + '
                         'docs/research/2026-09-26_multitrace_batch_proposal.md (array convention) + '
                         'docs/research/2026-09-26_batch_2d_spec_v1.md',
                 'n_cases': 6, 'n_exceptions': 0,
                 'derivation': 'mother_model_id = batch2d_slope_t2 frozen slope inputs; receiver-array '
                               'derivative only (title + #rx block), everything else byte-identical',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                  'convention': 'spec 4.2: any split is by group_id only, never by case_id, window or trace; '
                                'each slope variant is its own mother-model family with its own BG/NC',
                  'groups': groups_out}
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
        print('verify-only: all 6 MT inputs + cases/groups/budget/static_check byte-identical')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_slope_t2_mt.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}_mt'"),
                          encoding='utf-8', newline='\n')
    cmd_src = (ROOT / 'scripts/run_batch2d_mt_gpu.cmd').read_text(encoding='utf-8')
    assert 'run_approved_batch2d_mt.py' in cmd_src
    (ROOT / 'scripts/run_batch2d_slope_t2_mt_gpu.cmd').write_text(
        cmd_src.replace('run_approved_batch2d_mt.py', 'run_approved_batch2d_slope_t2_mt.py'),
        encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-SLOPE-T2-MT',
        'run_id': rid,
        'input_path': f'configs/research/batch2d_slope_t2_mt/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (per variant: BG, NC, TGT); any failure aborts the batch, '
                        'later attempts untouched. No retries.',
    } for rid in [m + '-MT33' for m in MT_ORDER]]

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
        'batch_id': f'{BATCH_ID}_mt',
        'approved_to_simulate': False,
        'approved_run_ids': [],
        'approved_run_ids_pending_agreement': [c['run_id'] for c in contracts],
        'approved_compute_budget': {
            'wall_minutes': 20, 'job_commit_GiB': 4, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': 6, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-27: "跑完的那批先给我看看bscan" — the finished batch2d_slope_t2 is single-trace '
            'by design and cannot yield B-scans; this MT derivative (33-rx array, same convention as '
            'batch2d_v1_mt) is the minimal way to visualize the slope tier-T2 variants as B-scans',
            'user 2026-09-27 five-point reply already authorized the first 2D slope batch ("首批做个2D的即可，'
            '具体你来指定，最大效率验证即可"); MT is the same six mothers with only the receiver array changed',
            'cost basis: slope base measured 17.8-23.6 s/case post vctip-warden; MT increment measured '
            'negligible on batch2d_v1_mt; 6-case conservative budget 20 min; cap wall 20 min/case',
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
        'note': 'batch2d_slope_t2_mt frozen 2026-09-27 (6 slope mothers x 33 rx, BASE, single-shot '
                'variable-offset array, anchor mt17 == single-trace Rx cell 666). approved_to_simulate '
                'stays false until the recorded basis is applied; G1/G3/G4/G5 untouched; no '
                'physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: 6 .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}_mt, approved_to_simulate=false, 6 contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
