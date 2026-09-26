"""Freeze the batch2d_v1 common-offset (CO) segment: 2 mothers x 11 traces (proposal 2026-09-26).

Generates 22 .in inputs (Tx/Rx pair translated along y, fixed 1.30 m offset; anchor t05 is
byte-identical to the mother input), cases/groups contracts, static_check.json, budget.json,
the chunked in-process launcher + worker, and rewrites the global execution gate with
approved_to_simulate=false. Execution requires the user's explicit agreement.

Read-only over archived contracts; no solver is launched here.
"""
import hashlib
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'configs/research/batch2d_v1'
OUT = ROOT / 'configs/research/batch2d_v1_co'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
N_TRACES = 11
DY_TRACE = 1.0
OFFSET = 1.30
Y_FIRST_RX = 12.65
RX_Z = 45.0
ANCHOR_INDEX = 4  # 0-based; t05 == mother geometry (Rx 16.65 / Tx 15.35)
TIME_WINDOW_S = 1200e-9
VARIANT = 'common_offset_11tr_1p0m'
BATCH_ID = 'batch2d_v1'
SEGMENT = 'co'
MOTHERS = ['B2D-C3m-BG', 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02']
CHUNK_SIZE = 6


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(group_id: str, case_id: str, variant: str) -> int:
    txt = f'{BATCH_ID}|{group_id}|{case_id}|{variant}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def build_co_input(mother_text: str, mother_id: str, run_id: str, k: int) -> str:
    """Replace only the Tx/Rx lines (and title, except the anchor trace which stays
    byte-identical to the mother input)."""
    rx_y = Y_FIRST_RX + k * DY_TRACE
    tx_y = rx_y - OFFSET
    lines = mother_text.splitlines()
    tx_idx = [i for i, l in enumerate(lines) if l.startswith('#hertzian_dipole:')]
    rx_idx = [i for i, l in enumerate(lines) if l.startswith('#rx:')]
    assert len(tx_idx) == 1 and len(rx_idx) == 1, run_id
    new_lines = list(lines)
    new_lines[tx_idx[0]] = f'#hertzian_dipole: x inf {tx_y:.2f} 45 impulse'
    new_lines[rx_idx[0]] = f'#rx: inf {rx_y:.2f} 45 measurement Ex'
    if k != ANCHOR_INDEX:
        title_idx = next(i for i, l in enumerate(lines) if l.startswith('#title:'))
        new_lines[title_idx] = f'#title: {run_id} batch2d_v1_co BASE common-offset b-scan'
    # preserve the mother's line endings exactly (anchor trace must stay byte-identical)
    nl = '\r\n' if '\r\n' in mother_text else '\n'
    return nl.join(new_lines) + nl


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base_cases = json.loads((BASE / 'cases.json').read_text(encoding='utf-8'))
    base_groups = json.loads((BASE / 'groups.json').read_text(encoding='utf-8'))
    by_id = {c['run_id']: c for c in base_cases['cases']}
    missing = [m for m in MOTHERS if m not in by_id]
    assert not missing, missing

    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    static = {
        'schema': 'batch2d_v1_co_static_check/1',
        'checked_date': '2026-09-26',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT,
                 'dt_formula_over_contract': dt_formula / DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'time_steps_contract': 20352,
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'common_offset': {
            'n_traces': N_TRACES, 'spacing_m': DY_TRACE, 'offset_m': OFFSET,
            'first_rx_y_m': Y_FIRST_RX, 'last_rx_y_m': Y_FIRST_RX + (N_TRACES - 1) * DY_TRACE,
            'first_tx_y_m': Y_FIRST_RX - OFFSET,
            'last_tx_y_m': Y_FIRST_RX + (N_TRACES - 1) * DY_TRACE - OFFSET,
            'z_m': RX_Z, 'anchor_trace_id': 't05',
            'anchor_rx_y_m': Y_FIRST_RX + ANCHOR_INDEX * DY_TRACE,
            'anchor_tx_y_m': Y_FIRST_RX + ANCHOR_INDEX * DY_TRACE - OFFSET,
            'trace_ids': [f't{k + 1:02d}' for k in range(N_TRACES)],
        },
        'trace_checks': [], 'anchor_input_byte_identical_to_mother': [],
        'materials_unchanged_vs_mother': [], 'source_unchanged_vs_mother': [],
    }
    assert cells == 2560000
    assert abs(dt_formula / DT_CONTRACT - 1.0) < 1e-12

    cases_out, groups_out, contracts = [], [], []
    runner_src = ROOT / 'scripts/run_approved_batch2d.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_co.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1'" in txt
    # The CO launcher differs structurally (chunked in-process execution), so it is
    # written from its own template below, not string-replaced from the serial launcher.
    assert runner_dst.exists(), 'write scripts/run_approved_batch2d_co.py before freezing'
    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    worker_path = ROOT / 'scripts/run_co_chunk_worker.py'
    assert worker_path.exists(), 'write scripts/run_co_chunk_worker.py before freezing'
    worker_sha = sha256_bytes(worker_path.read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    order = []  # execution order: all BG traces, then all A0 traces
    per_mother = {}
    for mother_id in MOTHERS:
        m = by_id[mother_id]
        mother_text = (BASE / m['file']).read_bytes().decode('utf-8')
        mother_sha = sha256_bytes((BASE / m['file']).read_bytes())
        entries = []
        for k in range(N_TRACES):
            run_id = f'{mother_id}-CO11-t{k + 1:02d}'
            co_text = build_co_input(mother_text, mother_id, run_id, k)
            in_path = OUT / (run_id + '.in')
            in_path.write_bytes(co_text.encode('utf-8'))
            in_bytes = in_path.read_bytes()
            in_sha = sha256_bytes(in_bytes)
            if k == ANCHOR_INDEX:
                assert in_bytes == (BASE / m['file']).read_bytes(), \
                    f'anchor {run_id} input is not byte-identical to mother'
                static['anchor_input_byte_identical_to_mother'].append(
                    {'run_id': run_id, 'mother': mother_id, 'input_sha256': in_sha})
            # static per-trace checks
            rx_y = Y_FIRST_RX + k * DY_TRACE
            tx_y = rx_y - OFFSET
            tx_c, rx_c = tx_y / DY, rx_y / DY
            z_c = RX_Z / DZ
            assert abs(tx_c - round(tx_c)) < 1e-9 and abs(rx_c - round(rx_c)) < 1e-9
            assert 1.0 < tx_y and rx_y < 31.0, f'tx/rx inside PML: {tx_y} {rx_y}'
            # everything except title/tx/rx byte-identical to mother
            def strip(s):
                return [l for l in s.splitlines()
                        if not l.startswith('#rx:') and not l.startswith('#title:')
                        and not l.startswith('#hertzian_dipole:')]
            assert strip(mother_text) == strip(co_text), run_id
            assert len([l for l in co_text.splitlines() if l.startswith('#rx:')]) == 1
            assert co_text.splitlines()[[i for i, l in enumerate(co_text.splitlines())
                                         if l.startswith('#rx:')][0]].endswith(' Ex')
            static['trace_checks'].append({
                'run_id': run_id, 'tx_y_m': tx_y, 'rx_y_m': rx_y,
                'tx_cell': int(round(tx_c)), 'rx_cell': int(round(rx_c)),
                'z_cell': int(round(z_c)), 'integer_cells': True, 'outside_pml': True})
            static['materials_unchanged_vs_mother'].append(run_id)
            static['source_unchanged_vs_mother'].append(run_id)

            entry = dict(m)
            entry.update(run_id=run_id, file=run_id + '.in', input_sha256=in_sha,
                         segment=SEGMENT, variant_tag=VARIANT,
                         seed=seed_for(m['group_id'], m['case_id'], f'{VARIANT}-t{k + 1:02d}'),
                         geometry=dict(m['geometry'], tx_m=[tx_y, RX_Z], rx_m=[rx_y, RX_Z]),
                         common_offset={
                             'n_traces': N_TRACES, 'trace_spacing_m': DY_TRACE,
                             'offset_m': OFFSET,
                             'rx_y_range_m': [Y_FIRST_RX, Y_FIRST_RX + (N_TRACES - 1) * DY_TRACE],
                             'tx_y_range_m': [Y_FIRST_RX - OFFSET,
                                              Y_FIRST_RX + (N_TRACES - 1) * DY_TRACE - OFFSET],
                             'anchor_trace_id': 't05',
                             'anchor_rx_y_m': 16.65, 'anchor_tx_y_m': 15.35,
                             'trace_id': f't{k + 1:02d}',
                             'trace_index': k,
                             'trace_order': 'construction',
                             'acquisition_geometry': 'common-offset b-scan (Tx-Rx pair '
                                                     'translated along y, fixed 1.30 m offset)',
                             'rx_output_components': ['Ex'],
                             'cross_geometry_pairing': 'forbidden: never pair with '
                                                       'variable-offset CSG traces',
                         })
            cases_out.append(entry)
            groups_out.append({'case_id': m['case_id'], 'group_id': m['group_id'],
                               'run_id': run_id,
                               'mother_model_hash': next(g['mother_model_hash']
                                                         for g in base_groups['groups']
                                                         if g['group_id'] == m['group_id']),
                               'seed': entry['seed'], 'variant_tag': VARIANT})
            entries.append({'run_id': run_id, 'input_sha256': in_sha})
        per_mother[mother_id] = entries
        order.extend(entries)

    # chunking: BG t01..t11 then A0 t01..t11, chunks of CHUNK_SIZE
    chunks = []
    for ci in range(0, len(order), CHUNK_SIZE):
        part = order[ci:ci + CHUNK_SIZE]
        chunks.append({'chunk_id': f'co_c{len(chunks) + 1:02d}', 'cases': part})
    for c in chunks:
        for case in c['cases']:
            case['chunk_id'] = c['chunk_id']

    for c in chunks:
        for case in c['cases']:
            run_id = case['run_id']
            contracts.append({
                'packet_id': 'BATCH2D-V1-CO',
                'input_path': f'configs/research/batch2d_v1_co/{run_id}.in',
                'input_sha256': case['input_sha256'],
                'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
                'runtime_identity_sha256': runtime_sha,
                'launcher_sha256': launcher_sha,
                'worker_sha256': worker_sha,
                'supervisor_sha256': supervisor_sha,
                'attempt_record': f'artifacts/research_checks/2026-09-26_{run_id}_attempt.json',
                'run_directory': f'artifacts/simulations/2026-09-26_{run_id}',
                'chunk_id': case['chunk_id'],
                'continuation': 'Chunks serial in frozen order (BG all, then A0); within a '
                                'chunk cases run in one process via gprMax.run() sequentially; '
                                'any failure aborts the batch, later chunks untouched. No retries.',
                'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
                'cuda_identity_sha256': cuda_sha,
                'run_id': run_id,
            })

    yfine_ref = {'wall_s': 1088.781, 'cells': 40960000, 'dt_s': 1.474e-11}
    wall_scaled = yfine_ref['wall_s'] * (cells / yfine_ref['cells']) * (yfine_ref['dt_s'] / DT_CONTRACT)
    budget = {
        'schema': 'batch2d_v1_co_budget/1',
        'date': '2026-09-26',
        'formulas': 'spec §3.1 (same as batch2d_v1_mt); execution via in-process multi-model '
                    'chunks (speedup probe 2026-09-26: ~21 s/case in-process vs 45-58 s/case '
                    'per-process; one-time ~30-35 s process cost per chunk)',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'dt_formula_s': dt_formula, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'wall_scaled_from_YFINE_s': wall_scaled,
                        'measured_inprocess_probe_s': [20.55, 21.27],
                        'planning_envelope_min': 2.5},
        'chunking': {'chunk_size': CHUNK_SIZE,
                     'n_chunks': len(chunks),
                     'chunk_case_counts': [len(c['cases']) for c in chunks],
                     'chunk_wall_cap_s': [300 + 150 * len(c['cases']) for c in chunks],
                     'measured_expectation_s': 35 + 21 * CHUNK_SIZE},
        'batch_totals': {
            'n_cases': 22,
            'planning_envelope_s': 22 * 150,
            'measured_estimate_s': 4 * 35 + 21 * 22,
            'conservative_budget_min': 30,
        },
        'hard_stops': {'wall_minutes_per_chunk': 20, 'job_commit_GiB': 4,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 22},
        'known_limitations': [
            'working scaling is not a guaranteed runtime',
            'main-field array excludes PML, update coefficients, CUDA context and host temporaries',
            'Job committed memory is not GPU VRAM',
            'estimates are for ranking and cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': f'{BATCH_ID}_co', 'segment': SEGMENT, 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-26_common_offset_subset_proposal.md '
                         '(pending user ratification at freeze time) + '
                         'docs/research/2026-09-26_batch_2d_spec_v1.md',
                 'n_cases': 22, 'n_exceptions': 0,
                 'derivation': 'mother_model_id identical to batch2d_v1 base cases (C3 BG + C3 A0); '
                               'Tx/Rx translated per trace, fixed 1.30 m offset; anchor t05 '
                               'byte-identical to mother input',
                 'cases': cases_out}
    groups_doc = {'batch_id': f'{BATCH_ID}_co', 'segment': SEGMENT, 'grid_tier': 'BASE',
                  'convention': 'spec §4.2: any split is by group_id only, never by case_id, window or trace',
                  'groups': groups_out}
    (OUT / 'cases.json').write_text(json.dumps(cases_doc, ensure_ascii=False, indent=1) + '\n',
                                    encoding='utf-8')
    (OUT / 'groups.json').write_text(json.dumps(groups_doc, ensure_ascii=False, indent=1) + '\n',
                                     encoding='utf-8')
    (OUT / 'budget.json').write_text(json.dumps(budget, ensure_ascii=False, indent=1) + '\n',
                                     encoding='utf-8')
    (OUT / 'static_check.json').write_text(json.dumps(static, ensure_ascii=False, indent=1) + '\n',
                                           encoding='utf-8')
    (OUT / 'chunks.json').write_text(json.dumps(
        {'schema': 'batch2d_v1_co_chunks/1', 'chunks': [
            {'chunk_id': c['chunk_id'],
             'run_ids': [x['run_id'] for x in c['cases']]} for c in chunks]},
        ensure_ascii=False, indent=1) + '\n', encoding='utf-8')

    old_gate = json.loads(GATE.read_text(encoding='utf-8'))
    new_gate = {
        'schema': old_gate['schema'],
        'updated_date': '2026-09-26',
        'authority': old_gate['authority'],
        'user_instruction': old_gate['user_instruction'],
        'target_solver': old_gate['target_solver'],
        'source_root_on_current_machine': old_gate['source_root_on_current_machine'],
        'batch_id': f'{BATCH_ID}_co',
        'approved_to_simulate': False,
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
            'common-offset subset proposal (docs/research/2026-09-26_common_offset_subset_proposal.md), '
            'ratified by user 2026-09-26 17:49 ("开始并执行"): answers whether trace-to-trace moveout on '
            'the variable-offset CSG changes operator-baseline (S4) conclusions',
            'parent proposal §2.4-B (2026-09-26_multitrace_batch_proposal.md) reserved this paired '
            'common-offset subset, contingent on the moveout-difference mattering',
            'speedup probe (docs/research/2026-09-26_speedup_probe_results.md): in-process multi-model '
            'execution measured 20.5-21.6 s/case vs 45-58 s/case per-process; chunk-of-6 wall ~3.5 min',
            'input package verified 2026-09-26: 22 .in derived from frozen C3 mothers; only '
            '#hertzian_dipole/#rx lines differ (title too, except anchor t05 which is byte-identical '
            'to the mother input); static_check.json + budget.json + chunks.json on disk',
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
        'note': 'batch2d_v1_co frozen 2026-09-26 (22 cases = 2 C3 mothers x 11 common-offset traces, '
                'BASE, double, chunked in-process execution). approved_to_simulate stays false until '
                'the user agrees to execute; results are stratified vs the MT33 batch and never '
                'cross-geometry paired; G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: 22 .in + cases/groups/budget/static_check/chunks')
    print(f'wrote {GATE.relative_to(ROOT)}: batch2d_v1_co, approved_to_simulate=false, '
          f'{len(contracts)} contracts, {len(chunks)} chunks')
    print('anchor t05 byte-identical to mothers:', static['anchor_input_byte_identical_to_mother'])
    print('NOTE: run this freeze ONLY after run_approved_batch2d_co.py and '
          'run_co_chunk_worker.py are written; contracts embed their SHA-256.')


if __name__ == '__main__':
    main()
