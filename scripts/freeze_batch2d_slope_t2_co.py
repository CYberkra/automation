"""Freeze the batch2d_slope_t2 common-offset (CO) segment: 3 BG mothers x 33 CO traces.

Basis: user 2026-09-27 — "共偏移+能否换成GSSI400MHz，跑一次至少确定一个标准正向样本".
The MT array (fixed Tx, variable-offset receivers) is a shot gather, not the classic
common-offset B-scan the user expects; this segment derives a true common-offset
profile (Tx = Rx - 1.30 m, both at z = 45, Rx y = 12.65..20.65 step 0.25, 33 traces,
one FDTD run per trace) for the BG-only models per the 2026-09-27 ablation decision
(target ablated; only cover clay + bedrock interface).

Mothers (BG only):
- configs/research/batch2d_slope_t2/B2D-C3mS2-BG.in   (slope 11.31 deg)
- configs/research/batch2d_slope_t2/B2D-C3mS2TZ-BG.in (slope + 1 m transition zone)
- configs/research/batch2d_v1/B2D-C3m-BG.in           (flat control)

Anchor: trace t17 (Rx 16.65, Tx 15.35) reproduces each mother's single-trace source/
receiver pair exactly, so its Ex must be bit-identical to the archived mother runs
(post-run check, same as mt17 anchor in the MT segments).

Mirrors scripts/freeze_batch2d_slope_t2_mt.py structure: generates .in inputs,
cases/groups contracts, static_check.json, budget.json, the runner clone, and rewrites
the global execution gate with approved_to_simulate=false. --verify-only rebuilds in
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
FLAT_BASE = ROOT / 'configs/research/batch2d_v1'
OUT = ROOT / 'configs/research/batch2d_slope_t2_co'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
N_TRACES = 33
DY_RX = 0.25
Y_FIRST = 12.65
OFFSET = 1.30
Z_LINE = 45.0
ANCHOR_INDEX = 16  # 0-based; t17 == mother single-trace pair (Rx 16.65, Tx 15.35)
TIME_WINDOW_S = 1200e-9
VARIANT = 'common_offset_33tx_0p25m_off1p30'
BATCH_ID = 'batch2d_slope_t2_co'
RUN_DATE = '2026-09-27'

MOTHERS = [
    ('B2D-C3mS2-BG', SLOPE_BASE, 'slope S2 (theta_eff 11.309932 deg)'),
    ('B2D-C3mS2TZ-BG', SLOPE_BASE, 'slope S2 + 1 m transition zone'),
    ('B2D-C3m-BG', FLAT_BASE, 'flat control (batch2d_v1 mother)'),
]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(mother_id: str, k: int) -> int:
    txt = f'{BATCH_ID}|{mother_id}|t{k + 1:02d}|{VARIANT}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def build_co_input(mother_text: str, run_id: str, mother_id: str, tx_y: float, rx_y: float, rx_id: str) -> str:
    lines = mother_text.splitlines()
    dip_idx = [i for i, l in enumerate(lines) if l.startswith('#hertzian_dipole:')]
    rx_idx = [i for i, l in enumerate(lines) if l.startswith('#rx:')]
    title_idx = [i for i, l in enumerate(lines) if l.startswith('#title:')]
    assert len(dip_idx) == 1 and len(rx_idx) == 1 and len(title_idx) == 1, mother_id
    assert lines[dip_idx[0]] == '#hertzian_dipole: x inf 15.35 45 impulse', mother_id
    new_lines = list(lines)
    new_lines[title_idx[0]] = f'#title: {run_id} {BATCH_ID} common-offset pair offset {OFFSET} m (mother {mother_id})'
    new_lines[dip_idx[0]] = f'#hertzian_dipole: x inf {tx_y:.2f} 45 impulse'
    new_lines[rx_idx[0]] = f'#rx: inf {rx_y:.2f} 45 {rx_id} Ex'
    return '\n'.join(new_lines) + '\n'


def collect():
    """Build all in-memory artifacts. Returns (inputs, cases_doc, groups_doc, budget, static)."""
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    static = {
        'schema': 'batch2d_slope_t2_co_static_check/1',
        'checked_date': RUN_DATE,
        'base_batches': ['batch2d_slope_t2 (slope BG mothers)', 'batch2d_v1 (flat BG mother)'],
        'ablation': 'BG only per user 2026-09-27 ablation decision (target ablated)',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT,
                 'dt_formula_over_contract': dt_formula / DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'time_steps_contract': 20352,
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'common_offset_array': {
            'n_traces': N_TRACES, 'spacing_m': DY_RX, 'offset_m': OFFSET,
            'rx_first_y_m': Y_FIRST, 'rx_last_y_m': Y_FIRST + (N_TRACES - 1) * DY_RX,
            'tx_first_y_m': Y_FIRST - OFFSET, 'tx_last_y_m': Y_FIRST + (N_TRACES - 1) * DY_RX - OFFSET,
            'z_m': Z_LINE, 'anchor_index': ANCHOR_INDEX,
            'anchor_rx_y_m': Y_FIRST + ANCHOR_INDEX * DY_RX,
            'anchor_tx_y_m': Y_FIRST + ANCHOR_INDEX * DY_RX - OFFSET,
            'trace_ids': [f't{k + 1:02d}' for k in range(N_TRACES)],
        },
        'mother_checks': [], 'trace_checks': [],
    }
    assert cells == 2560000
    assert abs(main_fields_gib - 0.3437627702951431) < 1e-9
    assert abs(static['common_offset_array']['anchor_rx_y_m'] - 16.65) < 1e-9
    assert abs(static['common_offset_array']['anchor_tx_y_m'] - 15.35) < 1e-9

    inputs, cases_out, groups_out = {}, [], []
    key_lines = {}
    for mother_id, base, desc in MOTHERS:
        mother_text = (base / (mother_id + '.in')).read_text(encoding='utf-8')
        keys = {k: [l for l in mother_text.splitlines() if l.startswith(k)]
                for k in ('#domain:', '#dx_dy_dz:', '#time_window:', '#pml_cells:', '#waveform:')}
        for k, v in keys.items():
            assert len(v) == 1, (mother_id, k)
            key_lines.setdefault(k, set()).add(v[0])
        static['mother_checks'].append({'mother': mother_id, 'base': str(base.relative_to(ROOT)),
                                        'desc': desc,
                                        'mother_sha256': hashlib.sha256(
                                            mother_text.encode('utf-8')).hexdigest()})
        for k in range(N_TRACES):
            rx_y = Y_FIRST + k * DY_RX
            tx_y = rx_y - OFFSET
            rx_cell, tx_cell = rx_y / DY, tx_y / DY
            assert abs(rx_cell - round(rx_cell)) < 1e-9 and abs(tx_cell - round(tx_cell)) < 1e-9
            assert 1.0 < tx_y and rx_y < 31.0, f'inside PML: tx {tx_y} rx {rx_y}'
            rx_id = f't{k + 1:02d}'
            run_id = f'{mother_id}-CO33-{rx_id}'
            text = build_co_input(mother_text, run_id, mother_id, tx_y, rx_y, rx_id)
            inputs[run_id] = text
            strip = lambda s: [l for l in s.splitlines()
                               if not l.startswith(('#rx:', '#title:', '#hertzian_dipole:'))]
            assert strip(mother_text) == strip(text), run_id
            if k == ANCHOR_INDEX:
                anchor_lines = {l.split(':', 1)[0] for l in text.splitlines()}
                assert '#rx: inf 16.65 45 t17 Ex' in text
                assert '#hertzian_dipole: x inf 15.35 45 impulse' in text
            static['trace_checks'].append({'run_id': run_id, 'rx_cell': int(round(rx_cell)),
                                           'tx_cell': int(round(tx_cell)), 'outside_pml': True})
            cases_out.append({
                'run_id': run_id, 'file': run_id + '.in',
                'input_sha256': sha256_bytes(text.encode('utf-8')),
                'mother_model_id': mother_id, 'case_id': f'{mother_id}|{rx_id}',
                'group_id': mother_id,
                'segment': 'co', 'variant_tag': VARIANT,
                'seed': seed_for(mother_id, k),
                'common_offset': {
                    'n_traces': N_TRACES, 'trace_spacing_m': DY_RX, 'offset_m': OFFSET,
                    'rx_y_m': rx_y, 'tx_y_m': tx_y, 'z_m': Z_LINE,
                    'trace_index': k, 'anchor_trace': k == ANCHOR_INDEX,
                    'acquisition_geometry': 'common-offset profile along y (Tx = Rx - 1.30 m, both z=45)',
                    'rx_output_components': ['Ex'],
                },
            })
        groups_out.append({'case_id': mother_id, 'group_id': mother_id,
                           'run_id_pattern': f'{mother_id}-CO33-t*',
                           'n_traces': N_TRACES, 'variant_tag': VARIANT})
    for k, v in key_lines.items():
        assert len(v) == 1, f'mothers disagree on {k}: {v}'
    static['mothers_share_grid_source_pml'] = {k: sorted(v) for k, v in key_lines.items()}

    measured = [17.812, 23.61]
    budget = {
        'schema': 'batch2d_slope_t2_co_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'dt_formula_s': dt_formula, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'measured_slope_mt_wall_range_s': measured,
                        'note': 'CO trace == mother geometry with moved Tx/Rx pair; per-trace cost '
                                'equals the measured slope single-trace/MT wall (17.8-23.6 s post '
                                'vctip-warden)'},
        'batch_totals': {'n_cases': 3 * N_TRACES,
                         'measured_estimate_s': 3 * N_TRACES * 25,
                         'conservative_budget_min': 90},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 3 * N_TRACES},
        'known_limitations': [
            'per-case supervised process startup included in measured range',
            'vctip.exe job-lingering can add ~15 min per case if the out-of-job singleton lapses',
            'estimates are for cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-27_slope_family_design_draft.md + '
                         'docs/research/2026-09-26_common_offset_subset_proposal.md (CO convention) + '
                         'decision_log 2026-09-27 ablation (BG only)',
                 'n_cases': 3 * N_TRACES, 'n_exceptions': 0,
                 'derivation': 'mother text byte-identical except #title/#hertzian_dipole/#rx lines; '
                               'anchor t17 reproduces each mother single-trace Tx/Rx pair',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                  'convention': 'spec 4.2: any split is by group_id (mother model) only, never by case, '
                                'window or trace; each mother is its own family',
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
        print(f'verify-only: all {len(inputs)} CO inputs + cases/groups/budget/static_check byte-identical')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_slope_t2_co.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-SLOPE-T2-CO',
        'run_id': rid,
        'input_path': f'configs/research/batch2d_slope_t2_co/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (S2 BG, S2TZ BG, flat C3 BG; t01..t33); any failure '
                        'aborts the batch, later attempts untouched. No retries.',
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
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': 3 * N_TRACES, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-27: "共偏移+能否换成GSSI400MHz，跑一次至少确定一个标准正向样本" — common-offset '
            'is the classic B-scan geometry the user expects; MT (fixed Tx variable-offset) cannot '
            'serve that role',
            'user 2026-09-27 ablation decision: BG only (cover clay + bedrock interface); TGT/NC '
            'variants suspended, archived results retained',
            'cost basis: measured 17.8-23.6 s/case post vctip-warden on identical grid; 99 cases '
            'conservative 90 min total; per-case wall cap 20 min',
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
        'note': 'batch2d_slope_t2_co frozen 2026-09-27 (3 BG mothers x 33 common-offset traces, '
                'offset 1.30 m, BASE, BG-only per ablation). approved_to_simulate stays false until '
                'the recorded basis is applied; G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
