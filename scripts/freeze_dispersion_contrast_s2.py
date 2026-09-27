"""Freeze the dispersion contrast segment: dispersive-material anchor traces.

Basis: user 2026-09-28 "做吧" — schedules the dispersion contrast batch proposed in
docs/research/2026-09-28_dispersion_parameters.md (step 1): quantify the effect of the
fitted Debye dispersion (configs/research/dispersion_materials_v0.1.json) on the
interface arrival time and amplitude before deciding a batch-wide material switch.

Design (anchor geometry Tx y=15.35 / Rx y=16.65, z=45, identical to the archived
mother single-trace runs, which are the non-dispersive baselines and are NOT rerun):

- B2D-C3mS2-BG-DISP-T17: mother B2D-C3mS2-BG; cover material replaced by
  dispersive Debye (eps_inf 18.017, dE 7.878, tau 6.4567e-9, sigma_dc 0.003);
  rock unchanged (non-dispersive per dispersion_materials_v0.1).
- B2D-C3mS2TZ-BG-DISP-T17: mother B2D-C3mS2TZ-BG; cover as above; tzone replaced by
  half-strength pole (eps_inf 12.5, dE 4.0, same tau, sigma_dc 0.003); rock unchanged.

Everything else byte-identical to the mothers (grid, window, PML, source, rx).

Baselines: artifacts/simulations/2026-09-27_B2D-C3mS2-BG (bit-identical geometry to
CO33-t17, certified 2026-09-28) and artifacts/simulations/2026-09-27_B2D-C3mS2TZ-BG.

Mirrors scripts/freeze_pml_sensitivity_s2bg.py structure; --verify-only rebuilds in
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
OUT = ROOT / 'configs/research/dispersion_contrast_s2'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
DISP = json.loads((ROOT / 'configs/research/dispersion_materials_v0.1.json').read_text(encoding='utf-8'))

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
TIME_WINDOW_S = 1200e-9
BATCH_ID = 'dispersion_contrast_s2_t17'
RUN_DATE = '2026-09-28'

COVER = DISP['materials']['cover_clay']
TZ = DISP['materials']['transition_zone']
assert not DISP['materials']['bedrock_sandstone']['dispersive']

COVER_MAT = f"#material: {COVER['epsilon_inf']} {COVER['sigma_dc']} 1 0 cover"
COVER_POLE = (f"#add_dispersion_debye: 1 {COVER['debye_poles'][0]['delta_epsilon']} "
              f"{COVER['debye_poles'][0]['tau_s']} cover")
TZ_MAT = f"#material: {TZ['epsilon_inf']} {TZ['sigma_dc']} 1 0 tzone"
TZ_POLE = (f"#add_dispersion_debye: 1 {TZ['debye_poles'][0]['delta_epsilon']} "
           f"{TZ['debye_poles'][0]['tau_s']} tzone")

VARIANTS = [
    # (run_id, mother_id, replacements [(old, new)], insertions [(anchor_line, new_line)], baseline dir)
    ('B2D-C3mS2-BG-DISP-T17', 'B2D-C3mS2-BG',
     [('#material: 16 0.01 1 0 cover', COVER_MAT)],
     [(COVER_MAT, COVER_POLE)],
     '2026-09-27_B2D-C3mS2-BG'),
    ('B2D-C3mS2TZ-BG-DISP-T17', 'B2D-C3mS2TZ-BG',
     [('#material: 16 0.01 1 0 cover', COVER_MAT),
      ('#material: 12 0.005 1 0 tzone', TZ_MAT)],
     [(COVER_MAT, COVER_POLE), (TZ_MAT, TZ_POLE)],
     '2026-09-27_B2D-C3mS2TZ-BG'),
]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(run_id: str) -> int:
    return int.from_bytes(hashlib.sha256(f'{BATCH_ID}|{run_id}'.encode('utf-8')).digest()[:8], 'big')


def build_input(mother_text: str, run_id: str, replacements, insertions) -> str:
    lines = mother_text.splitlines()
    assert '#hertzian_dipole: x inf 15.35 45 impulse' in lines
    assert '#rx: inf 16.65 45 measurement Ex' in lines
    out = []
    for l in lines:
        if l.startswith('#title:'):
            out.append(f'#title: {run_id} {BATCH_ID} dispersion contrast (dispersion_materials_v0.1)')
            continue
        new = l
        for old, rep in replacements:
            if l == old:
                new = rep
        out.append(new)
        for anchor, ins in insertions:
            if new == anchor:
                out.append(ins)
    text = '\n'.join(out) + '\n'
    for old, rep in replacements:
        assert old not in text and rep in text
    for anchor, ins in insertions:
        assert ins in text
    text.encode('ascii')  # gprMax reads inputs with the system codec (GBK here); ASCII-only
    return text


def collect():
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30

    static = {
        'schema': 'dispersion_contrast_s2_static_check/1',
        'checked_date': RUN_DATE,
        'question': 'effect of the fitted Debye dispersion (dispersion_materials_v0.1) on '
                    'interface arrival time and amplitude, anchor-trace geometry',
        'dispersion_config': 'configs/research/dispersion_materials_v0.1.json',
        'dispersion_config_sha256': sha256_bytes(
            (ROOT / 'configs/research/dispersion_materials_v0.1.json').read_bytes()),
        'anchor': {'tx_y_m': 15.35, 'rx_y_m': 16.65, 'z_m': 45.0},
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'variants': [{'run_id': r, 'mother': m, 'baseline_dir': f'artifacts/simulations/{b}'}
                     for r, m, _, _, b in VARIANTS],
        'controls': 'direct wave propagates in air only and must be unchanged; the surface '
                    'reflection may change slightly (ground reflection coefficient moves with '
                    'band-centre permittivity 16 -> ~18)',
    }
    assert cells == 2560000

    inputs, cases_out = {}, []
    for run_id, mother_id, replacements, insertions, baseline in VARIANTS:
        mother_text = (SLOPE_BASE / (mother_id + '.in')).read_text(encoding='utf-8')
        text = build_input(mother_text, run_id, replacements, insertions)
        inputs[run_id] = text
        cases_out.append({
            'run_id': run_id, 'file': run_id + '.in',
            'input_sha256': sha256_bytes(text.encode('utf-8')),
            'mother_model_id': mother_id, 'case_id': run_id, 'group_id': BATCH_ID,
            'segment': 'dispersion_contrast', 'seed': seed_for(run_id),
            'baseline_run_directory': f'artifacts/simulations/{baseline}',
            'anchor': {'tx_y_m': 15.35, 'rx_y_m': 16.65, 'z_m': 45.0},
        })

    budget = {
        'schema': 'dispersion_contrast_s2_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'measured_co_wall_s': 19.3,
                        'note': 'same grid/window as archived CO t17; Debye ADE adds one '
                                'polarisation accumulator per dispersive material (negligible)'},
        'batch_totals': {'n_cases': len(VARIANTS), 'measured_estimate_s': len(VARIANTS) * 25,
                         'conservative_budget_min': 5},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 1, 'retries': 0, 'max_fdtd_runs': len(VARIANTS)},
        'known_limitations': ['estimates are for cap-setting only, not ETA'],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': 'dispersion_contrast', 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-28_dispersion_parameters.md + '
                         'configs/research/dispersion_materials_v0.1.json',
                 'n_cases': len(VARIANTS), 'n_exceptions': 0,
                 'derivation': 'mother text byte-identical except #title, cover/tzone '
                               '#material lines and the inserted #add_dispersion_debye lines',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'dispersion_contrast', 'grid_tier': 'BASE',
                  'convention': 'single diagnostic group; not part of any dev/validation split',
                  'groups': [{'case_id': BATCH_ID, 'group_id': BATCH_ID,
                              'run_id_pattern': 'B2D-*-DISP-T17', 'n_traces': len(VARIANTS)}]}
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
        print(f'verify-only: all {len(inputs)} dispersion-contrast inputs + contracts byte-identical')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_dispersion_contrast_s2.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'DISP-CONTRAST-S2-T17',
        'run_id': rid,
        'input_path': f'configs/research/dispersion_contrast_s2/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (S2 BG, S2TZ BG); any failure aborts the batch, '
                        'later attempts untouched. No retries.',
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
            'user 2026-09-28: "做吧" — schedules the dispersion contrast batch proposed in '
            'docs/research/2026-09-28_dispersion_parameters.md (quantify dispersion effect on '
            'interface arrival/amplitude before any batch-wide material switch)',
            'parameters: configs/research/dispersion_materials_v0.1.json (fitted to literature '
            'anchors 2026-09-28; cover single-pole Debye, tzone half-strength pole, rock '
            'non-dispersive)',
            'cost basis: measured 16.8-23.3 s/case on identical grid 2026-09-28; 2 cases '
            'conservative 5 min total; per-case wall cap 20 min',
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
        'note': 'dispersion_contrast_s2_t17 frozen 2026-09-28 (2 anchor traces, dispersive '
                'cover/tzone vs archived non-dispersive baselines). approved_to_simulate stays '
                'false until the recorded basis is applied; G1/G3/G4/G5 untouched; no '
                'physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
