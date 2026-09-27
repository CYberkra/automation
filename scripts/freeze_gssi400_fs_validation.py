"""Freeze the GSSI 400 MHz free-space antenna pipeline-validation run (1 case).

Basis: user 2026-09-27 — "共偏移+能否换成GSSI400MHz，跑一次至少确定一个标准正向样本，
现在哪里有问题的话也好发现？在本地能跑得了吗". Purpose is a STANDARD PIPELINE SAMPLE:
run the official, published antenna model (Stadler et al. 2022 update, shipped in
gprMax v4 toolboxes/GPRAntennaModels) in free space exactly as the official example
specifies, on the local GPU install, so any local toolchain/solver defect shows up
against a known-good reference. It is NOT a project sample: the GSSI 400 band is far
above the project SFCW 20-170 MHz band; results must not be mixed with project data.

Generates: manifest.json (hashes of runner, official example, toolbox model, pulse
file, runtime/cuda identity) and rewrites the global execution gate with
approved_to_simulate=false. --verify-only rebuilds in memory and asserts equality.

No solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/gssi400_fs_validation_v1'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
V4ENV = ROOT / 'artifacts/local_checks/gprmax_v4_gpu_env'
TOOLBOX = V4ENV / 'Lib/site-packages/gprMax/toolboxes/GPRAntennaModels'
RUNNER = ROOT / 'scripts/run_gssi400_fs_gpu.py'
LAUNCHER = ROOT / 'scripts/run_gssi400_fs_supervised.py'
OFFICIAL_EXAMPLE = Path('E:/gprMax-v.4.0.0/gprMax-v.4.0.0/examples/gpr/antennas/antenna_like_GSSI_400_fs.py')

RUN_ID = 'GSSI400FS-VALID'
BATCH_ID = 'gssi400_fs_validation_v1'
RUN_DATE = '2026-09-27'

C = 299792458.0
DL = 0.002
DOMAIN = (0.340, 0.340, 0.318)
TIME_WINDOW_S = 15e-9


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def collect():
    dt_formula = 1.0 / (C * math.sqrt(3 * DL ** -2))
    nx, ny, nz = (round(v / DL) for v in DOMAIN)
    cells = nx * ny * nz
    manifest = {
        'schema': 'gssi400_fs_validation_v1_manifest/1',
        'date': RUN_DATE,
        'purpose': 'official antenna pipeline validation sample (NOT a project SFCW sample; '
                   'GSSI 400 band >> project 20-170 MHz)',
        'run_id': RUN_ID,
        'model': {
            'name': 'antenna_like_GSSI_400 (Stadler et al. 2022 update, official gprMax v4 toolbox)',
            'discretisation_m': DL,
            'domain_m': list(DOMAIN),
            'cells': {'nx': nx, 'ny': ny, 'nz': nz, 'total': cells},
            'time_window_s': TIME_WINDOW_S,
            'dt_formula_s': dt_formula,
            'time_steps_estimate': math.ceil(TIME_WINDOW_S / dt_formula),
            'antenna_external_dims_m': [0.3, 0.3, 0.178],
            'antenna_position_m': [0.170, 0.170, 0.100],
            'medium': 'free space',
            'precision': 'single (official example default; project batches use double)',
        },
        'hashes': {
            'runner_script': {'path': 'scripts/run_gssi400_fs_gpu.py', 'sha256': sha256_file(RUNNER)},
            'launcher_script': {'path': 'scripts/run_gssi400_fs_supervised.py', 'sha256': sha256_file(LAUNCHER)},
            'official_example_reference': {'path': str(OFFICIAL_EXAMPLE), 'sha256': sha256_file(OFFICIAL_EXAMPLE)},
            'toolbox_gssi_py': {'path': str(TOOLBOX / 'GSSI.py'), 'sha256': sha256_file(TOOLBOX / 'GSSI.py')},
            'toolbox_gssi_400_pulse': {'path': str(TOOLBOX / 'GSSI_400MHz_pulse.txt'),
                                       'sha256': sha256_file(TOOLBOX / 'GSSI_400MHz_pulse.txt')},
            'runtime_identity': {'path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
                                 'sha256': sha256_file(RUNTIME_ID)},
            'cuda_identity': {'path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
                              'sha256': sha256_file(CUDA_ID)},
        },
        'known_limitations': [
            'validation of install/toolbox/GPU chain only; no bearing on project 20-170 MHz samples',
            'single precision follows the official example default; double-precision project batches unchanged',
            'free-space only; no ground, no claim about antenna-ground coupling',
        ],
    }
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    manifest = collect()

    if args.verify_only:
        assert json.loads((OUT / 'manifest.json').read_text(encoding='utf-8')) == manifest
        print('verify-only: manifest byte-consistent')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=1) + '\n',
                                       encoding='utf-8')

    old_gate = json.loads(GATE.read_text(encoding='utf-8'))
    contract = {
        'packet_id': 'GSSI400-FS-VALIDATION',
        'run_id': RUN_ID,
        'input_path': 'scripts/run_gssi400_fs_gpu.py',
        'input_sha256': manifest['hashes']['runner_script']['sha256'],
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': manifest['hashes']['runtime_identity']['sha256'],
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': manifest['hashes']['cuda_identity']['sha256'],
        'launcher_sha256': manifest['hashes']['launcher_script']['sha256'],
        'supervisor_sha256': sha256_file(ROOT / 'scripts/bounded_windows_process.py'),
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{RUN_ID}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{RUN_ID}',
        'precision': 'single',
        'continuation': 'Single case; any failure aborts. No retries.',
    }
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
        'approved_run_ids_pending_agreement': [RUN_ID],
        'approved_compute_budget': {
            'wall_minutes': 30, 'job_commit_GiB': 8, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 1, 'threads': 8, 'backend': 'CUDA', 'precision': 'single',
            'retries': 0, 'max_fdtd_runs': 1, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 2,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-27: "能否换成GSSI400MHz，跑一次至少确定一个标准正向样本，现在哪里有问题的话'
            '也好发现？在本地能跑得了吗" — official published antenna model in free space as a standard '
            'pipeline-validation sample against local install/GPU chain',
            'cost basis: 170x170x159 cells (~4.6M), ~3898 steps, single precision — far below measured '
            'project 2D job peaks; conservative 30 min wall cap',
            'band caveat recorded in manifest: result is not a project sample and must not be mixed '
            'with 20-170 MHz SFCW evidence',
        ],
        'permitted_preparation': old_gate['permitted_preparation'],
        'requires_agreement_before_execution': old_gate['requires_agreement_before_execution'],
        'approved_execution_contracts': [contract],
        'execution_policy': old_gate['execution_policy'],
        'execution_outcome': {'status': 'pending_user_agreement'},
        'last_completed_execution_contract': {
            'batch_id': old_gate['batch_id'],
            'execution_outcome': old_gate.get('execution_outcome'),
            'note': old_gate.get('note'),
        },
        'note': 'gssi400_fs_validation_v1 frozen 2026-09-27 (official antenna_like_GSSI_400 free-space, '
                'single precision, 1 case). approved_to_simulate stays false until the recorded basis is '
                'applied; G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}/manifest.json')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, 1 contract')


if __name__ == '__main__':
    main()
