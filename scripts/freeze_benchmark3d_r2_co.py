"""Freeze the benchmark3d_r2_co batch: first 3D standard benchmark (rough bedrock-cover interface).

Basis: design draft docs/research/2026-10-01_3d_benchmark_design_draft.md (v0.1),
user sign-off "确认" 2026-10-01 22:23 covering S1-S5 (reduced domain 10x10x31 m physical,
13 traces CO, 600 ns, R2-grade interface RMS 0.16 m / CL 2.5 m isotropic, 5 cm grid with the
stated accuracy limit, optional flat-3D control deferred to after the first batch).

Composition (39 runs, frozen order):
  1) B3D5CM-C3mR2-BG-CO13-t01..t13  — 3D rough-interface benchmark, 5 cm, 38.0M cells
  2) B2D5CM-C3mR2-BG-CO13-t01..t13 — 2D TMx pair, same eta(y) slice, same 5 cm grid
  3) B2D-C3mR2-BG-CO13-t01..t13    — 2D TMx reference at the production 2.5 cm grid

Geometry (.in coords, 1 m PML offset on all sides):
  domain 12 x 12 x 33 m; rock z[1,9+eta], cover z[9+eta,12] (nominal 3 m), air z[12,33];
  antenna z=27 (15 m above surface); Tx y = 3.85+0.25k, Rx y = Tx+1.3, x=6.0; k=0..12.
  2D TMx versions share the same y/z layout with x collapsed to inf.

Interface field eta(x,y): zero-mean Gaussian-correlated random field,
  white noise -> Gaussian kernel sigma = CL/2 = 1.25 m (reflect-padded FFT filter, numpy-only)
  -> scaled to RMS 0.16 m -> clipped [-2.00, +1.95] -> quantized to 0.05 m.
  Seed 20261002; field stored in interface_field.npz, sha locked in sample_table.json.
  2D pairs use the slice eta(x=6.0, y) so 2D/3D differ ONLY by dimension.

Materials byte-identical to batch2d_b2_pilot_co (dispersion_materials_v0.1):
  rock (9, 0.001); cover (18.017, 0.003) + Debye (1 pole, dE 7.878, tau 6.4567e-09).
pml_cfs string copied byte-exact from configs/research/a0_3d_v1 (proven 3D 5cm on this machine)
  and batch2d_b2_pilot_co (proven 2D): quartic forward 0 0.21235349838321013.

Gate discipline: writes gate with approved_to_simulate=TRUE (user "确认" 2026-10-01 22:23, S5),
  batch_id=benchmark3d_r2_co, 39 contracts hash-locked (input/launcher/supervisor/runtime/cuda).
  Execution is additionally blocked on a pending OS reboot (NVML dead, CBS+WU reboot flags set,
  diagnosed 2026-10-01 ~22:30); the runner preflight re-verifies CUDA identity anyway.
Smoke protocol: --execute --max-runs 1 first (3D t01), inspect, then continue the batch.

--verify-only rebuilds all inputs in memory and asserts byte equality with disk; never touches
the gate. No solver is launched here.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/benchmark3d_r2_co'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
RUNNER = ROOT / 'scripts/run_approved_benchmark3d_r2.py'
SUPERVISOR = ROOT / 'scripts/bounded_windows_process.py'

BATCH_ID = 'benchmark3d_r2_co'
PACKET = 'BENCHMARK3D-R2-CO'
SEED = 20261002
DX = 0.05
N_TRACES = 13
RMS_M = 0.16
CL_M = 2.5
SIGMA_M = CL_M / 2.0          # convention: Gaussian kernel std = CL/2
ETA_CLIP = (-2.00, 1.95)      # z_if = 9 + eta in [7.00, 10.95], cover always >0, never touches surface
Z_IF0 = 9.0                   # nominal interface (.in coords, includes 1 m PML offset)
DOM = (12.0, 12.0, 33.0)      # .in domain incl. 1 m PML each side
PML3D = '20 20 20 20 20 20'
PML2D_5CM = '0 20 20 0 20 20'
PML2D_2P5CM = '0 40 40 0 40 40'
CFS = '#pml_cfs: constant forward 0 0 constant forward 1 1 quartic forward 0 0.21235349838321013'
TIME_WINDOW = '600e-9'
MATERIALS = ['#material: 9 0.001 1 0 rock',
             '#material: 18.017 0.003 1 0 cover',
             '#add_dispersion_debye: 1 7.878 6.4567e-09 cover']


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha_file(p: Path) -> str:
    return sha_bytes(p.read_bytes())


def make_field():
    """eta(x, y) on the full 12x12 m .in domain, 5 cm grid (240x240)."""
    n = int(round(DOM[0] / DX))
    rng = np.random.default_rng(SEED)
    white = rng.standard_normal((n, n))
    pad = int(round(4 * SIGMA_M / DX))
    w = np.pad(white, pad, mode='reflect')
    kx = np.fft.fftfreq(w.shape[0], d=DX)
    ky = np.fft.fftfreq(w.shape[1], d=DX)
    kx2, ky2 = np.meshgrid(kx, ky, indexing='ij')
    g = np.exp(-2.0 * np.pi ** 2 * SIGMA_M ** 2 * (kx2 ** 2 + ky2 ** 2))
    filt = np.real(np.fft.ifft2(np.fft.fft2(w) * g))[pad:pad + n, pad:pad + n]
    filt = filt - filt.mean()  # 零均值约定：强平滑后 DC 相对残余方差不可忽略，显式去除
    eta = filt / filt.std() * RMS_M
    eta = np.clip(eta, *ETA_CLIP)
    eta = np.round(eta / DX) * DX
    eta = np.clip(eta, *ETA_CLIP)
    assert eta.shape == (n, n) and abs(float(eta.mean())) < 0.02
    return eta


def trace_ys(k):
    tx = 3.85 + 0.25 * k
    return tx, tx + 1.3


def bin_z0(eta, i0, i1, j0, j1):
    """Mean interface height of a bin, re-quantized to the grid."""
    z = Z_IF0 + eta[i0:i1, j0:j1].mean()
    return round(round(z / DX) * DX, 2)


def fmt(v):
    s = f'{v:g}'
    return s


def header(title, domain_mode, domain, dx, pml):
    lines = [f'#title: {title}']
    if domain_mode:
        lines.append(f'#domain_mode: {domain_mode}')
    lines += [f'#domain: {domain}',
              f'#dx_dy_dz: {dx} {dx} {dx}',
              f'#time_window: {TIME_WINDOW}',
              '#omp_threads: 8',
              f'#pml_cells: {pml}',
              '#pml_formulation: HORIPML',
              '#waveform: impulse 1 1 impulse']
    return lines


def build_3d(eta, k, field_sha):
    tx, rx = trace_ys(k)
    tag = f't{k + 1:02d}'
    rid = f'B3D5CM-C3mR2-BG-CO13-{tag}'
    lines = header(f'{rid} benchmark3d_r2_co rough-interface 3D standard benchmark'
                   f' (RMS {RMS_M} CL {CL_M} seed {SEED} field {field_sha[:12]})',
                   '3D', '12 12 33', fmt(DX), PML3D)
    lines += [f'#hertzian_dipole: x 6 {fmt(tx)} 27 impulse',
              f'#rx: 6 {fmt(rx)} 27 {tag} Ex',
              CFS] + MATERIALS
    lines.append('#box: 1 1 1 11 11 12 rock')
    nb = int(round(10.0 / 0.25))
    for i in range(nb):
        for j in range(nb):
            x0 = 1.0 + 0.25 * i
            y0 = 1.0 + 0.25 * j
            z0 = bin_z0(eta, 20 + 5 * i, 20 + 5 * i + 5, 20 + 5 * j, 20 + 5 * j + 5)
            lines.append(f'#box: {fmt(x0)} {fmt(y0)} {fmt(z0)} {fmt(x0 + 0.25)} {fmt(y0 + 0.25)} 12 cover')
    return rid, ('\n'.join(lines) + '\n').encode('utf-8')


def build_2d(eta, k, dx, pml, label, field_sha):
    tx, rx = trace_ys(k)
    tag = f't{k + 1:02d}'
    rid = f'{label}-C3mR2-BG-CO13-{tag}'
    lines = header(f'{rid} benchmark3d_r2_co 2D pair of the 3D rough-interface benchmark'
                   f' (eta slice x=6.0, seed {SEED}, field {field_sha[:12]})',
                   'TM', 'inf 12 33', fmt(dx), pml)
    lines += [f'#hertzian_dipole: x inf {fmt(tx)} 27 impulse',
              f'#rx: inf {fmt(rx)} 27 {tag} Ex',
              CFS] + MATERIALS
    lines.append('#box: inf 1 1 inf 11 12 rock')
    ix = int(round(6.0 / DX))
    nb = int(round(10.0 / 0.25))
    for j in range(nb):
        y0 = 1.0 + 0.25 * j
        z0 = bin_z0(eta[ix:ix + 1, :], 0, 1, 20 + 5 * j, 20 + 5 * j + 5)
        lines.append(f'#box: inf {fmt(y0)} {fmt(z0)} inf {fmt(y0 + 0.25)} 12 cover')
    return rid, ('\n'.join(lines) + '\n').encode('utf-8')


def build_all(field_sha):
    eta = make_field()
    runs = []
    for k in range(N_TRACES):
        runs.append(build_3d(eta, k, field_sha))
    for k in range(N_TRACES):
        runs.append(build_2d(eta, k, DX, PML2D_5CM, 'B2D5CM', field_sha))
    for k in range(N_TRACES):
        runs.append(build_2d(eta, k, 0.025, PML2D_2P5CM, 'B2D', field_sha))
    return eta, runs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    eta = make_field()
    npz_path = OUT / 'interface_field.npz'
    if not a.verify_only:
        np.savez_compressed(npz_path, eta=eta)
    field_sha = sha_file(npz_path) if npz_path.exists() else None
    if field_sha is None and not a.verify_only:
        raise SystemExit('field npz missing after write')
    if a.verify_only:
        # rebuild expected field bytes for comparison via saved sha in sample table
        table = json.loads((OUT / 'sample_table.json').read_text(encoding='utf-8'))
        field_sha = table['field_sha256']
        assert sha_file(npz_path) == field_sha, 'interface field npz changed'

    eta2, runs = build_all(field_sha)
    assert np.array_equal(eta, eta2)

    if a.verify_only:
        for rid, content in runs:
            disk = (OUT / f'{rid}.in').read_bytes()
            assert disk == content, f'{rid} input changed'
        print(f'verify-only OK: {len(runs)} inputs byte-identical, field sha unchanged')
        return

    for rid, content in runs:
        (OUT / f'{rid}.in').write_bytes(content)

    table = {
        'schema': 'benchmark3d-r2-co-sample-table/1',
        'date': '2026-10-01',
        'batch_id': BATCH_ID,
        'seed': SEED,
        'field': {'kind': 'gaussian-correlated zero-mean random interface',
                  'rms_m': RMS_M, 'correlation_length_m': CL_M,
                  'kernel_sigma_m': SIGMA_M, 'sigma_convention': 'Gaussian kernel std = CL/2',
                  'clip_m': list(ETA_CLIP), 'quantize_m': DX,
                  'grid_m': DX, 'domain_m': [DOM[0], DOM[1]],
                  'filter': 'reflect-padded FFT Gaussian (numpy-only)',
                  'z_if_nominal_in_coords_m': Z_IF0},
        'field_file': 'interface_field.npz',
        'field_sha256': field_sha,
        'field_note': '2D pairs use the slice eta(x=6.0, y); 2D/3D differ only by dimension',
        'traces': [{'k': k, 'tx_y': trace_ys(k)[0], 'rx_y': trace_ys(k)[1],
                    'x': 6.0, 'z': 27.0} for k in range(N_TRACES)],
        'materials': MATERIALS,
        'pml_cfs': CFS,
        'time_window_s': TIME_WINDOW,
        'runs': [rid for rid, _ in runs],
        'basis': 'docs/research/2026-10-01_3d_benchmark_design_draft.md; user 确认 2026-10-01 22:23 (S1-S5)',
    }
    (OUT / 'sample_table.json').write_text(json.dumps(table, ensure_ascii=False, indent=1), encoding='utf-8')

    launcher_sha = sha_file(RUNNER)
    supervisor_sha = sha_file(SUPERVISOR)
    runtime_sha = sha_file(RUNTIME_ID)
    cuda_sha = sha_file(CUDA_ID)
    contracts = []
    for rid, content in runs:
        contracts.append({
            'packet_id': PACKET,
            'run_id': rid,
            'input_path': f'configs/research/{BATCH_ID}/{rid}.in',
            'input_sha256': sha_bytes(content),
            'runtime_identity_path': str(RUNTIME_ID.relative_to(ROOT)),
            'runtime_identity_sha256': runtime_sha,
            'cuda_identity_path': str(CUDA_ID.relative_to(ROOT)),
            'cuda_identity_sha256': cuda_sha,
            'launcher_sha256': launcher_sha,
            'supervisor_sha256': supervisor_sha,
            'attempt_record': f'artifacts/research_checks/2026-10-01_{rid}_attempt.json',
            'run_directory': f'artifacts/simulations/2026-10-01_{rid}',
            'continuation': 'Serial in frozen order (3D t01..t13, then 2D-5cm, then 2D-2.5cm); '
                            'any failure aborts the batch. No retries. Smoke: --max-runs 1 first.',
        })

    gate = json.loads(GATE.read_text(encoding='utf-8'))
    gate.update({
        'updated_date': '2026-10-01',
        'batch_id': BATCH_ID,
        'approved_to_simulate': True,
        'approved_run_ids': [c['run_id'] for c in contracts],
        'approved_run_ids_pending_agreement': [],
        'approved_compute_budget': {
            'wall_minutes': 30,
            'job_commit_GiB': 20,
            'minimum_available_RAM_GiB': 8,
            'output_GiB': 1,
            'threads': 8,
            'backend': 'CUDA',
            'precision': 'double',
            'retries': 0,
            'max_fdtd_runs': len(contracts),
            'official_postprocess_wall_minutes': 10,
            'device_id': 0,
            'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'approved_execution_contracts': contracts,
        'execution_outcome': {
            'status': 'approved_pending_execution',
            'note': 'user 确认 2026-10-01 22:23 (S1-S5); BLOCKED on OS reboot: NVML init fails and '
                    'CBS+WU RebootPending flags are both set (diagnosed ~22:30); runner preflight '
                    're-verifies CUDA identity and device name before any run',
        },
    })
    basis = ('benchmark3d_r2_co: user 2026-10-01 22:23 "确认" on design draft '
             'docs/research/2026-10-01_3d_benchmark_design_draft.md S1-S5 '
             '(first 3D standard benchmark; rough bedrock-cover interface RMS 0.16 m / CL 2.5 m; '
             '13+13+13 runs; flat-3D control deferred); premises signed 21:54 '
             '(flat ground surface; reduced domain accepted)')
    gate.setdefault('scope_expansion_basis', []).append(basis)
    GATE.write_text(json.dumps(gate, ensure_ascii=False, indent=1), encoding='utf-8')
    print(f'froze {BATCH_ID}: {len(runs)} runs, gate approved_to_simulate=TRUE '
          f'(execution blocked on reboot), field sha {field_sha[:16]}...')


if __name__ == '__main__':
    main()
