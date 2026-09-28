"""Freeze the batch3d_slope_t3_xcheck segment: 3D slope model vs 2D-5cm (2D/3D consistency).

Basis: user 2026-09-28 — "做3d". The t3 family became the batch-2D spec baseline
(decision_log 2026-09-28); whether cheap 2D-TM runs can stand in for 3D physics in
batch generation requires a direct 2D/3D consistency check on the SAME geometry,
SAME discretisation and SAME materials, so that the only remaining difference is
dimensionality (line source + invariant x vs point dipole + 3D spreading).

Design:
- 3D slope BG mother B3D5CM-C3mS2X-BG: domain 8 x 40 x 50 m @ 5 cm (128M cells),
  x-centred Hertzian dipole (4, 15.35, 45) / rx (4, 16.65, 45) = the t17 anchor pair;
  surface z=30; interface staircase follows the t3 smooth line z = 0.2y + 22.625
  clamped at the surface (natural outcrop), realised on the 5 cm grid with
  round-half-up per 0.25 m strip (max deviation from the line 0.025 m = half cell,
  documented in static_check); cover is the dispersive Debye of
  dispersion_materials_v0.1; rock non-dispersive; PML 20 cells (1.0 m, same physical
  thickness as the 2D t3 40-cell PML); HORIPML + same pml_cfs as archived 3D anchors.
- 2D-5cm slope BG mother B2D5CM-C3mS2X-BG + 33 CO traces: TM mode, domain inf x 40 x 50
  @ 5 cm, SAME strip realisation (same generator function => identical geometry),
  same materials, anchor t17 pair (Tx 15.35 / Rx 16.65, z=45) and the t3 CO array
  (Tx = Rx - 1.30 m, Rx y = 12.65..20.65 step 0.25).
- Comparisons enabled: (a) 3D-5cm vs 2D-5cm at the shared anchor isolates
  dimensionality; (b) 2D-5cm vs archived 2D-2.5cm t3 (batch2d_slope_t3_co) isolates
  grid/quantisation; together they bound how well batch 2D-2.5cm approximates 3D.

Cost basis (measured, archived): DEP3D_5CM_BG 32M cells -> 916 s wall, 8.57 GB peak
job commit (268 B/cell host). 128M cells => ~3660-4400 s wall (~61-73 min, Debye ADE
adds one polarisation accumulator), ~36 GB host commit. Machine: 63.2 GB RAM total
(~43 GB free at freeze), 15.8 GB VRAM free; field arrays ~6.1 GB + 1 Debye pole
accumulator ~1 GB. 2D-5cm CO traces ~5 s each.

Mirrors scripts/freeze_batch2d_slope_t3_co.py: generates mothers + inputs,
cases/groups contracts, static_check.json, budget.json, the runner clone, and
rewrites the global execution gate with approved_to_simulate=false. --verify-only
rebuilds in memory and asserts byte equality with disk without touching the gate.
No solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/batch3d_slope_t3_xcheck'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
DISP = json.loads((ROOT / 'configs/research/dispersion_materials_v0.1.json').read_text(encoding='utf-8'))

C = 299792458.0
C_M_NS = 0.299792458
D = 0.05                      # 5 cm grid (archived 3D convention)
SURF_Z = 30.0
STRIP_W = 0.25
Y1 = 40.0                     # slope domain y-extent (t3 outcrop geometry)
X1 = 8.0                      # 3D cross-line extent (archived 3D convention)
DOM_Z = 50.0
PML_CELLS_3D = 20             # 1.0 m at 5 cm (archived 3D convention)
PML_CELLS_2D = 20
TIME_WINDOW_S = 1200e-9
N_TRACES = 33
DY_RX = 0.25
Y_FIRST = 12.65
OFFSET = 1.30
Z_LINE = 45.0
ANCHOR_INDEX = 16
BATCH_ID = 'batch3d_slope_t3_xcheck'
RUN_DATE = '2026-09-28'

COVER = DISP['materials']['cover_clay']
assert not DISP['materials']['bedrock_sandstone']['dispersive']
COVER_MAT = f"#material: {COVER['epsilon_inf']} {COVER['sigma_dc']} 1 0 cover"
COVER_POLE = (f"#add_dispersion_debye: 1 {COVER['debye_poles'][0]['delta_epsilon']} "
              f"{COVER['debye_poles'][0]['tau_s']} cover")
ROCK_MAT = '#material: 9 0.001 1 0 rock'
PML_CFS = '#pml_cfs: constant forward 0 0 constant forward 1 1 quartic forward 0 0.21235349838321013'

M3 = 'B3D5CM-C3mS2X-BG'
M2 = 'B2D5CM-C3mS2X-BG'


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(run_id: str) -> int:
    return int.from_bytes(hashlib.sha256(f'{BATCH_ID}|{run_id}'.encode('utf-8')).digest()[:8], 'big')


def strips_5cm():
    """Interface staircase realised on the 5 cm grid from the t3 smooth line.

    zif_c = floor(z_line(ycen)/d + 0.5) (round-half-up, deterministic); strips whose
    interface reaches the surface are all-rock (natural outcrop). Returns
    (strips, max_dev_m, outcrop_interval)."""
    strips = []
    max_dev = 0.0
    n = round(Y1 / STRIP_W)
    ycrop = None
    for i in range(n):
        ya = i * STRIP_W
        yb = ya + STRIP_W
        ycen = ya + STRIP_W / 2
        z_line = min(0.2 * ycen + 22.625, SURF_Z)
        zif_c = math.floor(z_line / D + 0.5)
        zif = zif_c * D
        max_dev = max(max_dev, abs(zif - z_line))
        if zif >= SURF_Z - 1e-12:
            if ycrop is None:
                ycrop = (ya, yb)
            continue
        strips.append((ya, zif, yb))
    assert ycrop is not None
    # continuity + monotonicity
    for i in range(1, len(strips)):
        assert strips[i][0] == strips[i - 1][2]
        assert strips[i][1] >= strips[i - 1][1]
    assert abs(strips[-1][2] - ycrop[0]) < 1e-12
    return strips, max_dev, ycrop


def fermat_margin_assert():
    """Assert every CO trace's stationary interface-reflection point sits on the
    physical interface with >= 2 m margin to the PML start (39 m), at the dispersive
    band-centre velocity (er=18) and er=16."""
    ys_rx = np.array([Y_FIRST + k * DY_RX for k in range(N_TRACES)])
    yq_all = []
    for er in (16.0, 18.0):
        v2 = C_M_NS / math.sqrt(er)
        for rx_y in ys_rx:
            tx_y = rx_y - OFFSET
            y1g = np.linspace(rx_y - 8, rx_y + 14, 60)
            yqg = np.linspace(12.0, 36.9, 300)
            y2g = np.linspace(rx_y - 6, rx_y + 16, 60)
            Y1g, YQ, Y2g = np.meshgrid(y1g, yqg, y2g, indexing='ij')
            ZQ = np.minimum(0.2 * YQ + 22.625, SURF_Z)
            t = ((np.hypot(Y1g - tx_y, 15.0) + np.hypot(Y2g - rx_y, 15.0)) / C_M_NS
                 + (np.hypot(YQ - Y1g, SURF_Z - ZQ) + np.hypot(Y2g - YQ, SURF_Z - ZQ)) / v2)
            yq_all.append(float(YQ[np.unravel_index(int(np.argmin(t)), t.shape)]))
    yq_star = max(yq_all)
    pml_start = Y1 - PML_CELLS_2D * D
    assert yq_star + 2.0 <= pml_start, (yq_star, pml_start)
    return yq_star, pml_start


def grid_line_asserts(text: str):
    for l in text.splitlines():
        if l.startswith('#box:'):
            for v in l.split()[1:7]:
                if v != 'inf':
                    assert abs(float(v) / D - round(float(v) / D)) < 1e-9, l
    text.encode('ascii')  # gprMax reads inputs with the system codec (GBK here); ASCII-only


def build_mothers(strips):
    cover_3d = ['#box: 0 %.6g %.6g %.6g %.6g %.6g cover' % (ya, zif, X1, yb, SURF_Z)
                for ya, zif, yb in strips]
    cover_2d = ['#box: inf %.6g %.6g inf %.6g %.6g cover' % (ya, zif, yb, SURF_Z)
                for ya, zif, yb in strips]
    m3 = '\n'.join([
        f'#title: {M3} {BATCH_ID} 3D slope BG 5cm anchor t17 (t3 outcrop geometry, dispersive cover)',
        f'#domain: {X1:g} {Y1:g} {DOM_Z:g}',
        f'#dx_dy_dz: {D:g} {D:g} {D:g}',
        '#time_window: 1200e-9',
        '#omp_threads: 8',
        f'#pml_cells: {PML_CELLS_3D} {PML_CELLS_3D} {PML_CELLS_3D} {PML_CELLS_3D} {PML_CELLS_3D} {PML_CELLS_3D}',
        '#pml_formulation: HORIPML',
        '#waveform: impulse 1 1 impulse',
        '#hertzian_dipole: x 4 15.35 45 impulse',
        '#rx: 4 16.65 45 measurement Ex',
        PML_CFS,
        ROCK_MAT,
        f'#box: 0 0 0 {X1:g} {Y1:g} 30 rock',
        COVER_MAT,
        COVER_POLE,
        *cover_3d,
    ]) + '\n'
    m2 = '\n'.join([
        f'#title: {M2} {BATCH_ID} 2D-TM slope BG 5cm (same 5cm strips as {M3}, dispersive cover)',
        '#domain_mode: TM',
        f'#domain: inf {Y1:g} {DOM_Z:g}',
        f'#dx_dy_dz: {D:g} {D:g} {D:g}',
        '#time_window: 1200e-9',
        '#omp_threads: 8',
        f'#pml_cells: 0 {PML_CELLS_2D} {PML_CELLS_2D} 0 {PML_CELLS_2D} {PML_CELLS_2D}',
        '#pml_formulation: HORIPML',
        '#waveform: impulse 1 1 impulse',
        '#hertzian_dipole: x inf 15.35 45 impulse',
        '#rx: inf 16.65 45 measurement Ex',
        PML_CFS,
        ROCK_MAT,
        f'#box: inf 0 0 inf {Y1:g} 30 rock',
        COVER_MAT,
        COVER_POLE,
        *cover_2d,
    ]) + '\n'
    for t in (m3, m2):
        grid_line_asserts(t)
    return {M3: m3, M2: m2}


def build_co_input(mother_text: str, run_id: str, tx_y: float, rx_y: float, rx_id: str) -> str:
    lines = mother_text.splitlines()
    dip_idx = [i for i, l in enumerate(lines) if l.startswith('#hertzian_dipole:')]
    rx_idx = [i for i, l in enumerate(lines) if l.startswith('#rx:')]
    title_idx = [i for i, l in enumerate(lines) if l.startswith('#title:')]
    assert len(dip_idx) == 1 and len(rx_idx) == 1 and len(title_idx) == 1
    assert lines[dip_idx[0]] == '#hertzian_dipole: x inf 15.35 45 impulse'
    new_lines = list(lines)
    new_lines[title_idx[0]] = f'#title: {run_id} {BATCH_ID} common-offset pair offset {OFFSET} m (mother {M2})'
    new_lines[dip_idx[0]] = f'#hertzian_dipole: x inf {tx_y:.2f} 45 impulse'
    new_lines[rx_idx[0]] = f'#rx: inf {rx_y:.2f} 45 {rx_id} Ex'
    text = '\n'.join(new_lines) + '\n'
    text.encode('ascii')
    return text


def collect():
    strips, max_dev, ycrop = strips_5cm()
    yq_star, pml_start = fermat_margin_assert()

    nx3, ny3, nz3 = int(X1 / D), int(Y1 / D), int(DOM_Z / D)
    cells3 = nx3 * ny3 * nz3
    ny2 = int(Y1 / D)
    cells2 = ny2 * nz3
    assert cells3 == 128000000 and cells2 == 800000
    dt3 = 1.0 / (C * math.sqrt(3 * D ** -2))
    dt2 = 1.0 / (C * math.sqrt(2 * D ** -2))

    mothers = build_mothers(strips)
    inputs, cases_out, groups_out = {}, [], []

    # 2D-5cm CO traces first (fast, de-risks materials/geometry), 3D anchor last
    m2 = mothers[M2]
    for k in range(N_TRACES):
        rx_y = Y_FIRST + k * DY_RX
        tx_y = rx_y - OFFSET
        assert abs(rx_y / D - round(rx_y / D)) < 1e-9 and abs(tx_y / D - round(tx_y / D)) < 1e-9
        assert 1.0 < tx_y and rx_y < Y1 - 1.0
        rx_id = f't{k + 1:02d}'
        run_id = f'{M2}-CO33-{rx_id}'
        text = build_co_input(m2, run_id, tx_y, rx_y, rx_id)
        inputs[run_id] = text
        strip = lambda s: [l for l in s.splitlines()
                           if not l.startswith(('#rx:', '#title:', '#hertzian_dipole:'))]
        assert strip(m2) == strip(text), run_id
        if k == ANCHOR_INDEX:
            assert '#rx: inf 16.65 45 t17 Ex' in text
        cases_out.append({
            'run_id': run_id, 'file': run_id + '.in',
            'input_sha256': sha256_bytes(text.encode('utf-8')),
            'mother_model_id': M2, 'case_id': f'{M2}|{rx_id}', 'group_id': M2,
            'segment': 'co2d5cm', 'variant_tag': 'common_offset_33tx_0p25m_off1p30_5cm',
            'seed': seed_for(run_id),
            'common_offset': {'n_traces': N_TRACES, 'trace_spacing_m': DY_RX, 'offset_m': OFFSET,
                              'rx_y_m': rx_y, 'tx_y_m': tx_y, 'z_m': Z_LINE,
                              'trace_index': k, 'anchor_trace': k == ANCHOR_INDEX,
                              'rx_output_components': ['Ex']},
        })
    groups_out.append({'case_id': M2, 'group_id': M2, 'run_id_pattern': f'{M2}-CO33-t*',
                       'n_traces': N_TRACES, 'variant_tag': 'common_offset_33tx_0p25m_off1p30_5cm'})

    run3 = f'{M3}-T17'
    inputs[run3] = mothers[M3]
    cases_out.append({
        'run_id': run3, 'file': run3 + '.in',
        'input_sha256': sha256_bytes(mothers[M3].encode('utf-8')),
        'mother_model_id': M3, 'case_id': f'{M3}|t17', 'group_id': M3,
        'segment': 'anchor3d', 'variant_tag': 'true3d_5cm_anchor_t17',
        'seed': seed_for(run3),
        'anchor': {'tx': [4.0, 15.35, 45.0], 'rx': [4.0, 16.65, 45.0],
                   'note': 'same y/z pair as the 2D t17 anchor; x-centred point dipole'},
        'rx_output_components': ['Ex'],
    })
    groups_out.append({'case_id': M3, 'group_id': M3, 'run_id_pattern': f'{M3}-T*',
                       'n_traces': 1, 'variant_tag': 'true3d_5cm_anchor_t17'})

    static = {
        'schema': 'batch3d_slope_t3_xcheck_static_check/1',
        'checked_date': RUN_DATE,
        'question': 'does the batch-2D (TM, line source) slope model reproduce the true-3D '
                    '(point dipole, 3D spreading) interface response on identical geometry, '
                    'grid and materials — i.e. can t3 2D stand in for 3D in batch generation',
        'design': 'same smooth interface line (t3: z=0.2y+22.625 clamped at outcrop) realised '
                  'by one generator for both mothers; identical 5 cm grid, PML thickness 1.0 m, '
                  'HORIPML + archived pml_cfs, dispersive Debye cover (dispersion_materials_v0.1), '
                  'non-dispersive rock; only dimensionality differs',
        'grid_3d': {'nx': nx3, 'ny': ny3, 'nz': nz3, 'cells': cells3, 'd_m': D,
                    'dt_formula_s': dt3, 'time_steps_formula': round(TIME_WINDOW_S / dt3)},
        'grid_2d': {'ny': ny2, 'nz': nz3, 'cells': cells2, 'd_m': D,
                    'dt_formula_s': dt2, 'time_steps_formula': round(TIME_WINDOW_S / dt2)},
        'staircase': {'rule': 'zif_c = floor((0.2*ycen+22.625)/0.05 + 0.5), strips 0.25 m, '
                              'clamp at surface = natural outcrop',
                      'n_cover_strips': len(strips),
                      'max_deviation_from_line_m': round(max_dev, 6),
                      'outcrop_y_interval_m': list(ycrop)},
        'fermat_domain_sizing': {'yq_max_er16_er18_m': round(yq_star, 3),
                                 'margin_rule_m': 2.0, 'pml_start_m': pml_start,
                                 'assertion': 'max(yq) + 2.0 <= %.1f holds (%.3f)' % (pml_start, yq_star)},
        'dispersion_config': 'configs/research/dispersion_materials_v0.1.json',
        'dispersion_config_sha256': sha256_bytes(
            (ROOT / 'configs/research/dispersion_materials_v0.1.json').read_bytes()),
        'comparisons': {
            'a_dimensionality': f'{run3} vs {M2}-CO33-t17 (same grid, same materials, same strips)',
            'b_grid': f'{M2}-CO33-t* vs archived batch2d_slope_t3_co B2D-C3mS2X-BG-CO33-t* '
                      '(2.5 cm vs 5 cm, same line; staircase quantisation differs by design)'},
        'cost_basis': {
            'archived_dep3d_5cm': {'cells': 32000000, 'wall_s': 916.078,
                                   'peak_job_commit_GB': 8.57},
            'scaled_128M': {'wall_s_estimate': '3660-4400', 'host_commit_GB_estimate': '~36',
                            'note': 'linear in cells; Debye ADE adds one pole accumulator'},
            'machine': {'RAM_total_GB': 63.2, 'RAM_free_at_freeze_GB': 43.0,
                        'VRAM_free_at_freeze_GB': 15.8}},
        'mother_sha256': {m: sha256_bytes(t.encode('utf-8')) for m, t in mothers.items()},
    }

    budget = {
        'schema': 'batch3d_slope_t3_xcheck_budget/1',
        'date': RUN_DATE,
        'grid': {'d_m': D, 'cells_3d': cells3, 'cells_2d': cells2,
                 'time_window_ns': 1200},
        'single_case': {'wall_s_2d_estimate': 10, 'wall_s_3d_estimate': '3660-4400',
                        'host_commit_GB_3d_estimate': 36},
        'batch_totals': {'n_cases': N_TRACES + 1,
                         'measured_estimate_s': N_TRACES * 10 + 4400,
                         'conservative_budget_min': 150},
        'hard_stops': {'wall_minutes_per_case': 150, 'job_commit_GiB': 48,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': N_TRACES + 1},
        'known_limitations': [
            '3D estimate scaled from archived DEP3D_5CM_BG (32M cells, 916 s, 8.57 GB commit); '
            'Debye ADE overhead not included in the archived figure',
            'job_commit cap 48 GiB vs 63.2 GB total RAM; preflight RAM gate 40 GiB',
            'estimates are for cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': 'xcheck_2d3d', 'grid_tier': 'BASE5CM',
                 'spec': 'decision_log 2026-09-28 (t3 baseline + user "做3d"); geometry t3 line; '
                         'materials dispersion_materials_v0.1',
                 'n_cases': N_TRACES + 1, 'n_exceptions': 0,
                 'derivation': 'both mothers generated by one scripted strip realisation; '
                               '2D CO traces byte-identical to the 2D mother except '
                               '#title/#hertzian_dipole/#rx',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'xcheck_2d3d', 'grid_tier': 'BASE5CM',
                  'convention': 'spec 4.2: any split is by group_id (mother model) only',
                  'groups': groups_out}
    return mothers, inputs, cases_doc, groups_doc, budget, static


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    mothers, inputs, cases_doc, groups_doc, budget, static = collect()

    if args.verify_only:
        for run_id, text in inputs.items():
            assert (OUT / (run_id + '.in')).read_text(encoding='utf-8') == text, run_id
        assert json.loads((OUT / 'cases.json').read_text(encoding='utf-8')) == cases_doc
        assert json.loads((OUT / 'groups.json').read_text(encoding='utf-8')) == groups_doc
        assert json.loads((OUT / 'budget.json').read_text(encoding='utf-8')) == budget
        assert json.loads((OUT / 'static_check.json').read_text(encoding='utf-8')) == static
        print(f'verify-only: all {len(inputs)} xcheck inputs + cases/groups/budget/static_check byte-identical')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_batch3d_slope_t3_xcheck.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH3D-SLOPE-T3-XCHECK',
        'run_id': rid,
        'input_path': f'configs/research/batch3d_slope_t3_xcheck/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (2D-5cm CO t01..t33, then the 3D anchor); '
                        'any failure aborts the batch, later attempts untouched. No retries.',
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
            'wall_minutes': 150, 'job_commit_GiB': 48, 'minimum_available_RAM_GiB': 40,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': N_TRACES + 1, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 12,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-28: "做3d" — build the true-3D slope model for the 2D/3D consistency '
            'check that gates batch-scale 2D generation (t3 family is the batch-2D baseline, '
            'decision_log 2026-09-28)',
            'geometry: t3 smooth interface line realised at 5 cm by one scripted generator for '
            'both the 3D and 2D mothers (identical strips); 3D conventions mirror the archived '
            'a0_3d_v1/dep3d anchors (x extent 8 m, PML 20 cells, x-centred dipole)',
            'cost basis: archived DEP3D_5CM_BG 32M cells 916 s / 8.57 GB commit; 128M cells '
            '~61-73 min / ~36 GB; machine 63.2 GB RAM, 15.8 GB VRAM free; 34 cases conservative '
            '150 min total',
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
        'note': 'batch3d_slope_t3_xcheck frozen 2026-09-28 (33 2D-5cm CO traces + 1 true-3D 5cm '
                'anchor, identical strips/materials; 2D first, 3D last). approved_to_simulate '
                'stays false until the recorded basis is applied; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
