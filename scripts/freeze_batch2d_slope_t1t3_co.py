"""Freeze the batch2d_slope_t1t3 common-offset (CO) segment: T1/T3 slope tiers.

Basis: user 2026-09-29 "按建议开始" approving docs/research/2026-09-29_t1_t3_extension_proposal.md
decision points per recommendation: both T1 and T3 tiers in one batch (CO33 segment 1 only,
MT33 deferred), T1 domain widened y 40 -> 50 m for the natural outcrop (the 40 m domain
would leave 0.6 m cover at the high edge, repeating the t2 truncation invalidation), TGT
keeps the anchor target box at fixed z (single-variable principle), reference-window
construction NOT included in this batch.

Geometry (self-computed, asserted here):
- T3: tan=0.4, zc=22.6 m anchored at y=16 (t2 rule, half-span 16 m), domain y 40 m
  (unchanged), outcrop at y in [34.5, 34.75), 4.5 m before the PML inner edge (39 m).
- T1: tan=0.1, zc=27.0 m anchored at y=16, domain y 50 m, outcrop at y in [46.0, 46.25),
  2.75 m before the PML inner edge (49 m).
- Staircase rule unchanged: 0.25 m strips, grid-locked z_if = round to cell of
  zc + (ycen - 16) * tan; strips stop at the first strip whose z_if >= surface (30 m);
  beyond the outcrop bedrock reaches the surface (t3 machinery).

Materials: byte-identical to the t3 baseline (dispersion_materials_v0.1.json, hash-locked
below): dispersive Debye cover (18.017 / dE 7.878 / tau 6.4567e-9 / sigma 0.003), dispersive
tzone half-strength pole (12.5 / dE 4.0 / same tau / sigma 0.003), bedrock non-dispersive
(9 / 0.001). NOTE: the extension proposal's TZ text quoted the old draft's interpolation
assumption (12 / 0.005); the frozen dispersion config's tzone governs, as in batch2d_slope_t3_co.

Mothers (12 = 2 tiers x {no-TZ, TZ} x {BG, NC, TGT}), derived from the archived t2 mothers
by scripted line transforms (domain/rock-box resize, full staircase replacement, dispersive
material swap; NC nullcontrast and TGT target lines pass through unchanged):
  configs/research/batch2d_slope_t1t3/B2D-C3m{S1X,S1TZX,S3X,S3TZX}-{BG,NC,D10m-W4m-T0.5m-E20-S0.02}.in
CO array identical to batch2d_slope_t3_co (Tx = Rx - 1.30 m, z=45, Rx y 12.65..20.65 step
0.25, 33 traces, anchor t17 == mother single-trace pair). 12 x 33 = 396 FDTD runs.

NC semantics (unchanged from t2/v1 lineage): NC = BG geometry + target-shaped nullcontrast
box (eps 9 / sigma 0.001 = rock); expected bit-identical to its BG twin (verified post-run,
not assumed). TGT = BG + target box (20 / 0.02) at the anchor position (fixed z).

Mirrors scripts/freeze_batch2d_slope_t3_co.py. --verify-only rebuilds in memory and asserts
byte equality with disk without touching the gate. No solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'configs/research/batch2d_slope_t2'
OUT = ROOT / 'configs/research/batch2d_slope_t1t3_co'
MOTHER_DIR = ROOT / 'configs/research/batch2d_slope_t1t3'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
DISP_PATH = ROOT / 'configs/research/dispersion_materials_v0.1.json'
DISP = json.loads(DISP_PATH.read_text(encoding='utf-8'))

C = 299792458.0
C_M_NS = 0.299792458
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
N_TRACES = 33
DY_RX = 0.25
Y_FIRST = 12.65
OFFSET = 1.30
Z_LINE = 45.0
ANCHOR_INDEX = 16
TIME_WINDOW_S = 1200e-9
STRIP_W = 0.25
SURF_Z = 30.0
YMID = 16.0
VARIANT = 'common_offset_33tx_0p25m_off1p30'
BATCH_ID = 'batch2d_slope_t1t3_co'
RUN_DATE = '2026-09-29'

TIERS = {
    # tag: (tan_theta, zc_m, domain_y_m)
    'S1X': (0.1, 27.0, 50.0),
    'S1TZX': (0.1, 27.0, 50.0),
    'S3X': (0.4, 22.6, 40.0),
    'S3TZX': (0.4, 22.6, 40.0),
}
SCENES = {'BG': 'BG', 'NC': 'NC', 'TGT': 'D10m-W4m-T0.5m-E20-S0.02'}
SRC_TIER = {'S1X': 'S2', 'S1TZX': 'S2TZ', 'S3X': 'S2', 'S3TZX': 'S2TZ'}

COVER = DISP['materials']['cover_clay']
TZ = DISP['materials']['transition_zone']
ROCK = DISP['materials'].get('bedrock_sandstone') or DISP['materials'].get('bedstone_sandstone')
assert ROCK and not ROCK['dispersive']
COVER_MAT = f"#material: {COVER['epsilon_inf']} {COVER['sigma_dc']} 1 0 cover"
COVER_POLE = (f"#add_dispersion_debye: 1 {COVER['debye_poles'][0]['delta_epsilon']} "
              f"{COVER['debye_poles'][0]['tau_s']} cover")
TZ_MAT = f"#material: {TZ['epsilon_inf']} {TZ['sigma_dc']} 1 0 tzone"
TZ_POLE = (f"#add_dispersion_debye: 1 {TZ['debye_poles'][0]['delta_epsilon']} "
           f"{TZ['debye_poles'][0]['tau_s']} tzone")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(mother_id: str, k: int) -> int:
    txt = f'{BATCH_ID}|{mother_id}|t{k + 1:02d}|{VARIANT}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def staircase(tan: float, zc: float, y_max: float):
    """Grid-locked cover/tzone staircase strips from y=0 to the outcrop.

    Returns (cover_lines, tz_lines, n_strips, y_outcrop_lo, theta_eff_deg,
    max_quant_dev_m, zif_centers)."""
    zc_c = round(zc / DY)
    cover_lines, tz_lines, centers, zifs = [], [], [], []
    y_out_lo = None
    y = 0.0
    while y < y_max - 1e-9:
        ycen = y + STRIP_W / 2
        zif_c = zc_c + round((ycen - YMID) * tan / DY)
        zif = zif_c * DY
        if zif >= SURF_Z - 1e-12:
            if y_out_lo is None:
                y_out_lo = y
            y += STRIP_W
            continue
        assert zif > 5.0, f'interface too deep: {zif}'
        cover_lines.append('#box: inf %.6g %.6g inf %.6g %.6g cover' % (y, zif, y + STRIP_W, SURF_Z))
        tz_lines.append('#box: inf %.6g %.6g inf %.6g %.6g tzone' % (y, zif - 1.0, y + STRIP_W, zif))
        centers.append(ycen)
        zifs.append(zif)
        y += STRIP_W
    assert y_out_lo is not None, 'interface never outcrops inside the domain'
    # realised tilt by least squares over strip centres (report, never the nominal label)
    a, b = np.polyfit(centers, zifs, 1)
    theta_eff = math.degrees(math.atan(a))
    nominal = math.degrees(math.atan(tan))
    resid = np.array(zifs) - (a * np.array(centers) + b)
    return (cover_lines, tz_lines, len(cover_lines), y_out_lo, theta_eff,
            float(np.max(np.abs(resid))), {'slope': a, 'intercept': b, 'zifs': zifs})


def fermat_stationary(tan: float, zc: float, y_max: float, er: float):
    """Stationary interface-reflection point per CO trace (t3 method, per-tier line)."""
    v1 = C_M_NS
    v2 = C_M_NS / math.sqrt(er)
    b = zc - YMID * tan
    ys_rx = np.array([Y_FIRST + k * DY_RX for k in range(N_TRACES)])
    yq = []
    yq_hi = y_max - 40 * DY - 0.15
    for rx_y in ys_rx:
        tx_y = rx_y - OFFSET
        y1g = np.linspace(rx_y - 8, rx_y + 12, 70)
        yqg = np.linspace(8.0, yq_hi, 500)
        y2g = np.linspace(rx_y - 6, rx_y + 14, 70)
        Y1, YQ, Y2 = np.meshgrid(y1g, yqg, y2g, indexing='ij')
        ZQ = np.minimum(tan * YQ + b, SURF_Z)
        t = ((np.hypot(Y1 - tx_y, 15.0) + np.hypot(Y2 - rx_y, 15.0)) / v1
             + (np.hypot(YQ - Y1, SURF_Z - ZQ) + np.hypot(Y2 - YQ, SURF_Z - ZQ)) / v2)
        yq.append(float(YQ[np.unravel_index(int(np.argmin(t)), t.shape)]))
    return yq


def build_mother(new_id: str, tier: str, scene: str, cover_lines, tz_lines) -> str:
    tan, zc, y_max = TIERS[tier]
    with_tz = tier.endswith('TZX')
    src_id = f"B2D-C3m{SRC_TIER[tier]}-{SCENES[scene]}"
    src_text = (SRC / (src_id + '.in')).read_text(encoding='utf-8')
    lines = src_text.splitlines()
    assert lines.count('#domain: inf 32 50') == 1
    assert lines.count('#box: inf 0 0 inf 32 30 rock') == 1
    assert lines.count('#material: 16 0.01 1 0 cover') == 1
    n_old_cover = sum(1 for l in lines if l.startswith('#box:') and l.endswith(' cover'))
    n_old_tz = sum(1 for l in lines if l.startswith('#box:') and l.endswith(' tzone'))
    assert n_old_cover == round(32.0 / STRIP_W) and (n_old_tz == n_old_cover if with_tz else n_old_tz == 0)

    out = []
    for l in lines:
        if l.startswith('#title:'):
            out.append(f'#title: {new_id} {BATCH_ID} mother t1t3 extension '
                       f'(theta_eff realised, outcrop, dispersive; source {src_id})')
            continue
        if l == '#domain: inf 32 50':
            out.append(f'#domain: inf {y_max:g} 50')
            continue
        if l == '#box: inf 0 0 inf 32 30 rock':
            out.append(f'#box: inf 0 0 inf {y_max:g} 30 rock')
            continue
        if l == '#material: 16 0.01 1 0 cover':
            out.append(COVER_MAT)
            out.append(COVER_POLE)
            continue
        if with_tz and l == '#material: 12 0.005 1 0 tzone':
            out.append(TZ_MAT)
            out.append(TZ_POLE)
            continue
        if l.startswith('#box:') and (l.endswith(' cover') or l.endswith(' tzone')):
            continue  # old staircase dropped, new one inserted below
        out.append(l)

    # insert the new staircase after the cover material block (tz block right after)
    cov_mat_idx = [i for i, l in enumerate(out) if l == COVER_MAT]
    assert len(cov_mat_idx) == 1
    insert_at = cov_mat_idx[0] + 2
    new_block = list(cover_lines)
    if with_tz:
        tz_mat_idx = [i for i, l in enumerate(out) if l == TZ_MAT]
        assert len(tz_mat_idx) == 1
        new_block += list(tz_lines)
    out[insert_at:insert_at] = new_block

    text = '\n'.join(out) + '\n'
    assert '#domain: inf 32 50' not in text
    assert text.count(COVER_MAT) == 1 and text.count(COVER_POLE) == 1
    if with_tz:
        assert text.count(TZ_MAT) == 1 and text.count(TZ_POLE) == 1
    assert sum(1 for l in text.splitlines() if l.startswith('#box:') and l.endswith(' cover')) == len(cover_lines)
    text.encode('ascii')
    return text


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
    text = '\n'.join(new_lines) + '\n'
    text.encode('ascii')
    return text


def collect():
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    static = {
        'schema': 'batch2d_slope_t1t3_co_static_check/1',
        'checked_date': RUN_DATE,
        'basis': 'user 2026-09-29 "按建议开始" on docs/research/2026-09-29_t1_t3_extension_proposal.md',
        'base_batches': ['batch2d_slope_t2 (geometry/material-scene lineage source)',
                         'batch2d_slope_t3_co (t3 baseline: outcrop paradigm, dispersive materials, CO array)'],
        'dispersion_config': 'configs/research/dispersion_materials_v0.1.json',
        'dispersion_config_sha256': sha256_bytes(DISP_PATH.read_bytes()),
        'tz_material_note': 'frozen dispersion config tzone (12.5/dE 4.0/sigma 0.003) governs; '
                            'proposal text quoting the old draft interpolation (12/0.005) superseded',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'time_window_s': TIME_WINDOW_S,
                 'time_steps_contract': 20352, 'dt_formula_s': dt_formula,
                 'dt_contract_s': DT_CONTRACT},
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'common_offset_array': {
            'n_traces': N_TRACES, 'spacing_m': DY_RX, 'offset_m': OFFSET,
            'rx_first_y_m': Y_FIRST, 'z_m': Z_LINE, 'anchor_index': ANCHOR_INDEX,
            'trace_ids': [f't{k + 1:02d}' for k in range(N_TRACES)],
            'anchor_note': 't17 == mother single-trace pair (Rx 16.65, Tx 15.35); '
                           'no archived baseline for the new mothers, post-run Fermat arrival check'},
        'tiers': {}, 'mother_checks': [], 'trace_checks': [],
    }

    tier_info = {}
    for tier, (tan, zc, y_max) in TIERS.items():
        cover_lines, tz_lines, n_strips, y_out, theta_eff, qdev, fit = staircase(tan, zc, y_max)
        pml_start = y_max - 40 * DY
        yq16 = fermat_stationary(tan, zc, y_max, 16.0)
        yq18 = fermat_stationary(tan, zc, y_max, 18.0)
        yq_star = max(max(yq16), max(yq18))
        assert y_out + 2.0 <= pml_start, (tier, y_out, pml_start)
        assert yq_star + 2.0 <= pml_start, (tier, yq_star, pml_start)
        ny, nz = int(y_max / DY), int(50 / DZ)
        tier_info[tier] = (cover_lines, tz_lines)
        static['tiers'][tier] = {
            'tan_theta': tan, 'zc_m': zc, 'domain_y_m': y_max, 'n_strips': n_strips,
            'outcrop_y_lo_m': y_out, 'pml_inner_y_m': pml_start,
            'outcrop_to_pml_margin_m': round(pml_start - y_out, 3),
            'theta_eff_deg': round(theta_eff, 6),
            'nominal_deg': round(math.degrees(math.atan(tan)), 6),
            'theta_fit_residual_max_m': qdev,
            'fermat_yq_max_er16_m': round(max(yq16), 3),
            'fermat_yq_max_er18_m': round(max(yq18), 3),
            'fermat_margin_assertion': f'max(yq)+2.0 <= {pml_start} (yq_star {yq_star:.3f})',
            'cells': ny * nz, 'main_fields_ID_GiB': 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30,
        }

    # TZ pair interface-geometry identity (cover z rows must match across the pair)
    for base, tzv in (('S1X', 'S1TZX'), ('S3X', 'S3TZX')):
        zb = [l.split()[3] for l in tier_info[base][0]]
        zt = [l.split()[3] for l in tier_info[tzv][0]]
        assert zb == zt, f'TZ pair staircase mismatch {base}/{tzv}'

    mothers, inputs, cases_out, groups_out = {}, {}, [], []
    for tier in ('S3X', 'S3TZX', 'S1X', 'S1TZX'):  # execution order: cheaper 40 m domain first
        for scene in ('BG', 'NC', 'TGT'):
            new_id = f'B2D-C3m{tier}-{SCENES[scene]}'
            mother_text = build_mother(new_id, tier, scene, *tier_info[tier])
            mothers[new_id] = mother_text
            static['mother_checks'].append({
                'mother': new_id, 'tier': tier, 'scene': scene,
                'source_mother': f"B2D-C3m{SRC_TIER[tier]}-{SCENES[scene]}",
                'source_sha256': sha256_bytes((SRC / f"B2D-C3m{SRC_TIER[tier]}-{SCENES[scene]}.in").read_bytes()),
                'mother_sha256': sha256_bytes(mother_text.encode('utf-8')),
                'with_tz': tier.endswith('TZX'), 'dispersive_cover': True})
            for k in range(N_TRACES):
                rx_y = Y_FIRST + k * DY_RX
                tx_y = rx_y - OFFSET
                rx_cell, tx_cell = rx_y / DY, tx_y / DY
                assert abs(rx_cell - round(rx_cell)) < 1e-9 and abs(tx_cell - round(tx_cell)) < 1e-9
                lo = 1.0
                assert lo < tx_y and rx_y < min(TIERS[tier][2] - 1.0, 31.0) or rx_y < 31.0, 'inside PML'
                rx_id = f't{k + 1:02d}'
                run_id = f'{new_id}-CO33-{rx_id}'
                text = build_co_input(mother_text, run_id, new_id, tx_y, rx_y, rx_id)
                inputs[run_id] = text
                strip = lambda s: [l for l in s.splitlines()
                                   if not l.startswith(('#rx:', '#title:', '#hertzian_dipole:'))]
                assert strip(mother_text) == strip(text), run_id
                static['trace_checks'].append({'run_id': run_id, 'rx_cell': int(round(rx_cell)),
                                               'tx_cell': int(round(tx_cell)), 'outside_pml': True})
                cases_out.append({
                    'run_id': run_id, 'file': run_id + '.in',
                    'input_sha256': sha256_bytes(text.encode('utf-8')),
                    'mother_model_id': new_id, 'case_id': f'{new_id}|{rx_id}',
                    'group_id': new_id, 'segment': 'co', 'variant_tag': VARIANT,
                    'seed': seed_for(new_id, k),
                    'common_offset': {
                        'n_traces': N_TRACES, 'trace_spacing_m': DY_RX, 'offset_m': OFFSET,
                        'rx_y_m': rx_y, 'tx_y_m': tx_y, 'z_m': Z_LINE,
                        'trace_index': k, 'anchor_trace': k == ANCHOR_INDEX,
                        'acquisition_geometry': 'common-offset profile along y (Tx = Rx - 1.30 m, both z=45)',
                        'rx_output_components': ['Ex']}})
            groups_out.append({'case_id': new_id, 'group_id': new_id,
                               'run_id_pattern': f'{new_id}-CO33-t*',
                               'n_traces': N_TRACES, 'variant_tag': VARIANT})

    n = len(inputs)
    budget = {
        'schema': 'batch2d_slope_t1t3_co_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {
            'measured_t3_wall_range_s': [24.0, 29.0],
            'estimate_s3': 'same 40 m domain as t3 -> ~24-29 s',
            'estimate_s1': '50 m domain, x1.25 cells -> ~30-36 s (estimate, not commitment)',
            'job_commit_GiB_estimate': {'S3': 3.4, 'S1': 4.3},
            'note': 'Debye ADE accumulator negligible; vctip lingering can add ~15 min per case '
                    'if the singleton watcher lapses'},
        'batch_totals': {'n_cases': n, 'estimate_total_h': [2.9, 4.2],
                         'conservative_budget_min': 300},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 6,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': n},
        'known_limitations': ['estimates are for cap-setting only, not ETA',
                              'per-case supervised process startup included'],
    }
    cases_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-29_t1_t3_extension_proposal.md + '
                         'configs/research/dispersion_materials_v0.1.json (hash-locked)',
                 'n_cases': n, 'n_exceptions': 0,
                 'derivation': 't1t3 mothers derived from archived t2 mothers by scripted line '
                               'transforms (domain/rock-box resize to per-tier width, full '
                               'staircase replacement at the tier gradient, dispersive material '
                               'swap); NC nullcontrast / TGT target lines pass through; CO traces '
                               'byte-identical to their mother except #title/#hertzian_dipole/#rx',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                  'convention': 'spec 4.2: any split is by group_id (mother model) only; each '
                                'mother is its own family; tilt variants never share group_id',
                  'groups': groups_out}
    return mothers, inputs, cases_doc, groups_doc, budget, static


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    mothers, inputs, cases_doc, groups_doc, budget, static = collect()

    if args.verify_only:
        for mid, text in mothers.items():
            assert (MOTHER_DIR / (mid + '.in')).read_text(encoding='utf-8') == text, mid
        for run_id, text in inputs.items():
            assert (OUT / (run_id + '.in')).read_text(encoding='utf-8') == text, run_id
        for name, doc in (('cases.json', cases_doc), ('groups.json', groups_doc),
                          ('budget.json', budget), ('static_check.json', static)):
            assert json.loads((OUT / name).read_text(encoding='utf-8')) == doc, name
        print(f'verify-only: {len(mothers)} mothers + all {len(inputs)} CO inputs + '
              f'cases/groups/budget/static_check byte-identical')
        return

    MOTHER_DIR.mkdir(parents=True, exist_ok=True)
    for mid, text in mothers.items():
        (MOTHER_DIR / (mid + '.in')).write_text(text, encoding='utf-8', newline='\n')
    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_slope_t1t3_co.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')
    cmd_src = ROOT / 'scripts/run_batch2d_slope_t3_co_gpu.cmd'
    cmd_dst = ROOT / 'scripts/run_batch2d_slope_t1t3_co_gpu.cmd'
    cmd_dst.write_text(cmd_src.read_text(encoding='utf-8').replace(
        'run_approved_batch2d_slope_t3_co.py', 'run_approved_batch2d_slope_t1t3_co.py'),
        encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-SLOPE-T1T3-CO',
        'run_id': rid,
        'input_path': f'configs/research/batch2d_slope_t1t3_co/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (S3X, S3TZX, S1X, S1TZX; BG, NC, TGT; t01..t33); '
                        'any failure aborts the batch. No retries.',
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
        'approved_to_simulate': True,
        'approved_run_ids': [c['run_id'] for c in contracts],
        'approved_run_ids_pending_agreement': [],
        'approved_compute_budget': {
            'wall_minutes': 20, 'job_commit_GiB': 6, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': len(contracts), 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-29: "按建议开始" on the T1/T3 extension proposal — both tiers, T1 '
            'domain y 40->50 m (outcrop at y~46, margin 2.75 m to PML), TGT fixed anchor z, '
            'CO33 segment 1 only (396 runs), MT33 deferred, reference-window construction '
            'not in this batch',
            'geometry: t3 outcrop paradigm applied to T1/T3; zc anchored at y=16 (T3 22.6 m '
            'unchanged from t2 rule, T1 27.0 m); staircase rule identical (0.25 m strips, '
            'grid-locked); Fermat domain-sizing assertion embedded in static_check.json',
            'materials: byte-identical to batch2d_slope_t3_co (dispersion_materials_v0.1 '
            'hash-locked); proposal TZ text quoting the old draft interpolation superseded',
            'cost basis: t3 measured 24-29 s/case on the 40 m domain; S1 +25% cells; '
            'conservative 300 min total; per-case wall cap 20 min',
        ],
        'permitted_preparation': old_gate['permitted_preparation'],
        'requires_agreement_before_execution': old_gate['requires_agreement_before_execution'],
        'approved_execution_contracts': contracts,
        'execution_policy': old_gate['execution_policy'],
        'execution_outcome': {'status': 'approved_for_execution',
                              'basis': 'user 2026-09-29: "按建议开始" (T1/T3 proposal, per recommendation)'},
        'last_completed_execution_contract': {
            'batch_id': old_gate['batch_id'],
            'execution_outcome': old_gate.get('execution_outcome'),
            'note': old_gate.get('note'),
        },
        'note': 'batch2d_slope_t1t3_co frozen 2026-09-29 (12 mothers = T1/T3 x +/-TZ x BG/NC/TGT, '
                '33 CO traces each, BASE, dispersive cover/tzone per dispersion_materials_v0.1). '
                'New families S1X/S1TZX/S3X/S3TZX: mechanism/fidelity input only — no frozen '
                'reward contract applies until their own reference windows + tolerance '
                'recalibration exist (reward_weights v0.3 scope_limits). Test families '
                '{C5,C8} untouched; no labels; G4 maintained.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {MOTHER_DIR.relative_to(ROOT)}: {len(mothers)} mothers')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=true, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
