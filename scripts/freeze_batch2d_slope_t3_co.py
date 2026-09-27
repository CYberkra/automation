"""Freeze the batch2d_slope_t3 common-offset (CO) segment: revised slope family.

Basis: user 2026-09-28 — "1和2允许" (approve the batch-wide switch to the fitted
Debye dispersive materials of configs/research/dispersion_materials_v0.1.json, and
approve the slope-geometry revision that fixes the domain-truncation dominance found
2026-09-28: for the 15 m flight height + 11.31 deg dipping interface the Fermat
stationary point of the interface echo sits up-dip of every trace (yq ~ 27..35 m for
the CO array) while the t2 domain ends at y=32, so the archived t2 interface band was
dominated by the truncated domain edge, not by the true specular reflection).

Revision design (t3 family, "X" suffix):
- Slope mothers extend the domain y 32 -> 40 m; the interface staircase continues
  with the identical rule (zc 25.8 m grid-locked, +0.05 m per 0.25 m strip) until the
  natural outcrop at y ~ 37 m; beyond the outcrop bedrock reaches the surface. This
  replaces the fictitious truncation reflector with the true outcrop physics. A
  fictitious horizontal section was considered and rejected: it would imprint a
  constant-depth flat band that the geology does not have.
- All covers switch to the dispersive single-pole Debye cover (eps_inf 18.017,
  dE 7.878, tau 6.4567e-9 s, sigma_dc 0.003); tzone to the half-strength pole
  (eps_inf 12.5, dE 4.0, same tau, sigma_dc 0.003); rock stays non-dispersive
  (dispersion_materials_v0.1, literature-anchored 2026-09-28).
- Flat control keeps the batch2d_v1 geometry (no truncation issue for a flat
  interface; specular point sits under the trace) and only switches the cover to the
  dispersive material, so domain sizes differ inside the batch by design.

Domain sizing evidence (embedded static check): full 3-point Fermat (air legs from
the z=45 antennas to surface entry/exit points, clay legs to the specular point on
the extended interface, smooth-line stand-in for the staircase) gives stationary
points yq <= ~34.8 m for all 33 traces at both er=16 and er=18; the check asserts
max(yq) + 2.0 m <= PML start (39.0 m). The interface itself ends at the outcrop
(y ~ 37 m), 2 m before the PML.

Mothers (all BG per the 2026-09-27 ablation; generated deterministically here):
- configs/research/batch2d_slope_t3/B2D-C3mS2X-BG.in    (slope, extended, dispersive)
- configs/research/batch2d_slope_t3/B2D-C3mS2TZX-BG.in  (slope + 1 m dispersive TZ)
- configs/research/batch2d_slope_t3/B2D-C3mX-BG.in      (flat control, dispersive)

CO array identical to batch2d_slope_t2_co: Tx = Rx - 1.30 m, both z = 45, Rx
y = 12.65..20.65 step 0.25, 33 traces per mother, anchor t17 == mother single-trace
pair (no archived baseline exists for the t3 mothers; t17 validation is the
post-run Fermat-arrival check, not bit identity).

Mirrors scripts/freeze_batch2d_slope_t2_co.py / freeze_dispersion_contrast_s2.py:
generates mothers + 99 .in inputs, cases/groups contracts, static_check.json,
budget.json, the runner clone, and rewrites the global execution gate with
approved_to_simulate=false. --verify-only rebuilds in memory and asserts byte
equality with disk without touching the gate. No solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SLOPE_T2 = ROOT / 'configs/research/batch2d_slope_t2'
FLAT_BASE = ROOT / 'configs/research/batch2d_v1'
OUT = ROOT / 'configs/research/batch2d_slope_t3_co'
MOTHER_DIR = ROOT / 'configs/research/batch2d_slope_t3'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'
DISP = json.loads((ROOT / 'configs/research/dispersion_materials_v0.1.json').read_text(encoding='utf-8'))

C = 299792458.0
C_M_NS = 0.299792458
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
BATCH_ID = 'batch2d_slope_t3_co'
RUN_DATE = '2026-09-28'

# t2 mother geometry constants (grid-locked staircase, gen_slope_family_geometry.py)
SURF_Z = 30.0
ZC_C = round(25.8 / DY)           # 1032
YMID_T2 = 16.0
TAN_EFF = 0.2                     # realised theta_eff 11.309932 deg
STRIP_W = 0.25
Y2_OLD, Y2_NEW = 32.0, 40.0       # domain y-extent old -> new
OUTCROP_MARGIN_M = 2.0            # asserted margin: max stationary yq + 2 <= PML start

COVER = DISP['materials']['cover_clay']
TZ = DISP['materials']['transition_zone']
assert not DISP['materials']['bedrock_sandstone']['dispersive']
COVER_MAT = f"#material: {COVER['epsilon_inf']} {COVER['sigma_dc']} 1 0 cover"
COVER_POLE = (f"#add_dispersion_debye: 1 {COVER['debye_poles'][0]['delta_epsilon']} "
              f"{COVER['debye_poles'][0]['tau_s']} cover")
TZ_MAT = f"#material: {TZ['epsilon_inf']} {TZ['sigma_dc']} 1 0 tzone"
TZ_POLE = (f"#add_dispersion_debye: 1 {TZ['debye_poles'][0]['delta_epsilon']} "
           f"{TZ['debye_poles'][0]['tau_s']} tzone")

MOTHERS = [
    # (new_mother_id, source_dir, source_mother_id, extend_domain, with_tz)
    ('B2D-C3mS2X-BG', SLOPE_T2, 'B2D-C3mS2-BG', True, False),
    ('B2D-C3mS2TZX-BG', SLOPE_T2, 'B2D-C3mS2TZ-BG', True, True),
    ('B2D-C3mX-BG', FLAT_BASE, 'B2D-C3m-BG', False, False),
]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(mother_id: str, k: int) -> int:
    txt = f'{BATCH_ID}|{mother_id}|t{k + 1:02d}|{VARIANT}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def strip_zif_c(ycen: float) -> int:
    """Grid-cell interface level of the strip centred at ycen (t2 rule, unchanged)."""
    return ZC_C + round((ycen - YMID_T2) * TAN_EFF / DY)


def extension_strips():
    """Cover/tz box lines continuing the staircase from y=32 to the outcrop (~37 m).

    Strips whose interface would sit at/above the surface are all-rock (no boxes).
    Returns (cover_lines, tz_lines, y_outcrop_lo, y_outcrop_hi)."""
    cover_lines, tz_lines = [], []
    n = round((Y2_NEW - Y2_OLD) / STRIP_W)
    ycrop_lo = ycrop_hi = None
    for i in range(n):
        ya = Y2_OLD + i * STRIP_W
        yb = ya + STRIP_W
        ycen = ya + STRIP_W / 2
        zif = strip_zif_c(ycen) * DY
        if zif >= SURF_Z - 1e-12:
            if ycrop_lo is None:
                ycrop_lo, ycrop_hi = ya, yb
            continue
        cover_lines.append("#box: inf %.6g %.6g inf %.6g %.6g cover" % (ya, zif, yb, SURF_Z))
        tz_lines.append("#box: inf %.6g %.6g inf %.6g %.6g tzone" % (ya, zif - 1.0, yb, zif))
    assert ycrop_lo is not None, 'interface never outcrops inside the extended domain'
    return cover_lines, tz_lines, ycrop_lo, ycrop_hi


def fermat_stationary(er: float):
    """Stationary interface-reflection point per CO trace on the extended geometry.

    Full 3-point Fermat: air legs from the z=45 antennas to surface entry/exit
    points, clay legs to the specular point Q on the smooth interface line
    z = 0.2*y + 22.625 clamped at the surface (outcrop). Smooth line is a stand-in
    for the staircase (deviation <= 0.025 m + quantization, immaterial for the
    domain-margin assertion)."""
    v1 = C_M_NS
    v2 = C_M_NS / math.sqrt(er)
    ys_rx = np.array([Y_FIRST + k * DY_RX for k in range(N_TRACES)])
    yq_max = []
    for rx_y in ys_rx:
        tx_y = rx_y - OFFSET
        y1g = np.linspace(rx_y - 8, rx_y + 14, 80)
        yqg = np.linspace(12.0, 36.9, 400)
        y2g = np.linspace(rx_y - 6, rx_y + 16, 80)
        Y1, YQ, Y2 = np.meshgrid(y1g, yqg, y2g, indexing='ij')
        ZQ = np.minimum(0.2 * YQ + 22.625, SURF_Z)
        t = ((np.hypot(Y1 - tx_y, 15.0) + np.hypot(Y2 - rx_y, 15.0)) / v1
             + (np.hypot(YQ - Y1, SURF_Z - ZQ) + np.hypot(Y2 - YQ, SURF_Z - ZQ)) / v2)
        yq_max.append(float(YQ[np.unravel_index(int(np.argmin(t)), t.shape)]))
    return yq_max


def build_mother(new_id: str, src_dir: Path, src_id: str, extend: bool, with_tz: bool,
                 ext_cover, ext_tz, ycrop) -> str:
    src_text = (src_dir / (src_id + '.in')).read_text(encoding='utf-8')
    lines = src_text.splitlines()
    assert lines.count('#domain: inf 32 50') == 1
    assert lines.count('#box: inf 0 0 inf 32 30 rock') == 1
    assert lines.count('#material: 16 0.01 1 0 cover') == 1
    if with_tz:
        assert lines.count('#material: 12 0.005 1 0 tzone') == 1

    out = []
    for l in lines:
        if l.startswith('#title:'):
            out.append(f'#title: {new_id} {BATCH_ID} mother t3 revised slope family '
                       f'(dispersive cover, domain y<={Y2_NEW:g} with outcrop, source {src_id})')
            continue
        if extend and l == '#domain: inf 32 50':
            out.append(f'#domain: inf {Y2_NEW:g} 50')
            continue
        if extend and l == '#box: inf 0 0 inf 32 30 rock':
            out.append(f'#box: inf 0 0 inf {Y2_NEW:g} 30 rock')
            continue
        if l == '#material: 16 0.01 1 0 cover':
            out.append(COVER_MAT)
            out.append(COVER_POLE)
            continue
        if with_tz and l == '#material: 12 0.005 1 0 tzone':
            out.append(TZ_MAT)
            out.append(TZ_POLE)
            continue
        out.append(l)

    if extend:
        # insert extension cover strips after the last existing cover box,
        # extension tz strips after the last existing tz box (S2TZX only)
        cover_idx = [i for i, l in enumerate(out) if l.startswith('#box:') and l.endswith(' cover')]
        assert len(cover_idx) == round(Y2_OLD / STRIP_W), f'{src_id}: unexpected cover strip count'
        first_ext = ext_cover[0]
        assert abs(float(first_ext.split()[2]) - Y2_OLD) < 1e-9
        prev_last = out[cover_idx[-1]].split()  # #box: x0 y0 z0 x1 y1 z1 id
        assert abs(float(prev_last[5]) - Y2_OLD) < 1e-9, 'last cover strip must end at the old edge'
        # staircase continuity across the seam: z rises by exactly one 0.05 m step
        assert abs(float(first_ext.split()[3]) - float(prev_last[3]) - 0.05) < 1e-9
        out[cover_idx[-1] + 1:cover_idx[-1] + 1] = ext_cover
        if with_tz:
            tz_idx = [i for i, l in enumerate(out) if l.startswith('#box:') and l.endswith(' tzone')]
            assert len(tz_idx) == round(Y2_OLD / STRIP_W)
            out[tz_idx[-1] + 1:tz_idx[-1] + 1] = ext_tz

    text = '\n'.join(out) + '\n'
    assert '#material: 16 0.01 1 0 cover' not in text and COVER_MAT in text and COVER_POLE in text
    if with_tz:
        assert '#material: 12 0.005 1 0 tzone' not in text and TZ_MAT in text and TZ_POLE in text
    if extend:
        assert f'#domain: inf {Y2_NEW:g} 50' in text and '#domain: inf 32 50' not in text
    text.encode('ascii')  # gprMax reads inputs with the system codec (GBK here); ASCII-only
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
    ny_new, nz = int(Y2_NEW / DY), int(50 / DZ)
    cells_new = ny_new * nz
    main_fields_gib_new = 72 * 2 * (ny_new + 1) * (nz + 1) / 2 ** 30
    ny_old = int(Y2_OLD / DY)
    cells_old = ny_old * nz
    main_fields_gib_old = 72 * 2 * (ny_old + 1) * (nz + 1) / 2 ** 30

    ext_cover, ext_tz, ycrop_lo, ycrop_hi = extension_strips()
    assert len(ext_cover) == 20 and len(ext_tz) == 20
    assert 36.5 <= ycrop_lo and ycrop_hi <= 37.5, (ycrop_lo, ycrop_hi)

    # domain-sizing evidence: stationary points of the interface echo (er 16 and 18)
    yq16 = fermat_stationary(16.0)
    yq18 = fermat_stationary(18.0)
    yq_star = max(max(yq16), max(yq18))
    pml_start = Y2_NEW - 40 * DY
    assert yq_star + OUTCROP_MARGIN_M <= pml_start, (yq_star, pml_start)

    mothers, inputs, cases_out, groups_out = {}, {}, [], []
    static = {
        'schema': 'batch2d_slope_t3_co_static_check/1',
        'checked_date': RUN_DATE,
        'revision': 't3: dispersive cover/tzone (dispersion_materials_v0.1); slope domain '
                    'y 32->40 m with natural outcrop at y ~ %.3g m replacing the fictitious '
                    'truncation reflector (2026-09-28 Fermat finding)' % (0.5 * (ycrop_lo + ycrop_hi)),
        'base_batches': ['batch2d_slope_t2 (slope BG mothers, geometry source)',
                         'batch2d_v1 (flat BG mother, geometry source)'],
        'ablation': 'BG only per user 2026-09-27 ablation decision (target ablated)',
        'dispersion_config': 'configs/research/dispersion_materials_v0.1.json',
        'dispersion_config_sha256': sha256_bytes(
            (ROOT / 'configs/research/dispersion_materials_v0.1.json').read_bytes()),
        'grid_slope': {'dy_m': DY, 'dz_m': DZ, 'ny': ny_new, 'nz': nz, 'cells': cells_new,
                       'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT,
                       'dt_formula_over_contract': dt_formula / DT_CONTRACT},
        'grid_flat': {'dy_m': DY, 'dz_m': DZ, 'ny': ny_old, 'nz': nz, 'cells': cells_old},
        'time_window_s': TIME_WINDOW_S,
        'time_steps_contract': 20352,
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'order': 'x0 y0 z0 xmax ymax zmax (V4 source '
                'fdtd_grid.set_pml_thickness)', 'physical_thickness_m': 1.0},
        'outcrop': {'y_interval_m': [ycrop_lo, ycrop_hi],
                    'rule': 'staircase strips stop when z_if >= surface; beyond, bedrock '
                            'reaches the surface',
                    'extension_cover_strips': len(ext_cover), 'extension_tz_strips': len(ext_tz)},
        'fermat_domain_sizing': {
            'method': 'full 3-point Fermat (air legs + clay legs, smooth-line interface '
                      'clamped at outcrop), stationary reflection point per trace',
            'yq_max_er16_m': round(max(yq16), 3), 'yq_max_er18_m': round(max(yq18), 3),
            'margin_rule_m': OUTCROP_MARGIN_M, 'pml_start_m': pml_start,
            'assertion': 'max(yq) + %.1f <= %.1f holds (%.3f)' % (OUTCROP_MARGIN_M, pml_start, yq_star)},
        'common_offset_array': {
            'n_traces': N_TRACES, 'spacing_m': DY_RX, 'offset_m': OFFSET,
            'rx_first_y_m': Y_FIRST, 'rx_last_y_m': Y_FIRST + (N_TRACES - 1) * DY_RX,
            'tx_first_y_m': Y_FIRST - OFFSET, 'tx_last_y_m': Y_FIRST + (N_TRACES - 1) * DY_RX - OFFSET,
            'z_m': Z_LINE, 'anchor_index': ANCHOR_INDEX,
            'anchor_rx_y_m': Y_FIRST + ANCHOR_INDEX * DY_RX,
            'anchor_tx_y_m': Y_FIRST + ANCHOR_INDEX * DY_RX - OFFSET,
            'trace_ids': [f't{k + 1:02d}' for k in range(N_TRACES)],
            'anchor_note': 't17 == mother single-trace pair; t3 mothers have no archived run, '
                           'so the anchor validates by post-run Fermat arrival, not bit identity',
        },
        'mother_checks': [], 'trace_checks': [],
    }
    assert cells_new == 3200000 and cells_old == 2560000

    for new_id, src_dir, src_id, extend, with_tz in MOTHERS:
        mother_text = build_mother(new_id, src_dir, src_id, extend, with_tz,
                                   ext_cover, ext_tz, (ycrop_lo, ycrop_hi))
        mothers[new_id] = mother_text
        keys = {k: [l for l in mother_text.splitlines() if l.startswith(k)]
                for k in ('#dx_dy_dz:', '#time_window:', '#pml_cells:', '#waveform:')}
        for k, v in keys.items():
            assert len(v) == 1, (new_id, k)
        static['mother_checks'].append({
            'mother': new_id, 'source_mother': src_id,
            'source_sha256': sha256_bytes((src_dir / (src_id + '.in')).read_bytes()),
            'mother_sha256': sha256_bytes(mother_text.encode('utf-8')),
            'extended_domain': extend, 'dispersive_cover': True, 'dispersive_tzone': with_tz})
        for k in range(N_TRACES):
            rx_y = Y_FIRST + k * DY_RX
            tx_y = rx_y - OFFSET
            rx_cell, tx_cell = rx_y / DY, tx_y / DY
            assert abs(rx_cell - round(rx_cell)) < 1e-9 and abs(tx_cell - round(tx_cell)) < 1e-9
            assert 1.0 < tx_y and rx_y < 31.0, f'inside PML: tx {tx_y} rx {rx_y}'
            rx_id = f't{k + 1:02d}'
            run_id = f'{new_id}-CO33-{rx_id}'
            text = build_co_input(mother_text, run_id, new_id, tx_y, rx_y, rx_id)
            inputs[run_id] = text
            strip = lambda s: [l for l in s.splitlines()
                               if not l.startswith(('#rx:', '#title:', '#hertzian_dipole:'))]
            assert strip(mother_text) == strip(text), run_id
            if k == ANCHOR_INDEX:
                assert '#rx: inf 16.65 45 t17 Ex' in text
                assert '#hertzian_dipole: x inf 15.35 45 impulse' in text
            static['trace_checks'].append({'run_id': run_id, 'rx_cell': int(round(rx_cell)),
                                           'tx_cell': int(round(tx_cell)), 'outside_pml': True})
            cases_out.append({
                'run_id': run_id, 'file': run_id + '.in',
                'input_sha256': sha256_bytes(text.encode('utf-8')),
                'mother_model_id': new_id, 'case_id': f'{new_id}|{rx_id}',
                'group_id': new_id,
                'segment': 'co', 'variant_tag': VARIANT,
                'seed': seed_for(new_id, k),
                'common_offset': {
                    'n_traces': N_TRACES, 'trace_spacing_m': DY_RX, 'offset_m': OFFSET,
                    'rx_y_m': rx_y, 'tx_y_m': tx_y, 'z_m': Z_LINE,
                    'trace_index': k, 'anchor_trace': k == ANCHOR_INDEX,
                    'acquisition_geometry': 'common-offset profile along y (Tx = Rx - 1.30 m, both z=45)',
                    'rx_output_components': ['Ex'],
                },
            })
        groups_out.append({'case_id': new_id, 'group_id': new_id,
                           'run_id_pattern': f'{new_id}-CO33-t*',
                           'n_traces': N_TRACES, 'variant_tag': VARIANT})

    budget = {
        'schema': 'batch2d_slope_t3_co_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells_slope': cells_new, 'cells_flat': cells_old,
                 'dt_s': DT_CONTRACT, 'dt_formula_s': dt_formula,
                 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB_slope': main_fields_gib_new,
                        'main_fields_ID_GiB_flat': main_fields_gib_old,
                        'measured_t2_wall_range_s': [17.812, 29.0],
                        'note': 'slope traces scale ~x1.25 cells vs the measured t2 grid; Debye '
                                'ADE adds one polarisation accumulator per dispersive material '
                                '(negligible); flat traces match the measured t2 grid'},
        'batch_totals': {'n_cases': 3 * N_TRACES,
                         'measured_estimate_s': 3 * N_TRACES * 38,
                         'conservative_budget_min': 150},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 6,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 3 * N_TRACES},
        'known_limitations': [
            'per-case supervised process startup included in the measured range',
            'vctip.exe job-lingering can add ~15 min per case if the out-of-job singleton lapses',
            'job_commit raised 4->6 GiB: slope domain cells x1.25 vs the 3.16 GB measured peak',
            'estimates are for cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                 'spec': 'decision_log 2026-09-28 (Fermat truncation finding + batch-wide '
                         'dispersive switch) + docs/research/2026-09-28_dispersion_parameters.md + '
                         'configs/research/dispersion_materials_v0.1.json',
                 'n_cases': 3 * N_TRACES, 'n_exceptions': 0,
                 'derivation': 't3 mothers derived from t2/v1 mothers by scripted line transforms '
                               '(domain/rock-box extension, staircase continuation to outcrop, '
                               'dispersive material swap); CO traces byte-identical to their '
                               'mother except #title/#hertzian_dipole/#rx',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                  'convention': 'spec 4.2: any split is by group_id (mother model) only, never by case, '
                                'window or trace; each mother is its own family',
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
        assert json.loads((OUT / 'cases.json').read_text(encoding='utf-8')) == cases_doc
        assert json.loads((OUT / 'groups.json').read_text(encoding='utf-8')) == groups_doc
        assert json.loads((OUT / 'budget.json').read_text(encoding='utf-8')) == budget
        assert json.loads((OUT / 'static_check.json').read_text(encoding='utf-8')) == static
        print(f'verify-only: {len(mothers)} t3 mothers + all {len(inputs)} CO inputs + '
              f'cases/groups/budget/static_check byte-identical')
        return

    MOTHER_DIR.mkdir(parents=True, exist_ok=True)
    for mid, text in mothers.items():
        (MOTHER_DIR / (mid + '.in')).write_text(text, encoding='utf-8', newline='\n')
    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d_mt.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_slope_t3_co.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-SLOPE-T3-CO',
        'run_id': rid,
        'input_path': f'configs/research/batch2d_slope_t3_co/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (S2X BG, S2TZX BG, flat C3mX BG; t01..t33); '
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
            'wall_minutes': 20, 'job_commit_GiB': 6, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': 3 * N_TRACES, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'user 2026-09-28: "1和2允许" — batch-wide switch to the fitted Debye dispersive '
            'materials (dispersion_materials_v0.1; measured effect on the interface echo '
            '+2.5 ns / -3.3..-6.6 dB in dispersion_contrast_s2_t17) AND slope-geometry revision '
            'fixing the domain-truncation dominance (Fermat stationary point up-dip of every '
            'trace, t2 domain ended at y=32); remaining S2TZ/flat CO cases decided: full 99-case '
            'refreeze on the revised family, t2 CO results archived as superseded for slope '
            'interpretation',
            'geometry revision: domain y 32->40 m, interface staircase continued with the '
            'identical grid-locked rule to the natural outcrop at y~37 m (a fictitious '
            'horizontal section was rejected: it would imprint a constant-depth band); '
            'Fermat domain-sizing assertion embedded in static_check.json',
            'cost basis: measured 17.8-29.0 s/case on the t2 grid; slope traces ~x1.25 cells; '
            '99 cases conservative 150 min total; per-case wall cap 20 min',
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
        'note': 'batch2d_slope_t3_co frozen 2026-09-28 (3 revised BG mothers x 33 common-offset '
                'traces, offset 1.30 m, BASE, BG-only per ablation; dispersive cover/tzone; '
                'slope domains extended to the outcrop). approved_to_simulate stays false until '
                'the recorded basis is applied; t2 CO runs archived, not overwritten; no '
                'physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {MOTHER_DIR.relative_to(ROOT)}: {len(mothers)} t3 mothers')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, {len(contracts)} contracts')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
