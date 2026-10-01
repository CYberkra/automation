"""Freeze the batch2d_b2_pilot_co batch: B2 pilot training-data sampling (12 mothers, 396 CO33 runs).

Basis: user 2026-10-01 "可以，开始" on the pretraining prerequisite ledger P1-1, following the
B2 sampling design draft docs/research/2026-10-01_b2_sampling_design_draft.md with all §6
decision points approved per recommendation ("都认", 2026-10-01).

Composition (A层 8 + B层 4 = 12 mothers; A层 BG archives reused, NOT re-run):
  A层 (reward-evaluable, dev families C1/C3, tolerance v0.2 + weights v0.3 apply):
    flat C1 trio:  B2D-C1mX-{BG, NC, D10m-W4m-T0.5m-E20-S0.02}   (from batch2d_v1 C1m sources)
    flat C3 pair:  B2D-C3mX-{NC, D10m-W4m-T0.5m-E20-S0.02}       (BG archived 2026-09-28)
    slope T2 pair: B2D-C3mS2X-{NC, D10m-W4m-T0.5m-E20-S0.02}     (BG archived 2026-09-28)
    slope T2 TZ:   B2D-C3mS2TZX-NC                               (BG archived; TGT deferred to main batch)
  B层 (new-family probes, NO frozen reward applies until own reference windows + recalibration):
    B2D-C1p5mS1X-BG (tan 0.1, 50 m domain), B2D-C1p5mS3X-BG (tan 0.4, 40 m),
    B2D-C2mS3X-BG, B2D-C2mS3TZX-BG (tan 0.4, 40 m, TZ probe)

zc rule disclosure: A层 slope mothers keep the archived t3 geometry (zc=25.8 legacy t2 rule,
anchor cover 4.2 m) byte-exact so NC/TGT pair with the archived BG. B层 new tiers use
zc = 30 - h (vertical cover at anchor y=16 equals the label h); the two rules are NOT nested
and B层 labels are new families anyway. T1 at h>=2.0 and T3 at h=2.5 fail the outcrop/PML
margin rule inside existing domains, so B层 covers h in {1.5, 2.0} only; h=2.5 deferred.

Materials byte-identical to the t3 baseline (dispersion_materials_v0.1.json, hash-locked).
NC semantics: BG geometry + target-shaped nullcontrast box (rock params); NC-BG bit-identity
is a post-run requirement, re-verified, never assumed. CO array identical to t3/t1t3.

Gate discipline: this freeze writes the gate with approved_to_simulate=FALSE and all 396
run_ids in approved_run_ids_pending_agreement. Execution requires explicit user sign-off.

Mirrors scripts/freeze_batch2d_slope_t1t3_co.py. --verify-only rebuilds in memory and asserts
byte equality with disk without touching the gate. No solver is launched here.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FLAT_SRC = ROOT / 'configs/research/batch2d_v1'
SLOPE_SRC = ROOT / 'configs/research/batch2d_slope_t2'
T3_ARCHIVE = ROOT / 'configs/research/batch2d_slope_t3'
OUT = ROOT / 'configs/research/batch2d_b2_pilot_co'
MOTHER_DIR = ROOT / 'configs/research/batch2d_b2_pilot_mothers'
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
BATCH_ID = 'batch2d_b2_pilot_co'
RUN_DATE = '2026-10-01'

TGT_SUFFIX = 'D10m-W4m-T0.5m-E20-S0.02'

# slope tiers: tag -> (tan_theta, zc_m, domain_y_m, with_tz)
SLOPE_TIERS = {
    'S2X': (0.2, 25.8, 40.0, False),      # A层: archived t3 geometry, zc legacy rule
    'S2TZX': (0.2, 25.8, 40.0, True),
    'C1p5mS1X': (0.1, 28.5, 50.0, False),  # B层: zc = 30 - h (anchor cover == label)
    'C1p5mS3X': (0.4, 28.5, 40.0, False),
    'C2mS3X': (0.4, 28.0, 40.0, False),
    'C2mS3TZX': (0.4, 28.0, 40.0, True),
}
SLOPE_SRC_ID = {False: 'B2D-C3mS2', True: 'B2D-C3mS2TZ'}  # by with_tz

# (new_mother_id, layer, kind, source info)
MOTHERS = [
    # A层 flat (from batch2d_v1; dispersive swap only)
    ('B2D-C1mX-BG', 'A', 'flat', 'B2D-C1m-BG'),
    ('B2D-C1mX-NC', 'A', 'flat', 'B2D-C1m-D10m-NC'),
    ('B2D-C1mX-' + TGT_SUFFIX, 'A', 'flat', 'B2D-C1m-' + TGT_SUFFIX),
    ('B2D-C3mX-NC', 'A', 'flat', 'B2D-C3m-D10m-NC'),
    ('B2D-C3mX-' + TGT_SUFFIX, 'A', 'flat', 'B2D-C3m-' + TGT_SUFFIX),
    # A层 slope T2 (from t2 sources; full staircase rebuild at archived t3 geometry)
    ('B2D-C3mS2X-NC', 'A', 'slope', ('S2X', 'NC')),
    ('B2D-C3mS2X-' + TGT_SUFFIX, 'A', 'slope', ('S2X', 'TGT')),
    ('B2D-C3mS2TZX-NC', 'A', 'slope', ('S2TZX', 'NC')),
    # B层 new-family probes (BG only)
    ('B2D-C1p5mS1X-BG', 'B', 'slope', ('C1p5mS1X', 'BG')),
    ('B2D-C1p5mS3X-BG', 'B', 'slope', ('C1p5mS3X', 'BG')),
    ('B2D-C2mS3X-BG', 'B', 'slope', ('C2mS3X', 'BG')),
    ('B2D-C2mS3TZX-BG', 'B', 'slope', ('C2mS3TZX', 'BG')),
]

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
    """Grid-locked cover/tzone staircase strips from y=0 to the outcrop (t1t3 machinery)."""
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
    a, b = np.polyfit(centers, zifs, 1)
    theta_eff = math.degrees(math.atan(a))
    resid = np.array(zifs) - (a * np.array(centers) + b)
    return (cover_lines, tz_lines, len(cover_lines), y_out_lo, theta_eff,
            float(np.max(np.abs(resid))), {'slope': a, 'intercept': b})


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


def build_flat_mother(new_id: str, src_id: str) -> str:
    """v1 flat mother -> dispersive cover swap + retitle; boxes (incl. NC/target) unchanged."""
    src_text = (FLAT_SRC / (src_id + '.in')).read_text(encoding='utf-8')
    lines = src_text.splitlines()
    assert lines.count('#domain: inf 32 50') == 1
    assert lines.count('#material: 16 0.01 1 0 cover') == 1
    n_cover_box = sum(1 for l in lines if l.startswith('#box:') and l.endswith(' cover'))
    assert n_cover_box == 1
    cover_z = 29.0 if '-C1m-' in src_id else 27.0
    assert f'#box: inf 0 {cover_z:g} inf 32 30 cover' in lines
    if '-NC' in src_id:
        assert lines.count('#material: 9 0.001 1 0 nullcontrast') == 1
        assert lines.count('#box: inf 14 19.75 inf 18 20.25 nullcontrast') == 1
    if TGT_SUFFIX in src_id:
        assert lines.count('#material: 20 0.02 1 0 target') == 1
        assert lines.count('#box: inf 14 19.75 inf 18 20.25 target') == 1

    out = []
    for l in lines:
        if l.startswith('#title:'):
            out.append(f'#title: {new_id} {BATCH_ID} mother flat dispersive '
                       f'(source batch2d_v1 {src_id})')
            continue
        if l == '#material: 16 0.01 1 0 cover':
            out.append(COVER_MAT)
            out.append(COVER_POLE)
            continue
        out.append(l)
    text = '\n'.join(out) + '\n'
    assert '#material: 16 0.01 1 0 cover' not in text
    assert text.count(COVER_MAT) == 1 and text.count(COVER_POLE) == 1
    text.encode('ascii')
    return text


def build_slope_mother(new_id: str, tier: str, scene: str, cover_lines, tz_lines) -> str:
    """t2 slope source -> per-tier staircase rebuild + dispersive swap (t1t3 machinery).

    NC nullcontrast / TGT target lines pass through unchanged from the t2 source."""
    tan, zc, y_max, with_tz = SLOPE_TIERS[tier]
    src_scene = {'BG': 'BG', 'NC': 'NC', 'TGT': TGT_SUFFIX}[scene]
    src_id = f'{SLOPE_SRC_ID[with_tz]}-{src_scene}'
    src_text = (SLOPE_SRC / (src_id + '.in')).read_text(encoding='utf-8')
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
            out.append(f'#title: {new_id} {BATCH_ID} mother tier {tier} scene {scene} '
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
        'schema': 'batch2d_b2_pilot_co_static_check/1',
        'checked_date': RUN_DATE,
        'basis': 'B2 sampling design draft docs/research/2026-10-01_b2_sampling_design_draft.md; '
                 'user 2026-10-01 "都认" on draft section 6 decision points + "可以，开始" on the '
                 'prerequisite ledger P1-1',
        'composition': {
            'A_layer': 'reward-evaluable dev-family cells (C1 flat trio; C3 flat NC/TGT; '
                       'T2 slope S2X NC/TGT; S2TZX NC); BG archives 2026-09-28 reused, not re-run',
            'B_layer': 'new-family probes (cover 1.5/2.0 m x tilt T1/T3, one TZ probe); no frozen '
                       'reward applies until own reference windows + tolerance recalibration '
                       '(reward_weights v0.3 scope_limits)',
            'zc_rule_disclosure': 'A slope mothers keep archived t3 geometry (zc=25.8, legacy t2 '
                                  'rule, anchor cover 4.2 m) for byte-exact BG pairing; B tiers '
                                  'use zc=30-h (anchor cover == label); rules not nested; '
                                  'h=2.5 and T1 h>=2.0 fail outcrop/PML margin in existing '
                                  'domains and are deferred to the main batch',
        },
        'dispersion_config': 'configs/research/dispersion_materials_v0.1.json',
        'dispersion_config_sha256': sha256_bytes(DISP_PATH.read_bytes()),
        'grid': {'dy_m': DY, 'dz_m': DZ, 'time_window_s': TIME_WINDOW_S,
                 'time_steps_contract': 20352, 'dt_formula_s': dt_formula,
                 'dt_contract_s': DT_CONTRACT},
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'common_offset_array': {
            'n_traces': N_TRACES, 'spacing_m': DY_RX, 'offset_m': OFFSET,
            'rx_first_y_m': Y_FIRST, 'z_m': Z_LINE, 'anchor_index': ANCHOR_INDEX,
            'trace_ids': [f't{t + 1:02d}' for t in range(N_TRACES)],
            'anchor_note': 't17 == mother single-trace pair (Rx 16.65, Tx 15.35)'},
        'nc_bit_identity': 'NC-BG paired difference must be bitwise zero post-run; re-verified '
                           'on the new geometries, never assumed inherited',
        'tiers': {}, 'mother_checks': [], 'trace_checks': [],
        'archive_reuse': ['B2D-C3mX-BG', 'B2D-C3mS2X-BG', 'B2D-C3mS2TZX-BG'],
    }

    # tier geometry + Fermat margins
    tier_geom = {}
    for tier, (tan, zc, y_max, with_tz) in SLOPE_TIERS.items():
        cover_lines, tz_lines, n_strips, y_out, theta_eff, qdev, fit = staircase(tan, zc, y_max)
        pml_start = y_max - 40 * DY
        yq16 = fermat_stationary(tan, zc, y_max, 16.0)
        yq18 = fermat_stationary(tan, zc, y_max, 18.0)
        yq_star = max(max(yq16), max(yq18))
        assert y_out + 2.0 <= pml_start, (tier, y_out, pml_start)
        assert yq_star + 2.0 <= pml_start, (tier, yq_star, pml_start)
        tier_geom[tier] = (cover_lines, tz_lines)
        ny, nz = int(y_max / DY), int(50 / DZ)
        static['tiers'][tier] = {
            'tan_theta': tan, 'zc_m': zc, 'domain_y_m': y_max, 'with_tz': with_tz,
            'anchor_cover_m': round(SURF_Z - zc, 3), 'n_strips': n_strips,
            'outcrop_y_lo_m': y_out, 'pml_inner_y_m': pml_start,
            'outcrop_to_pml_margin_m': round(pml_start - y_out, 3),
            'theta_eff_deg': round(theta_eff, 6),
            'theta_fit_residual_max_m': qdev,
            'fermat_yq_max_er16_m': round(max(yq16), 3),
            'fermat_yq_max_er18_m': round(max(yq18), 3),
            'fermat_margin_assertion': f'max(yq)+2.0 <= {pml_start} (yq_star {yq_star:.3f})',
            'cells': ny * nz, 'main_fields_ID_GiB': 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30,
        }

    # cross-check: rebuilt S2X/S2TZX staircase must byte-match the archived t3 BG mothers
    for tier, arch_id in (('S2X', 'B2D-C3mS2X-BG'), ('S2TZX', 'B2D-C3mS2TZX-BG')):
        arch = (T3_ARCHIVE / (arch_id + '.in')).read_text(encoding='utf-8').splitlines()
        arch_cover = [l for l in arch if l.startswith('#box:') and l.endswith(' cover')]
        arch_tz = [l for l in arch if l.startswith('#box:') and l.endswith(' tzone')]
        new_cover, new_tz = tier_geom[tier]
        assert arch_cover == new_cover, f'{tier}: rebuilt staircase != archived t3 mother'
        assert arch_tz == (new_tz if tier.endswith('TZX') else []), f'{tier}: tz mismatch'
        static.setdefault('staircase_cross_check', []).append(
            {'tier': tier, 'archived_mother': arch_id,
             'cover_strips': len(new_cover), 'byte_identical_box_lines': True})

    mothers, inputs, cases_out, groups_out = {}, {}, [], []
    # execution order: 32 m flat, then 40 m slope, then 50 m S1 probe
    order_key = {'B2D-C1mX-BG': 0, 'B2D-C1mX-NC': 1, 'B2D-C1mX-' + TGT_SUFFIX: 2,
                 'B2D-C3mX-NC': 3, 'B2D-C3mX-' + TGT_SUFFIX: 4,
                 'B2D-C3mS2X-NC': 5, 'B2D-C3mS2X-' + TGT_SUFFIX: 6, 'B2D-C3mS2TZX-NC': 7,
                 'B2D-C1p5mS3X-BG': 8, 'B2D-C2mS3X-BG': 9, 'B2D-C2mS3TZX-BG': 10,
                 'B2D-C1p5mS1X-BG': 11}
    for new_id, layer, kind, src in sorted(MOTHERS, key=lambda m: order_key[m[0]]):
        if kind == 'flat':
            mother_text = build_flat_mother(new_id, src)
            src_sha = sha256_bytes((FLAT_SRC / (src + '.in')).read_bytes())
            tier_tag = 'flat'
        else:
            tier, scene = src
            mother_text = build_slope_mother(new_id, tier, scene, *tier_geom[tier])
            src_scene = {'BG': 'BG', 'NC': 'NC', 'TGT': TGT_SUFFIX}[scene]
            src_id = f"{SLOPE_SRC_ID[SLOPE_TIERS[tier][3]]}-{src_scene}"
            src_sha = sha256_bytes((SLOPE_SRC / (src_id + '.in')).read_bytes())
            tier_tag = tier
        mothers[new_id] = mother_text
        keys = {k: [l for l in mother_text.splitlines() if l.startswith(k)]
                for k in ('#domain:', '#dx_dy_dz:', '#time_window:', '#pml_cells:', '#waveform:')}
        for k, v in keys.items():
            assert len(v) == 1, (new_id, k)
        static['mother_checks'].append({
            'mother': new_id, 'layer': layer, 'tier': tier_tag,
            'source_sha256': src_sha,
            'mother_sha256': sha256_bytes(mother_text.encode('utf-8')),
            'reward_evaluable': layer == 'A'})
        for k in range(N_TRACES):
            rx_y = Y_FIRST + k * DY_RX
            tx_y = rx_y - OFFSET
            rx_cell, tx_cell = rx_y / DY, tx_y / DY
            assert abs(rx_cell - round(rx_cell)) < 1e-9 and abs(tx_cell - round(tx_cell)) < 1e-9
            y_max = 32.0 if kind == 'flat' else SLOPE_TIERS[tier][2]
            assert rx_y < min(y_max - 1.0, 31.0) and tx_y > 1.0, 'inside PML'
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
                'layer': layer, 'seed': seed_for(new_id, k),
                'common_offset': {
                    'n_traces': N_TRACES, 'trace_spacing_m': DY_RX, 'offset_m': OFFSET,
                    'rx_y_m': rx_y, 'tx_y_m': tx_y, 'z_m': Z_LINE,
                    'trace_index': k, 'anchor_trace': k == ANCHOR_INDEX,
                    'acquisition_geometry': 'common-offset profile along y (Tx = Rx - 1.30 m, both z=45)',
                    'rx_output_components': ['Ex']}})
        groups_out.append({'case_id': new_id, 'group_id': new_id, 'layer': layer,
                           'run_id_pattern': f'{new_id}-CO33-t*',
                           'n_traces': N_TRACES, 'variant_tag': VARIANT})

    n = len(inputs)
    assert n == 396, n
    budget = {
        'schema': 'batch2d_b2_pilot_co_budget/1',
        'date': RUN_DATE,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {
            'measured_t3_wall_range_s': [24.0, 29.0],
            'estimate_flat_32m_s': [18.0, 24.0],
            'estimate_slope_40m_s': [24.0, 29.0],
            'estimate_s1_50m_s': [30.0, 36.0],
            'note': 'estimates for cap-setting only, not ETA; vctip lingering can add '
                    '~15 min per case if the singleton watcher lapses'},
        'batch_totals': {'n_cases': n,
                         'mothers_by_domain': {'32m_flat': 5, '40m_slope': 6, '50m_s1': 1},
                         'estimate_total_h': [2.4, 3.6],
                         'conservative_budget_min': 300},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 6,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': n},
        'known_limitations': ['serial fail-abort, retries=0, CUDA double, no CPU fallback',
                              'GPU blocked at freeze time: NVML init failure pending OS reboot '
                              '(CBS_RebootPending=True, 2026-10-01); preflight must re-verify'],
    }
    cases_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-10-01_b2_sampling_design_draft.md + '
                         'configs/research/dispersion_materials_v0.1.json (hash-locked)',
                 'n_cases': n, 'n_exceptions': 0,
                 'derivation': 'A-layer flat mothers from batch2d_v1 (dispersive swap only); '
                               'A-layer slope NC/TGT and B-layer probes from t2 sources with '
                               'full staircase rebuild (t1t3 machinery); S2X/S2TZX staircase '
                               'byte-verified against archived t3 BG mothers; NC/TGT lines pass '
                               'through; CO traces byte-identical to mother except '
                               '#title/#hertzian_dipole/#rx',
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': 'co', 'grid_tier': 'BASE',
                  'convention': 'any split is by group_id (mother model) only; each mother is '
                                'its own family; tilt variants never share group_id; test '
                                'families {C5,C8} untouched',
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
    runner_dst = ROOT / 'scripts/run_approved_batch2d_b2_pilot_co.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1_mt'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1_mt'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')
    cmd_src = ROOT / 'scripts/run_batch2d_slope_t3_co_gpu.cmd'
    cmd_dst = ROOT / 'scripts/run_batch2d_b2_pilot_co_gpu.cmd'
    cmd_dst.write_text(cmd_src.read_text(encoding='utf-8').replace(
        'run_approved_batch2d_slope_t3_co.py', 'run_approved_batch2d_b2_pilot_co.py'),
        encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-B2-PILOT-CO',
        'run_id': rid,
        'input_path': f'configs/research/batch2d_b2_pilot_co/{rid}.in',
        'input_sha256': sha256_bytes(inputs[rid].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{rid}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{rid}',
        'continuation': 'Serial in frozen order (flat C1 trio, flat C3 pair, S2X pair, S2TZX NC, '
                        'then B-layer S3 probes, S1 probe last; t01..t33 within each mother); '
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
        'approved_to_simulate': False,
        'approved_run_ids': [],
        'approved_run_ids_pending_agreement': [c['run_id'] for c in contracts],
        'approved_compute_budget': {
            'wall_minutes': 20, 'job_commit_GiB': 6, 'minimum_available_RAM_GiB': 8,
            'output_GiB': 2, 'threads': 8, 'backend': 'CUDA', 'precision': 'double',
            'retries': 0, 'max_fdtd_runs': len(contracts), 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'B2 sampling design draft docs/research/2026-10-01_b2_sampling_design_draft.md; '
            'user 2026-10-01 "都认" on draft section 6 (two-layer design, 12-mother pilot, '
            'B-layer attached, T1 50 m domain, CO33 only) and "可以，开始" on ledger P1-1',
            'A-layer: flat C1 trio + flat C3 NC/TGT + S2X NC/TGT + S2TZX NC (8 mothers); '
            'archived BG (2026-09-28 t3_co) reused per draft section 5; B-layer: 4 new-family '
            'probes (cover 1.5/2.0 m, T1/T3, one TZ); B-layer reward NOT applicable until own '
            'reference windows + tolerance recalibration',
            'zc rule disclosure: A slope mothers byte-match archived t3 geometry (zc=25.8 '
            'legacy); B tiers use zc=30-h; rules not nested; h=2.5 / T1 h>=2 deferred',
            'GPU prerequisite: NVML init failure diagnosed 2026-10-01 '
            '(CBS_RebootPending=True); execution preflight must re-verify nvidia-smi after reboot',
        ],
        'permitted_preparation': old_gate['permitted_preparation'],
        'requires_agreement_before_execution': old_gate['requires_agreement_before_execution'],
        'approved_execution_contracts': contracts,
        'execution_policy': old_gate['execution_policy'],
        'execution_outcome': {'status': 'frozen_pending_user_agreement',
                              'basis': 'gate written with approved_to_simulate=false per B2 '
                                       'draft section 7; awaiting explicit user sign-off'},
        'last_completed_execution_contract': {
            'batch_id': old_gate['batch_id'],
            'execution_outcome': old_gate.get('execution_outcome'),
            'note': old_gate.get('note'),
        },
        'note': 'batch2d_b2_pilot_co frozen 2026-10-01 (12 mothers = A8+B4, 33 CO traces each, '
                'BASE, dispersive cover/tzone per dispersion_materials_v0.1). B-layer families '
                'C1p5mS1X/C1p5mS3X/C2mS3X/C2mS3TZX: mechanism/fidelity input only. Test families '
                '{C5,C8} untouched; no labels; G4 maintained; NOT approved to simulate.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {MOTHER_DIR.relative_to(ROOT)}: {len(mothers)} mothers')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(inputs)} .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=FALSE, '
          f'{len(contracts)} contracts pending agreement')
    print('launcher_sha256', launcher_sha)


if __name__ == '__main__':
    main()
