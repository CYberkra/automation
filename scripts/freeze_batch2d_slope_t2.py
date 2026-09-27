"""Freeze the first slope_family_v0 2D batch: tier T2 (tan_theta=0.2) +/- transition zone.

User authorization: 2026-09-27 five-point reply to the slope_family_v0 design draft
(docs/research/2026-09-27_slope_family_design_draft.md section 9), point 1:
"首批做个2D的即可，具体你来指定，最大效率验证即可" and point 5: "排期" (schedule approved,
scheme A first). kimi designates the first batch as: C3 family x T2 x {noTZ, TZ} x {BG, NC, TGT}
= 6 cases, coarse BASE grid (dy=dz=0.025), everything except the cover-bedrock geometry
byte-identical to the frozen anchor B2D-C3m-D10m-W4m-T0.5m-E20-S0.02 (batch2d_v1).

6 cases in execution order (family BG precedes NC precedes TGT):
  B2D-C3mS2-BG, B2D-C3mS2-NC, B2D-C3mS2-D10m-W4m-T0.5m-E20-S0.02,
  B2D-C3mS2TZ-BG, B2D-C3mS2TZ-NC, B2D-C3mS2TZ-D10m-W4m-T0.5m-E20-S0.02

Naming caveat (recorded, not hidden): the "C3m" family tag is inherited from the flat-anchor
naming convention; the T2 tier shifts the cover center from 3.0 m to 4.2 m
(zc 27.0 -> 25.8, cover_thickness_range [1.025, 7.375]) to keep min cover >= 1 m
(geometry.json constraints). Cross-tier comparisons must not attribute cover-center
differences to slope angle alone.

--verify-only rebuilds everything in memory and asserts byte equality with disk
(two-run byte determinism gate); it does not touch the execution gate.

Read-only over archived contracts; no solver is launched here.
"""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANCHOR = ROOT / 'configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in'
ANCHOR_NC = ROOT / 'configs/research/batch2d_v1/B2D-C3m-D10m-NC.in'
GEOMETRY = ROOT / 'configs/research/slope_family_v0/geometry.json'
GEOMETRY_SHA256 = '760f1b71146edc4475dc7357e1ba8ebace7bf04f774385ff0b5123c92ff4b3de'
OUT = ROOT / 'configs/research/batch2d_slope_t2'
GATE = ROOT / 'configs/research/gprmax_v4_execution_gate.json'
RUNTIME_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json'
CUDA_ID = ROOT / 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json'

C = 299792458.0
DY = DZ = 0.025
DT_CONTRACT = 5.896635841874211e-11
TIME_WINDOW_S = 1200e-9
BATCH_ID = 'batch2d_slope_t2'
SEGMENT = 'base'
TZ_MATERIAL_LINE = '#material: 12 0.005 1 0 tzone'
RUN_DATE = '2026-09-27'

ROLES = ('BG', 'NC', 'TGT')
TIERS = (('S2', '2d_T2', 'slope_t2', False), ('S2TZ', '2d_T2TZ', 'slope_t2_tz', True))


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def seed_for(group_id: str, case_id: str, variant: str) -> int:
    txt = f'{BATCH_ID}|{group_id}|{case_id}|{variant}'
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def parse_box(line: str):
    m = re.fullmatch(r'#box: inf (\S+) (\S+) inf (\S+) (\S+) (\S+)', line)
    assert m, line
    return float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)), m.group(5)


def build_input(anchor_head, cover_lines, tz_lines, role, run_id, title_tag):
    """Assemble one .in text. anchor_head = mother lines up to and including the cover
    material line (i.e. everything except title/cover-box/target lines)."""
    lines = [f'#title: {run_id} {BATCH_ID} {title_tag}']
    lines += anchor_head
    lines += cover_lines
    if tz_lines:
        lines.append(TZ_MATERIAL_LINE)
        lines += tz_lines
    if role == 'NC':
        lines.append('#material: 9 0.001 1 0 nullcontrast')
        lines.append('#box: inf 14 19.75 inf 18 20.25 nullcontrast')
    elif role == 'TGT':
        lines.append('#material: 20 0.02 1 0 target')
        lines.append('#box: inf 14 19.75 inf 18 20.25 target')
    return '\n'.join(lines) + '\n'


def static_check_input(run_id, text, variant, role, mother_lines):
    """Per-case static checks; returns the check record. Raises on any violation."""
    lines = text.splitlines()
    assert lines[0].startswith(f'#title: {run_id} ')
    # Everything outside the geometry block must be byte-identical to the mother
    # (mother lines 1..12 header/source/rx/pml, then rock material/box, cover material).
    head = [l for l in lines if not l.startswith('#box:') and not l.startswith('#material:')
            and not l.startswith('#title:')]
    mother_head = [l for l in mother_lines if not l.startswith('#box:') and not l.startswith('#material:')
                   and not l.startswith('#title:')]
    assert head == mother_head, f'{run_id}: non-geometry header differs from mother'
    boxes = [parse_box(l) for l in lines if l.startswith('#box:')]
    mats = {l.split()[-1]: l for l in lines if l.startswith('#material:')}
    assert mats['rock'] == '#material: 9 0.001 1 0 rock'
    assert mats['cover'] == '#material: 16 0.01 1 0 cover'
    assert ('tzone' in mats) == bool(variant['transition_zone'])
    if variant['transition_zone']:
        assert mats['tzone'] == TZ_MATERIAL_LINE
    # grid alignment of every box coordinate
    for y0, z0, y1, z1, _name in boxes:
        for v in (y0, z0, y1, z1):
            c = v / DY
            assert abs(c - round(c)) < 1e-9, f'{run_id}: off-grid coordinate {v}'
    rock_boxes = [b for b in boxes if b[4] == 'rock']
    assert rock_boxes == [(0.0, 0.0, 32.0, 30.0, 'rock')]
    cover_boxes = [b for b in boxes if b[4] == 'cover']
    tz_boxes = [b for b in boxes if b[4] == 'tzone']
    assert len(cover_boxes) == variant['n_strips']
    strips = variant['strips']
    for i, b in enumerate(cover_boxes):
        s = strips[i]
        assert b == (s['y0'], s['z_if'], s['y1'], 30.0, 'cover'), f'{run_id} cover strip {i}'
    if variant['transition_zone']:
        assert len(tz_boxes) == variant['n_strips']
        for i, b in enumerate(tz_boxes):
            s = strips[i]
            assert b == (s['y0'], s['z_if'] - 1.0, s['y1'], s['z_if'], 'tzone'), f'{run_id} tz strip {i}'
    else:
        assert not tz_boxes
    # y partition of [0, 32]
    assert cover_boxes[0][0] == 0.0 and cover_boxes[-1][2] == 32.0
    for a, b in zip(cover_boxes, cover_boxes[1:]):
        assert a[2] == b[0]
    # target / nullcontrast slot and clearance vs interface and TZ bottom
    z_if_min = min(s['z_if'] for s in strips)
    slot = [b for b in boxes if b[4] in ('target', 'nullcontrast')]
    if role == 'BG':
        assert not slot
        clearance = None
    else:
        assert slot == [(14.0, 19.75, 18.0, 20.25, slot[0][4])]
        clearance = round(z_if_min - (1.0 if variant['transition_zone'] else 0.0) - 20.25, 6)
        assert clearance > 0
    return {
        'run_id': run_id,
        'header_byte_identical_to_mother': True,
        'n_cover_strips': len(cover_boxes),
        'n_tz_strips': len(tz_boxes),
        'all_box_coordinates_grid_aligned': True,
        'y_partition_0_32': True,
        'strips_match_geometry_json': True,
        'z_interface_min_m': z_if_min,
        'target_slot_clearance_m': clearance,
        'materials': sorted(mats),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()

    assert sha256_bytes(GEOMETRY.read_bytes()) == GEOMETRY_SHA256, 'geometry.json changed'
    geometry = json.loads(GEOMETRY.read_text(encoding='utf-8'))
    assert geometry['contract_id'] == 'slope_family_v0'
    assert geometry['required_material_line'] == TZ_MATERIAL_LINE

    mother_text = ANCHOR.read_text(encoding='utf-8')
    mother_lines = mother_text.splitlines()
    cover_box_i = mother_lines.index('#box: inf 0 27 inf 32 30 cover')
    assert mother_lines[cover_box_i + 1:] == ['#material: 20 0.02 1 0 target',
                                              '#box: inf 14 19.75 inf 18 20.25 target']
    anchor_head = mother_lines[1:cover_box_i]  # header + rock + cover material, no title
    nc_text = ANCHOR_NC.read_text(encoding='utf-8')
    assert '#material: 9 0.001 1 0 nullcontrast' in nc_text

    # ---------- grid/time invariants (same domain as batch2d_v1 base) ----------
    dt_formula = 1.0 / (C * math.sqrt(DY ** -2 + DZ ** -2))
    ny, nz = int(32 / DY), int(50 / DZ)
    cells = ny * nz
    main_fields_gib = 72 * 2 * (ny + 1) * (nz + 1) / 2 ** 30
    assert cells == 2560000
    assert abs(main_fields_gib - 0.3437627702951431) < 1e-9

    static = {
        'schema': 'batch2d_slope_t2_static_check/1',
        'checked_date': RUN_DATE,
        'anchor_mother': 'configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in',
        'anchor_mother_sha256': sha256_bytes(ANCHOR.read_bytes()),
        'geometry_contract': 'configs/research/slope_family_v0/geometry.json',
        'geometry_sha256': GEOMETRY_SHA256,
        'grid': {'dy_m': DY, 'dz_m': DZ, 'ny': ny, 'nz': nz, 'cells': cells,
                 'dt_formula_s': dt_formula, 'dt_contract_s': DT_CONTRACT,
                 'dt_formula_over_contract': dt_formula / DT_CONTRACT},
        'time_window_s': TIME_WINDOW_S,
        'time_steps_contract': 20352,
        'pml': {'cells': [0, 40, 40, 0, 40, 40], 'physical_thickness_m': 1.0},
        'naming_caveat': ('C3m family tag inherited from the flat anchor; T2 tier cover center is '
                          '4.2 m (zc 25.8, range [1.025, 7.375]), not 3.0 m. Do not attribute '
                          'cross-tier cover-center differences to slope angle alone.'),
        'case_checks': [],
    }

    order = []  # (run_id, tier_tag, variant_key, variant_tag, has_tz, role)
    for tier_tag, variant_key, variant_tag, has_tz in TIERS:
        for role in ROLES:
            suffix = {'BG': 'BG', 'NC': 'NC',
                      'TGT': 'D10m-W4m-T0.5m-E20-S0.02'}[role]
            order.append((f'B2D-C3m{tier_tag}-{suffix}', tier_tag, variant_key, variant_tag, has_tz, role))
    assert len(order) == 6 and len({o[0] for o in order}) == 6

    cases_out, groups_out, inputs = [], [], {}
    for run_id, tier_tag, variant_key, variant_tag, has_tz, role in order:
        variant = geometry['variants'][variant_key]
        text = build_input(anchor_head, variant['cover_box_lines'],
                           variant['tz_box_lines'] if has_tz else [], role, run_id,
                           f'slope family v0 tier {tier_tag} {"BG" if role == "BG" else role} '
                           f'(theta_eff {variant["theta_eff_deg"]} deg, zc {variant["zc_m"]} m)')
        check = static_check_input(run_id, text, variant, role, mother_lines)
        static['case_checks'].append(check)
        inputs[run_id] = text
        in_sha = sha256_bytes(text.encode('utf-8'))

        group_id = f'B2D-C3m{tier_tag}'
        case_id = run_id
        seed = seed_for(group_id, case_id, variant_tag)
        skeleton = {
            'domain_m': [None, 32.0, 50.0], 'dy_dz_m': [DY, DZ],
            'pml_cells': [0, 40, 40, 0, 40, 40], 'time_window_s': TIME_WINDOW_S,
            'tx_m': [15.35, 45.0], 'rx_m': [16.65, 45.0], 'component': 'Ex',
            'surface_z_m': 30.0, 'materials': {'rock': [9.0, 0.001], 'cover': [16.0, 0.01]},
            'slope': {'tier': 'T2', 'tan_theta_nominal': 0.2,
                      'theta_eff_deg': variant['theta_eff_deg'], 'zc_m': variant['zc_m'],
                      'z_interface_range_m': variant['z_interface_range'],
                      'transition_zone': has_tz,
                      'tz_material': [12.0, 0.005] if has_tz else None},
        }
        mother_hash = sha256_bytes(json.dumps(skeleton, sort_keys=True).encode('utf-8'))
        cases_out.append({
            'case_id': case_id, 'run_id': run_id, 'file': run_id + '.in',
            'input_sha256': in_sha, 'batch_id': BATCH_ID, 'segment': SEGMENT,
            'grid_tier': 'BASE', 'family': f'C3{tier_tag}',
            'cover_thickness_center_m': round(sum(variant['cover_thickness_range']) / 2, 6),
            'role': role.lower(), 'mother_model_id': f'B2D-C3m{tier_tag}-BG',
            'group_id': group_id, 'seed': seed, 'variant_tag': variant_tag,
            'scan_axis': {
                'slope_tier': 'T2', 'tan_theta_nominal': 0.2,
                'theta_eff_deg': variant['theta_eff_deg'],
                'transition_zone': has_tz,
                'depth_m': 10.0 if role == 'TGT' else None,
                'width_m': 4.0 if role == 'TGT' else None,
                'thickness_m': 0.5 if role == 'TGT' else None,
                'eps_r': 20.0 if role == 'TGT' else None,
                'sigma_s_m': 0.02 if role == 'TGT' else None,
            },
            'materials': {'rock': [9.0, 0.001], 'cover': [16.0, 0.01],
                          'tzone': [12.0, 0.005] if has_tz else None,
                          'target': [20.0, 0.02] if role == 'TGT' else None,
                          'nullcontrast': [9.0, 0.001] if role == 'NC' else None},
            'geometry': {
                'domain_m': [None, 32.0, 50.0], 'dy_dz_m': [DY, DZ],
                'pml_cells': [0, 40, 40, 0, 40, 40], 'pml_physical_thickness_m': 1.0,
                'surface_z_m': 30.0,
                'cover_strips_n': variant['n_strips'],
                'cover_strip_width_m': variant['strip_width_m'],
                'zc_m': variant['zc_m'],
                'z_interface_range_m': variant['z_interface_range'],
                'cover_thickness_range_m': variant['cover_thickness_range'],
                'max_quantization_deviation_m': variant['max_quantization_deviation_m'],
                'target_box_m': [14.0, 19.75, 18.0, 20.25] if role != 'BG' else None,
                'target_z_fixed_horizontal': True,
                'tx_m': [15.35, 45.0], 'rx_m': [16.65, 45.0], 'component': 'Ex',
            },
            'cells': cells, 'dt_s': DT_CONTRACT, 'time_steps': 20352,
            'main_fields_ID_GiB': main_fields_gib,
            'static_check': {'grid_aligned': True,
                             'target_slot_clearance_m': check['target_slot_clearance_m']},
            'exceptions': None,
        })
        groups_out.append({'case_id': case_id, 'group_id': group_id, 'run_id': run_id,
                           'mother_model_hash': mother_hash, 'seed': seed,
                           'variant_tag': variant_tag})

    # ---------- budget.json (spec 3.1 formulas; measured base-batch anchors) ----------
    yfine_ref = {'wall_s': 1088.781, 'cells': 40960000, 'dt_s': 1.474e-11}
    wall_scaled = yfine_ref['wall_s'] * (cells / yfine_ref['cells']) * (yfine_ref['dt_s'] / DT_CONTRACT)
    budget = {
        'schema': 'batch2d_slope_t2_budget/1',
        'date': RUN_DATE,
        'formulas': 'spec 3.1: cells=(32/dy)*(50/dz); dt=1/(c*sqrt(dy^-2+dz^-2)); '
                    'main_fields_ID_GiB=72B*2*(ny+1)(nz+1)/2^30; wall=ref*(cells/cells_ref)*(dt_ref/dt)',
        'grid': {'dy_m': DY, 'dz_m': DZ, 'cells': cells, 'dt_s': DT_CONTRACT,
                 'dt_formula_s': dt_formula, 'time_steps': 20352, 'time_window_ns': 1200},
        'single_case': {'main_fields_ID_GiB': main_fields_gib,
                        'wall_scaled_from_YFINE_s': wall_scaled,
                        'measured_base_batch_mean_s': 48.3,
                        'measured_base_batch_range_s_excl_first': [29.859, 59.218],
                        'strip_count_note': '129-257 #box commands per case vs 2-3 in the flat anchor; '
                                            'scene-build overhead expected small vs solver wall, unmeasured'},
        'batch_totals': {'n_cases': 6, 'planning_envelope_s': 6 * 150,
                         'measured_mean_estimate_s': 6 * 48.3,
                         'conservative_budget_min': 20},
        'hard_stops': {'wall_minutes_per_case': 20, 'job_commit_GiB': 4,
                       'output_GiB': 2, 'retries': 0, 'max_fdtd_runs': 6},
        'known_limitations': [
            'working scaling is not a guaranteed runtime',
            'main-field array excludes PML, update coefficients, CUDA context and host temporaries',
            'Job committed memory is not GPU VRAM',
            'estimates are for ranking and cap-setting only, not ETA',
        ],
    }

    cases_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                 'spec': 'docs/research/2026-09-27_slope_family_design_draft.md (kimi-accepted draft; '
                         'first-batch scope designated by kimi under user 2026-09-27 five-point reply) + '
                         'docs/research/2026-09-26_batch_2d_spec_v1.md',
                 'n_cases': 6, 'n_exceptions': 0,
                 'derivation': 'cover-bedrock geometry from slope_family_v0 geometry.json variants '
                               '2d_T2/2d_T2TZ; all other lines byte-identical to the frozen batch2d_v1 '
                               'C3 anchor mother',
                 'naming_caveat': static['naming_caveat'],
                 'cases': cases_out}
    groups_doc = {'batch_id': BATCH_ID, 'segment': SEGMENT, 'grid_tier': 'BASE',
                  'convention': 'spec 4.2: any split is by group_id only, never by case_id, window or trace; '
                                'each slope variant is its own mother-model family with its own BG/NC '
                                '(never paired against the flat batch2d_v1 backgrounds)',
                  'groups': groups_out}

    if args.verify_only:
        for run_id, text in inputs.items():
            disk = (OUT / (run_id + '.in')).read_text(encoding='utf-8')
            assert disk == text, f'{run_id}: disk differs from regeneration'
        assert json.loads((OUT / 'cases.json').read_text(encoding='utf-8')) == cases_doc
        assert json.loads((OUT / 'groups.json').read_text(encoding='utf-8')) == groups_doc
        assert json.loads((OUT / 'budget.json').read_text(encoding='utf-8')) == budget
        old_static = json.loads((OUT / 'static_check.json').read_text(encoding='utf-8'))
        assert old_static == static
        print('verify-only: all 6 inputs + cases/groups/budget/static_check byte-identical to regeneration')
        return

    OUT.mkdir(parents=True, exist_ok=True)
    for run_id, text in inputs.items():
        (OUT / (run_id + '.in')).write_text(text, encoding='utf-8', newline='\n')

    runner_src = ROOT / 'scripts/run_approved_batch2d.py'
    runner_dst = ROOT / 'scripts/run_approved_batch2d_slope_t2.py'
    txt = runner_src.read_text(encoding='utf-8')
    assert "BATCH = 'batch2d_v1'" in txt
    runner_dst.write_text(txt.replace("BATCH = 'batch2d_v1'", f"BATCH = '{BATCH_ID}'"),
                          encoding='utf-8', newline='\n')
    cmd_src = ROOT / 'scripts/run_batch2d_gpu.cmd'
    cmd_dst = ROOT / 'scripts/run_batch2d_slope_t2_gpu.cmd'
    cmd_txt = cmd_src.read_text(encoding='utf-8')
    assert 'run_approved_batch2d.py' in cmd_txt
    cmd_dst.write_text(cmd_txt.replace('run_approved_batch2d.py', 'run_approved_batch2d_slope_t2.py'),
                       encoding='utf-8', newline='\n')

    launcher_sha = sha256_bytes(runner_dst.read_bytes())
    supervisor_sha = sha256_bytes((ROOT / 'scripts/bounded_windows_process.py').read_bytes())
    runtime_sha = sha256_bytes(RUNTIME_ID.read_bytes())
    cuda_sha = sha256_bytes(CUDA_ID.read_bytes())

    contracts = [{
        'packet_id': 'BATCH2D-SLOPE-T2',
        'run_id': run_id,
        'input_path': f'configs/research/batch2d_slope_t2/{run_id}.in',
        'input_sha256': sha256_bytes(inputs[run_id].encode('utf-8')),
        'runtime_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json',
        'runtime_identity_sha256': runtime_sha,
        'cuda_identity_path': 'artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json',
        'cuda_identity_sha256': cuda_sha,
        'launcher_sha256': launcher_sha,
        'supervisor_sha256': supervisor_sha,
        'attempt_record': f'artifacts/research_checks/{RUN_DATE}_{run_id}_attempt.json',
        'run_directory': f'artifacts/simulations/{RUN_DATE}_{run_id}',
        'continuation': 'Serial in frozen order (per variant: BG, NC, TGT); any failure aborts the batch, '
                        'later attempts untouched. No retries.',
    } for run_id, *_ in order]

    (OUT / 'cases.json').write_text(json.dumps(cases_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'groups.json').write_text(json.dumps(groups_doc, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'budget.json').write_text(json.dumps(budget, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    (OUT / 'static_check.json').write_text(json.dumps(static, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')

    # ---------- gate rewrite: approved_to_simulate=false, order pending agreement ----------
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
            'retries': 0, 'max_fdtd_runs': 6, 'official_postprocess_wall_minutes': 10,
            'device_id': 0, 'minimum_available_VRAM_GiB': 4,
            'VRAM_limit_semantics': 'Preflight availability threshold; not a hard device-memory quota',
        },
        'scope_expansion_basis': [
            'slope_family_v0 design draft (docs/research/2026-09-27_slope_family_design_draft.md), '
            'kimi-accepted 2026-09-27; user five-point reply authorizes first 2D batch '
            '("首批做个2D的即可，具体你来指定，最大效率验证即可") and schedule ("排期", scheme A first)',
            'first-batch designation by kimi: C3 family x tier T2 (tan_theta=0.2, theta_eff=11.309932 deg) '
            'x {noTZ, +1m tzone} x {BG, NC, TGT} = 6 cases, BASE coarse grid; T2 is the middle tier, '
            'maximally informative per case for interface-tilt response and TZ effect screening',
            'geometry contract configs/research/slope_family_v0/geometry.json sha256 '
            + GEOMETRY_SHA256 + ' (two-run byte-identical, registry slope_family_v0)',
            'cost basis: batch2d_v1 coarse measured 29.9-59.2 s/case excl. first, job peak ~2.72 GiB; '
            'same domain/grid/time window; 6-case conservative budget 20 min; cap wall 20 min/case',
            'naming caveat recorded in cases.json/static_check.json: C3m tag inherited, T2 cover center '
            'is 4.2 m not 3.0 m; cross-tier cover-center differences must not be attributed to slope alone',
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
        'note': 'batch2d_slope_t2 frozen 2026-09-27 (6 cases: slope_family_v0 tier T2 +/- 1m transition '
                'zone x BG/NC/TGT, BASE coarse, target fixed horizontal z 19.75-20.25 with >=1.375 m '
                'clearance to TZ bottom). approved_to_simulate stays false until the agreement recorded '
                'in user_agreement is applied; G1/G3/G4/G5 untouched; no physical/training labels.',
    }
    GATE.write_text(json.dumps(new_gate, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: 6 .in + cases/groups/budget/static_check')
    print(f'wrote {GATE.relative_to(ROOT)}: {BATCH_ID}, approved_to_simulate=false, 6 contracts')
    print('launcher_sha256', launcher_sha)
    for run_id, *_ in order:
        print(' ', run_id, sha256_bytes(inputs[run_id].encode('utf-8'))[:16])


if __name__ == '__main__':
    main()
