"""Freeze the batch2d_v1 BASE event table v0.1 (docs/research/2026-09-26_event_table_proposal.md §7.1 steps 4-6).

Pure arithmetic + JSON output: no solver, no FDTD, no network. Window centres are self-computed
two-way travel times (proposal §3.1, t_basis=two_way_selfcomputed_v1, c=299792458 m/s);
processing results are never used to place or resize windows.

Time axis (verified 2026-09-26 against archived h5 + gprMax SFCW toolbox): t_i = i*dt with
receiver TimeSampleOffset=0 (gprMax ``rx.times`` = time_offset + dt*arange(n)); BASE
dt = 5.896635841874211e-11 s, N_analysis = min(time_steps, floor(1200ns/dt)+1) = 20351,
t_last = 1199.9653938214021 ns (NOT a T/(N-1) endpoint convention).

Mask indices use inward rounding: i_lo = ceil(t_lo/dt), i_hi = floor(t_hi/dt) (proposal §3.6(b)).
"""

import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'configs/research/batch2d_v1/cases.json'
GROUPS = ROOT / 'configs/research/batch2d_v1/groups.json'
OUT = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'

C_M_S = 299792458.0
HALF_WIDTH_NS = 80.0
TIME_WINDOW_NS = 1200.0
SCHEMA = 'event_window_v0.1'
T_BASIS = 'two_way_selfcomputed_v1'
CONV = 'ti_i_mul_dt_rxoffset0'  # verified: gprMax rx.times, t_i = i*dt, rx offset 0

FAMILIES = {'C1': 1, 'C3': 3, 'C5': 5, 'C8': 8}  # family -> cover thickness h (m)

# Proposal nominal BASE indices (§3.3 / §4) -- used as a cross-check of this script's rounding.
NOMINAL = {
    ('C1', 'COV'): (793, 3506), ('C1', 2): (1133, 3845), ('C1', 5): (2151, 4863), ('C1', 10): (3848, 6561),
    ('C3', 'COV'): (1698, 4411), ('C3', 5): (2377, 5090), ('C3', 10): (4074, 6787), ('C3', 20): (7468, 10181),
    ('C5', 'COV'): (2604, 5316), ('C5', 10): (4301, 7013), ('C5', 20): (7695, 10407),
    ('C8', 'COV'): (3961, 6674), ('C8', 10): (4640, 7352), ('C8', 20): (8034, 10747),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def t_center_ns(h: float, d) -> float:
    """Self-computed two-way travel time in ns. d='COV' -> cover bottom interface."""
    if d == 'COV':
        return (30.0 + 8.0 * h) / C_M_S * 1e9
    return (30.0 + 2.0 * h + 6.0 * d) / C_M_S * 1e9


def window_ns(h: float, d):
    tc = t_center_ns(h, d)
    return (tc - HALF_WIDTH_NS, tc + HALF_WIDTH_NS)


def fmt6(x: float) -> str:
    return f'{x:.6f}'


def window_hash(event_id, role, family, mother, group, t_lo, t_hi) -> str:
    canon = '|'.join([
        SCHEMA, event_id, role, family.lower(), mother, group, 'sample', 'ns',
        T_BASIS, '299792458', fmt6(t_lo), fmt6(t_hi), fmt6(HALF_WIDTH_NS), 'all', CONV,
    ])
    return hashlib.sha256(canon.encode('utf-8')).hexdigest()[:16]


def mask_hash(event_id, grid_tier, dt_s, i_lo, i_hi, n_samples) -> str:
    canon = '|'.join([event_id, grid_tier, f'{dt_s:.17e}', str(i_lo), str(i_hi), str(n_samples), CONV])
    return hashlib.sha256(canon.encode('utf-8')).hexdigest()[:16]


def main():
    cases_doc = json.loads(CASES.read_text(encoding='utf-8'))
    groups_doc = json.loads(GROUPS.read_text(encoding='utf-8'))
    cases = cases_doc['cases']
    dts = {c['dt_s'] for c in cases}
    assert len(dts) == 1, dts
    dt_s = dts.pop()
    n_raw = cases[0]['time_steps']
    n_samples = min(n_raw, int(math.floor(TIME_WINDOW_NS * 1e-9 / dt_s)) + 1)
    assert n_samples == 20351 and n_raw == 20352, (n_samples, n_raw)

    def inward(lo_ns, hi_ns):
        return math.ceil(lo_ns * 1e-9 / dt_s), math.floor(hi_ns * 1e-9 / dt_s)

    # role lookups from the executed contract (proposal §8.1: 29 cases = 4 BG + 4 NC + 1 OFF + 20 target)
    by_role = {}
    for c in cases:
        by_role.setdefault(c['role'], []).append(c)
    assert {k: len(v) for k, v in by_role.items()} == {'bg': 4, 'nc': 4, 'off': 1, 'target': 20}
    bg_by_fam = {c['family']: c for c in by_role['bg']}
    nc_by_fam = {c['family']: c for c in by_role['nc']}
    off_case = by_role['off'][0]
    assert off_case['family'] == 'C3'

    # actual depth ladder per family from executed targets (§8.2: no windows for removed cases)
    depths_by_fam = {}
    for c in by_role['target']:
        d = int(c['run_id'].split('-D')[1].split('m')[0])
        depths_by_fam.setdefault(c['family'], set()).add(d)

    entries = []

    def add(event_id, role, family, mother, group, t_lo, t_hi, extra=None):
        i_lo, i_hi = inward(t_lo, t_hi)
        assert 0 <= i_lo <= i_hi < n_samples, (event_id, i_lo, i_hi)
        e = dict(
            event_id=event_id, role=role, family=family.lower(), mother_model_id=mother, group_id=group,
            reference_type='paired_contrast' if role in ('cover_interface', 'target') else 'negative_control',
            reference_state='numerically_unresolved',
            domain='time_real', axes=['sample', 'trace'], time_unit='ns', event_axis='sample',
            t_basis=T_BASIS, c_m_s=int(C_M_S), half_width_ns=HALF_WIDTH_NS,
            t_lo_ns=t_lo, t_hi_ns=t_hi, trace_mask='all',
            grid_tier='BASE', dt_s=dt_s, n_samples=n_samples,
            index_lo=i_lo, index_hi=i_hi, n_window_samples=i_hi - i_lo + 1,
            time_axis_convention_version=CONV,
            window_hash=window_hash(event_id, role, family, mother, group, t_lo, t_hi),
            mask_hash=mask_hash(event_id, 'BASE', dt_s, i_lo, i_hi, n_samples),
            isolated_event_eligible=False,
            input_window='entire_predeclared_window',
        )
        if extra:
            e.update(extra)
        entries.append(e)

    # Table A: cover-interface events (one per family, applies to every case of the family)
    for fam, h in FAMILIES.items():
        lo, hi = window_ns(h, 'COV')
        add(f'EV-{fam}m-COV', 'cover_interface', fam, bg_by_fam[fam]['run_id'], f'B2D-{fam}m', lo, hi,
            extra=dict(applies_to='all_cases_of_family'))

    # Table B: target events, 1:1 with executed target cases
    for c in sorted(by_role['target'], key=lambda x: x['run_id']):
        fam = c['family']
        d = int(c['run_id'].split('-D')[1].split('m')[0])
        lo, hi = window_ns(FAMILIES[fam], d)
        add('EV-' + c['run_id'].removeprefix('B2D-'), 'target', fam, c['run_id'], c['group_id'], lo, hi)

    # Table C part 1: NC zero-contrast controls, W1 (geometry-matched D10m window) + W2 (full window)
    for fam in FAMILIES:
        lo, hi = window_ns(FAMILIES[fam], 10)
        add(f'EV-{fam}m-D10m-NC-W1', 'nc_zero', fam, nc_by_fam[fam]['run_id'], f'B2D-{fam}m', lo, hi,
            extra=dict(window_kind='geometry_matched_d10m', zero_reference_policy='no_epsilon_ratio; absolute_residual_only'))
        add(f'EV-{fam}m-D10m-NC-W2', 'nc_zero', fam, nc_by_fam[fam]['run_id'], f'B2D-{fam}m', 0.0, TIME_WINDOW_NS,
            extra=dict(window_kind='full_predeclared_window', mask='all_ones',
                       zero_reference_policy='no_epsilon_ratio; absolute_residual_only',
                       index_lo=0, index_hi=n_samples - 1))

    # Table C part 2: OFF far-offset control (C3, on-path D10m target window; NOT a zero-difference control)
    lo, hi = window_ns(FAMILIES['C3'], 10)
    add('EV-C3m-D10m-OFF-W1', 'off_path', 'C3', off_case['run_id'], 'B2D-C3m', lo, hi,
        extra=dict(window_kind='on_path_target_window', side_pml_margin_below_13m=True,
                   listed_separately=True, zero_difference_control=False))

    # Table C part 3: BG absent-target controls, one per (family, executed depth)
    for fam in FAMILIES:
        for d in sorted(depths_by_fam[fam]):
            lo, hi = window_ns(FAMILIES[fam], d)
            add(f'EV-{fam}m-BG-NEG-D{d}m', 'bg_absent', fam, bg_by_fam[fam]['run_id'], f'B2D-{fam}m', lo, hi,
                extra=dict(window_kind='target_window_on_family_bg'))

    # W2 rows use the full window: fix their (already overwritten) indices to [0, N-1]
    for e in entries:
        if e.get('window_kind') == 'full_predeclared_window':
            e['index_lo'], e['index_hi'] = 0, n_samples - 1
            e['n_window_samples'] = n_samples
            e['mask_hash'] = mask_hash(e['event_id'], 'BASE', dt_s, 0, n_samples - 1, n_samples)

    # Cross-checks against the proposal tables
    assert len(entries) == 43, len(entries)
    win_by_key = {}
    for e in entries:
        fam_key = e['family'].capitalize()
        if e['role'] == 'cover_interface':
            win_by_key[(fam_key, 'COV')] = (e['t_lo_ns'], e['t_hi_ns'], e['index_lo'], e['index_hi'])
        elif e['role'] == 'target':
            d = int(e['event_id'].split('-D')[1].split('m')[0])
            win_by_key[(fam_key, d)] = (e['t_lo_ns'], e['t_hi_ns'], e['index_lo'], e['index_hi'])
    for (fam, d), (i_lo, i_hi) in NOMINAL.items():
        got = win_by_key[(fam, d)]
        assert (got[2], got[3]) == (i_lo, i_hi), (fam, d, got, (i_lo, i_hi))
        # ns window values must match the proposal §3.3 table at 6 decimals
        h = FAMILIES[fam]
        lo, hi = window_ns(h, d)
        assert (fmt6(lo), fmt6(hi)) == (fmt6(got[0]), fmt6(got[1])), (fam, d)

    # Overlap matrix (centre distance < 160 ns => overlap), recorded per proposal §4.7
    overlaps = []
    for fam in FAMILIES:
        items = [('COV', t_center_ns(FAMILIES[fam], 'COV'))]
        for d in sorted(depths_by_fam[fam]):
            items.append((d, t_center_ns(FAMILIES[fam], d)))
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                dist = abs(items[i][1] - items[j][1])
                if dist < 2 * HALF_WIDTH_NS:
                    ov = 2 * HALF_WIDTH_NS - dist
                    overlaps.append(dict(family=fam, a=str(items[i][0]), b=str(items[j][0]),
                                         centre_distance_ns=dist, overlap_ns=ov))
    assert len(overlaps) == 9, len(overlaps)  # proposal §4.7: C1 4 + C3 3 + C5 1 + C8 1 = 9 overlapping pairs

    table = dict(
        schema_version='event_table_v0.1', batch_id='batch2d_v1', segment='base', grid_tier='BASE',
        status='frozen', frozen_date='2026-09-26',
        proposal='docs/research/2026-09-26_event_table_proposal.md',
        approved_by='user_2026-09-26', acceptance_note='kimi accepted; user approved with proposal defaults (roles, +/-80 ns, BG-NEG included, no narrow auxiliary window)',
        time_axis=dict(convention_version=CONV, formula='t_i = i*dt (i from 0), receiver TimeSampleOffset=0',
                       endpoint_convention='i_mul_dt_not_T_over_N_minus_1', dt_s=dt_s, n_samples=n_samples,
                       t_last_ns=(n_samples - 1) * dt_s * 1e9,
                       evidence='archived h5 SampleInterval/TimeSampleOffset + gprMax toolboxes/SFCW/processing.py rx.times; analyze_batch2d.py uses rx.times[:n]'),
        window_rule=dict(half_width_ns=HALF_WIDTH_NS, rounding='inward: i_lo=ceil(t_lo/dt), i_hi=floor(t_hi/dt)',
                         centre_rule='self-computed two-way travel time; processing peaks never used'),
        declarations=dict(
            R_c='not_computed_no_legitimate_pure_clutter_window',
            physical_thresholds_null=True, tolerances_null=True,
            training_labels_generated=False, training_eligible=False,
            physical_label_eligible=False, clean_truth_generated=False,
            reference_state='numerically_unresolved',
            isolated_event_granted=False,
            cross_family_differencing=False, cross_tier_differencing=False,
            negative_controls_merged_into_targets=False,
            data_split_unit='group_id_only',
            evaluation_before_freeze='undetermined',
        ),
        overlap_matrix=overlaps,
        n_entries=len(entries),
        inputs={p.as_posix(): sha256_file(p) for p in (CASES, GROUPS)},
        entries=entries,
    )
    OUT.write_text(json.dumps(table, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'wrote {OUT.relative_to(ROOT)}: {len(entries)} entries, sha256={sha256_file(OUT)}')
    roles = {}
    for e in entries:
        roles[e['role']] = roles.get(e['role'], 0) + 1
    print('role counts:', roles)
    print('nominal-index cross-check: 14/14 windows OK; overlap pairs: 9')


if __name__ == '__main__':
    main()
