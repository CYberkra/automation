"""One-shot S5 confirmation on frozen test families C5/C8.

Thresholds are read only from the frozen G4 program and decision log. This
runner creates an exclusive attempt marker before evaluating test outcomes;
failed formal attempts remain consumed. --selftest uses only synthetic rows.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_damage_ladder import (  # noqa: E402
    L_GAINS, L_LAMBDAS, L_WIDTHS, R_THRESHOLDS, catalogue,
)
from run_g4_mission_capability import (  # noqa: E402
    EVENT_TABLE, MISSION, MISSION_SHA256,
    TEST_CHUNK_GLOB, TEST_FAMILIES, load_chunks, validate_records,
)
from run_g4_sensitivity_fliprate import (  # noqa: E402
    JITTER_REPLICATES, level_class_map, p80, pair_list,
)
from run_g4_threshold_derivation import (  # noqa: E402
    DERIVATION, IDENTITY_ANCHOR, REF_BUDGET, STABILITY_RULE, S4_SPLIT,
    agg_level_axis, agg_value_axis, agg_window_axis, build_metric_index,
    pair_level_axis, pair_value_axis, pair_window_axis, window_query,
)

INPUTS = ROOT / 'configs/research/g4_s5_test_confirmation_inputs_v1.0.json'
ATTEMPT = ROOT / 'artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_attempt.json'
N_THRESHOLD = 1.0
MISSION_CLASS = 'mission_relevant'
WEAK_CLASS = 'weak_event_band'
DAMAGE_TYPES = ('amplitude_scale', 'polarity_flip', 'sample_shift', 'trace_deletion')
COMBOS = (('mt', 'c5'), ('mt', 'c8'))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def verify_inputs():
    frozen = read_json(INPUTS)
    for rel, expected in frozen['inputs_sha256'].items():
        actual = sha256(ROOT / rel)
        assert actual == expected, f'input SHA mismatch: {rel}: {actual} != {expected}'
    assert frozen['script'] == 'scripts/run_g4_test_confirmation.py'
    assert frozen['schema'] == 'g4-s5-test-confirmation-inputs/1'
    return frozen


def candidate_ids():
    ids = [c['id'] for c in catalogue()]
    ids += [f'L_lam{lam}_w{width}_q{gain}' for lam in L_LAMBDAS
            for width in L_WIDTHS for gain in L_GAINS]
    ids += [f'R_tau{tau}_q1' for tau in R_THRESHOLDS]
    assert len(ids) == 57 and len(set(ids)) == 57
    return sorted(ids)


def epsilon_index(fliprate):
    index = {}
    for row in fliprate['strata']:
        if row['quantity'] != 'D':
            continue
        key = (row['damage_type'], row['mission_class'], row['geometry'], row['family'])
        assert key not in index, f'duplicate dev epsilon stratum: {key}'
        index[key] = row['value_jitter']['epsilon']
    return index


def test_epsilon(eps_index, damage_type, mission_class, family):
    sources = {}
    for dev_family in ('c1', 'c3'):
        key = (damage_type, mission_class, 'mt', dev_family)
        assert key in eps_index, f'missing frozen dev epsilon source: {key}'
        sources[dev_family] = eps_index[key]
    return min(sources.values()), sources


def index_export(rows, key_fn, label):
    out = {}
    for row in rows:
        key = key_fn(row)
        assert key not in out, f'duplicate {label} row: {key}'
        out[key] = row
    return out


def split_levels(mapping):
    out = {MISSION_CLASS: {}, WEAK_CLASS: {}}
    for mission_class in (MISSION_CLASS, WEAK_CLASS):
        for damage_type, levels in mapping[mission_class].items():
            if isinstance(levels, list):
                out[mission_class][damage_type] = [1.0 if v == 'all' else v for v in levels]
    return out


def epsilon_rows(epsilon_map):
    return [
        {'damage_type': dt, 'mission_class': cls, 'geometry': geometry,
         'family': family, **value}
        for (dt, cls, geometry, family), value in sorted(epsilon_map.items())
    ]


def capability_a80(candidate_list, family, damage_type, levels, cap, epsilon):
    per_candidate = []
    for cid in candidate_list:
        per_level = []
        values = []
        for level in levels[damage_type]:
            row = cap.get((cid, damage_type, repr(level), 'mt', family, MISSION_CLASS))
            available = row is not None and row['n'] > 0 and row['D_p80'] is not None
            if not available:
                per_level.append({'level': level, 'available': False,
                                  'D_p80': None, 'within_epsilon': None, 'n': None})
                continue
            value = float(row['D_p80'])
            values.append(value)
            per_level.append({'level': level, 'available': True,
                              'D_p80': value, 'within_epsilon': value <= epsilon,
                              'n': row['n']})
        per_candidate.append({
            'candidate_id': cid,
            'a80_upper_bound_D_p80': max(values) if values else None,
            'a80_statement_complete': len(values) == len(levels[damage_type]),
            'per_level': per_level,
        })
    return per_candidate


def make_test_strata(records, class_map):
    """Same row shape as the dev flip-rate strata, with test-family filtering."""
    strata = {}
    for row in records:
        if (row['family'] not in TEST_FAMILIES or row['availability'] != 'ran'
                or row['row_type'] != 'event'):
            continue
        damage_type = row['damage']['type']
        level = row['damage']['level']
        mission_class = class_map.get((damage_type, level))
        if mission_class is None:
            mission_class = class_map.get((damage_type, 'all'))
        if mission_class is None:
            continue
        key = (damage_type, mission_class, row['geometry'], row['family'])
        strata.setdefault(key, {}).setdefault(row['candidate_id'], []).append(
            (level, row['event_id'], row['metrics']['waveform']['metrics']['nrmse']))
    return strata


def admissibility(candidate_list, family, cap, nc, levels, epsilon_map):
    d_rows = []
    n_rows = []
    accepted = []
    for cid in candidate_list:
        per_level = []
        d_pass = True
        for damage_type in DAMAGE_TYPES:
            eps = epsilon_map[(damage_type, MISSION_CLASS, 'mt', family)]['epsilon']
            for level in levels[damage_type]:
                key = (cid, damage_type, repr(level), 'mt', family, MISSION_CLASS)
                row = cap.get(key)
                available = row is not None and row['n'] > 0 and row['D_p80'] is not None
                value = None if not available else float(row['D_p80'])
                passed = bool(available and value <= eps)
                d_pass = d_pass and passed
                per_level.append({
                    'damage_type': damage_type, 'level': level,
                    'available': available, 'D_p80': value,
                    'epsilon': eps, 'pass': passed,
                })
        d_rows.append({'candidate_id': cid, 'pass': d_pass, 'per_level': per_level})

        q_rows = []
        n_pass = True
        for q in (2, 4):
            row = nc.get((cid, 'mt', family, q))
            available = row is not None and row['n'] > 0 and row['N_b_energy_ratio_p80'] is not None
            value = None if not available else float(row['N_b_energy_ratio_p80'])
            passed = bool(available and value <= N_THRESHOLD)
            n_pass = n_pass and passed
            q_rows.append({'q': q, 'available': available,
                           'N_b_energy_ratio_p80': value, 'threshold': N_THRESHOLD,
                           'pass': passed})
        n_rows.append({'candidate_id': cid, 'pass': n_pass, 'per_q': q_rows})
        if d_pass and n_pass:
            accepted.append(cid)
    return d_rows, n_rows, accepted


def stability_rows(strata, eps_index, admissible_by_combo):
    rows = []
    for key in sorted(strata):
        damage_type, mission_class, geometry, family = key
        if damage_type == 'nc_amplify':
            continue
        by_candidate = strata[key]
        cids = sorted(by_candidate)
        full_values = {cid: p80([v for _level, _event, v in by_candidate[cid]]) for cid in cids}
        pairs = pair_list(cids)
        levels, level_pairs = pair_level_axis(by_candidate, pairs)
        events, window_pairs, window_reps = pair_window_axis(by_candidate, full_values, pairs)
        eps, sources = test_epsilon(eps_index, damage_type, mission_class, family)
        value_pairs, value_reps = pair_value_axis(full_values, pairs, cids, repr(key), eps)
        # Aggregates mirror the frozen development implementation; verdicts use
        # only the zero-flip clauses from the frozen stability rule.
        level_summary = agg_level_axis(levels, level_pairs)
        window_summary = agg_window_axis(len(events), window_pairs, window_reps)
        value_summary = agg_value_axis(value_pairs, eps, {})
        admissible_here = sorted(set(admissible_by_combo[(geometry, family)]) & set(cids))
        relations = {}
        pair_rows = []
        for a, b in pair_list(admissible_here):
            la, wa, va = level_pairs[(a, b)], window_pairs[(a, b)], value_pairs[(a, b)]
            diff = full_values[a] - full_values[b]
            better = None if diff == 0 else (a if diff < 0 else b)
            determined = bool(diff != 0 and la['flips'] == 0 and wa['flips'] == 0 and va['flips'] == 0)
            relations[(a, b)] = better if determined else None
            pair_rows.append({
                'a': a, 'b': b, 'base_D_p80_a': full_values[a], 'base_D_p80_b': full_values[b],
                'better': better, 'level_axis_flips': la['flips'], 'level_axis_ties': la['ties'],
                'window_jackknife_flips': wa['flips'], 'window_jackknife_comparisons': wa['comparisons'],
                'window_jackknife_ties': wa['ties'], 'value_jitter_flips': va['flips'],
                'value_jitter_comparisons': va['comparisons'], 'value_jitter_ties': va['ties'],
                'determined': determined,
            })
        tops = [cid for cid in admissible_here if all(
            relations.get((min(cid, other), max(cid, other))) == cid
            for other in admissible_here if other != cid)]
        if not admissible_here:
            verdict = 'no_admissible_candidate'
        elif tops:
            assert len(tops) == 1, f'multiple stable tops in {key}: {tops}'
            verdict = 'unique_stable_top'
        elif any(row['determined'] for row in pair_rows):
            verdict = 'partial_only'
        else:
            verdict = 'set'
        rows.append({
            'stratum_key': ['D', damage_type, mission_class, geometry, family],
            'quantity': 'D', 'damage_type': damage_type, 'mission_class': mission_class,
            'geometry': geometry, 'family': family,
            'epsilon_value_jitter': eps,
            'epsilon_source_dev_strata': {
                dev_family: ['D', damage_type, mission_class, 'mt', dev_family]
                for dev_family in ('c1', 'c3')},
            'epsilon_source_values': sources,
            'n_candidates_in_stratum': len(cids),
            'admissible_present': admissible_here,
            'n_admissible_present': len(admissible_here),
            'n_determined_pairs': sum(x['determined'] for x in pair_rows),
            'n_pairs_admissible': len(pair_rows),
            'stability_verdict': verdict,
            'unique_stable_top_candidate': tops[0] if tops else None,
            'singleton_admissible_set_vacuously_unique': bool(tops) and len(admissible_here) == 1,
            'stability_applied': True,
            'cross_scale_caveat': False,
            'axis_summary': {'level_axis': level_summary,
                             'window_jackknife': window_summary,
                             'value_jitter': value_summary},
            'pair_rows': pair_rows,
        })
    return rows


def selftest():
    # Synthetic only: missing mission level fails; both N q rows are required.
    cids = ['A', 'B']
    cap = {('A', 'amplitude_scale', repr(0.5), 'mt', 'c5', MISSION_CLASS):
           {'n': 2, 'D_p80': 0.01}}
    nc = {('A', 'mt', 'c5', 2): {'n': 1, 'N_b_energy_ratio_p80': 1.0},
          ('A', 'mt', 'c5', 4): {'n': 1, 'N_b_energy_ratio_p80': 1.0}}
    eps = {(dt, MISSION_CLASS, 'mt', 'c5'): {'epsilon': 0.02}
           for dt in DAMAGE_TYPES}
    d, n, accepted = admissibility(cids, 'c5', cap, nc,
                                   {'amplitude_scale': [0.5, 0.1],
                                    'polarity_flip': [], 'sample_shift': [],
                                    'trace_deletion': []}, eps)
    assert d[0]['pass'] is False and d[0]['per_level'][1]['available'] is False
    assert n[0]['pass'] is True and n[1]['pass'] is False and accepted == []
    eps_dev = {
        ('amplitude_scale', MISSION_CLASS, 'mt', 'c1'): 0.03,
        ('amplitude_scale', MISSION_CLASS, 'mt', 'c3'): 0.02,
    }
    borrowed, sources = test_epsilon(eps_dev, 'amplitude_scale', MISSION_CLASS, 'c5')
    assert borrowed == 0.02 and sources == {'c1': 0.03, 'c3': 0.02}
    singleton = stability_rows({
        ('amplitude_scale', MISSION_CLASS, 'mt', 'c5'):
            {'A': [(0.5, 'e1', 0.1), (0.5, 'e2', 0.1)]}},
        {**eps_dev,
         ('amplitude_scale', MISSION_CLASS, 'mt', 'c5'): borrowed},
        {('mt', 'c5'): ['A']})[0]
    assert singleton['stability_verdict'] == 'unique_stable_top'
    assert singleton['singleton_admissible_set_vacuously_unique'] is True
    rows = epsilon_rows({('amplitude_scale', MISSION_CLASS, 'mt', 'c5'):
                         {'epsilon': borrowed, 'sources': sources}})
    json.dumps({'epsilon_source_rows': rows}, sort_keys=True)
    a80 = capability_a80(['A'], 'c5', 'amplitude_scale',
                         {'amplitude_scale': [0.5, 0.1]},
                         {('A', 'amplitude_scale', repr(0.5), 'mt', 'c5', MISSION_CLASS):
                          {'n': 0, 'D_p80': None}}, borrowed)
    assert a80[0]['a80_statement_complete'] is False
    assert a80[0]['per_level'][0]['available'] is False
    levels = {'amplitude_scale': [0.5, 0.1], 'polarity_flip': [1.0],
              'sample_shift': [-16], 'trace_deletion': [4]}
    complete_cap = {}
    for dt, damage_levels in levels.items():
        for level in damage_levels:
            complete_cap[('X', dt, repr(level), 'mt', 'c5', MISSION_CLASS)] = {
                'n': 3, 'D_p80': 0.02 if (dt, level) == ('amplitude_scale', 0.5) else 0.01}
    complete_nc = {('X', 'mt', 'c5', 2): {'n': 2, 'N_b_energy_ratio_p80': 0.9},
                   ('X', 'mt', 'c5', 4): {'n': 2, 'N_b_energy_ratio_p80': 1.0}}
    complete_eps = {(dt, MISSION_CLASS, 'mt', 'c5'): {'epsilon': 0.02}
                    for dt in DAMAGE_TYPES}
    d, n, accepted = admissibility(['X'], 'c5', complete_cap, complete_nc, levels, complete_eps)
    assert d[0]['pass'] and n[0]['pass'] and accepted == ['X']
    missing_q4 = dict(complete_nc)
    del missing_q4[('X', 'mt', 'c5', 4)]
    assert admissibility(['X'], 'c5', complete_cap, missing_q4, levels, complete_eps)[2] == []
    high_q4 = dict(complete_nc)
    high_q4[('X', 'mt', 'c5', 4)] = {'n': 2, 'N_b_energy_ratio_p80': 1.000001}
    assert admissibility(['X'], 'c5', complete_cap, high_q4, levels, complete_eps)[2] == []
    high_d = dict(complete_cap)
    high_d[('X', 'amplitude_scale', repr(0.5), 'mt', 'c5', MISSION_CLASS)] = {
        'n': 3, 'D_p80': 0.020001}
    assert admissibility(['X'], 'c5', high_d, complete_nc, levels, complete_eps)[2] == []
    print('synthetic selftest passed: complete pass at D=epsilon/N=1, just-over failures, missing D/q4 failures, epsilon min/provenance, singleton vacuity, serialization')


def preflight():
    frozen = verify_inputs()
    derivation = read_json(DERIVATION)
    mission = read_json(MISSION)
    stability_rule = read_json(STABILITY_RULE)
    split = read_json(S4_SPLIT)
    fliprate = read_json(ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json')
    capability = read_json(ROOT / 'artifacts/research_checks/2026-09-27_g4_mission_capability_test_r1/results.json')
    assert derivation['schema_version'] == 'g4_threshold_derivation_v1.0'
    assert derivation['status'] == 'frozen_program_only_no_threshold_values_derived'
    assert mission['status'] == stability_rule['status'] == split['status'] == 'frozen'
    assert derivation['derivation_rules']['D_constraint']['formula']
    assert derivation['derivation_rules']['N_constraint']['semantics'].startswith('N_th = 1.0')
    assert split['groups']['test'] == ['B2D-C5m', 'B2D-C8m']
    assert split['groups']['dev'] == ['B2D-C1m', 'B2D-C3m']
    assert fliprate['n_strata'] == len(fliprate['strata']) == 24
    assert capability['split'] == 'test' and capability['families'] == ['c5', 'c8']
    assert capability['ladder_records_merged'] == 12768
    candidates = candidate_ids()
    ladder = load_chunks(TEST_CHUNK_GLOB, run='r1')
    validate_records(ladder, TEST_FAMILIES)
    assert len(ladder) == 12768
    assert {row['candidate_id'] for row in ladder} == set(candidates)
    assert {row['geometry'] for row in capability['capability']} == {'mt'}
    assert {row['family'] for row in capability['capability']} == {'c5', 'c8'}
    combos = {(row['geometry'], row['family']) for row in capability['capability']}
    assert combos == set(COMBOS), f'test export combinations changed: {sorted(combos)}'
    id_cap = [row for row in capability['capability']
              if row['candidate_id'] == IDENTITY_ANCHOR and row['mission_class'] == MISSION_CLASS]
    id_nc = [row for row in capability['negative_control']
             if row['candidate_id'] == IDENTITY_ANCHOR]
    assert len(id_cap) == 18 and all(row['D_p80'] == 0.0 for row in id_cap)
    assert len(id_nc) == 4 and all(row['N_b_energy_ratio_p80'] == 1.0 for row in id_nc)
    mapping = split_levels(mission['ladder_level_mapping'])
    eps_dev = epsilon_index(fliprate)
    class_map = level_class_map(mission['ladder_level_mapping'])
    strata = make_test_strata(ladder, class_map)
    expected = {(dt, cls, 'mt', family)
                for family in ('c5', 'c8')
                for cls, by_type in mapping.items()
                for dt in by_type}
    assert len(expected) == 14 and set(strata) == expected
    for dt, cls, _geometry, family in expected:
        test_epsilon(eps_dev, dt, cls, family)
    return frozen, derivation, mission, stability_rule, split, fliprate, capability, ladder, candidates, eps_dev, mapping, strata


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', help='new output directory; one-shot formal run only')
    ap.add_argument('--selftest', action='store_true', help='synthetic logic check; does not read test outcomes')
    ap.add_argument('--preflight', action='store_true', help='read-only hash/status/coverage checks; no threshold evaluation')
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if args.preflight:
        prepared = preflight()
        print('static preflight passed: input hashes/statuses, 12768 records, 57 candidates, '
              '14 D strata, identity and epsilon sources; no thresholds applied')
        return
    assert args.output_dir, '--output-dir required for formal S5 execution'
    output = Path(args.output_dir)
    if not output.is_absolute():
        output = ROOT / output
    output.relative_to(ROOT)
    assert not output.exists(), f'output directory already exists: {output}'
    prepared = preflight()
    assert not ATTEMPT.exists(), f'S5 attempt already exists; never rerun: {ATTEMPT}'
    (frozen, derivation, mission, stability_rule, split, fliprate,
     capability, ladder, candidates, eps_dev, levels_by_class, strata) = prepared
    ATTEMPT.parent.mkdir(parents=True, exist_ok=True)
    attempt = {
        'schema': 'g4-s5-one-shot-attempt/1', 'attempt': 1,
        'script': 'scripts/run_g4_test_confirmation.py',
        'script_sha256': sha256(__file__), 'input_manifest_sha256': sha256(INPUTS),
        'inputs_sha256': frozen['inputs_sha256'],
        'command': [sys.executable, str(Path(__file__).relative_to(ROOT)),
                    '--output-dir', str(output.relative_to(ROOT))],
        'threshold_application_started': True,
        'failed_attempts_are_consumed': True,
    }
    with ATTEMPT.open('x', encoding='utf-8') as handle:
        json.dump(attempt, handle, indent=2, sort_keys=True)
        handle.write('\n')

    started = time.perf_counter()
    budget = read_json(REF_BUDGET)

    capability_index = index_export(capability['capability'],
        lambda r: (r['candidate_id'], r['damage_type'], repr(r['damage_level']),
                   r['geometry'], r['family'], r['mission_class']), 'capability')
    nc_index = index_export(capability['negative_control'],
        lambda r: (r['candidate_id'], r['geometry'], r['family'], int(r['nc_amplify_q'])),
        'negative-control')
    if {r['candidate_id'] for r in capability['capability']} | {
            r['candidate_id'] for r in capability['negative_control']} != set(candidates):
        raise ValueError('test capability export does not cover the frozen 57-candidate set')

    mission_levels = levels_by_class[MISSION_CLASS]
    epsilon_map = {}
    for family in ('c5', 'c8'):
        for damage_type in DAMAGE_TYPES:
            eps, sources = test_epsilon(eps_dev, damage_type, MISSION_CLASS, family)
            epsilon_map[(damage_type, MISSION_CLASS, 'mt', family)] = {
                'epsilon': eps, 'sources': sources}
    combo_results = []
    admissible_by_combo = {}
    for _geometry, family in COMBOS:
        d_rows, n_rows, accepted = admissibility(
            candidates, family, capability_index, nc_index, mission_levels, epsilon_map)
        assert len(d_rows) == len(n_rows) == len(candidates) == 57
        identity_d = next(row for row in d_rows if row['candidate_id'] == IDENTITY_ANCHOR)
        identity_n = next(row for row in n_rows if row['candidate_id'] == IDENTITY_ANCHOR)
        assert identity_d['pass'] and all(cell['available'] and cell['D_p80'] == 0.0
                                           for cell in identity_d['per_level'])
        assert identity_n['pass'] and all(cell['available'] and
                                          cell['N_b_energy_ratio_p80'] == 1.0
                                          for cell in identity_n['per_q'])
        admissible_by_combo[('mt', family)] = accepted
        combo_results.append({
            'geometry': 'mt', 'family': family,
            'n_candidates': 57, 'n_D_pass': sum(r['pass'] for r in d_rows),
            'n_N_pass': sum(r['pass'] for r in n_rows),
            'n_admissible': len(accepted), 'admissible_candidates': accepted,
            'epsilon_by_damage_type': {
                dt: epsilon_map[(dt, MISSION_CLASS, 'mt', family)] for dt in DAMAGE_TYPES},
            'D_constraint': {'mission_class_only': MISSION_CLASS,
                             'mission_levels': mission_levels, 'per_candidate': d_rows},
            'N_constraint': {'threshold': N_THRESHOLD,
                             'joint_reading': 'both q=2 and q=4 available and <= 1.0',
                             'per_candidate': n_rows},
        })

    stable = stability_rows(strata, eps_dev, admissible_by_combo)
    assert len(stable) == 14 and all(r['quantity'] == 'D' for r in stable)

    arrival = build_metric_index([{
        'event_id': row['event_id'], 'family': row['family'],
        'value': row['arrival_time_class']['direct_evidence']['grid_tier_shift_ns'],
        'direction': row['arrival_time_class']['direct_evidence']['direction']}
        for row in budget['per_event_window']], 'value', 'direction')
    capability_statements = []
    for _geometry, family in COMBOS:
        event_ids = sorted({r['event_id'] for r in ladder
                            if r['family'] == family and r['row_type'] == 'event'})
        arrival_summary = window_query(arrival, event_ids, family)
        for damage_type in DAMAGE_TYPES:
            eps = epsilon_map[(damage_type, MISSION_CLASS, 'mt', family)]['epsilon']
            per_candidate = capability_a80(candidates, family, damage_type,
                                            mission_levels, capability_index, eps)
            capability_statements.append({
                'geometry': 'mt', 'family': family, 'damage_type': damage_type,
                'mission_class': MISSION_CLASS, 'headline_excludes_weak_event_band': True,
                'statistical_term': 'a80 upper bound of erasure from constructed reference',
                'mission_tolerance_epsilon': eps,
                'systematic_bias_term': {
                    'grid_chain_direction_median_ns': arrival_summary['median'],
                    'grid_chain_direction_median_abs_ns': (None if arrival_summary['median'] is None
                                                           else abs(arrival_summary['median'])),
                    'directions_observed': arrival_summary['directions'],
                    'n_event_windows': arrival_summary['n_event_windows'],
                    'unit': 'ns', 'reported_separately': True,
                    'mixed_into_statistical_term': False, 'depth_conversion_performed': False,
                },
                'per_candidate': per_candidate,
            })

    result = {
        'schema': 'g4-s5-test-confirmation/1',
        'program': 'configs/research/g4_threshold_derivation_v1.0.json',
        'program_sha256': sha256(DERIVATION),
        'one_shot_attempt': str(ATTEMPT.relative_to(ROOT)),
        'g4_relieved': False, 'solver_invoked': False, 'training_eligible': False,
        'thresholds_are_simulation_domain_only': True,
        'constructed_reference_is_physical_truth': False,
        'no_ranking_no_selection_no_training': True,
        'test_split': {'families': ['c5', 'c8'], 'geometries_present': ['mt'],
                       'co_rows_absent_by_frozen_event_coverage': True},
        'candidate_set': {'n_candidates': 57, 'candidate_ids': candidates,
                          'unavailable_candidates_retained_and_fail_closed': True},
        'epsilon_rule': ('per test (damage_type, mission_class, geometry=mt, family) stratum: '
                         'min of frozen dev mt/c1 and mt/c3 value_jitter epsilon; conservative '
                         'development calibration reuse, not test recalibration'),
        'epsilon_source_rows': epsilon_rows(epsilon_map),
        'N_threshold': N_THRESHOLD,
        'frozen_derivation_rules': derivation['derivation_rules'],
        'stability_rule': stability_rule['determined_order_rule'],
        'combo_summary': combo_results,
        'strata_stability': stable,
        'N_b_stability': {'applied': False, 'reason': 'cross-scale caveat; energy ratio vs NRMSE epsilon'},
        'stability_interpretation': ('A singleton admissible set is reported as vacuously unique under the '
                                     'frozen pairwise rule; this does not establish an algorithm winner or '
                                     'superiority.'),
        'capability_statements': capability_statements,
        'hard_limits': derivation['hard_limits'],
        'disclosures': derivation['disclosures'],
        'formal_rules_applied_once': True,
    }
    assert len(capability_statements) == 8
    assert all(len(block['per_candidate']) == 57 for block in capability_statements)
    output.mkdir(parents=True)
    result_path = output / 'results.json'
    result_path.write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + '\n',
                           encoding='utf-8')
    run_manifest = {
        'attempt_marker_sha256': sha256(ATTEMPT),
        'script_sha256': sha256(__file__), 'input_manifest_sha256': sha256(INPUTS),
        'inputs_sha256': frozen['inputs_sha256'], 'results_sha256': sha256(result_path),
        'one_shot': True, 'solver_invoked': False, 'training_eligible': False,
        'g4_relieved': False, 'n_records': len(ladder),
        'n_candidates': len(candidates), 'n_combos': len(COMBOS),
        'n_d_stability_strata': len(stable), 'n_n_b_stability_strata': 0,
        'admissible_counts': {r['family']: r['n_admissible'] for r in combo_results},
        'stability_verdict_histogram': {
            verdict: sum(r['stability_verdict'] == verdict for r in stable)
            for verdict in sorted({r['stability_verdict'] for r in stable})},
        'total_wall_s': time.perf_counter() - started,
    }
    (output / 'run_manifest.json').write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True, ensure_ascii=False) + '\n', encoding='utf-8')
    print('wrote', result_path)
    print('sha256:', sha256(result_path))
    print('admissible per combo:', run_manifest['admissible_counts'])
    print('stability verdicts:', run_manifest['stability_verdict_histogram'])


if __name__ == '__main__':
    main()
