"""Independent structural and arithmetic acceptance for one S5 test run."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / 'configs/research/g4_s5_test_confirmation_inputs_v1.0.json'
CAPABILITY = ROOT / 'artifacts/research_checks/2026-09-27_g4_mission_capability_test_r1/results.json'
FLIPRATE = ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json'
DAMAGE_TYPES = ('amplitude_scale', 'polarity_flip', 'sample_shift', 'trace_deletion')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def check(output_dir):
    output = Path(output_dir)
    if not output.is_absolute():
        output = ROOT / output
    result_path = output / 'results.json'
    manifest_path = output / 'run_manifest.json'
    result = read_json(result_path)
    manifest = read_json(manifest_path)
    attempt_path = ROOT / result['one_shot_attempt']
    attempt = read_json(attempt_path)
    frozen = read_json(INPUTS)
    assert attempt['attempt'] == 1 and attempt['threshold_application_started'] is True
    assert attempt['script_sha256'] == sha256(ROOT / 'scripts/run_g4_test_confirmation.py')
    assert attempt['input_manifest_sha256'] == sha256(INPUTS)
    assert manifest['attempt_marker_sha256'] == sha256(attempt_path)
    assert manifest['script_sha256'] == attempt['script_sha256']
    assert manifest['input_manifest_sha256'] == attempt['input_manifest_sha256']
    assert manifest['results_sha256'] == sha256(result_path)
    for rel, expected in frozen['inputs_sha256'].items():
        assert sha256(ROOT / rel) == expected, f'input SHA changed: {rel}'
        assert attempt['inputs_sha256'][rel] == expected

    assert result['schema'] == 'g4-s5-test-confirmation/1'
    assert result['formal_rules_applied_once'] is True
    assert result['g4_relieved'] is False and result['solver_invoked'] is False
    assert result['training_eligible'] is False and result['no_ranking_no_selection_no_training'] is True
    candidates = result['candidate_set']['candidate_ids']
    assert len(candidates) == 57 and len(set(candidates)) == 57
    cap = read_json(CAPABILITY)
    flip = read_json(FLIPRATE)
    epsilon_dev = {}
    for row in flip['strata']:
        if row['quantity'] == 'D':
            key = (row['damage_type'], row['mission_class'], row['geometry'], row['family'])
            assert key not in epsilon_dev
            epsilon_dev[key] = row['value_jitter']['epsilon']
    cap_index = {}
    for row in cap['capability']:
        key = (row['candidate_id'], row['damage_type'], repr(row['damage_level']),
               row['geometry'], row['family'], row['mission_class'])
        assert key not in cap_index
        cap_index[key] = row
    nc_index = {}
    for row in cap['negative_control']:
        key = (row['candidate_id'], row['geometry'], row['family'], int(row['nc_amplify_q']))
        assert key not in nc_index
        nc_index[key] = row

    assert {r['family'] for r in cap['capability']} == {'c5', 'c8'}
    assert {r['geometry'] for r in cap['capability']} == {'mt'}
    assert len(result['combo_summary']) == 2
    for combo in result['combo_summary']:
        geometry, family = combo['geometry'], combo['family']
        assert geometry == 'mt' and family in ('c5', 'c8')
        by_d = {row['candidate_id']: row for row in combo['D_constraint']['per_candidate']}
        by_n = {row['candidate_id']: row for row in combo['N_constraint']['per_candidate']}
        assert set(by_d) == set(by_n) == set(candidates)
        levels = combo['D_constraint']['mission_levels']
        expected_count = sum(len(levels[dt]) for dt in DAMAGE_TYPES)
        accepted = []
        for cid in candidates:
            all_d = True
            seen_cells = []
            for dt in DAMAGE_TYPES:
                dev_sources = {d: epsilon_dev[(dt, 'mission_relevant', 'mt', d)]
                               for d in ('c1', 'c3')}
                epsilon = min(dev_sources.values())
                declared = combo['epsilon_by_damage_type'][dt]
                assert declared['epsilon'] == epsilon
                assert declared['sources'] == dev_sources
                for level in levels[dt]:
                    row = cap_index.get((cid, dt, repr(level), 'mt', family, 'mission_relevant'))
                    available = row is not None and row['n'] > 0 and row['D_p80'] is not None
                    passed = bool(available and row['D_p80'] <= epsilon)
                    seen_cells.append((dt, level, available, passed))
                    all_d = all_d and passed
            assert len(seen_cells) == expected_count
            out_d = by_d[cid]
            assert out_d['pass'] == all_d
            assert len(out_d['per_level']) == expected_count
            assert [(x['damage_type'], x['level'], x['available'], x['pass'])
                    for x in out_d['per_level']] == seen_cells

            q_pass = []
            for q in (2, 4):
                row = nc_index.get((cid, 'mt', family, q))
                available = row is not None and row['n'] > 0 and row['N_b_energy_ratio_p80'] is not None
                q_pass.append(bool(available and row['N_b_energy_ratio_p80'] <= 1.0))
            all_n = all(q_pass)
            assert by_n[cid]['pass'] == all_n
            assert [x['pass'] for x in by_n[cid]['per_q']] == q_pass
            if all_d and all_n:
                accepted.append(cid)
        assert combo['admissible_candidates'] == accepted
        assert combo['n_admissible'] == len(accepted)

    stable = result['strata_stability']
    assert len(stable) == 14 and all(row['quantity'] == 'D' for row in stable)
    seen_keys = set()
    for row in stable:
        key = (row['quantity'], row['damage_type'], row['mission_class'],
               row['geometry'], row['family'])
        assert key not in seen_keys
        seen_keys.add(key)
        assert row['epsilon_source_values']['c1'] == epsilon_dev[
            (row['damage_type'], row['mission_class'], 'mt', 'c1')]
        assert row['epsilon_source_values']['c3'] == epsilon_dev[
            (row['damage_type'], row['mission_class'], 'mt', 'c3')]
        assert row['epsilon_value_jitter'] == min(row['epsilon_source_values'].values())
        assert row['n_pairs_admissible'] == len(row['pair_rows'])
        for pair in row['pair_rows']:
            determined = (pair['base_D_p80_a'] != pair['base_D_p80_b']
                          and pair['level_axis_flips'] == 0
                          and pair['window_jackknife_flips'] == 0
                          and pair['value_jitter_flips'] == 0)
            assert pair['determined'] == determined
        count = row['n_admissible_present']
        assert count == len(row['admissible_present'])
        assert row['singleton_admissible_set_vacuously_unique'] is (count == 1 and
                                                                      row['stability_verdict'] == 'unique_stable_top')
        if count == 0:
            assert row['stability_verdict'] == 'no_admissible_candidate'
        elif count == 1:
            assert row['stability_verdict'] == 'unique_stable_top'
            assert row['unique_stable_top_candidate'] == row['admissible_present'][0]
            assert row['pair_rows'] == []
    assert result['N_b_stability']['applied'] is False
    assert len(result['capability_statements']) == 8
    for block in result['capability_statements']:
        assert block['headline_excludes_weak_event_band'] is True
        assert block['systematic_bias_term']['reported_separately'] is True
        assert block['systematic_bias_term']['mixed_into_statistical_term'] is False
        assert block['systematic_bias_term']['depth_conversion_performed'] is False
        assert len(block['per_candidate']) == 57
    print('S5 acceptance passed: frozen one-shot record, 57 x 2 D/N arithmetic, epsilon sources, 14 stability strata, a80/systematic separation')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', required=True)
    args = ap.parse_args()
    check(args.output_dir)


if __name__ == '__main__':
    main()
