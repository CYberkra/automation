"""Root acceptance: recompute S5 D/N from raw archived ladder, not runner helpers."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).parent


def read(rel):
    return json.loads((ROOT / rel).read_text(encoding='utf-8'))


def p80(values):
    ordered = sorted(values)
    position = (len(ordered) - 1) * 0.8
    lower = int(position)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[min(lower + 1, len(ordered) - 1)] * fraction


def main():
    result_rel = 'artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_r1/results.json'
    result = read(result_rel)
    mission = read('configs/research/g4_mission_tolerance_v0.2.json')
    levels = {dt: [1.0 if v == 'all' else v for v in vals]
              for dt, vals in mission['ladder_level_mapping']['mission_relevant'].items()
              if isinstance(vals, list)}
    flip = read('artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json')
    eps = {dt: min(row['value_jitter']['epsilon'] for row in flip['strata']
                   if row['quantity'] == 'D' and row['damage_type'] == dt
                   and row['mission_class'] == 'mission_relevant'
                   and row['geometry'] == 'mt' and row['family'] in ('c1', 'c3'))
           for dt in levels}
    d, n = defaultdict(list), defaultdict(list)
    records = []
    for chunk in range(8):
        records.extend(read(f'artifacts/research_checks/2026-09-27_damage_ladder_test_r1_c{chunk:02d}/records.json')['records'])
    assert len(records) == 12768
    candidates = sorted({r['candidate_id'] for r in records})
    assert len(candidates) == 57 and candidates == result['candidate_set']['candidate_ids']
    for row in records:
        if row['availability'] != 'ran':
            continue
        fam, cid = row['family'], row['candidate_id']
        dt, level = row['damage']['type'], row['damage']['level']
        assert row['geometry'] == 'mt'
        if row['row_type'] == 'event' and dt in levels and level in levels[dt]:
            d[fam, cid, dt, level].append(row['metrics']['waveform']['metrics']['nrmse'])
        elif row['row_type'] == 'negative_control' and dt == 'nc_amplify':
            n[fam, cid, level].append(row['metrics']['N_b']['energy_ratio'])
    counts, cells = {}, 0
    for combo in result['combo_summary']:
        fam = combo['family']
        declared_d = {row['candidate_id']: row for row in combo['D_constraint']['per_candidate']}
        declared_n = {row['candidate_id']: row for row in combo['N_constraint']['per_candidate']}
        assert combo['D_constraint']['mission_levels'] == levels
        accepted = []
        for cid in candidates:
            d_flags, n_flags = [], []
            for cell in declared_d[cid]['per_level']:
                dt, level = cell['damage_type'], cell['level']
                values = d[fam, cid, dt, level]
                value = p80(values) if values else None
                assert cell['available'] == bool(values)
                assert cell['D_p80'] == value and cell['epsilon'] == eps[dt]
                passed = bool(values) and value <= eps[dt]
                assert cell['pass'] == passed
                d_flags.append(passed)
                cells += 1
            assert len(d_flags) == sum(map(len, levels.values())) == 9
            for cell in declared_n[cid]['per_q']:
                values = n[fam, cid, cell['q']]
                value = p80(values) if values else None
                assert cell['available'] == bool(values)
                assert cell['N_b_energy_ratio_p80'] == value and cell['threshold'] == 1.0
                passed = bool(values) and value <= 1.0
                assert cell['pass'] == passed
                n_flags.append(passed)
                cells += 1
            assert len(n_flags) == 2
            assert declared_d[cid]['pass'] == all(d_flags)
            assert declared_n[cid]['pass'] == all(n_flags)
            if all(d_flags + n_flags):
                accepted.append(cid)
        assert accepted == combo['admissible_candidates']
        counts[fam] = {'admissible': accepted, 'n_D_pass': combo['n_D_pass'], 'n_N_pass': combo['n_N_pass']}
    assert set(counts) == {'c5', 'c8'} and cells == 1254
    summary = {'status': 'pass', 'raw_records': len(records), 'candidate_family_pairs': 114,
               'independently_recomputed_D_N_cells': cells, 'counts': counts,
               'result_sha256': hashlib.sha256((ROOT / result_rel).read_bytes()).hexdigest(),
               'formal_runner_reexecuted': False, 'thresholds_changed': False}
    (OUT / 'results.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
