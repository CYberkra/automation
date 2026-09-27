"""G4 step 5 S2: execution of the frozen threshold-derivation program on the
development group data.

Standard of proof (口径来源): the whole procedure is pinned by
`configs/research/g4_threshold_derivation_v1.0.json`
(SHA-256 d114fd193259008a4882ba6a3e3af3ab9ed507fd608fe3cf72a2432d99283eca).
That file is the ONLY specification source and is hash-gated at start-up; no
rule, margin or threshold is re-decided here.

What this script does, and does not do:

  * applies the frozen D constraint (mission-relevant ladder levels vs the
    flip-rate epsilon of the same stratum) and N constraint (conservative joint
    reading of q in {2,4} against N_th = 1.0, the identity anchor constant);
  * reports admissibility per (geometry x family), never relaxing D_th when a
    stratum ends up with no admissible candidate;
  * applies the frozen stability rule ONCE, on the flip-rate v2 methodology,
    and ONLY inside the 21 D strata; the 3 N_b strata get a reported
    cross-scale caveat line instead of a per-pair verdict;
  * reports stratum structure (`unique_stable_top` / `partial_only` / `set`)
    and a80-style capability upper bounds with an explicitly SEPARATED
    systematic bias term.

NO ranking, NO selection, NO training (不排名、不选优、不训练), NO labels, NO solver
invocation, NO relaxation. Every input is read-only and sha256-gated; any
mutation of a gated input aborts the run.

This is the S2 draft implemented per docs/research/task_codebuddy_g4_s2_derivation.md
and awaits kimi acceptance (independent recomputation + regression).
"""
import argparse
import hashlib
import json
import math
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

# --- frozen methodology is reused verbatim from the step-4 flip-rate runner ---
from run_g4_sensitivity_fliprate import (  # noqa: E402
    CHUNK_GLOB, DEV_FAMILIES, EVENT_TABLE, EVENT_TABLE_SHA256, IDENTITY_ANCHOR,
    JITTER_REPLICATES, MISSION, MISSION_SHA256, N_CHUNKS, N_RECORDS, NC_TYPE,
    REF_BUDGET, REF_BUDGET_SHA256, build_strata, level_class_map, level_axis,
    p80, pair_list, quantile, sha256_file, value_axis, window_axis,
)

ROOT = Path(__file__).resolve().parents[1]

DERIVATION = ROOT / 'configs/research/g4_threshold_derivation_v1.0.json'
DERIVATION_SHA256 = 'd114fd193259008a4882ba6a3e3af3ab9ed507fd608fe3cf72a2432d99283eca'
MISSION_TOLERANCE = MISSION
MISSION_TOLERANCE_SHA256 = MISSION_SHA256
STABILITY_RULE = ROOT / 'configs/research/g4_fliprate_stability_rule_v0.1.json'
STABILITY_RULE_SHA256 = '8ee61d0b7777f188d84fd754e3302ade9fe678ed0eae1b7f3f0ade4c451b65f6'
S4_SPLIT = ROOT / 'configs/research/s4_baseline_group_v0.1.json'
S4_SPLIT_SHA256 = '8676f37d6aee603b7f2481779ffe361d2878735316863da0bed0db82e16b5593'
CAPABILITY = ROOT / 'artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r1/results.json'
CAPABILITY_SHA256 = '2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8'
FLIPRATE = ROOT / 'artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json'
FLIPRATE_SHA256 = 'f96d317a3e046dfcc9efb01da55d49f78ab09c6fdb562c37e4ce93c6480e558f'
LADDER_MERGED_STRIPPED_SHA256 = '6bc13f63c5985fa0b38b4724c2ffa0cf61716053751746520766ef6fa0e4d0d9'

DAMAGE_TYPES = ('amplitude_scale', 'polarity_flip', 'sample_shift', 'trace_deletion')
MISSION_RELEVANT = 'mission_relevant'
WEAK_EVENT_BAND = 'weak_event_band'
NC_QS = (2, 4)
N_THRESHOLD = 1.0
LEVEL_ALIAS = {'all': 1.0}
PLAN = 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md'
STEP = ('G4 step 5 S2 (frozen threshold-derivation program v1.0 executed once on the '
        'development group)')


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def assert_sha(path, expected, label):
    got = sha256_file(path)
    assert got == expected, f'{label} gate: sha256 mismatch ({got} != {expected})'
    return got


def load_ladder_chunks():
    """Merge logic copied verbatim from scripts/check_damage_ladder_acceptance.py merged()."""
    recs = []
    for k in range(N_CHUNKS):
        d = load_json(ROOT / CHUNK_GLOB.format(k) / 'records.json')
        for r in d['records']:
            r.pop('resource', None)
        recs.extend(d['records'])
    stripped = json.dumps(recs, sort_keys=True, ensure_ascii=False).encode('utf-8')
    got = hashlib.sha256(stripped).hexdigest()
    assert got == LADDER_MERGED_STRIPPED_SHA256, (
        f'damage-ladder merged-resource-stripped gate: sha256 mismatch ({got})')
    return recs


def canonical_level(level):
    if isinstance(level, str):
        return LEVEL_ALIAS[level] if level in LEVEL_ALIAS else float(level)
    return float(level)


def frozen_levels(mapping):
    return {dt: [canonical_level(v) for v in levels] for dt, levels in mapping.items()}


def index_capability(rows):
    idx = {}
    dupes = []
    for r in rows:
        key = (r['candidate_id'], r['damage_type'], canonical_level(r['damage_level']),
               r['geometry'], r['family'], r['mission_class'])
        if key in idx:
            dupes.append(key)
        idx[key] = r
    return idx, dupes


def index_negative_control(rows):
    idx = {}
    dupes = []
    for r in rows:
        key = (r['candidate_id'], r['geometry'], r['family'], int(r['nc_amplify_q']))
        if key in idx:
            dupes.append(key)
        idx[key] = r
    return idx, dupes


# --------------------------------------------------------------------------
# per-pair re-derivation of the three flip axes
#
# These return PERT-PAIR flip/tie counts. Summing them over the pairs (and, for
# the window/jitter axes, over replicates) reproduces the v2 aggregate counts
# bit-for-bit; both checks are executed against v2 in main(). The Pair ordering
# convention is the v2 one: pairs come from the sorted candidate list, and the
# tie convention is exact float equality, excluded from the flip denominators.
# --------------------------------------------------------------------------
def pair_diff_state(base_diff, other_values, a, b):
    """Return 'flip' / 'no_flip' / 'tie' / 'skip' for one pair vs one replicate."""
    if base_diff == 0.0:
        return 'tie'
    oa = other_values.get(a)
    ob = other_values.get(b)
    if oa is None or ob is None:
        return 'skip'
    other_diff = oa - ob
    if other_diff == 0.0:
        return 'tie'
    return 'flip' if (base_diff > 0) != (other_diff > 0) else 'no_flip'


def pair_level_axis(records_by_cid, pairs):
    by_level = {}
    for cid, rows in records_by_cid.items():
        per_level = {}
        for lvl, _event_id, val in rows:
            per_level.setdefault(lvl, []).append(val)
        for lvl, vals in per_level.items():
            by_level.setdefault(lvl, {})[cid] = p80(vals)
    levels = sorted(by_level, key=repr)
    out = {}
    for a, b in pairs:
        signs = []
        ties = 0
        for lvl in levels:
            va = by_level[lvl].get(a)
            vb = by_level[lvl].get(b)
            if va is None or vb is None:
                continue
            diff = va - vb
            if diff == 0.0:
                ties += 1
                continue
            signs.append(diff > 0)
        comparable = len(signs) >= 2
        out[(a, b)] = {
            'comparable': comparable,
            'flips': 1 if (comparable and len(set(signs)) > 1) else 0,
            'ties': ties,
        }
    return levels, out


def agg_level_axis(levels, per_pair):
    flipped = sum(p['flips'] for p in per_pair.values())
    comparable = sum(1 for p in per_pair.values() if p['comparable'])
    ties = sum(p['ties'] for p in per_pair.values())
    if len(levels) < 2 or comparable == 0:
        return {'applicable': False, 'comparisons': 0, 'flip_rate': None,
                'flipped_pairs': 0, 'n_levels': len(levels), 'ties': 0}
    return {'applicable': True, 'comparisons': comparable,
            'flip_rate': flipped / comparable, 'flipped_pairs': flipped,
            'n_levels': len(levels), 'ties': ties}


def pair_window_axis(records_by_cid, full_values, pairs):
    events = sorted({event_id for rows in records_by_cid.values()
                     for _lvl, event_id, _v in rows})
    per_pair = {p: {'flips': 0, 'comparisons': 0, 'ties': 0} for p in pairs}
    reps = []
    for event in events:
        jack = {}
        for cid, rows in records_by_cid.items():
            vals = [v for _lvl, e, v in rows if e != event]
            if vals:
                jack[cid] = p80(vals)
        r_flips = 0
        r_comparisons = 0
        r_ties = 0
        for a, b in pairs:
            state = pair_diff_state(full_values[a] - full_values[b], jack, a, b)
            if state == 'tie':
                per_pair[(a, b)]['ties'] += 1
                r_ties += 1
            elif state == 'skip':
                continue
            else:
                per_pair[(a, b)]['comparisons'] += 1
                r_comparisons += 1
                if state == 'flip':
                    per_pair[(a, b)]['flips'] += 1
                    r_flips += 1
        reps.append((r_flips, r_comparisons, r_ties))
    return events, per_pair, reps


def agg_window_axis(n_events, per_pair, reps):
    flips = sum(p['flips'] for p in per_pair.values())
    comparisons = sum(p['comparisons'] for p in per_pair.values())
    ties = sum(p['ties'] for p in per_pair.values())
    if n_events < 2 or comparisons == 0:
        return {'applicable': False, 'comparisons': 0, 'flip_rate': None, 'flips': 0,
                'max_replicate_flip_rate': None, 'n_replicates': n_events, 'ties': 0}
    max_rate = None
    for r_flips, r_comparisons, _r_ties in reps:
        if r_comparisons > 0:
            rate = r_flips / r_comparisons
            if max_rate is None or rate > max_rate:
                max_rate = rate
    return {'applicable': True, 'comparisons': comparisons,
            'flip_rate': flips / comparisons, 'flips': flips,
            'max_replicate_flip_rate': max_rate, 'n_replicates': n_events, 'ties': ties}


def pair_value_axis(full_values, pairs, cids, stratum_key, epsilon):
    per_pair = {p: {'flips': 0, 'comparisons': 0, 'ties': 0} for p in pairs}
    reps = []
    for k in range(JITTER_REPLICATES):
        jittered = {}
        for cid in cids:
            rng = random.Random(f'{stratum_key}|{k}|{cid}')
            jittered[cid] = max(0.0, full_values[cid] + rng.uniform(-epsilon, epsilon))
        r_flips = 0
        r_comparisons = 0
        r_ties = 0
        for a, b in pairs:
            state = pair_diff_state(full_values[a] - full_values[b], jittered, a, b)
            if state == 'tie':
                per_pair[(a, b)]['ties'] += 1
                r_ties += 1
            elif state == 'skip':
                continue
            else:
                per_pair[(a, b)]['comparisons'] += 1
                r_comparisons += 1
                if state == 'flip':
                    per_pair[(a, b)]['flips'] += 1
                    r_flips += 1
        reps.append((r_flips, r_comparisons, r_ties))
    return per_pair, reps


def agg_value_axis(per_pair, epsilon, jitter_meta):
    flips = sum(p['flips'] for p in per_pair.values())
    comparisons = sum(p['comparisons'] for p in per_pair.values())
    ties = sum(p['ties'] for p in per_pair.values())
    out = {
        'comparisons': comparisons,
        'epsilon': epsilon,
        'flip_rate': (flips / comparisons) if comparisons else None,
        'flips': flips,
        'replicates': JITTER_REPLICATES,
        'ties': ties,
    }
    out.update(jitter_meta)
    return out


def diff_dict(label, mine, theirs, keys, mismatches):
    for key in keys:
        if mine.get(key) != theirs.get(key):
            mismatches.append(f'{label}.{key}: rederived={mine.get(key)!r} v2={theirs.get(key)!r}')


def window_query(budget_metric, event_ids, family):
    """(event_id, family) resolution with the v2 epsilon fallback discipline."""
    by_pair, by_family, all_vals = budget_metric
    values = []
    directions = []
    n_fallback_family = 0
    n_fallback_global = 0
    for event_id in event_ids:
        entry = by_pair.get((event_id, family))
        if entry is None:
            entry = by_family.get(family)
            n_fallback_family += 1
        if entry is None:
            entry = all_vals['global_median']
            n_fallback_global += 1
        values.append(entry['value'])
        if entry.get('direction') is not None:
            directions.append(entry['direction'])
    return {
        'median': quantile(sorted(values), 0.5) if values else None,
        'directions': sorted(set(directions)),
        'n_event_windows': len(values),
        'n_fallback_family_median': n_fallback_family,
        'n_fallback_global_median': n_fallback_global,
    }


def build_metric_index(per_windows, value_key, direction_key=None):
    """{(event_id, family): {'value':..., 'direction':...}} plus family/global medians."""
    by_pair = {}
    by_family = {}
    all_values = []
    for w in per_windows:
        value = w[value_key]
        node = {'value': value}
        if direction_key is not None:
            node['direction'] = w[direction_key]
        key = (w['event_id'], w['family'])
        by_pair[key] = node
        by_family.setdefault(w['family'], []).append(value)
        all_values.append(value)
    family_median = {f: {'value': quantile(sorted(v), 0.5)} for f, v in by_family.items()}
    global_value = quantile(sorted(all_values), 0.5)
    global_median = {'value': global_value, 'direction': None}
    return by_pair, family_median, {'global_median': global_median}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', required=True, help='must not exist: evidence is never overwritten')
    ap.add_argument('--run-tag', default='r1')
    args = ap.parse_args()
    out = Path(args.output_dir)
    if not out.is_absolute():
        out = ROOT / out
    assert not out.exists(), f'output dir must not exist: {out}'
    started = time.perf_counter()

    # ---------------- input gates (all nine inputs) ------------------------
    derivation_sha = assert_sha(DERIVATION, DERIVATION_SHA256, 'derivation program')
    mission_sha = assert_sha(MISSION_TOLERANCE, MISSION_TOLERANCE_SHA256, 'mission tolerance v0.2')
    rule_sha = assert_sha(STABILITY_RULE, STABILITY_RULE_SHA256, 'stability rule v0.1')
    event_sha = assert_sha(EVENT_TABLE, EVENT_TABLE_SHA256, 'event table')
    s4_sha = assert_sha(S4_SPLIT, S4_SPLIT_SHA256, 'S4 baseline split')
    cap_sha = assert_sha(CAPABILITY, CAPABILITY_SHA256, 'capability export v3')
    flip_sha = assert_sha(FLIPRATE, FLIPRATE_SHA256, 'fliprate v2')
    budget_sha = assert_sha(REF_BUDGET, REF_BUDGET_SHA256, 'reference budget')

    derivation = load_json(DERIVATION)
    assert derivation['schema_version'] == 'g4_threshold_derivation_v1.0', 'derivation schema mismatch'
    mission = load_json(MISSION_TOLERANCE)
    assert mission['status'] == 'frozen', 'mission tolerance config is not frozen'
    rule = load_json(STABILITY_RULE)
    assert rule['status'] == 'frozen', 'stability rule is not frozen'
    s4 = load_json(S4_SPLIT)
    assert s4['status'] == 'frozen', 'S4 baseline split is not frozen'
    capability = load_json(CAPABILITY)
    fliprate = load_json(FLIPRATE)
    budget = load_json(REF_BUDGET)

    assert derivation['inputs_sha256']['capability_v3_r1'] == cap_sha, 'capability hash not the frozen one'
    assert derivation['inputs_sha256']['fliprate_v2_r1'] == flip_sha, 'fliprate hash not the frozen one'
    assert derivation['damage_ladder_merged_stripped_sha256'] == LADDER_MERGED_STRIPPED_SHA256

    # S4 dev families must be the ones the ladder reports
    dev_families = sorted(g.split('-')[1].lower().rstrip('m') for g in s4['groups']['dev'])
    assert tuple(sorted(DEV_FAMILIES)) == tuple(dev_families), f'dev family mismatch: {dev_families}'

    levels_map = frozen_levels(mission['ladder_level_mapping'][MISSION_RELEVANT])
    assert set(levels_map) == set(DAMAGE_TYPES), 'mission-relevant damage types changed'

    # ---------------- ladder re-derivation ---------------------------------
    ladder = load_ladder_chunks()
    assert len(ladder) == N_RECORDS, f'unexpected ladder record count: {len(ladder)}'
    idents = [r for r in ladder if r['candidate_id'] == IDENTITY_ANCHOR and r['availability'] == 'ran']
    assert idents and all(r['metrics']['waveform']['metrics']['nrmse'] == 0.0 for r in idents), \
        'identity anchor D!=0 found'

    level_class = level_class_map(mission['ladder_level_mapping'])
    strata_records = build_strata(ladder, level_class)

    # ---------------- reference budget metric indices ----------------------
    deviation = build_metric_index(
        [{'event_id': w['event_id'], 'family': w['family'],
          'value': w['shape_class']['direct_evidence']['deviation_from_unity']}
         for w in budget['per_event_window']], 'value')
    arrival = build_metric_index(
        [{'event_id': w['event_id'], 'family': w['family'],
          'value': w['arrival_time_class']['direct_evidence']['grid_tier_shift_ns'],
          'direction': w['arrival_time_class']['direct_evidence']['direction']}
         for w in budget['per_event_window']], 'value', 'direction')

    # ---------------- per-pair re-derivation + ACCEPTANCE GATE -------------
    v2_index = {}
    for row in fliprate['strata']:
        key = (row['quantity'], row['damage_type'], row['mission_class'],
               row['geometry'], row['family'])
        assert key not in v2_index, f'duplicate fliprate v2 stratum {key}'
        v2_index[key] = row
    assert len(v2_index) == fliprate['n_strata'] == 24, 'fliprate v2 stratum count changed'

    mismatches = []
    per_pair_tables = {}
    strata_info = {}
    for key in sorted(strata_records):
        records_by_cid = strata_records[key]
        dt_, cls, geo, fam = key
        quantity = 'N_b' if dt_ == NC_TYPE else 'D'
        v2_key = (quantity, dt_, cls, geo, fam)
        assert v2_key in v2_index, f'stratum {v2_key} missing from fliprate v2'
        v2 = v2_index[v2_key]

        cids = sorted(records_by_cid)
        pairs = pair_list(cids)
        full_values = {cid: p80([v for _lvl, _e, v in records_by_cid[cid]]) for cid in cids}
        events = sorted({event_id for rows in records_by_cid.values()
                         for _lvl, event_id, _v in rows})

        assert cids == sorted(v2['full_values']), f'{v2_key}: candidate set differs from v2'
        for cid in cids:
            assert full_values[cid] == v2['full_values'][cid], f'{v2_key}/{cid}: full value differs'

        dev_window = window_query(deviation, events, fam)
        median_dev = dev_window['median']
        epsilon = math.sqrt(2.0 * median_dev)
        jitter_meta = {
            'median_deviation_from_unity': median_dev,
            'n_event_windows': dev_window['n_event_windows'],
            'n_fallback_family_median': dev_window['n_fallback_family_median'],
            'n_fallback_global_median': dev_window['n_fallback_global_median'],
        }

        # re-derivation A: the verbatim v2 functions (imported)
        levels_a, agg_a_level = level_axis(records_by_cid, pairs)
        agg_a_window = window_axis(records_by_cid, full_values, pairs)
        agg_a_value = value_axis(full_values, pairs, cids, repr(key), epsilon)
        agg_a_value.update(jitter_meta)

        # re-derivation B: per-pair tables (what the stability rule consumes)
        levels_b, tab_level = pair_level_axis(records_by_cid, pairs)
        events_b, tab_window, reps_window = pair_window_axis(records_by_cid, full_values, pairs)
        tab_value, _reps_value = pair_value_axis(full_values, pairs, cids, repr(key), epsilon)
        assert levels_a == levels_b == v2['levels'], f'{v2_key}: level grid differs from v2'
        assert events_b == events, f'{v2_key}: event set differs'

        agg_level = agg_level_axis(levels_b, tab_level)
        agg_window = agg_window_axis(len(events_b), tab_window, reps_window)
        agg_value = agg_value_axis(tab_value, epsilon, jitter_meta)

        diff_dict(f'{v2_key}/level_axis', agg_level, v2['level_axis'],
                  ('applicable', 'comparisons', 'flip_rate', 'flipped_pairs', 'n_levels', 'ties'), mismatches)
        diff_dict(f'{v2_key}/window_jackknife', agg_window, v2['window_jackknife'],
                  ('applicable', 'comparisons', 'flip_rate', 'flips', 'max_replicate_flip_rate',
                   'n_replicates', 'ties'), mismatches)
        diff_dict(f'{v2_key}/value_jitter', agg_value, v2['value_jitter'],
                  ('comparisons', 'epsilon', 'flip_rate', 'flips', 'replicates', 'ties',
                   'median_deviation_from_unity', 'n_event_windows',
                   'n_fallback_family_median', 'n_fallback_global_median'), mismatches)
        # the two independent re-derivations must also agree with each other
        diff_dict(f'{v2_key}/level_axis(cross)', agg_level, agg_a_level,
                  ('applicable', 'comparisons', 'flipped_pairs', 'n_levels', 'ties'), mismatches)
        diff_dict(f'{v2_key}/window_jackknife(cross)', agg_window, agg_a_window,
                  ('applicable', 'comparisons', 'flips', 'n_replicates', 'ties'), mismatches)
        diff_dict(f'{v2_key}/value_jitter(cross)', agg_value, agg_a_value,
                  ('comparisons', 'flips', 'ties', 'epsilon'), mismatches)

        arrival_window = window_query(arrival, events, fam)
        per_pair_tables[key] = {
            'level': tab_level,
            'window': tab_window,
            'value': tab_value,
        }
        strata_info[key] = {
            'quantity': quantity,
            'damage_type': dt_,
            'mission_class': cls,
            'geometry': geo,
            'family': fam,
            'cids': cids,
            'pairs': pairs,
            'full_values': full_values,
            'levels': levels_b,
            'events': events,
            'epsilon': v2['value_jitter']['epsilon'],
            'jitter_meta': jitter_meta,
            'arrival_window': arrival_window,
        }

    per_pair_reproduction_ok = not mismatches
    if mismatches:
        print('ACCEPTANCE GATE FAILED: per-pair re-derivation does not reproduce fliprate v2',
              file=sys.stderr)
        for line in mismatches[:50]:
            print('  ' + line, file=sys.stderr)
        raise SystemExit(2)

    # ---------------- capability / negative-control indices ----------------
    cap_rows = capability['capability']
    nc_rows = capability['negative_control']
    assert len(cap_rows) == capability['n_capability_rows'] == 2065, 'capability row count changed'
    assert len(nc_rows) == capability['n_negative_control_rows'] == 324, 'nc row count changed'
    cap_index, cap_dupes = index_capability(cap_rows)
    nc_index, nc_dupes = index_negative_control(nc_rows)
    assert not cap_dupes, f'duplicate capability rows: {cap_dupes[:5]}'
    assert not nc_dupes, f'duplicate negative-control rows: {nc_dupes[:5]}'

    candidates = sorted({r['candidate_id'] for r in cap_rows} | {r['candidate_id'] for r in nc_rows})
    combos = sorted({(r['geometry'], r['family']) for r in cap_rows
                     if r['family'] in DEV_FAMILIES})
    # every mission-relevant stratum of every combo must exist (needed for epsilon)
    for geo, fam in combos:
        for dt_ in DAMAGE_TYPES:
            key = (dt_, MISSION_RELEVANT, geo, fam)
            assert key in strata_info, f'mission-relevant stratum missing: {key}'
    # cross-scale caveat disclosure: located by prefix, not by list position
    cross_scale_note = next((d for d in derivation['disclosures']
                             if d.startswith('cross-scale caveat')), None)
    assert cross_scale_note is not None, 'cross-scale caveat disclosure not found in the frozen program'

    # ---------------- D and N constraints ----------------------------------
    combo_summary = []
    admissible_by_combo = {}
    for geo, fam in combos:
        epsilons = {dt_: strata_info[(dt_, MISSION_RELEVANT, geo, fam)]['epsilon']
                    for dt_ in DAMAGE_TYPES}
        admissible = []
        per_candidate = []
        for cid in candidates:
            d_blocks = []
            d_pass = True
            for dt_ in DAMAGE_TYPES:
                rows_ok = True
                any_exceed = False
                per_level = []
                n_missing = 0
                biggest = None
                for lvl in levels_map[dt_]:
                    key = (cid, dt_, lvl, geo, fam, MISSION_RELEVANT)
                    row = cap_index.get(key)
                    if row is None:
                        n_missing += 1
                        rows_ok = False
                        per_level.append({'level': lvl, 'available': False,
                                          'D_p80': None, 'within_epsilon': None,
                                          'n': None})
                        continue
                    dp80 = float(row['D_p80'])
                    within = dp80 <= epsilons[dt_]
                    any_exceed = any_exceed or (not within)
                    biggest = dp80 if biggest is None else max(biggest, dp80)
                    per_level.append({'level': lvl, 'available': True, 'D_p80': dp80,
                                      'within_epsilon': within, 'n': row['n']})
                reasons = []
                if not rows_ok:
                    reasons.append('missing_or_unavailable_capability_row')
                if any_exceed:
                    reasons.append('D_p80_above_mission_tolerance_epsilon')
                all_ok = not reasons
                block = {
                    'damage_type': dt_,
                    'epsilon': epsilons[dt_],
                    'levels_frozen': list(levels_map[dt_]),
                    'n_levels_missing': n_missing,
                    'max_D_p80_over_available_levels': biggest,
                    'per_level': per_level,
                    'all_levels_available': rows_ok,
                    'all_levels_within_epsilon': rows_ok and not any_exceed,
                    'pass': all_ok,
                    'verdict_reason': (None if all_ok else ' + '.join(reasons)),
                }
                d_blocks.append(block)
                d_pass = d_pass and block['pass']

            nc_blocks = []
            n_pass = True
            for q in NC_QS:
                row = nc_index.get((cid, geo, fam, q))
                if row is None:
                    nc_blocks.append({'nc_amplify_q': q, 'available': False,
                                      'N_b_energy_ratio_p80': None, 'within_n_th': None,
                                      'n': None})
                    n_pass = False
                    continue
                val = float(row['N_b_energy_ratio_p80'])
                within = val <= N_THRESHOLD
                n_pass = n_pass and within
                nc_blocks.append({'nc_amplify_q': q, 'available': True,
                                  'N_b_energy_ratio_p80': val, 'within_n_th': within,
                                  'n': row['n']})

            admissible_flag = bool(d_pass and n_pass)
            if admissible_flag:
                admissible.append(cid)
            per_candidate.append({
                'candidate_id': cid,
                'D_constraint': d_blocks,
                'D_pass': d_pass,
                'N_constraint': {'n_th': N_THRESHOLD,
                                 'joint_reading': 'CONSERVATIVE JOINT (both q rows must pass)',
                                 'per_q': nc_blocks,
                                 'pass': n_pass},
                'admissible': admissible_flag,
            })
        admissible_by_combo[(geo, fam)] = admissible
        n_d_pass = sum(1 for row in per_candidate if row['D_pass'])
        n_n_pass = sum(1 for row in per_candidate if row['N_constraint']['pass'])
        combo_summary.append({
            'geometry': geo,
            'family': fam,
            'n_candidates_total': len(candidates),
            'n_D_pass': n_d_pass,
            'n_N_pass': n_n_pass,
            'n_admissible': len(admissible),
            'admissible_candidates': list(admissible),
            'epsilon_by_damage_type': epsilons,
            'per_candidate': per_candidate,
        })

    # identity anchor acceptance check: it must clear the D constraint everywhere
    anchor_ok = True
    anchor_detail = []
    for geo, fam in combos:
        present = IDENTITY_ANCHOR in admissible_by_combo[(geo, fam)]
        zero = True
        for dt_ in DAMAGE_TYPES:
            for lvl in levels_map[dt_]:
                row = cap_index.get((IDENTITY_ANCHOR, dt_, lvl, geo, fam, MISSION_RELEVANT))
                if row is None or float(row['D_p80']) != 0.0:
                    zero = False
        anchor_ok = anchor_ok and present and zero
        anchor_detail.append({'geometry': geo, 'family': fam,
                              'admissible': present, 'D_p80_zero_on_all_mission_levels': zero})

    # ---------------- stability rule applied once (21 D strata) ------------
    strata_verdicts = []
    nb_caveat_rows = []
    for key in sorted(strata_info):
        info = strata_info[key]
        quantity = info['quantity']
        admissible_combo = sorted(admissible_by_combo.get((info['geometry'], info['family']), []))
        cids_set = set(info['cids'])
        present_cids = [cid for cid in admissible_combo if cid in cids_set]
        tables = per_pair_tables[key]
        pair_rows = []
        relation = {}
        for a, b in [(x, y) for x in present_cids for y in present_cids if x < y]:
            lv = tables['level'][(a, b)]
            wd = tables['window'][(a, b)]
            vl = tables['value'][(a, b)]
            base_diff = info['full_values'][a] - info['full_values'][b]
            better = None if base_diff == 0.0 else (a if base_diff < 0.0 else b)
            determined = (base_diff != 0.0 and lv['flips'] == 0 and wd['flips'] == 0
                          and vl['flips'] == 0)
            relation[(a, b)] = better if determined else None
            pair_rows.append({
                'a': a, 'b': b,
                'base_D_p80_a': info['full_values'][a],
                'base_D_p80_b': info['full_values'][b],
                'better': better,
                'level_axis_flips': lv['flips'],
                'level_axis_ties': lv['ties'],
                'window_jackknife_flips': wd['flips'],
                'window_jackknife_comparisons': wd['comparisons'],
                'window_jackknife_ties': wd['ties'],
                'value_jitter_flips': vl['flips'],
                'value_jitter_comparisons': vl['comparisons'],
                'value_jitter_ties': vl['ties'],
                'determined': determined,
            })
        tops = []
        for cid in present_cids:
            others = [o for o in present_cids if o != cid]
            if all(relation.get((min(cid, o), max(cid, o))) == cid for o in others):
                tops.append(cid)
        if not present_cids:
            verdict = 'no_admissible_candidate'
        elif tops:
            assert len(tops) == 1, f'{key}: more than one unique stable top {tops}'
            verdict = 'unique_stable_top'
        elif any(r['determined'] for r in pair_rows):
            verdict = 'partial_only'
        else:
            verdict = 'set'

        row = {
            'stratum_key': (info['quantity'], info['damage_type'], info['mission_class'],
                            info['geometry'], info['family']),
            'quantity': quantity,
            'damage_type': info['damage_type'],
            'mission_class': info['mission_class'],
            'geometry': info['geometry'],
            'family': info['family'],
            'n_candidates_in_stratum': len(info['cids']),
            'n_admissible_candidates': len(admissible_combo),
            'n_admissible_present_in_stratum': len(present_cids),
            'admissible_present': present_cids,
            'admissible_not_present_in_stratum': [c for c in admissible_combo if c not in cids_set],
            'n_determined_pairs': sum(1 for r in pair_rows if r['determined']),
            'n_pairs_admissible': len(pair_rows),
            'stability_verdict': verdict,
            'unique_stable_top_candidate': tops[0] if tops else None,
            'singleton_admissible_set_vacuously_unique': bool(tops) and len(present_cids) == 1,
            'stability_applied': True,
            'cross_scale_caveat': False,
            'pair_rows': pair_rows,
        }
        if quantity == 'N_b':
            # the 3 N_b strata get NO per-pair stability verdict: epsilon is
            # calibrated on the shape-class (NRMSE) scale and applied to a
            # quadratic energy-ratio quantity there.
            row['pair_rows'] = []
            row['stability_applied'] = False
            row['stability_verdict'] = 'not_applied_cross_scale_caveat'
            row['cross_scale_caveat'] = True
            row['note'] = cross_scale_note
            row['unique_stable_top_candidate'] = None
            row['singleton_admissible_set_vacuously_unique'] = False
            nb_caveat_rows.append(row)
        strata_verdicts.append(row)

    d_verdict_rows = [r for r in strata_verdicts if r['quantity'] == 'D']
    assert len(d_verdict_rows) == 21, f'expected 21 D strata, got {len(d_verdict_rows)}'
    assert len(nb_caveat_rows) == 3, f'expected 3 N_b strata, got {len(nb_caveat_rows)}'

    # ---------------- a80-style capability statements ----------------------
    capability_statements = []
    for geo, fam in combos:
        for dt_ in DAMAGE_TYPES:
            key = (dt_, MISSION_RELEVANT, geo, fam)
            info = strata_info[key]
            arr = info['arrival_window']
            per_candidate_rows = []
            for cid in candidates:
                per_level = []
                values = []
                missing = 0
                for lvl in levels_map[dt_]:
                    row = cap_index.get((cid, dt_, lvl, geo, fam, MISSION_RELEVANT))
                    if row is None:
                        missing += 1
                        per_level.append({'level': lvl, 'available': False, 'D_p80': None,
                                          'within_epsilon': None, 'n': None})
                        continue
                    dp80 = float(row['D_p80'])
                    values.append(dp80)
                    per_level.append({'level': lvl, 'available': True, 'D_p80': dp80,
                                      'within_epsilon': dp80 <= info['epsilon'], 'n': row['n']})
                per_candidate_rows.append({
                    'candidate_id': cid,
                    'mission_class': MISSION_RELEVANT,
                    'levels_included': list(levels_map[dt_]),
                    'n_levels_missing': missing,
                    'a80_upper_bound_D_p80': max(values) if values else None,
                    'a80_statement_complete': missing == 0 and bool(values),
                    'per_level': per_level,
                })
            capability_statements.append({
                'geometry': geo,
                'family': fam,
                'damage_type': dt_,
                'mission_class': MISSION_RELEVANT,
                'headline_excludes_weak_event_band': True,
                'weak_event_band_levels_excluded': mission['ladder_level_mapping']
                                                          [WEAK_EVENT_BAND].get(dt_, []),
                'mission_tolerance_epsilon': info['epsilon'],
                'statistical_term': {
                    'kind': 'a80 upper bound of erasure from the constructed reference',
                    'definition': ('D_p80 = p80 of per-record nrmse of the operator output vs '
                                   'identity(damaged input); D->1 means the damage is erased'),
                    'how_to_read': ('upper bound over the frozen mission-relevant levels; lower is '
                                    'smaller erasure; constructed reference, not physical truth'),
                },
                'systematic_bias_term': {
                    'kind': 'time-of-arrival systematic bias under velocity deviation',
                    'source': ('reference_uncertainty_budget_r1 per_event_window '
                               'arrival_time_class.direct_evidence.grid_tier_shift_ns'),
                    'quantity': 'envelope-peak BASE<->FINE2 inter-tier offset of the stratum event windows',
                    'grid_chain_direction_median_ns': arr['median'],
                    'grid_chain_direction_median_abs_ns': None if arr['median'] is None
                                                         else abs(arr['median']),
                    'unit': 'ns',
                    'sign_convention': 'negative = FINE2 envelope peak earlier (budget direction field)',
                    'directions_observed': arr['directions'],
                    'n_event_windows': arr['n_event_windows'],
                    'n_fallback_family_median': arr['n_fallback_family_median'],
                    'n_fallback_global_median': arr['n_fallback_global_median'],
                    'window_selection': ('same (event_id, family) resolution and family/global median '
                                         'fallback discipline as the stratum epsilon'),
                    'event_ids': list(info['events']),
                    'reported_separately': True,
                    'depth_conversion_performed': False,
                    'mixed_into_statistical_term': False,
                },
                'per_candidate': per_candidate_rows,
            })

    # ---------------- result ----------------------------------------------
    result = {
        'schema': 'g4-threshold-derivation/1',
        'plan': PLAN,
        'plan_step': STEP,
        'derivation_program': str(DERIVATION.relative_to(ROOT)),
        'derivation_program_sha256': derivation_sha,
        'derivation_rules': derivation['derivation_rules'],
        'hard_limits': derivation['hard_limits'],
        'disclosures': derivation['disclosures'],
        'g4_relieved': False,
        'thresholds_are_simulation_domain_only': True,
        'constructed_reference_declaration': (
            'D is measured against identity(damaged input): a constructed reference, exact by '
            'construction, NOT physical clean truth; every statement here is simulation-domain and '
            'must not be extrapolated to field performance'),
        'no_ranking_no_selection_no_training': ('structure only: admissibility and pairwise '
                                                'determination are reported, no configuration is '
                                                'selected for deployment and no labels are generated'),
        'frozen_inputs_sha256': {
            str(DERIVATION.relative_to(ROOT)): derivation_sha,
            str(MISSION_TOLERANCE.relative_to(ROOT)): mission_sha,
            str(STABILITY_RULE.relative_to(ROOT)): rule_sha,
            str(EVENT_TABLE.relative_to(ROOT)): event_sha,
            str(S4_SPLIT.relative_to(ROOT)): s4_sha,
            str(CAPABILITY.relative_to(ROOT)): cap_sha,
            str(FLIPRATE.relative_to(ROOT)): flip_sha,
            str(REF_BUDGET.relative_to(ROOT)): budget_sha,
            'artifacts/research_checks/2026-09-26_damage_ladder_r1_c00..07/records.json'
            ' (resource-stripped merge)': LADDER_MERGED_STRIPPED_SHA256,
        },
        'methodology': {
            'scope': 'development group only (S4 frozen: {B2D-C1m, B2D-C3m}); MT and CO layered '
                     'separately, never merged or cross-paired; families reported separately',
            'candidate_set_source': 'all configurations present in capability export v3 '
                                    '(union of capability and negative_control rows); '
                                    'unavailable configurations recorded as unavailable, never scored zero',
            'D_constraint': derivation['derivation_rules']['D_constraint'],
            'N_constraint': derivation['derivation_rules']['N_constraint'],
            'epsilon_source': 'fliprate v2 stratum value_jitter.epsilon of the same '
                              '(damage_type, mission_relevant, geometry, family) stratum',
            'frozen_mission_levels': {dt_: list(map(canonical_level, lv))
                                      for dt_, lv in mission['ladder_level_mapping'][MISSION_RELEVANT].items()},
            'n_constraint_joint_reading': 'conservative joint: BOTH q=2 and q=4 rows must have '
                                          'N_b_energy_ratio_p80 <= 1.0 (max over q <= N_th)',
            'n_threshold': N_THRESHOLD,
            'per_pair_rederivation': 'from the frozen v2 methodology (verbatim functions of '
                                     'scripts/run_g4_sensitivity_fliprate.py plus per-pair tables); '
                                     'acceptance gate below compares per stratum per axis',
            'determined_pair_definition': rule['determined_order_rule']['statement_zh'],
            'determined_pair_operationalization': ('a pair counts as determined only if the base '
                                                   'ordering itself is strict (base difference != 0) '
                                                   'AND zero flips are observed on all three axes; '
                                                   'ties are excluded and counted separately'),
            'stratum_verdict_operationalization': ('no admissible candidate present -> '
                                                   'no_admissible_candidate; exactly one candidate '
                                                   'determined-better than every other admissible -> '
                                                   'unique_stable_top; some determined pairs but no '
                                                   'unique top -> partial_only; no determined pair at '
                                                   'all -> set (protocol v0.2 section 7 semantics)'),
            'stability_rule_applied_once': True,
            'n_b_strata_per_pair_verdict': 'NOT computed (cross-scale caveat reported instead)',
            'tie_policy': 'exact float equality = tie; ties excluded from flip denominators and '
                          'counted separately',
            'weak_event_band': 'never enters the capability headline (protection clause)',
        },
        'checks': {
            'per_pair_rederivation_reproduces_fliprate_v2_aggregates': per_pair_reproduction_ok,
            'n_strata_verified': len(strata_info),
            'n_d_strata': len(d_verdict_rows),
            'n_n_b_strata': len(nb_caveat_rows),
            'identity_anchor_admissible_and_D_zero': anchor_ok,
            'identity_anchor_detail': anchor_detail,
            'identity_anchor_candidate_id': IDENTITY_ANCHOR,
        },
        'candidate_set': {
            'n_candidates': len(candidates),
            'candidate_ids': candidates,
            'source': str(CAPABILITY.relative_to(ROOT)),
            'n_candidates_in_ladder': len({r['candidate_id'] for r in ladder}),
            'candidates_in_ladder_not_in_capability': sorted(
                {r['candidate_id'] for r in ladder} - set(candidates)),
        },
        'combo_summary': combo_summary,
        'strata_stability': strata_verdicts,
        'n_b_strata_cross_scale_caveat': nb_caveat_rows,
        'capability_statements': capability_statements,
    }
    assert anchor_ok, 'identity anchor failed the D constraint in at least one stratum'

    out.mkdir(parents=True)
    text = json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False)
    (out / 'results.json').write_bytes(text.encode('utf-8'))

    manifest = {
        'run_tag': args.run_tag,
        'runner': str(Path(__file__).relative_to(ROOT)),
        'stop_reason': None,
        'inputs_gated': 9,
        'n_candidates': len(candidates),
        'n_combos': len(combos),
        'per_combo_admissible': {f'{g}/{f}': len(admissible_by_combo[(g, f)])
                                 for g, f in combos},
        'n_d_strata': len(d_verdict_rows),
        'n_n_b_strata': len(nb_caveat_rows),
        'stability_verdict_histogram': {v: sum(1 for r in strata_verdicts
                                               if r['stability_verdict'] == v)
                                        for v in sorted({r['stability_verdict']
                                                         for r in strata_verdicts})},
        'checks': {
            'per_pair_rederivation_reproduces_fliprate_v2_aggregates': per_pair_reproduction_ok,
            'identity_anchor_admissible_and_D_zero': anchor_ok,
        },
        'hard_limits_source': str(DERIVATION.relative_to(ROOT)),
        'g4_relieved': False,
        'solver_invoked': False,
        'training_eligible': False,
        'total_wall_s': time.perf_counter() - started,
    }
    (out / 'run_manifest.json').write_bytes(
        (json.dumps(manifest, indent=1, sort_keys=True, ensure_ascii=False) + '\n').encode('utf-8'))
    print('wrote', out / 'results.json')
    print('sha256:', sha256_file(out / 'results.json'))
    print('admissible per combo:', manifest['per_combo_admissible'])
    print('stability verdicts:', manifest['stability_verdict_histogram'])


if __name__ == '__main__':
    main()
