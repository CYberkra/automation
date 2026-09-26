"""G4 step 4: preregistered sensitivity check of pairwise candidate ordering.

Deterministic, CPU read-only. Gates on the frozen mission tolerance config, the
frozen event table and the archived reference-uncertainty budget, then merges
the r1 damage-ladder records (resource-stripped, identical across r1/r2 per
acceptance evidence) and reports, per stratum, how often the pairwise order of
two candidates flips under three perturbation axes:

  level_axis       order recomputed at each frozen ladder level,
  window_jackknife leave-one-event-out replicates as the executable proxy of a
                   window-boundary +/-delta perturbation,
  value_jitter     deterministic uniform jitter of magnitude
                   epsilon = sqrt(2 * median deviation_from_unity) taken from
                   the reference-uncertainty budget shape class.

Stratification: MT/CO geometries are layered separately and never merged or
cross-paired; families are reported separately on dev families {C1, C3} only.

Flip rates are diagnostic distributions only. No ranking, no selection, no
pass/fail thresholding, no capability labels: this step fixes no "frequent
flip" cutoff (reserved for step 5) and does not relieve G4.
"""
import argparse
import hashlib
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.1.json'
MISSION_SHA256 = '2d630c01cc95f3c13febbc19c4bf8f3b5a53b5fb82bbabcb2f68e7537ef78b28'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
REF_BUDGET = ROOT / 'artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json'
REF_BUDGET_SHA256 = '38cf98bab723358dc8d75486ef6824057d42c32a9030b56f7101400fc87a8cc1'
CHUNK_GLOB = 'artifacts/research_checks/2026-09-26_damage_ladder_r1_c{:02d}'
N_CHUNKS = 8
N_RECORDS = 34086
DEV_FAMILIES = ('c1', 'c3')
Q = 0.8
JITTER_REPLICATES = 200
MISSION_CLASSES = ('mission_relevant', 'weak_event_band')
NC_TYPE = 'nc_amplify'
NC_CLASS = 'negative_control'
IDENTITY_ANCHOR = 'B0_G1_BG'


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_bytes_concat(paths):
    h = hashlib.sha256()
    for p in paths:
        h.update(Path(p).read_bytes())
    return h.hexdigest()


def load_chunks():
    recs = []
    for k in range(N_CHUNKS):
        d = json.loads((ROOT / CHUNK_GLOB.format(k) / 'records.json').read_text(encoding='utf-8'))
        for r in d['records']:
            r.pop('resource', None)
        recs.extend(d['records'])
    return recs


def quantile(sorted_vals, q):
    """Deterministic linear-interpolation quantile (numpy default method)."""
    if not sorted_vals:
        return None
    n = len(sorted_vals)
    pos = (n - 1) * q
    lo = int(pos)
    hi = min(lo + 1, n - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def p80(vals):
    return quantile(sorted(vals), Q)


def ordering(values):
    """Ascending value = better (lower D / lower residual ratio); cid tie-break."""
    return sorted(values, key=lambda c: (values[c], c))


def level_class_map(mapping):
    out = {}
    for cls in MISSION_CLASSES:
        block = mapping.get(cls, {})
        for damage_type, levels in block.items():
            if not isinstance(levels, list):
                continue
            for lvl in levels:
                out[(damage_type, lvl)] = cls
    return out


def build_strata(recs, level_class):
    """Return {stratum_key: {(level, event_id, value) per candidate}}.

    stratum_key = (damage_type, mission_class, geometry, family); the level is
    kept per record so the level axis can recompute the order per level.
    """
    strata = {}
    for r in recs:
        if r['family'] not in DEV_FAMILIES or r['availability'] != 'ran':
            continue
        cid = r['candidate_id']
        geo = r['geometry']
        fam = r['family']
        if r['row_type'] == 'event':
            dt_ = r['damage']['type']
            lvl = r['damage']['level']
            cls = level_class.get((dt_, lvl))
            if cls is None:
                # polarity_flip is frozen as level 'all' -> every event row.
                cls = level_class.get((dt_, 'all'))
            if cls is None:
                continue
            d = r['metrics']['waveform']['metrics']['nrmse']
            key = (dt_, cls, geo, fam)
        elif r['row_type'] == 'negative_control' and r['damage']['type'] == NC_TYPE:
            lvl = r['damage']['level']
            d = r['metrics']['N_b']['energy_ratio']
            key = (NC_TYPE, NC_CLASS, geo, fam)
        else:
            continue
        strata.setdefault(key, {}).setdefault(cid, []).append((lvl, r['event_id'], d))
    return strata


def pair_list(cids):
    return [(cids[i], cids[j]) for i in range(len(cids)) for j in range(i + 1, len(cids))]


def sign_compare(base_values, other_values, pairs):
    """Count order flips of `other_values` against `base_values` over pairs.

    Exact float equality on the difference is a tie: excluded from the
    comparison denominator and counted separately.
    """
    flips = 0
    comparisons = 0
    ties = 0
    for a, b in pairs:
        base_diff = base_values[a] - base_values[b]
        if base_diff == 0.0:
            ties += 1
            continue
        oa = other_values.get(a)
        ob = other_values.get(b)
        if oa is None or ob is None:
            continue
        other_diff = oa - ob
        if other_diff == 0.0:
            ties += 1
            continue
        comparisons += 1
        if (base_diff > 0) != (other_diff > 0):
            flips += 1
    return flips, comparisons, ties


def level_axis(records_by_cid, pairs):
    """Flip rate across the frozen ladder levels present in the stratum."""
    by_level = {}
    for cid, rows in records_by_cid.items():
        per_level = {}
        for lvl, _event_id, val in rows:
            per_level.setdefault(lvl, []).append(val)
        for lvl, vals in per_level.items():
            by_level.setdefault(lvl, {})[cid] = p80(vals)
    levels = sorted(by_level, key=repr)
    empty = {
        'applicable': False,
        'n_levels': len(levels),
        'flip_rate': None,
        'flipped_pairs': 0,
        'comparisons': 0,
        'ties': 0,
    }
    if len(levels) < 2:
        return levels, empty
    flipped_pairs = 0
    comparable_pairs = 0
    ties = 0
    for a, b in pairs:
        signs = []
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
        if len(signs) >= 2:
            comparable_pairs += 1
            if len(set(signs)) > 1:
                flipped_pairs += 1
    if comparable_pairs == 0:
        return levels, empty
    return levels, {
        'applicable': True,
        'n_levels': len(levels),
        'flip_rate': flipped_pairs / comparable_pairs,
        'flipped_pairs': flipped_pairs,
        'comparisons': comparable_pairs,
        'ties': ties,
    }


def window_axis(records_by_cid, full_values, pairs):
    """Leave-one-event-out jackknife replicates vs the full-stratum ordering."""
    events = sorted({event_id for rows in records_by_cid.values() for _lvl, event_id, _v in rows})
    empty = {
        'applicable': False,
        'n_replicates': len(events),
        'flip_rate': None,
        'max_replicate_flip_rate': None,
        'flips': 0,
        'comparisons': 0,
        'ties': 0,
    }
    if len(events) < 2:
        return empty
    flips = 0
    comparisons = 0
    ties = 0
    max_rate = None
    for e in events:
        jack = {}
        for cid, rows in records_by_cid.items():
            vals = [v for _lvl, event_id, v in rows if event_id != e]
            if vals:
                jack[cid] = p80(vals)
        f, c, t = sign_compare(full_values, jack, pairs)
        flips += f
        comparisons += c
        ties += t
        if c > 0:
            rate = f / c
            if max_rate is None or rate > max_rate:
                max_rate = rate
    if comparisons == 0:
        return empty
    return {
        'applicable': True,
        'n_replicates': len(events),
        'flip_rate': flips / comparisons,
        'max_replicate_flip_rate': max_rate,
        'flips': flips,
        'comparisons': comparisons,
        'ties': ties,
    }


def value_axis(full_values, pairs, cids, stratum_key, epsilon):
    """Deterministic uniform jitter of the ranking quantity, K replicates."""
    flips = 0
    comparisons = 0
    ties = 0
    for k in range(JITTER_REPLICATES):
        jittered = {}
        for cid in cids:
            rng = random.Random(f'{stratum_key}|{k}|{cid}')
            jittered[cid] = max(0.0, full_values[cid] + rng.uniform(-epsilon, epsilon))
        f, c, t = sign_compare(full_values, jittered, pairs)
        flips += f
        comparisons += c
        ties += t
    return {
        'replicates': JITTER_REPLICATES,
        'epsilon': epsilon,
        'flip_rate': (flips / comparisons) if comparisons else None,
        'flips': flips,
        'comparisons': comparisons,
        'ties': ties,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-dir', required=True)
    args = ap.parse_args()
    out = Path(args.output_dir)
    assert not out.exists(), f'output dir must not exist: {out}'

    assert sha256_file(MISSION) == MISSION_SHA256, 'mission tolerance gate mismatch'
    assert sha256_file(EVENT_TABLE) == EVENT_TABLE_SHA256, 'event table gate mismatch'
    assert sha256_file(REF_BUDGET) == REF_BUDGET_SHA256, 'reference budget gate mismatch'

    mission = json.loads(MISSION.read_text(encoding='utf-8'))
    assert mission['status'] == 'frozen', 'mission tolerance config is not frozen'
    level_class = level_class_map(mission['ladder_level_mapping'])

    recs = load_chunks()
    assert len(recs) == N_RECORDS, f'unexpected ladder record count: {len(recs)}'

    # identity anchor self-consistency (constructed reference)
    idents = [r for r in recs if r['candidate_id'] == IDENTITY_ANCHOR and r['availability'] == 'ran']
    anchor_ok = all(r['metrics']['waveform']['metrics']['nrmse'] == 0.0 for r in idents)
    assert anchor_ok, 'identity anchor D!=0 found'

    budget = json.loads(REF_BUDGET.read_text(encoding='utf-8'))
    windows = [(w['event_id'], w['family'],
                w['shape_class']['direct_evidence']['deviation_from_unity'])
               for w in budget['per_event_window']]
    win_by_pair = {(e, f): d for e, f, d in windows}
    win_by_family = {}
    all_dev = []
    for e, f, d in windows:
        win_by_family.setdefault(f, []).append(d)
        all_dev.append(d)
    family_median = {f: quantile(sorted(v), 0.5) for f, v in win_by_family.items()}
    global_median = quantile(sorted(all_dev), 0.5)

    strata = build_strata(recs, level_class)
    rows = []
    for key in sorted(strata):
        records_by_cid = strata[key]
        dt_, cls, geo, fam = key
        cids = sorted(records_by_cid)
        full_values = {cid: p80([v for _lvl, _e, v in records_by_cid[cid]]) for cid in cids}
        pairs = pair_list(cids)
        events = sorted({event_id for rows_ in records_by_cid.values()
                         for _lvl, event_id, _v in rows_})

        devs = []
        n_fallback_family = 0
        n_fallback_global = 0
        for event_id in events:
            dev = win_by_pair.get((event_id, fam))
            if dev is None:
                dev = family_median.get(fam)
                n_fallback_family += 1
            if dev is None:
                dev = global_median
                n_fallback_global += 1
            devs.append(dev)
        median_dev = quantile(sorted(devs), 0.5)
        epsilon = math.sqrt(2.0 * median_dev)

        levels, lvl = level_axis(records_by_cid, pairs)
        rows.append({
            'quantity': 'N_b' if dt_ == NC_TYPE else 'D',
            'damage_type': dt_,
            'mission_class': cls,
            'geometry': geo,
            'family': fam,
            'n_candidates': len(cids),
            'n_pairs': len(pairs),
            'levels': levels,
            'n_events': len(events),
            'full_values': {cid: full_values[cid] for cid in cids},
            'full_ordering': ordering(full_values),
            'level_axis': lvl,
            'window_jackknife': window_axis(records_by_cid, full_values, pairs),
            'value_jitter': value_axis(full_values, pairs, cids, repr(key), epsilon),
        })
        rows[-1]['value_jitter']['median_deviation_from_unity'] = median_dev
        rows[-1]['value_jitter']['n_event_windows'] = len(devs)
        rows[-1]['value_jitter']['n_fallback_family_median'] = n_fallback_family
        rows[-1]['value_jitter']['n_fallback_global_median'] = n_fallback_global

    chunk_paths = [ROOT / CHUNK_GLOB.format(k) / 'records.json' for k in range(N_CHUNKS)]
    result = {
        'schema': 'g4-sensitivity-fliprate/1',
        'plan': 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md',
        'plan_section': 'G4 step 4 (section 5 preregistration: pairwise order flip rate under perturbation)',
        'mission_tolerance_config': str(MISSION.relative_to(ROOT)),
        'mission_tolerance_sha256': MISSION_SHA256,
        'event_table_sha256': EVENT_TABLE_SHA256,
        'reference_budget': str(REF_BUDGET.relative_to(ROOT)),
        'reference_budget_sha256': REF_BUDGET_SHA256,
        'ladder_source': CHUNK_GLOB.replace('{:02d}', '00..07') + ' (r1; r2 byte-identical after resource strip per acceptance)',
        'ladder_records_merged': len(recs),
        'dev_families': list(DEV_FAMILIES),
        'layering': 'MT and CO layered separately, never merged or cross-paired; families reported separately',
        'methodology': {
            'ranking_quantity': 'D_p80 = p80 of per-record nrmse vs identity(damaged input); lower = damage preserved; constructed reference, not physical truth',
            'nb_quantity': 'N_b energy_ratio p80 on nc_amplify rows; diagnostic residual ratio, not SNR/threshold',
            'level_axis': 'order recomputed at each frozen ladder level (operationalization of frozen tolerance grid); flip = sign change of pairwise difference between any two levels',
            'window_axis': 'leave-one-event-out jackknife as executable proxy of window-boundary +/-delta perturbation',
            'value_axis': 'deterministic uniform jitter of magnitude epsilon = sqrt(2*median deviation_from_unity over the stratum\'s event windows, from reference budget shape class); conversion assumes zero-mean equal-energy signals',
            'jitter_replicates': JITTER_REPLICATES,
            'tie_policy': 'exact float equality = tie; ties excluded from comparison denominators and counted separately',
            'preregistration_mapping': 'plan section 5: tau perturbation -> level_axis + value_axis; window +/-delta -> window_axis; no frequent-flip cutoff is fixed in this step (reserved for step 5)',
            'discipline': 'flip rates are diagnostic distributions only; no ranking, no selection, no pass/fail thresholding, no labels',
            'stratum_key': '(damage_type, mission_class, geometry, family); level is not part of the key, it is the level_axis grid; N_b strata carry damage_type=nc_amplify and mission_class=negative_control',
            'levels_ordering': 'level labels are the frozen-grid levels present in the stratum, ordered by repr for determinism; polarity_flip has one level ("all" -> 1.0) so its level axis is not applicable',
            'value_axis_reference': 'each jitter replicate is compared with the full-stratum ordering, identically to the jackknife replicates',
            'epsilon_fallback': 'deviation_from_unity looked up per (event_id, family) of the stratum; unmatched windows fall back to the family median of all budget windows, then to the global median (counts reported per stratum under value_jitter)',
            'value_axis_clamp': 'jittered D and N_b values are clamped at 0 (both quantities are non-negative)',
            'jitter_seed': 'random.Random("<repr(stratum_key)>|<replicate>|<candidate_id>").uniform(-epsilon, +epsilon); no other randomness, no wall clock',
        },
        'strata': rows,
        'n_strata': len(rows),
        'hard_limits': {
            'absolute_accuracy_claimed': False,
            'clean_truth_generated': False,
            'constructed_reference_not_physical_truth': True,
            'grid_convergence_certified': False,
            'no_ranking_no_selection_no_thresholding': True,
            'physical_acceptance_threshold': None,
            'reference_state': 'numerically_unresolved',
            'simulation_domain_only': True,
            'solver_invoked': False,
            'training_eligible': False,
            'training_labels_generated': False,
            'G4': 'NOT relieved by this analysis; candidates stay undetermined',
        },
        'provenance': {
            str(MISSION.relative_to(ROOT)): MISSION_SHA256,
            str(EVENT_TABLE.relative_to(ROOT)): EVENT_TABLE_SHA256,
            str(REF_BUDGET.relative_to(ROOT)): REF_BUDGET_SHA256,
            'artifacts/research_checks/2026-09-26_damage_ladder_r1_c00..07/records.json': sha256_bytes_concat(chunk_paths),
        },
    }

    out.mkdir(parents=True)
    text = json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False)
    (out / 'results.json').write_bytes(text.encode('utf-8'))
    print('wrote', out / 'results.json')
    print('sha256:', sha256_file(out / 'results.json'))
    print('strata:', len(rows))


if __name__ == '__main__':
    main()
