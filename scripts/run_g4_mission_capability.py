"""G4 step 3: export mission capability values from the damage ladder.

Deterministic, CPU read-only. Gates on the frozen mission tolerance config
(freeze_g4_mission_tolerance.py output), the frozen event table, and the 16
archived damage-ladder chunk directories. Merges the r1 ladder records
(resource-stripped, identical across r1/r2 per acceptance evidence), splits
damage instances into mission-relevant vs weak-event band per the frozen
mapping, and exports per-operator D distributions (median + a80-style p80
upper bound), stratified by geometry (MT/CO layered separately) and family.

No ranking, no selection, no pass/fail thresholding: undetermined
determination is G4 step 4 business. Run twice (r1/r2) and compare bytes.

Split mode (--split {dev,test}, default 'dev'): dev re-exports the frozen
capability values from the dev-family ladder {c1, c3} (unchanged, byte-for-byte
identical to the v3_r1 export); test exports the same aggregation over the
test families {c5, c8} from the test-family ladder. Both are pure data export:
no threshold is applied here, that happens once in S5.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.2.json'
MISSION_SHA256 = 'ee039fb1fab3f1147e8a5fe53809bbb8ae9d6f51aafc2feee0fb90071a6c57ea'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
CHUNK_GLOB = 'artifacts/research_checks/2026-09-26_damage_ladder_r1_c{:02d}'
TEST_CHUNK_GLOB = 'artifacts/research_checks/2026-09-27_damage_ladder_test_r1_c{:02d}'
N_CHUNKS = 8
DEV_FAMILIES = ('c1', 'c3')
TEST_FAMILIES = ('c5', 'c8')
DEV_N_RECORDS = 34086


def build_mission_levels(mission_json):
    """Derive {(damage_type, level): band} from ladder_level_mapping.

    Only list-valued keys inside each band are damage types; list elements are
    levels kept with the type JSON gave them. polarity_flip has no numeric grid,
    its 'all' entry maps to the float 1.0 used by the ladder records.
    """
    levels = {}
    for band in ('mission_relevant', 'weak_event_band'):
        block = mission_json['ladder_level_mapping'][band]
        for damage_type, value in block.items():
            if not isinstance(value, list):
                continue
            for level in value:
                if level == 'all':
                    level = 1.0
                levels[(damage_type, level)] = band
    return levels


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_chunks(glob=CHUNK_GLOB):
    recs = []
    for k in range(N_CHUNKS):
        d = json.loads((ROOT / glob.format(k) / 'records.json').read_text(encoding='utf-8'))
        for r in d['records']:
            r.pop('resource', None)
        recs.extend(d['records'])
    return recs


def expected_record_layout(families):
    """Build the required record set from frozen metadata; no gathers or solve."""
    from run_damage_ladder import (MT_ANALYSIS, expand_dev_rows, damage_instances,
                                  catalogue, L_LAMBDAS, L_WIDTHS, L_GAINS, R_THRESHOLDS)
    table = json.loads(EVENT_TABLE.read_text(encoding='utf-8'))
    mt = json.loads(MT_ANALYSIS.read_text(encoding='utf-8'))
    mothers = {c['mother']: c['run_id'] for c in mt['cases']}
    candidates = [c['id'] for c in catalogue()]
    candidates += [f'L_lam{lam}_w{w}_q{q}' for lam in L_LAMBDAS
                   for w in L_WIDTHS for q in L_GAINS]
    candidates += [f'R_tau{tau}_q1' for tau in R_THRESHOLDS]
    expected = {}
    for e, geom, case in expand_dev_rows(table, mothers, families):
        role = 'negative_control' if e['role'] in ('nc_zero', 'bg_absent', 'off_path') else 'event'
        for dtype, level in damage_instances(role):
            prefix = f"{e['event_id']}::{geom}:{case}::{dtype}:{level}"
            for cid in candidates:
                key = f'{prefix}::{cid}'
                if key in expected:
                    raise ValueError(f'duplicate expected record: {key}')
                expected[key] = (e['family'], geom, case, cid, role, dtype, level)
    if not expected or {v[0] for v in expected.values()} != set(families):
        raise ValueError('frozen metadata does not cover requested families')
    return expected


def validate_records(recs, families):
    """Unavailable candidates stay present; missing/duplicate rows never vanish."""
    expected = expected_record_layout(families)
    seen = set()
    for r in recs:
        key = r['record_id']
        actual = (r['family'], r['geometry'], r['case'], r['candidate_id'],
                  r['row_type'], r['damage']['type'], r['damage']['level'])
        if key in seen or expected.get(key) != actual:
            raise ValueError(f'duplicate, unexpected or inconsistent record: {key}')
        seen.add(key)
        if r['availability'] not in ('ran', 'unavailable'):
            raise ValueError(f'invalid availability: {key}')
        if r['candidate_id'] == 'B0_G1_BG':
            if (r['availability'] != 'ran' or
                    r['metrics']['waveform']['metrics']['nrmse'] != 0.0):
                raise ValueError(f'identity anchor unavailable or D!=0: {key}')
    if seen != expected.keys():
        raise ValueError(f'incomplete ladder: missing {len(expected.keys() - seen)} records')


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


def describe(vals):
    s = sorted(vals)
    return {
        'n': len(s),
        'D_median': quantile(s, 0.5),
        'D_p80': quantile(s, 0.8),
        'D_min': s[0],
        'D_max': s[-1],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--split', choices=['dev', 'test'], default='dev',
                    help="dev={c1,c3} frozen export (unchanged); "
                         "test={c5,c8} from the test-family ladder")
    args = ap.parse_args()
    split = args.split
    out = Path(args.output_dir)
    assert not out.exists(), f'output dir must not exist: {out}'

    if split == 'dev':
        chunk_glob = CHUNK_GLOB
        source_note = ' (r1; r2 byte-identical after resource strip per acceptance)'
        families = DEV_FAMILIES
    else:
        chunk_glob = TEST_CHUNK_GLOB
        source_note = (' (r1; r2 byte-identity for the test split is to be '
                       'demonstrated by acceptance, not assumed here)')
        families = TEST_FAMILIES

    assert sha256_file(MISSION) == MISSION_SHA256, 'mission tolerance gate mismatch'
    assert sha256_file(EVENT_TABLE) == EVENT_TABLE_SHA256, 'event table gate mismatch'
    recs = load_chunks(chunk_glob)
    if split == 'dev':
        assert len(recs) == DEV_N_RECORDS, f'unexpected ladder record count: {len(recs)}'
    validate_records(recs, families)

    mission = json.loads(MISSION.read_text(encoding='utf-8'))
    mission_levels = build_mission_levels(mission)

    cap = {}
    nc = {}
    for r in recs:
        if r['family'] not in families or r['availability'] != 'ran':
            continue
        cid = r['candidate_id']
        geo = r['geometry']
        fam = r['family']
        base = r['baseline']
        if r['row_type'] == 'event':
            dt_ = r['damage']['type']
            lvl = r['damage']['level']
            cls = mission_levels.get((dt_, lvl))
            if cls is None:
                continue
            d = r['metrics']['waveform']['metrics']['nrmse']
            key = (cid, base, dt_, repr(lvl), cls, geo, fam)
            cap.setdefault(key, []).append(d)
        elif r['row_type'] == 'negative_control' and r['damage']['type'] == 'nc_amplify':
            q = r['damage']['level']
            nb = r['metrics']['N_b']['energy_ratio']
            key = (cid, base, repr(q), geo, fam)
            nc.setdefault(key, []).append(nb)

    capability = []
    for (cid, base, dt_, lvl, cls, geo, fam), vals in sorted(cap.items()):
        row = {
            'candidate_id': cid,
            'baseline': base,
            'damage_type': dt_,
            'damage_level': json.loads(lvl),
            'mission_class': cls,
            'geometry': geo,
            'family': fam,
            'D_definition': 'nrmse of operator output vs identity(damaged input); constructed reference, not physical truth',
        }
        row.update(describe(vals))
        capability.append(row)

    negative_control = []
    for (cid, base, q, geo, fam), vals in sorted(nc.items()):
        s = sorted(vals)
        negative_control.append({
            'candidate_id': cid,
            'baseline': base,
            'nc_amplify_q': json.loads(q),
            'geometry': geo,
            'family': fam,
            'N_b_energy_ratio_median': quantile(s, 0.5),
            'N_b_energy_ratio_p80': quantile(s, 0.8),
            'n': len(s),
            'note': 'N_b is a diagnostic residual ratio on amplified negative controls, not SNR/threshold',
        })

    result = {
        'schema': 'g4-mission-capability/1',
        'plan': 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md',
        'mission_tolerance_config': str(MISSION.relative_to(ROOT)),
        'mission_tolerance_sha256': MISSION_SHA256,
        'event_table_sha256': EVENT_TABLE_SHA256,
        'ladder_source': chunk_glob.replace('{:02d}', '00..07') + source_note,
        'ladder_records_merged': len(recs),
        'identity_anchor': 'B0_G1_BG D==0 on all ran instances',
        'a80_note': 'D_p80 is the a80-style conservative upper bound of erasure (D->1 = damage erased); capability statements must add the explicit systematic bias term per mission_tolerance.conservativeness',
        'dev_families': list(DEV_FAMILIES),
        'layering': 'MT and CO layered separately, never merged or cross-paired; families reported separately',
        'n_capability_rows': len(capability),
        'n_negative_control_rows': len(negative_control),
        'capability': capability,
        'negative_control': negative_control,
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
            'G4': 'NOT relieved by this export; candidates stay undetermined',
        },
        'provenance': {
            str(MISSION.relative_to(ROOT)): MISSION_SHA256,
            str(EVENT_TABLE.relative_to(ROOT)): EVENT_TABLE_SHA256,
        },
    }
    if split == 'test':
        result['split'] = 'test'
        result['families'] = list(families)
        del result['dev_families']
    assert mission['status'] == 'frozen'

    out.mkdir(parents=True)
    text = json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False)
    (out / 'results.json').write_bytes(text.encode('utf-8'))
    print('wrote', out / 'results.json')
    print('sha256:', sha256_file(out / 'results.json'))
    print('capability rows:', len(capability), 'nc rows:', len(negative_control))


if __name__ == '__main__':
    main()
