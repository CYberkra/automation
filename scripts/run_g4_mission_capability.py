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
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MISSION = ROOT / 'configs/research/g4_mission_tolerance_v0.1.json'
MISSION_SHA256 = '2d630c01cc95f3c13febbc19c4bf8f3b5a53b5fb82bbabcb2f68e7537ef78b28'
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
CHUNK_GLOB = 'artifacts/research_checks/2026-09-26_damage_ladder_r1_c{:02d}'
N_CHUNKS = 8
DEV_FAMILIES = ('c1', 'c3')

MISSION_LEVELS = {
    ('amplitude_scale', 0.5): 'mission_relevant',
    ('amplitude_scale', 0.1): 'mission_relevant',
    ('amplitude_scale', 0.9): 'weak_event_band',
    ('polarity_flip', 1.0): 'mission_relevant',
    ('sample_shift', 4): 'mission_relevant',
    ('sample_shift', 16): 'mission_relevant',
    ('sample_shift', -4): 'mission_relevant',
    ('sample_shift', -16): 'mission_relevant',
    ('sample_shift', 1): 'weak_event_band',
    ('sample_shift', -1): 'weak_event_band',
    ('trace_deletion', 4): 'mission_relevant',
    ('trace_deletion', 8): 'mission_relevant',
    ('trace_deletion', 1): 'weak_event_band',
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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
    args = ap.parse_args()
    out = Path(args.output_dir)
    assert not out.exists(), f'output dir must not exist: {out}'

    assert sha256_file(MISSION) == MISSION_SHA256, 'mission tolerance gate mismatch'
    assert sha256_file(EVENT_TABLE) == EVENT_TABLE_SHA256, 'event table gate mismatch'
    recs = load_chunks()
    assert len(recs) == 34086, f'unexpected ladder record count: {len(recs)}'

    # identity anchor self-consistency (constructed reference)
    idents = [r for r in recs if r['candidate_id'] == 'B0_G1_BG' and r['availability'] == 'ran']
    anchor_ok = all(r['metrics']['waveform']['metrics']['nrmse'] == 0.0 for r in idents)
    assert anchor_ok, 'identity anchor D!=0 found'

    cap = {}
    nc = {}
    for r in recs:
        if r['family'] not in DEV_FAMILIES or r['availability'] != 'ran':
            continue
        cid = r['candidate_id']
        geo = r['geometry']
        fam = r['family']
        base = r['baseline']
        if r['row_type'] == 'event':
            dt_ = r['damage']['type']
            lvl = r['damage']['level']
            cls = MISSION_LEVELS.get((dt_, lvl))
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

    mission = json.loads(MISSION.read_text(encoding='utf-8'))
    result = {
        'schema': 'g4-mission-capability/1',
        'plan': 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md',
        'mission_tolerance_config': str(MISSION.relative_to(ROOT)),
        'mission_tolerance_sha256': MISSION_SHA256,
        'event_table_sha256': EVENT_TABLE_SHA256,
        'ladder_source': CHUNK_GLOB.replace('{:02d}', '00..07') + ' (r1; r2 byte-identical after resource strip per acceptance)',
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
    assert mission['status'] == 'frozen'

    out.mkdir(parents=True)
    text = json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False)
    (out / 'results.json').write_bytes(text.encode('utf-8'))
    print('wrote', out / 'results.json')
    print('sha256:', sha256_file(out / 'results.json'))
    print('capability rows:', len(capability), 'nc rows:', len(negative_control))


if __name__ == '__main__':
    main()
