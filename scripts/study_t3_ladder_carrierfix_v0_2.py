"""Dev-side t3 damage-ladder recompute under the official carrier (v0.2).

Re-runs the exact 2026-09-28 damage ladder (same mothers, same damage
instances, same seed path, same metric battery) twice per trace load:
once with the legacy 95 MHz carrier (rebuilt from the same complex
envelope) and once with the official real_bandpass (20 MHz carrier, 2x
amplitude). Output is a paired old/new/delta record set plus the archived
r1 values for direct comparison.

Frozen ladder outputs are NOT overwritten; this writes new files only.
No solver runs; no test-family data.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import sfcw_official_loader_v0_2 as L
from study_t3_damage_ladder import (
    MOTHERS, AMP_DB, SHIFT_SMP, NC_AMP_DB, SEED, NC_WIN, N_TRACES,
    fermat_times, event_mask, metrics, apply_damage)

ARCH_R1 = ROOT / 'artifacts/research_checks/2026-09-28_t3_damage_ladder_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-10-02_t3_ladder_carrierfix_v0_2_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-10-02_t3_ladder_carrierfix_v0_2_r2.json'

KEYS = ('D', 'a', 'A', 'H', 'rho', 'Nb_ratio',
        'arrival_drift_ns_median', 'arrival_drift_ns_max_abs')


def build_records():
    records = []
    for fam, mother, geo in MOTHERS:
        ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
            else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
        t, s_off, s_leg = L.load_bscan_both(mother, date='2026-09-28',
                                            n_traces=N_TRACES)
        t_ev = fermat_times(ifz)
        m_ev = event_mask(t, t_ev)
        m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s_off.shape)
        damages = [('identity', 0.0)]
        damages += [('amp_db', v) for v in AMP_DB]
        damages += [('polarity', 0.0)]
        damages += [('shift_smp', v) for v in SHIFT_SMP]
        damages += [('delete', 0.0)]
        damages += [('nc_amp', v) for v in NC_AMP_DB]
        for rep, s in (('legacy95', s_leg), ('official20', s_off)):
            rng = np.random.default_rng(SEED)
            for kind, level in damages:
                z = apply_damage(s, m_ev, m_nc, t, kind, level, rng)
                rec = {'family': fam, 'representation': rep,
                       'damage': kind, 'level': level,
                       **metrics(z, s, m_ev, m_nc, t, t_ev)}
                records.append(rec)
        print(fam, 'done:', len(damages), 'damage instances x 2 representations')
    return records


def pair_delta(records):
    """Join legacy/official rows of the same family+damage into deltas."""
    by_key = {}
    for r in records:
        by_key.setdefault((r['family'], r['damage'], r['level']), {})[
            r['representation']] = r
    out = []
    for (fam, kind, level), pair in sorted(by_key.items()):
        row = {'family': fam, 'damage': kind, 'level': level}
        for k in KEYS:
            a, b = pair['legacy95'][k], pair['official20'][k]
            row[f'{k}_legacy'] = a
            row[f'{k}_official'] = b
            row[f'{k}_delta'] = (None if a is None or b is None
                                 else round(b - a, 6))
        out.append(row)
    return out


def main():
    records = build_records()
    arch = json.loads(ARCH_R1.read_text(encoding='utf-8'))
    doc = {
        'schema': 't3_damage_ladder_carrierfix/1',
        'date': '2026-10-02',
        'basis': 'model design review 2026-10-02 §2 (75 MHz carrier bug); '
                 'dev-side recompute, frozen outputs untouched',
        'representations': {
            'legacy95': 'Re(env*exp(2j*pi*95MHz*t)) — buggy, kept for diff',
            'official20': 'tr.real_bandpass = 2*Re(env*exp(2j*pi*20MHz*t)) '
                          '— official gprMax chain'},
        'loader_module': 'scripts/sfcw_official_loader_v0_2.py',
        'archived_r1_sha256': hashlib.sha256(ARCH_R1.read_bytes()).hexdigest(),
        'note': 'envelope arrays and 501-pt complex frequency responses are '
                'unaffected by the carrier fix; only signed-waveform metrics '
                'are recomputed. No solver runs; no test-family data.',
        'records': records,
        'paired_delta': pair_delta(records),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    rec2 = build_records()
    doc2 = dict(doc, records=rec2, paired_delta=pair_delta(rec2))
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('records:', len(records),
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
