"""Dev-side S1/S3 reference-window recompute under the official carrier (v0.2).

Recomputes the frozen measured-peak registration (freeze_s1s3_reference_
window_v0_1) under both signed representations from ONE load per trace:
legacy 95 MHz carrier vs official real_bandpass (20 MHz, 2x amplitude).
Same geometry (STATIC sha-asserted), same exact full-line Fermat tables
(imported from the freeze script), same +-15 ns search band, same h5
suffix search. Frozen contract and caches are NOT touched.

No solver runs; BG development data only; no test family.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import sfcw_official_loader_v0_2 as L
from freeze_s1s3_reference_window_v0_1 import (
    STATIC, STATIC_SHA256, TIERS, N_TRACES, ER_FERMAT, SEARCH_HALF_NS,
    NC_WIN, PRE_EVENT_BAND, RUN_DATE, find_h5, fermat_times)

CONTRACT = ROOT / 'configs/research/s1s3_reference_window_contract_v0.1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-10-02_s1s3_refwindow_carrierfix_v0_2_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-10-02_s1s3_refwindow_carrierfix_v0_2_r2.json'


def measured_peak(t_f, sig, t):
    m = (t >= t_f - SEARCH_HALF_NS) & (t <= t_f + SEARCH_HALF_NS)
    assert m.any()
    idx = np.where(m)[0]
    j = idx[int(np.argmax(np.abs(sig[idx])))]
    return float(t[j])


def build():
    assert hashlib.sha256(STATIC.read_bytes()).hexdigest() == STATIC_SHA256
    static = json.loads(STATIC.read_text(encoding='utf-8'))
    contract = json.loads(CONTRACT.read_text(encoding='utf-8'))
    families = {}
    for fam in TIERS:
        s = static['tiers'][fam]
        t_f18 = fermat_times(s['tan_theta'], s['zc_m'], s['domain_y_m'],
                             ER_FERMAT)
        mother = f'B2D-C3m{fam}-BG'
        rows = []
        for k in range(N_TRACES):
            rid = f'{mother}-CO33-t{k + 1:02d}'
            h5 = find_h5(rid)
            assert h5 is not None, rid
            tr = L.trace_time_response(h5)
            t = L.time_ns(tr)
            sig_off = L.signed_official(tr)
            sig_leg = L.signed_legacy(tr)
            m_nc = (t >= NC_WIN[0]) & (t <= NC_WIN[1])
            m_fl = (t >= PRE_EVENT_BAND[0]) & (t <= PRE_EVENT_BAND[1])
            rows.append({
                'trace': k + 1,
                't_measured_legacy_ns': round(
                    measured_peak(t_f18[k], sig_leg, t), 4),
                't_measured_official_ns': round(
                    measured_peak(t_f18[k], sig_off, t), 4),
                'nc_rms_legacy': round(float(np.sqrt(np.mean(sig_leg[m_nc] ** 2))), 12),
                'nc_rms_official': round(float(np.sqrt(np.mean(sig_off[m_nc] ** 2))), 12),
                'floor_rms_legacy': round(float(np.sqrt(np.mean(sig_leg[m_fl] ** 2))), 12),
                'floor_rms_official': round(float(np.sqrt(np.mean(sig_off[m_fl] ** 2))), 12),
            })
        frozen_t = contract['families'][fam]['t_measured_ns']
        for k, row in enumerate(rows):
            row['t_measured_frozen_ns'] = frozen_t[k]
            row['peak_shift_samples'] = round(
                (row['t_measured_official_ns']
                 - row['t_measured_legacy_ns']) / (t[1] - t[0]), 3)
        n_moved = sum(1 for r in rows if r['peak_shift_samples'] != 0)
        legacy_matches_frozen = all(
            abs(r['t_measured_legacy_ns'] - r['t_measured_frozen_ns']) < 1e-9
            for r in rows)
        families[fam] = {
            'rows': rows,
            'n_traces_peak_moved': n_moved,
            'legacy_matches_frozen_contract': legacy_matches_frozen,
            'max_abs_peak_shift_ns': max(
                abs(r['t_measured_official_ns'] - r['t_measured_legacy_ns'])
                for r in rows),
        }
        print(fam, 'done; peaks moved:', n_moved, '/ 33;',
              'legacy==frozen:', legacy_matches_frozen)
    return {
        'schema': 's1s3_refwindow_carrierfix/1',
        'date': '2026-10-02',
        'basis': 'model design review 2026-10-02 §2; dev-side recompute, '
                 'frozen contract/caches untouched',
        'search_half_ns': SEARCH_HALF_NS,
        'run_date': RUN_DATE,
        'note': 'official real_bandpass = 2x legacy amplitude, so absolute '
                'RMS readouts double; peak TIMES may move by integer '
                'reconstruction samples (0.832 ns grid). No solver runs; '
                'no test family.',
        'families': families,
    }


def main():
    doc = build()
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    doc2 = build()
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
