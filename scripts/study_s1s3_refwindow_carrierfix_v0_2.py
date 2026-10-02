"""Dev-side S1/S3 reference-window recompute under the official carrier (v0.2).

Recomputes the frozen measured-peak registration (freeze_s1s3_reference_
window_v0_1) under both signed representations from ONE load per trace:
legacy 95 MHz carrier vs official real_bandpass (20 MHz, 2x amplitude).
Same geometry (STATIC sha-asserted), same exact full-line Fermat tables
(imported from the freeze script), same +-15 ns search band, same h5
frozen recovery directories and per-trace H5 hashes. Frozen contract and caches
are NOT touched; accepted output uses new v0.3 names.

No solver runs; BG development data only; no test family.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import sfcw_official_loader_v0_2 as L
from sfcw_carrierfix_acceptance import (accept_s1s3, provenance, require,
                                       sha256, verify_input, write_new_pair)
from freeze_s1s3_reference_window_v0_1 import (
    STATIC, STATIC_SHA256, TIERS, N_TRACES, ER_FERMAT, SEARCH_HALF_NS,
    NC_WIN, PRE_EVENT_BAND, RUN_DATE, fermat_times)

CONTRACT = ROOT / 'configs/research/s1s3_reference_window_contract_v0.1.json'
CONTRACT_SHA256 = '831326d83c37c3f2a17feeba29df694f7fde71c2f8abe45be120f86dc6ad4848'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-10-02_s1s3_refwindow_carrierfix_v0_3_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-10-02_s1s3_refwindow_carrierfix_v0_3_r2.json'


def measured_peak(t_f, sig, t):
    m = (t >= t_f - SEARCH_HALF_NS) & (t <= t_f + SEARCH_HALF_NS)
    assert m.any()
    idx = np.where(m)[0]
    j = idx[int(np.argmax(np.abs(sig[idx])))]
    return float(t[j])


def build():
    require(sha256(STATIC) == STATIC_SHA256, 'S1/S3 static geometry SHA mismatch')
    require(sha256(CONTRACT) == CONTRACT_SHA256, 'S1/S3 frozen contract SHA mismatch')
    L.verify_official_runtime()
    static = json.loads(STATIC.read_text(encoding='utf-8'))
    contract = json.loads(CONTRACT.read_text(encoding='utf-8'))
    families = {}
    audited_inputs = []
    for fam in TIERS:
        s = static['tiers'][fam]
        t_f18 = fermat_times(s['tan_theta'], s['zc_m'], s['domain_y_m'],
                             ER_FERMAT)
        frozen_family = contract['families'][fam]
        require([round(v, 4) for v in t_f18] == frozen_family['t_fermat_er18_ns'],
                'S1/S3 computed Fermat table differs from frozen contract')
        expected_inputs = {r['run_id']: r for r in frozen_family['inputs']}
        require(len(expected_inputs) == N_TRACES == len(frozen_family['inputs']),
                'S1/S3 frozen input coverage mismatch')
        mother = f'B2D-C3m{fam}-BG'
        rows = []
        t_axis = None
        for k in range(N_TRACES):
            rid = f'{mother}-CO33-t{k + 1:02d}'
            expected = expected_inputs[rid]
            h5 = L.SIMS / expected['source_dir'] / f'{rid}.h5'
            audited_inputs.append(verify_input(h5, expected))
            tr = L.trace_time_response(h5)
            verify_input(h5, expected)
            t = L.time_ns(tr)
            require(t_axis is None or np.array_equal(t_axis, t), 'S1/S3 reconstructed axes differ')
            t_axis = t
            sig_off = L.signed_official(tr)
            sig_leg = L.signed_legacy(tr)
            m_nc = (t >= NC_WIN[0]) & (t <= NC_WIN[1])
            m_fl = (t >= PRE_EVENT_BAND[0]) & (t <= PRE_EVENT_BAND[1])
            require(m_nc.any() and m_fl.any(), 'empty NC/floor RMS window')
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
        acceptance = accept_s1s3(rows, frozen_family, contract['constants'])
        families[fam] = {
            'rows': rows,
            'n_traces_peak_moved': n_moved,
            'legacy_matches_frozen_contract': legacy_matches_frozen,
            'acceptance': acceptance,
            'max_abs_peak_shift_ns': max(
                abs(r['t_measured_official_ns'] - r['t_measured_legacy_ns'])
                for r in rows),
        }
        print(fam, 'done; peaks moved:', n_moved, '/ 33;',
              'legacy==frozen:', legacy_matches_frozen)
    return {
        'schema': 's1s3_refwindow_carrierfix/2',
        'date': '2026-10-02',
        'basis': 'model design review 2026-10-02 §2; dev-side recompute, '
                 'frozen contract/caches untouched',
        'search_half_ns': SEARCH_HALF_NS,
        'run_date': RUN_DATE,
        'note': 'official real_bandpass uses a 2x real-part convention and '
                'a corrected carrier; RMS ratios need not equal two. Peak TIMES may move by integer '
                'reconstruction samples (0.832 ns grid). No solver runs; '
                'no test family.',
        'families': families,
        'provenance': provenance(__file__, audited_inputs, [CONTRACT, STATIC,
            ROOT / 'scripts/freeze_s1s3_reference_window_v0_1.py',
            ROOT / 'scripts/sfcw_carrierfix_acceptance.py']),
    }


def main():
    doc = build()
    doc2 = build()
    write_new_pair(OUT_R1, OUT_R2, doc, doc2)
    print('accepted S1/S3 families:', len(doc['families']))
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
