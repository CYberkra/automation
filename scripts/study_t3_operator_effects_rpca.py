"""RPCA addendum to the t3 operator effects table: draft §4 step 4 completion.

The 2026-09-28 ladder-results plan (docs/research/2026-09-28_t3_damage_ladder_
results.md §3) listed RPCA among the operators for the effects table, but the
executed table (study_t3_operator_effects.py) did not include it and the
catalogue (operator_catalogue_v0.1) has no RPCA entry. This addendum closes
that gap as an OUT-OF-CATALOGUE mechanism reference (same status as
local_mean_control), on the identical loading chain (batch2d_slope_t3_co
archived h5, official SFCW chain, Hann window, signed radargram at 95 MHz).

Three lambda runs: standard PCP choice lam0 = 1/sqrt(max(shape)) plus half and
double, to bracket the sparsity-weight sensitivity. Same battery as the main
table: D_e/a/A/H/rho on the output (no end-gain involved), a_raw_event,
contrast delta vs the NC floor (250-400 ns), Nb_ratio.

Deterministic; r1/r2 byte-identity asserted. Development families only; no
solver runs; no test-family data; no frozen contract touched.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from research_operator_contract import rpca_control
from study_t3_operator_effects import battery, contrast_db
from study_t3_damage_ladder import (MOTHERS, NC_WIN, load_bscan, fermat_times,
                                    event_mask)

ROOT = Path(__file__).resolve().parents[1]
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_t3_operator_effects_rpca_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_t3_operator_effects_rpca_r2.json'


def run_family(fam, mother, geo):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)          # [trace, sample]
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev)         # [trace, sample]
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s.shape)
    x = np.ascontiguousarray(s.T)      # [sample, trace]
    m_ev_t = m_ev.T
    m_nc_t = m_nc.T
    c_in = contrast_db(x, m_ev_t, m_nc_t)
    nb_in = float(np.mean(x[m_nc_t] ** 2))
    lam0 = 1.0 / (max(x.shape) ** 0.5)
    rows = []
    for lam_factor in (0.5, 1.0, 2.0):
        lam = lam0 * lam_factor
        y, diag = rpca_control(x, lam)
        b = battery(y, x, m_ev_t)
        rows.append({'family': fam, 'mother': mother,
                     'config': f'rpca_lam{lam_factor:g}x_lam0(out-of-catalogue)',
                     'available': True,
                     **b,
                     'a_raw_event': round(battery(y, x, m_ev_t)['a'], 6),
                     'contrast_db_in': round(c_in, 3),
                     'contrast_db_out': round(contrast_db(y, m_ev_t, m_nc_t), 3),
                     'contrast_db_delta': round(contrast_db(y, m_ev_t, m_nc_t) - c_in, 3),
                     'Nb_ratio': round(float(np.mean(y[m_nc_t] ** 2)) / nb_in, 6),
                     'background': 'rpca', 'parameter': round(lam, 9),
                     'end_gain': 1, 'order': 'B',
                     'diag': {k: diag[k] for k in
                              ('lam_over_standard', 'status', 'rank_L',
                               'iterations', 'rel_residual', 'sparse_fraction')}})
    print(fam, 'done:', len(rows), 'rpca rows')
    return rows


def build_records():
    records = []
    for fam, mother, geo in MOTHERS:
        records.extend(run_family(fam, mother, geo))
    return records


def main():
    doc = {
        'schema': 't3_operator_effects_rpca/1',
        'date': '2026-09-29',
        'basis': 'addendum to docs/research/2026-09-28_t3_operator_effects_and_'
                 'sensitivity.md; closes the RPCA item of the ladder-results '
                 'plan §3 (RPCA was planned but absent from the executed table '
                 'and from operator_catalogue_v0.1)',
        'data': 'batch2d_slope_t3_co archived h5 (BG only); official SFCW chain, '
                'Hann window, signed radargram at 95 MHz band centre; loading '
                'chain identical to study_t3_operator_effects / '
                'study_t3_damage_ladder',
        'operators': {'out_of_catalogue': 'rpca_control (inexact ALM principal '
                      'component pursuit, Lin-Chen-Ma 2009), lam0 = '
                      '1/sqrt(max(sample,trace)); runs at 0.5x / 1.0x / 2.0x; '
                      'no end-gain; mechanism reference, not a catalogue '
                      'candidate (catalogue freeze untouched)'},
        'metric_notes': ['same battery and windows as the 2026-09-28 effects '
                         'table; D_e computed directly on the rpca output (no '
                         'gain compensation involved)',
                         'contrast = 20log10(event peak / NC-window RMS); NC '
                         'window 250-400 ns; event window = Fermat line '
                         '(er=18) +/- 10 ns'],
        'note': 'development family; no solver runs; no test-family data; feeds '
                'the user-confirmation step of reward weights v0.1',
        'records': build_records(),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    rec2 = build_records()
    doc2 = dict(doc, records=rec2)
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('records:', len(doc['records']),
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
