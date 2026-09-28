"""Operator effects table on the t3 BG CO B-scans: reward draft v0.1 §4 step 4.

Applies catalogue background/gain configurations (research_operator_contract,
operator_catalogue_v0.1) plus the out-of-catalogue local-mean control to the
archived batch2d_slope_t3_co signed radargrams (official SFCW chain, Hann
window; same loading chain as study_t3_damage_ladder) and records, per
family x configuration:

  D_e/a/A/H/rho   event-window preservation battery (v0.2 style) computed on
                  the gain-compensated output (audit_before_gain) vs the input,
                  so the fixed shared gain itself is not counted as damage
  a_raw           event-window amplitude factor of the raw output (shows the
                  gain the configuration actually applied to the event)
  contrast_db_*   event-peak to NC-floor contrast before/after and its delta
                  (the background-suppression positive indicator; slope family
                  reference value ~80 dB for mean removal, reward draft §5)
  Nb_ratio        NC-window (250-400 ns) mean-square ratio out/in (negative
                  indicator; interference amplification)

Expected key row: mean removal (lambda=1.0) on the flat C3mX family drives
D_e toward ~1 (the horizontal interface is absorbed by the mean trace) while
the same configuration on the slope families yields the ~80 dB contrast
gain — the measured positive/negative duality the reward design must encode.

Deterministic; run writes r1 then reruns in-process and asserts byte-identity.
Development family only; no solver runs; no test-family data.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from research_operator_contract import (ConfigUnavailable, apply_configuration,
                                        local_mean_control)
from study_t3_damage_ladder import (MOTHERS, NC_WIN, load_bscan, fermat_times,
                                    event_mask)

ROOT = Path(__file__).resolve().parents[1]
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r2.json'

CONFIGS = ['B0_G1_BG', 'B2_G1_BG', 'B3_G1_BG', 'B4_G1_BG', 'B5_G1_BG',
           'B0_G2_BG', 'B0_G4_BG', 'B3_G2_BG', 'B4_G2_BG']
LOCAL_MEAN_WIDTH = 11  # odd, <= 33 traces; out-of-catalogue adaptive reference


def battery(z, s, m_ev):
    """v0.2-style preservation battery of z vs reference s on the event window."""
    zs, ss = z[m_ev], s[m_ev]
    denom = float(np.linalg.norm(ss))
    assert denom > 0
    a = float(np.dot(zs, ss) / np.dot(ss, ss))
    zn = float(np.linalg.norm(zs))
    rho = float(np.dot(zs, ss) / (zn * denom)) if zn > 0 else None
    return {
        'D_e': round(float(np.linalg.norm(zs - ss) / denom), 6),
        'a': round(a, 6),
        'A': round(abs(a - 1.0), 6),
        'H': round(float(np.linalg.norm(zs - a * ss) / denom), 6),
        'rho': (round(rho, 6) if rho is not None else None),
    }


def contrast_db(z, m_ev, m_nc):
    peak = float(np.max(np.abs(z[m_ev])))
    floor = float(np.sqrt(np.mean(z[m_nc] ** 2)))
    assert floor > 0
    return 20.0 * np.log10(peak / floor)


def run_family(fam, mother, geo):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)          # [trace, sample]
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev)         # [trace, sample]
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s.shape)
    x = np.ascontiguousarray(s.T)      # contract axes: [sample, trace]
    m_ev_t = m_ev.T
    m_nc_t = m_nc.T
    c_in = contrast_db(x, m_ev_t, m_nc_t)
    nb_in = float(np.mean(x[m_nc_t] ** 2))
    rows = []
    for cfg in CONFIGS:
        try:
            r = apply_configuration(x, cfg)
        except ConfigUnavailable as exc:
            # Contract guard: configuration legitimately unavailable for this
            # input (e.g. svd_cutoff_gap_unresolved). Recorded, not patched.
            rows.append({'family': fam, 'mother': mother, 'config': cfg,
                         'available': False, 'failure_reason': str(exc)})
            continue
        y_ag = r['audit_before_gain']  # gain-compensated: non-gain damage only
        y = r['output']
        row = {'family': fam, 'mother': mother, 'config': cfg, 'available': True,
               **battery(y_ag, x, m_ev_t),
               'a_raw_event': round(battery(y, x, m_ev_t)['a'], 6),
               'contrast_db_in': round(c_in, 3),
               'contrast_db_out': round(contrast_db(y, m_ev_t, m_nc_t), 3),
               'contrast_db_delta': round(contrast_db(y, m_ev_t, m_nc_t) - c_in, 3),
               'Nb_ratio': round(float(np.mean(y[m_nc_t] ** 2)) / nb_in, 6),
               'background': r['config']['background'],
               'parameter': r['config']['parameter'],
               'end_gain': r['config']['end_gain'],
               'order': ''.join(r['config']['order'])}
        rows.append(row)
    y_lm = local_mean_control(x, LOCAL_MEAN_WIDTH, strength=1.0)
    rows.append({'family': fam, 'mother': mother,
                 'config': f'local_mean_w{LOCAL_MEAN_WIDTH}_s1.0(out-of-catalogue)',
                 'available': True,
                 **battery(y_lm, x, m_ev_t),
                 'a_raw_event': round(battery(y_lm, x, m_ev_t)['a'], 6),
                 'contrast_db_in': round(c_in, 3),
                 'contrast_db_out': round(contrast_db(y_lm, m_ev_t, m_nc_t), 3),
                 'contrast_db_delta': round(contrast_db(y_lm, m_ev_t, m_nc_t) - c_in, 3),
                 'Nb_ratio': round(float(np.mean(y_lm[m_nc_t] ** 2)) / nb_in, 6),
                 'background': 'local_mean', 'parameter': LOCAL_MEAN_WIDTH,
                 'end_gain': 1, 'order': 'B'})
    print(fam, 'done:', len(rows), 'operator rows')
    return rows


def build_records():
    records = []
    for fam, mother, geo in MOTHERS:
        records.extend(run_family(fam, mother, geo))
    return records


def main():
    doc = {
        'schema': 't3_operator_effects/1',
        'date': '2026-09-28',
        'basis': 'reward metrics draft v0.1 §4 step 4 (operator effects table); '
                 'damage ladder + tolerance proposal committed at 6501783',
        'data': 'batch2d_slope_t3_co archived h5 (BG only); official SFCW chain, '
                'Hann window, signed radargram at 95 MHz band centre; loading chain '
                'identical to study_t3_damage_ladder',
        'operators': {'catalogue': 'operator_catalogue_v0.1.json subset '
                      '(background identity/mean/svd x end_gain 1/2/4, BG order); '
                      'axes [sample, trace]; mask=None (dense valid window)',
                      'out_of_catalogue': 'local_mean_control width=%d strength=1.0 '
                      '(mechanism reference, not a catalogue candidate)' % LOCAL_MEAN_WIDTH},
        'metric_notes': ['D_e battery computed on audit_before_gain (gain-compensated) '
                         'so the fixed shared gain is not counted as damage',
                         'a_raw_event shows the amplitude factor the configuration '
                         'actually applied inside the event window',
                         'contrast = 20log10(event peak / NC-window RMS); NC window '
                         '250-400 ns; event window = Fermat line (er=18) +/- 10 ns'],
        'note': 'development family; no solver runs; no test-family data; effects '
                'table feeds tolerance sensitivity check (draft §4 step 4 cont.)',
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
