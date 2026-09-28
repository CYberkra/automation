"""Scalar reward weight derivation v0.1 (reward draft §4 step 5, background class).

Draft §2.1 proposed R_bg = w_c*log(contrast gain) - preservation penalties
- w_b*max(log Nb, 0), weights TBD. Step 5 derives the weights from dB-unit
identities and the frozen tolerance contract — NOT from candidate outcomes
(Cawley-Talbot clause):

  R_bg = Δcontrast_dB                      (measured positive term)
         - 20*log10(1/(1-D_e))             (preservation: D_e → dB-equivalent by
                                            the 20log10 amplitude identity; the
                                            -1 dB ladder rung maps to exactly
                                            1.0 dB, so the coefficient is 1 by
                                            construction, not fitted)
         - max(10*log10(Nb_ratio), 0)      (noise amplification in dB; only the
                                            positive part counts — suppression
                                            of NC energy is already credited in
                                            Δcontrast)

All three terms are in dB, so w_c = w_d = w_b = 1 by the equal-dB-equal-scale
principle: reward = NET event-to-interference improvement. This also encodes
draft §2.2 ("brightness is not SNR") by construction: a pure gain changes
Δcontrast by ~0 and pays the Nb penalty, scoring below identity.

Feasibility uses the FROZEN contract reward_tolerance_contract_v0.1.json
(SHA-256 asserted). Gain-class quality labels remain undetermined (v0.2 §6:
only known-multiplicative-attenuation controls carry quality labels; those
cases do not exist yet) — gain rows here are diagnostics within the shared
table, not a gain-class reward.

Sensitivity: the two derived coefficients perturbed x0.5/x1.5 one at a time;
pairwise ranking inversions per family reported. r1/r2 byte-identity asserted.
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
CONTRACT_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'
OPERATORS = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-28_reward_weights_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-28_reward_weights_r2.json'


def reward(row, c_d=1.0, c_b=1.0):
    pres = -20.0 * np_log10(1.0 - row['D_e']) if row['D_e'] < 1.0 else float('inf')
    nb = max(10.0 * np_log10(row['Nb_ratio']), 0.0)
    total = row['contrast_db_delta'] - c_d * pres - c_b * nb
    return {'preservation_db': round(pres, 6), 'nb_penalty_db': round(nb, 6),
            'R_bg_db': round(total, 6)}


def np_log10(x):
    import math
    return math.log10(x)


def evaluate(tol):
    rows = json.loads(OPERATORS.read_text(encoding='utf-8'))['records']
    out = []
    for r in rows:
        base = {'family': r['family'], 'config': r['config']}
        if not r.get('available', True):
            out.append({**base, 'feasible': False, 'reason': r['failure_reason']})
            continue
        gates = {'A': r['A'] <= tol['tau_A'], 'D_e': r['D_e'] <= tol['tau_D'],
                 'Nb': r['Nb_ratio'] <= tol['eps_Nb']}
        if not all(gates.values()):
            out.append({**base, 'feasible': False,
                        'reason': 'constraint:' + ','.join(k for k, v in gates.items() if not v)})
            continue
        out.append({**base, 'feasible': True, 'D_e': r['D_e'], 'A': r['A'],
                    'Nb_ratio': r['Nb_ratio'], 'contrast_db_delta': r['contrast_db_delta'],
                    **reward(r)})
    return out


def rank_flips(rows_a, rows_b):
    """Pairwise inversions of R among jointly feasible rows, per family."""
    flips = []
    fams = sorted({r['family'] for r in rows_a})
    for fam in fams:
        ra = {r['config']: r['R_bg_db'] for r in rows_a
              if r['family'] == fam and r['feasible']}
        rb = {r['config']: r['R_bg_db'] for r in rows_b
              if r['family'] == fam and r['feasible']}
        common = sorted(set(ra) & set(rb))
        for i in range(len(common)):
            for j in range(i + 1, len(common)):
                ca, cb = common[i], common[j]
                if (ra[ca] - ra[cb]) * (rb[ca] - rb[cb]) < 0:
                    flips.append(f'{fam}|{ca}|{cb}')
    return flips


def build():
    assert hashlib.sha256(CONTRACT.read_bytes()).hexdigest() == CONTRACT_SHA256, \
        'frozen tolerance contract hash mismatch'
    tol = {k: v['value'] for k, v in
           json.loads(CONTRACT.read_text(encoding='utf-8'))['tolerances'].items()
           if k in ('tau_A', 'tau_D', 'eps_Nb')}
    central = evaluate(tol)
    axes = []
    for coef in ('c_d', 'c_b'):
        for f in (0.5, 1.5):
            rows = json.loads(OPERATORS.read_text(encoding='utf-8'))['records']
            kw = {coef: f}
            varied = []
            for r in rows:
                base = {'family': r['family'], 'config': r['config']}
                match = next(c for c in central
                             if c['family'] == r['family'] and c['config'] == r['config'])
                if not match['feasible']:
                    varied.append({**base, 'feasible': False, 'reason': match['reason']})
                    continue
                varied.append({**base, 'feasible': True, **reward(r, **kw)})
            flips = rank_flips(central, varied)
            axes.append({'coefficient': coef, 'factor': f,
                         'pairwise_inversions': flips, 'n_inversions': len(flips)})
    return {'tolerances_frozen': tol, 'formula': 'R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) '
            '- max(10*log10(Nb_ratio), 0); coefficients 1 by dB-unit construction',
            'rows': central, 'sensitivity': axes}


def main():
    report = build()
    doc = {
        'schema': 'reward_weights/0.1',
        'date': '2026-09-28',
        'basis': 'reward draft v0.1 §4 step 5; tolerances frozen in '
                 'reward_tolerance_contract_v0.1.json (user "确认", commit 619170f)',
        'weight_derivation': 'coefficients = 1 by dB-unit identities (equal dB = equal '
                             'scale; net event-to-interference improvement); NOT fitted '
                             'to candidate outcomes (Cawley-Talbot)',
        'note': 'background-suppression class scalar proposal (unfrozen, user '
                'confirmation required before contract freeze); gain-class quality '
                'labels undetermined (v0.2 §6: no controlled-attenuation cases yet)',
        **report,
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = dict(doc, **build())
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    for fam in ('S2X', 'S2TZX', 'C3mX'):
        feas = sorted((r for r in report['rows'] if r['family'] == fam and r['feasible']),
                      key=lambda r: -r['R_bg_db'])
        print(fam, 'feasible ranking:')
        for r in feas:
            print('  %-44s R=%+8.2f dB (contrast %+6.1f, pres %5.2f, nb %4.2f)'
                  % (r['config'], r['R_bg_db'], r['contrast_db_delta'],
                     r['preservation_db'], r['nb_penalty_db']))
        infeas = [r['config'] for r in report['rows'] if r['family'] == fam and not r['feasible']]
        print('  infeasible:', infeas)
    for a in report['sensitivity']:
        print('sensitivity', a['coefficient'], 'x', a['factor'],
              '| inversions:', a['n_inversions'])
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
