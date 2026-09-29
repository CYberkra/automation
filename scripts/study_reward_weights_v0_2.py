"""Scalar reward weights v0.2: re-issue of the §4 step-5 ranking with RPCA rows.

Identical formula, tolerances (SHA-asserted frozen contract) and sensitivity
protocol as v0.1 (study_reward_weights_v0_1.py, 2026-09-28); the only change is
the candidate table: the 2026-09-29 RPCA addendum rows join the 2026-09-28
effects rows. RPCA rows are reported under their catalogue v0.2 ids
(B7/B8/B9_G1_BG = lambda factor 0.5/1.0/2.0 x lambda0); local_mean stays an
out-of-catalogue mechanism reference. Still UNFROZEN: user confirmation
required before any contract freeze. r1/r2 byte-identity asserted.
"""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from study_reward_weights_v0_1 import CONTRACT, CONTRACT_SHA256, reward, rank_flips

EFFECTS_V1 = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r1.json'
RPCA_R1 = ROOT / 'artifacts/research_checks/2026-09-29_t3_operator_effects_rpca_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_reward_weights_v0_2_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_reward_weights_v0_2_r2.json'

RPCA_ID = {'rpca_lam0.5x_lam0(out-of-catalogue)': 'B7_G1_BG',
           'rpca_lam1x_lam0(out-of-catalogue)': 'B8_G1_BG',
           'rpca_lam2x_lam0(out-of-catalogue)': 'B9_G1_BG'}


def combined_records():
    rows = json.loads(EFFECTS_V1.read_text(encoding='utf-8'))['records']
    rows = rows + json.loads(RPCA_R1.read_text(encoding='utf-8'))['records']
    for r in rows:
        r['config'] = RPCA_ID.get(r['config'], r['config'])
    return rows


def evaluate(tol):
    out = []
    for r in combined_records():
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


def build():
    assert hashlib.sha256(CONTRACT.read_bytes()).hexdigest() == CONTRACT_SHA256, \
        'frozen tolerance contract hash mismatch'
    tol = {k: v['value'] for k, v in
           json.loads(CONTRACT.read_text(encoding='utf-8'))['tolerances'].items()
           if k in ('tau_A', 'tau_D', 'eps_Nb')}
    central = evaluate(tol)
    axes = []
    recs = {r['config']: r for r in combined_records()}
    for coef in ('c_d', 'c_b'):
        for f in (0.5, 1.5):
            kw = {coef: f}
            varied = []
            for r in combined_records():
                base = {'family': r['family'], 'config': r['config']}
                match = next(c for c in central
                             if c['family'] == r['family'] and c['config'] == r['config'])
                if not match['feasible']:
                    varied.append({**base, 'feasible': False, 'reason': match['reason']})
                    continue
                varied.append({**base, 'feasible': True, **reward(r, **kw)})
            axes.append({'coefficient': coef, 'factor': f,
                         'pairwise_inversions': rank_flips(central, varied),
                         'n_inversions': len(rank_flips(central, varied))})
    return {'tolerances_frozen': tol,
            'formula': 'R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) '
                       '- max(10*log10(Nb_ratio), 0); coefficients 1 by dB-unit '
                       'construction (unchanged from v0.1)',
            'rows': central, 'sensitivity': axes}


def main():
    report = build()
    doc = {
        'schema': 'reward_weights/0.2',
        'date': '2026-09-29',
        'basis': 'v0.1 (2026-09-28, unfrozen) re-issued with the RPCA addendum '
                 'rows under catalogue v0.2 ids; tolerances still frozen in '
                 'reward_tolerance_contract_v0.1.json (SHA asserted)',
        'weight_derivation': 'coefficients = 1 by dB-unit identities (unchanged, '
                             'Cawley-Talbot); RPCA rows added as candidates only',
        'note': 'background-suppression class scalar proposal (UNFROZEN, user '
                'confirmation required); gain-class quality labels undetermined',
        **report,
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = dict(doc, **build())
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
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
              '| inversions:', a['n_inversions'], a['pairwise_inversions'][:4])


if __name__ == '__main__':
    main()
