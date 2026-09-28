"""Tolerance sensitivity (±50%) flip-rate report: reward draft v0.1 §4 step 4 cont.

Reads the committed damage-ladder and operator-effects JSON outputs and applies
the proposed tolerances (τ_A=0.20, τ_D=0.95, ε(Nb ratio)=1.5 — proposal, NOT
frozen) to every record, then perturbs each tolerance axis independently by
±50% and counts judgment flips versus the central values.

Judgments:
  ladder record  pass  = A<=τ_A and D<=τ_D and (Nb_ratio is None or Nb_ratio<=ε)
  operator row   admissible = available and A<=τ_A and D_e<=τ_D and Nb_ratio<=ε
  operator best  per family: max contrast_db_delta among admissible rows

A flip is any change of pass/admissible/best under one-axis perturbation.
This is a sensitivity measurement of the tolerance PROPOSAL on a coarse ladder;
flips here mean the classification of ladder rungs/operator rows moves with the
tolerance, not that any tolerance value is wrong. Thresholds are not fitted to
candidate outcomes (Cawley–Talbot): the proposal comes from the ladder scale.

Deterministic pure-array computation; r1/r2 byte-identity asserted in-process.
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LADDER = ROOT / 'artifacts/research_checks/2026-09-28_t3_damage_ladder_r1.json'
OPERATORS = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-28_tolerance_sensitivity_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-28_tolerance_sensitivity_r2.json'

CENTRAL = {'tau_A': 0.20, 'tau_D': 0.95, 'eps_Nb': 1.5}
AXES = ['tau_A', 'tau_D', 'eps_Nb']
FACTORS = [0.5, 1.5]


def ladder_pass(r, tol):
    nb = r['Nb_ratio']
    return (r['A'] <= tol['tau_A'] and r['D'] <= tol['tau_D']
            and (nb is None or nb <= tol['eps_Nb']))


def operator_admissible(r, tol):
    return bool(r.get('available', True)) and \
        r['A'] <= tol['tau_A'] and r['D_e'] <= tol['tau_D'] and r['Nb_ratio'] <= tol['eps_Nb']


def judge(tol):
    ladder = json.loads(LADDER.read_text(encoding='utf-8'))['records']
    operators = json.loads(OPERATORS.read_text(encoding='utf-8'))['records']
    lad = {f"{r['family']}|{r['damage']}|{r['level']}": ladder_pass(r, tol) for r in ladder}
    ope = {f"{r['family']}|{r['config']}": operator_admissible(r, tol) for r in operators}
    best = {}
    fams = sorted({r['family'] for r in operators})
    for fam in fams:
        adm = [r for r in operators
               if r['family'] == fam and ope[f"{r['family']}|{r['config']}"]]
        best[fam] = (max(adm, key=lambda r: r['contrast_db_delta'])['config'] if adm else None)
    return lad, ope, best


def build():
    lad0, ope0, best0 = judge(CENTRAL)
    report = {'central_tolerances': CENTRAL, 'perturbation': 'one axis at a time, x0.5 / x1.5',
              'ladder_judgments_central': lad0, 'operator_admissible_central': ope0,
              'operator_best_central': best0, 'axes': []}
    for axis in AXES:
        for f in FACTORS:
            tol = dict(CENTRAL)
            tol[axis] = round(CENTRAL[axis] * f, 6)
            lad, ope, best = judge(tol)
            lad_flip = sorted(k for k in lad0 if lad[k] != lad0[k])
            ope_flip = sorted(k for k in ope0 if ope[k] != ope0[k])
            best_flip = sorted(k for k in best0 if best[k] != best0[k])
            report['axes'].append({
                'axis': axis, 'factor': f, 'value': tol[axis],
                'ladder_flips': lad_flip, 'ladder_flip_rate': round(len(lad_flip) / len(lad0), 6),
                'operator_flips': ope_flip, 'operator_flip_rate': round(len(ope_flip) / len(ope0), 6),
                'best_flips': best_flip})
    return report


def main():
    report = build()
    doc = {
        'schema': 'tolerance_sensitivity/1',
        'date': '2026-09-28',
        'basis': 'reward metrics draft v0.1 §4 step 4 (tolerance ±50% sensitivity); '
                 'tolerances are the unfrozen proposal from the damage ladder '
                 '(commit 6501783), derived from ladder scale, not from candidate outcomes',
        'inputs': {'ladder': LADDER.name, 'operator_effects': OPERATORS.name},
        'note': 'development family; flips measure proposal sensitivity on a coarse '
                'ladder; no tolerance is fitted or frozen here',
        **report,
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = dict(doc, **build())
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    n_lad = len(report['ladder_judgments_central'])
    n_ope = len(report['operator_admissible_central'])
    print('central judgments: ladder', n_lad, '| operator', n_ope)
    for a in report['axes']:
        print(f"{a['axis']} x{a['factor']}: ladder flips {len(a['ladder_flips'])} "
              f"| operator flips {len(a['operator_flips'])} | best flips {a['best_flips']}")
    print('records-axes:', len(report['axes']),
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
