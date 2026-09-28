"""Freeze the reward tolerance contract v0.1 (reward draft §4: user-confirmed proposal).

User confirmed the tolerance proposal ("确认", 2026-09-28) after the damage
ladder (6501783) and operator-effects + sensitivity (99785ff) units. This
script writes configs/research/reward_tolerance_contract_v0.1.json with the
full evidence hashes embedded and asserts every anchor claim against the
committed evidence files before writing, so the contract cannot drift from
its basis (freeze-script assertion discipline).

Frozen tolerances (development-family scope, NOT physical field thresholds):
  tau_A = 0.20   event amplitude-factor deviation (ladder: -1 dB pass / -3 dB fail)
  tau_D = 0.95   event-window total distortion (1-sample shift pass / delete, >=2-sample fail)
  eps_Nb = 1.5   negative-control window energy ratio limit (structural floor 1.0:
                 eps < 1 would forbid identity, measured anchor in sensitivity r1)
  arrival drift  diagnostic only, not a hard constraint
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
LADDER = ROOT / 'artifacts/research_checks/2026-09-28_t3_damage_ladder_r1.json'
OPERATORS = ROOT / 'artifacts/research_checks/2026-09-28_t3_operator_effects_r1.json'
SENSITIVITY = ROOT / 'artifacts/research_checks/2026-09-28_tolerance_sensitivity_r1.json'

TOL = {'tau_A': 0.20, 'tau_D': 0.95, 'eps_Nb': 1.5}


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_anchors():
    lad = json.loads(LADDER.read_text(encoding='utf-8'))['records']
    sens = json.loads(SENSITIVITY.read_text(encoding='utf-8'))

    def lrow(fam, dmg, lvl):
        return next(r for r in lad
                    if r['family'] == fam and r['damage'] == dmg and r['level'] == lvl)

    # τ_A anchor: -1 dB rung below, -3 dB rung above, all three families.
    for fam in ('S2X', 'S2TZX', 'C3mX'):
        assert lrow(fam, 'amp_db', -1.0)['A'] < TOL['tau_A'] < lrow(fam, 'amp_db', -3.0)['A']
    # τ_D anchor: 1-sample shift below, delete above; identity exact zero.
    for fam in ('S2X', 'S2TZX', 'C3mX'):
        assert lrow(fam, 'identity', 0.0)['D'] == 0.0
        assert lrow(fam, 'shift_smp', 1)['D'] < TOL['tau_D'] < lrow(fam, 'delete', 0.0)['D']
    # ε anchor: pure gain inside, nc +6 dB outside; identity exactly 1.0.
    ope = json.loads(OPERATORS.read_text(encoding='utf-8'))['records']
    gain_nb = [r['Nb_ratio'] for r in ope if r['config'] == 'B0_G4_BG']
    assert gain_nb and max(gain_nb) < TOL['eps_Nb'] < lrow('S2X', 'nc_amp', 6.0)['Nb_ratio']
    assert lrow('S2X', 'identity', 0.0)['Nb_ratio'] == 1.0
    # Sensitivity anchors: τ_D zero flips both directions; ε×0.5 forbids identity.
    axes = {(a['axis'], a['factor']): a for a in sens['axes']}
    assert axes[('tau_D', 0.5)]['ladder_flips'] == axes[('tau_D', 1.5)]['ladder_flips'] == []
    assert axes[('tau_D', 0.5)]['operator_flips'] == axes[('tau_D', 1.5)]['operator_flips'] == []
    eps_half = axes[('eps_Nb', 0.5)]
    assert 'S2X|identity|0.0' in eps_half['ladder_flips']
    assert 'S2X|B0_G1_BG' in eps_half['operator_flips']
    assert axes[('eps_Nb', 1.5)]['ladder_flips'] == axes[('eps_Nb', 1.5)]['operator_flips'] == []
    return sens


def main():
    sens = verify_anchors()
    doc = {
        'schema': 'reward-tolerance-contract/0.1',
        'status': 'frozen',
        'date': '2026-09-28',
        'frozen_by': 'user confirmation "确认" (2026-09-28) of the proposal in '
                     'docs/research/2026-09-28_t3_damage_ladder_results.md §2, after '
                     'operator-effects + sensitivity evidence (docs/research/'
                     '2026-09-28_t3_operator_effects_and_sensitivity.md)',
        'tolerances': {
            'tau_A': {'value': TOL['tau_A'],
                      'meaning': 'event amplitude-factor deviation |a-1| limit',
                      'anchor': 'damage ladder: -1 dB rung passes (A=0.109), -3 dB rung '
                                'fails (A=0.292); resolution limited by ladder spacing '
                                '(+/-50% perturbation moves the rung boundary, recorded '
                                'as ladder property, not refitted)'},
            'tau_D': {'value': TOL['tau_D'],
                      'meaning': 'event-window total distortion D limit',
                      'anchor': 'ladder: 1-sample (0.832 ns) shift passes (D~0.86-0.89), '
                                'delete fails (D=1.0); +/-50% sensitivity: zero flips '
                                'both directions (robust on current evidence set)'},
            'eps_Nb': {'value': TOL['eps_Nb'],
                       'meaning': 'negative-control window (250-400 ns) mean-square '
                                  'ratio out/in limit',
                       'anchor': 'structural floor 1.0 measured: eps x0.5 = 0.75 makes '
                                 'identity inadmissible (sensitivity r1); 1.5 margin '
                                 'admits pure gain (Nb=1.273) and excludes nc_amp +6 dB '
                                 '(Nb=4.96)'},
            'arrival_drift': {'value': None,
                              'meaning': 'diagnostic only; not a hard constraint at v0.1 '
                                         '(0.832 ns reconstruction-axis quantisation)'},
        },
        'judgement_rules': {
            'implementation': 'scripts/study_tolerance_sensitivity.py (judge)',
            'ladder_pass': 'A<=tau_A and D<=tau_D and (Nb_ratio is None or Nb_ratio<=eps_Nb)',
            'operator_admissible': 'available and A<=tau_A and D_e<=tau_D and Nb_ratio<=eps_Nb',
            'hard_constraints_first': True,
            'scalar_weights': 'TBD (draft §4 step 5, after this freeze; not derived from '
                              'candidate outcomes, Cawley-Talbot clause)',
        },
        'evidence': {
            'damage_ladder_r1': {'path': LADDER.name, 'sha256': sha256(LADDER)},
            'operator_effects_r1': {'path': OPERATORS.name, 'sha256': sha256(OPERATORS)},
            'tolerance_sensitivity_r1': {'path': SENSITIVITY.name, 'sha256': sha256(SENSITIVITY)},
        },
        'scope_limits': [
            'values are properties of the t3 development families plus the stated '
            'dispersion material assumption; not field-data thresholds',
            'do not tighten beyond the material-assumption sensitivity '
            '(+/-3.3..6.6 dB) before the dispersion field-data study',
            'all peak-based readouts quantised to the 0.832 ns reconstruction axis',
            'test families {C5,C8} untouched; no solver runs; no training',
            'gain-compensated outputs (audit_before_gain) used for the D_e battery so '
            'the fixed shared gain is not counted as damage',
        ],
    }
    text = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('anchors verified | wrote', OUT.relative_to(ROOT))
    print('sha256', hashlib.sha256(text.encode('utf-8')).hexdigest())


if __name__ == '__main__':
    main()
