"""Freeze the reward weights v0.2 (reward draft §4 step 5, user-confirmed).

User confirmed the v0.2 ranking ("确认这版，继续推进吧", 2026-09-29) after the
RPCA addendum and catalogue v0.2 freeze (d9fa423). Freezes the scalar reward
formula and coefficient-1 derivation into an independent contract; the
tolerance contract (reward_tolerance_contract_v0.1.json) stays as-is and is
referenced, not modified.

Freeze-script assertion discipline: every anchor claim is asserted against
committed evidence files before writing; build runs twice in-process and must
be byte-identical.
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/reward_weights_contract_v0.2.json'
TOL = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
CAT2 = ROOT / 'configs/research/operator_catalogue_v0.2.json'
W_R1 = ROOT / 'artifacts/research_checks/2026-09-29_reward_weights_v0_2_r1.json'
W_R2 = ROOT / 'artifacts/research_checks/2026-09-29_reward_weights_v0_2_r2.json'
RPCA = ROOT / 'artifacts/research_checks/2026-09-29_t3_operator_effects_rpca_r1.json'
DOC = ROOT / 'docs/research/2026-09-29_reward_weights_v0.2.md'

TOL_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build():
    # Anchor: tolerance contract untouched and hash-locked.
    assert sha256(TOL) == TOL_SHA256, 'frozen tolerance contract hash mismatch'
    # Anchor: r1/r2 byte-identical and schema/lineage correct.
    t1, t2 = W_R1.read_text(encoding='utf-8'), W_R2.read_text(encoding='utf-8')
    assert t1 == t2, 'weights r1/r2 mismatch'
    w = json.loads(t1)
    assert w['schema'] == 'reward_weights/0.2'
    assert w['tolerances_frozen'] == {'tau_A': 0.2, 'tau_D': 0.95, 'eps_Nb': 1.5}
    rows = w['rows']
    assert len(rows) == 13 * 3
    # Anchor: top-1 claims per family (the ranking the user confirmed).
    def top1(fam):
        feas = [r for r in rows if r['family'] == fam and r['feasible']]
        return max(feas, key=lambda r: r['R_bg_db'])
    assert top1('S2X')['config'] == 'B9_G1_BG' and abs(top1('S2X')['R_bg_db'] - 27.02) < 0.01
    assert top1('S2TZX')['config'] == 'B9_G1_BG' and abs(top1('S2TZX')['R_bg_db'] - 25.31) < 0.01
    c3 = [r for r in rows if r['family'] == 'C3mX' and r['feasible']]
    assert sorted(c3, key=lambda r: -r['R_bg_db'])[0]['config'] == 'B7_G1_BG'
    # Anchor: sensitivity structure (8 inversions all local_mean <-> mean/svd).
    axes = {(a['coefficient'], a['factor']): a['n_inversions'] for a in w['sensitivity']}
    assert axes == {('c_d', 0.5): 0, ('c_d', 1.5): 8, ('c_b', 0.5): 0, ('c_b', 1.5): 0}
    # Anchor: catalogue v0.2 carries the RPCA ids the ranking reports.
    cat = json.loads(CAT2.read_text(encoding='utf-8'))
    ids = {c['id'] for c in cat['configurations']}
    assert {'B7_G1_BG', 'B8_G1_BG', 'B9_G1_BG'} <= ids and cat['candidate_count'] == 30

    doc = {
        'schema': 'reward-weights-contract/0.2',
        'date': '2026-09-29',
        'status': 'frozen',
        'frozen_by': 'user "确认这版，继续推进吧" (2026-09-29)',
        'reward_background_class': {
            'formula': 'R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) '
                       '- max(10*log10(Nb_ratio), 0)',
            'coefficients': {'w_c': 1, 'w_d': 1, 'w_b': 1},
            'derivation': 'coefficients = 1 by dB-unit identities (equal dB = '
                          'equal scale; net event-to-interference improvement); '
                          'NOT fitted to candidate outcomes (Cawley-Talbot)',
            'domain': 'feasible rows only (hard gates of the tolerance contract '
                      'apply first); D_e < 1 required by the log identity',
            'inputs': 'per-family effects rows: contrast_db_delta, D_e, '
                      'Nb_ratio from the v0.2 battery (event window = Fermat '
                      'line +/-10 ns; NC window 250-400 ns)',
        },
        'tolerance_contract': {'path': TOL.name, 'sha256': TOL_SHA256,
                               'relation': 'referenced, not modified'},
        'catalogue': {'path': CAT2.name, 'sha256': sha256(CAT2),
                      'note': 'RPCA ids B7/B8/B9_G1_BG reported by this ranking'},
        'confirmed_ranking_top': {
            'S2X': 'B9_G1_BG +27.02 dB',
            'S2TZX': 'B9_G1_BG +25.31 dB',
            'C3mX': 'B7_G1_BG +30.02 dB (only non-identity feasible)',
        },
        'scope_limits': [
            'values are properties of the t3 development families plus the '
            'dispersion material assumption; not field-data thresholds',
            'background-suppression class only; gain-class quality labels '
            'remain undetermined (no controlled multiplicative-attenuation '
            'cases yet)',
            'new scenes require their own reference/window construction and '
            'recalibration of tolerances before the reward is applied',
            'test families {C5,C8} untouched; no solver runs; no training',
        ],
        'evidence': {
            'reward_weights_v0_2_r1': {'path': W_R1.name, 'sha256': sha256(W_R1)},
            'reward_weights_v0_2_r2': {'path': W_R2.name, 'sha256': sha256(W_R2)},
            'rpca_effects_r1': {'path': RPCA.name, 'sha256': sha256(RPCA)},
            'doc': {'path': 'docs/research/2026-09-29_reward_weights_v0.2.md'},
        },
    }
    return doc


def main():
    text = json.dumps(build(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    again = json.dumps(build(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text == again, 'two build runs disagree'
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('anchors verified | wrote', OUT.relative_to(ROOT))
    print('sha256', hashlib.sha256(text.encode('utf-8')).hexdigest())


if __name__ == '__main__':
    main()
