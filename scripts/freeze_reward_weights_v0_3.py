"""Freeze the reward weights v0.3 (gain class, design step 5 closure).

User approved the R_gain construction ("你先推进着做吧", 2026-09-29) after the
v0.2 tolerance freeze (SHA 2b2b6aff...). Adds the gain-class reward to the
weights contract as an independent document; the v0.2 background-class section
and the tolerance contract are referenced, not modified.

Freeze-script assertion discipline: every anchor claim is asserted against
committed evidence files before writing; build runs twice in-process and must
be byte-identical.
"""

import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/reward_weights_contract_v0.3.json'
W2 = ROOT / 'configs/research/reward_weights_contract_v0.2.json'
TOL2 = ROOT / 'configs/research/reward_tolerance_contract_v0.2.json'
G_R1 = ROOT / 'artifacts/research_checks/2026-09-29_gain_weights_r1.json'
G_R2 = ROOT / 'artifacts/research_checks/2026-09-29_gain_weights_r2.json'
DOC = ROOT / 'docs/research/2026-09-29_gain_weights_proposal.md'

W2_SHA256 = '0427b6ed58328e5ae01d2e0bfbda90460ef1a68beea318808a3f2ae2ed69bfb0'
TOL2_SHA256 = '2b2b6aff6a6bb3710f6d361a61ca710946a52e6d1cfeb4dc8828a64d4c9a5fee'

# Sensitivity anchor: zero inversions inside feasible cells (degenerate-robust,
# honestly flagged); ungated diagnostic inversion counts per axis (recorded,
# NOT refit against).
SENS_UNGATED = {
    ('c_rec', 0.5): 209, ('c_rec', 1.5): 138,
    ('c_under', 0.5): 209, ('c_under', 1.5): 138,
    ('c_clip', 0.5): 234, ('c_clip', 1.5): 215,
    ('c_nb', 0.5): 328, ('c_nb', 1.5): 237,
}


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build():
    # Anchor: v0.2 contracts untouched and hash-locked.
    assert sha256(W2) == W2_SHA256, 'frozen weights v0.2 hash mismatch'
    assert sha256(TOL2) == TOL2_SHA256, 'frozen tolerance v0.2 hash mismatch'
    # Anchor: gain weights r1/r2 byte-identical and schema/lineage correct.
    g1, g2 = G_R1.read_text(encoding='utf-8'), G_R2.read_text(encoding='utf-8')
    assert g1 == g2, 'gain weights r1/r2 mismatch'
    g = json.loads(g1)
    assert g['schema'] == 'gain_weights/0.1'
    assert len(g['records']) == 864
    # Anchor: structural identities from the proposal.
    sc = g['structural_checks']
    assert sc['oracle_scores_exactly_plus_due'] is True
    assert sc['identity_diagnostic_minus_due'] is True
    feas = [r for r in g['records'] if r['feasible_v0_2']]
    assert len(feas) == 72 and sc['feasible_rows_v0_2'] == 72
    # Anchor: every feasible row is the contract-derived oracle and reads
    # exactly +|due| per attenuation level (dB identity, 5 levels x cells).
    assert all(r['candidate'].startswith('M1') for r in feas)
    by_level = defaultdict(set)
    for r in feas:
        by_level[r['level_db']].add(round(r['R_gain_dB'], 6))
    assert by_level == {-3.0: {3.0}, -6.0: {6.0}, -10.0: {10.0},
                        -20.0: {20.0}, -40.0: {40.0}}, 'oracle identity broken'
    # Anchor: sensitivity counts (feasible cells carry all discriminative power).
    sens = g['coefficient_sensitivity']
    assert set(sens) == {f'{c}x{f}' for c, f in SENS_UNGATED}
    for (coef, factor), n_ungated in SENS_UNGATED.items():
        a = sens[f'{coef}x{factor}']
        assert a['value'] == factor
        assert a['inversions_feasible_cells'] == 0
        assert a['inversions_all_ungated_diagnostic'] == n_ungated

    doc = {
        'schema': 'reward-weights-contract/0.3',
        'date': '2026-09-29',
        'status': 'frozen',
        'frozen_by': 'user "你先推进着做吧" (2026-09-29), gain class step 5',
        'supersedes': {'path': W2.name, 'sha256': W2_SHA256,
                       'relation': 'v0.2 background-class section unchanged; '
                                   'this contract adds the gain class'},
        'reward_gain_class': {
            'formula': 'R_gain = c_rec*min(drec,due) - c_under*(due-drec)+ '
                       '- c_clip*10log10(1/(1-clip_ratio)) '
                       '- c_nb*max(10log10(nc_ratio),0)',
            'coefficients': {'c_rec': 1, 'c_under': 1, 'c_clip': 1, 'c_nb': 1},
            'derivation': 'all coefficients = 1 by dB-unit identities (benefit '
                          'capped at known loss; 1 dB short = 1 dB penalty; '
                          'clip = destroyed-energy dB equivalent; nc term '
                          'isomorphic to R_bg); NOT fitted to candidate '
                          'outcomes (Cawley-Talbot)',
            'domain': 'feasible rows only (hard gates of tolerance contract '
                      'v0.2 apply first, incl. tau_clip=0.002 on clip_frac)',
            'identity_anchor': 'oracle reads exactly +|due| at every level '
                               '(3/6/10/20/40 dB); identity diagnostic = '
                               '-due; feasible set on the effects table = '
                               'oracle only (72/864)',
            'sensitivity_note': 'zero inversions inside feasible cells '
                                '(degenerate-robust, honestly flagged); '
                                'ungated diagnostic inversions 138-328 '
                                'recorded, not refit against',
        },
        'tolerance_contract': {'path': TOL2.name, 'sha256': TOL2_SHA256,
                               'relation': 'referenced, not modified'},
        'scope_limits': [
            'values are properties of the t3 gain effects table (3 families x '
            '24 cells x 12 candidates); not field-data thresholds',
            'gain-class quality labels require controlled multiplicative-'
            'attenuation cases (known-decay control examples) before the '
            'reward is applied beyond the contract-derived oracle',
            'new scenes require their own reference/window construction and '
            'recalibration of tolerances before the reward is applied',
            'test families {C5,C8} untouched; no solver runs; no training '
            '(G4 maintained)',
        ],
        'evidence': {
            'gain_weights_r1': {'path': G_R1.name, 'sha256': sha256(G_R1)},
            'gain_weights_r2': {'path': G_R2.name, 'sha256': sha256(G_R2)},
            'doc': {'path': 'docs/research/2026-09-29_gain_weights_proposal.md'},
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
