"""Freeze the operator catalogue v0.2 (RPCA lambda axis joins the catalogue).

v0.1 (27 entries) is preserved untouched so every v0.1-era frozen hash and
check stays valid. v0.2 adds three RPCA background entries — lambda axis
0.5/1.0/2.0 x standard PCP lambda0 = 1/sqrt(max(sample,trace)), end_gain 1,
BG order only (ids B7/B8/B9_G1_BG) — on the evidence of the 2026-09-29 RPCA
effects addendum (user "按你的建议最佳实践"). RPCA damage-ladder coverage and
gain-combination rows are deferred (recorded in scope_limits).

Freeze-script assertion discipline: every anchor claim is asserted against
committed files before writing; build is run twice in-process and must be
byte-identical (rerun determinism), so the contract cannot drift from its
basis.

Usage: python scripts/freeze_operator_catalogue_v0_2.py
"""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_operator_contract import VERSION, VERSION_02, RPCA_LAM_FACTORS, catalogue

OUT = ROOT / 'configs/research/operator_catalogue_v0.2.json'
V1 = ROOT / 'configs/research/operator_catalogue_v0.1.json'
RPCA_R1 = ROOT / 'artifacts/research_checks/2026-09-29_t3_operator_effects_rpca_r1.json'
MODULE = ROOT / 'scripts/research_operator_contract.py'


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build():
    v1 = json.loads(V1.read_text(encoding='utf-8'))
    v2_configs = catalogue('0.2')
    v1_configs = catalogue('0.1')

    # Anchor assertions against committed files.
    assert v1['configurations'] == v1_configs, 'v0.1 file disagrees with module catalogue()'
    assert len(v1_configs) == 27 and v1['candidate_count'] == 27
    assert len(v2_configs) == 30 and len({c['id'] for c in v2_configs}) == 30
    by_id = {c['id']: c for c in v2_configs}
    # v0.1 subset preserved exactly (identical id/background/parameter/end_gain/order).
    for c in v1_configs:
        assert by_id[c['id']] == c, f'v0.1 entry drifted: {c["id"]}'
    # New RPCA axis: ids, order, lambda factors.
    rpca_ids = [f'B{7 + i}_G1_BG' for i in range(3)]
    for cid, f in zip(rpca_ids, RPCA_LAM_FACTORS):
        c = by_id[cid]
        assert c['background'] == 'rpca' and c['parameter'] == f
        assert c['end_gain'] == 1 and c['order'] == 'BG'
    # Evidence file records exactly these three lambda runs on the three t3 families.
    rec = json.loads(RPCA_R1.read_text(encoding='utf-8'))['records']
    assert len(rec) == 9
    for r in rec:
        assert r['background'] == 'rpca' and r['available']
    fams = {r['family'] for r in rec}
    assert fams == {'S2X', 'S2TZX', 'C3mX'}

    doc = {
        'schema': 'research-operator-catalogue/0.2',
        'operator_version': VERSION_02,
        'status': 'array_contract_checked_not_physics_validated',
        'supersedes': 'operator_catalogue_v0.1.json (v0.1 preserved; default '
                      'catalogue() remains v0.1 for backward compatibility)',
        'date': '2026-09-29',
        'frozen_by': 'user "按你的建议最佳实践" (2026-09-29), on '
                     'docs/research/2026-09-29_rpca_effects_addendum.md',
        'fdtd_execution_enabled': False,
        'training_enabled': False,
        'domain': 'time_real',
        'axes': ['sample', 'trace'],
        'input_window': 'entire_predeclared_window',
        'gain_coordinate': 'u=i/(n_samples-1)',
        'gap_rtol': 1e-08,
        'numerical_rank_floor': 'max(n_samples,n_traces)*float64_eps',
        'candidate_count': 30,
        'legal_count_is_input_dependent': True,
        'local_mean_control_in_catalogue': False,
        'rpca_lambda_axis': {
            'lambda0': '1/sqrt(max(n_samples,n_traces)) (standard PCP choice)',
            'factors': list(RPCA_LAM_FACTORS),
            'solver': 'rpca_control: inexact ALM principal component pursuit '
                      '(Lin-Chen-Ma 2009), deterministic float64, max_iter=500, tol=1e-7',
            'measured_basis': 't3 operator effects: 2.0x tops slope-family R_bg '
                              '(+27.0/+25.3 dB), 0.5x is the only admissible background '
                              'removal on the flat family (D_e=0.014); no single factor '
                              'wins both families',
        },
        'configurations': v2_configs,
        'evidence': {
            'catalogue_v0.1': {'path': V1.name, 'sha256': sha256(V1)},
            'rpca_effects_r1': {'path': RPCA_R1.name, 'sha256': sha256(RPCA_R1)},
            'module': {'path': 'scripts/research_operator_contract.py',
                       'sha256': sha256(MODULE)},
        },
        'scope_limits': [
            'RPCA entries have no damage-ladder coverage yet (deferred); tolerance '
            'contract anchors do not depend on operators, Cawley-Talbot clause intact',
            'RPCA x end_gain>1 and GB order combinations not in catalogue v0.2',
            'RPCA compute cost is orders of magnitude above mean/SVD; cost accounting '
            'is diagnostic-only in the 2026-09-29 addendum, not a reward cost term',
            'same validity class as v0.1: array contract checked, not physics validated',
        ],
    }
    return doc


def main():
    doc = build()
    text = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    again = json.dumps(build(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text == again, 'two build runs disagree'
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('anchors verified | wrote', OUT.relative_to(ROOT))
    print('sha256', hashlib.sha256(text.encode('utf-8')).hexdigest())


if __name__ == '__main__':
    main()
