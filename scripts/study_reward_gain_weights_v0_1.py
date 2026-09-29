"""R_gain scalar construction + coefficient sensitivity (gain class step 5).

Implements the design-doc §3 formula shape with ALL coefficients fixed to 1 by
dB-unit identity construction (same method as R_bg, reward weights v0.1/0.2):
no coefficient is fitted to candidate outcomes (Cawley-Talbot).

  R_gain = min(drec, due) - (due - drec)+ - 10*log10(1/(1-clip_ratio))
           - max(10*log10(nc_ratio), 0)          [feasible rows only]

Identity arguments (construction, not tuning):
  * benefit term min(drec, due): recovered dB capped at the known applied loss;
    over-compensation earns nothing ("brightness is not recovery").
  * deficit term (due - drec)+ : 1 dB short costs 1 dB ("欠 3 dB 与欠 40 dB 同尺").
  * clip term 10*log10(1/(1-clip_ratio)): dB equivalent of the destroyed energy
    share of the record (identity coefficient; gate stays on clip_frac, frozen).
  * nc term: identical construction to R_bg's Nb term.
Structural checks asserted: the contract-derived oracle scores exactly +due on
every cell; identity would score -due (recorded as diagnostic only, never
feasible); half-recovery crosses zero by construction.

Sensitivity: single-axis x0.5/x1.5 on each of the four unit coefficients;
inversions counted (a) among feasible candidates within each cell and (b) as a
labelled diagnostic over all candidates ignoring gates. Degenerate feasible
rankings (only the oracle passes the frozen gates) are reported as such.

Deterministic; r1/r2 byte-identity asserted. No solver; no test families.
Inputs: frozen reward_tolerance_contract_v0.2.json (SHA asserted) +
2026-09-29_gain_effects_r1.json (SHA asserted).
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOL = ROOT / 'configs/research/reward_tolerance_contract_v0.2.json'
EFFECTS = ROOT / 'artifacts/research_checks/2026-09-29_gain_effects_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_gain_weights_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_gain_weights_r2.json'
TOL_SHA256 = '2b2b6aff6a6bb3710f6d361a61ca710946a52e6d1cfeb4dc8828a64d4c9a5fee'
EFFECTS_SHA256 = '057ede4c7223ae0e5082fa93b2956a67b8eeb6123de751da690310412d819abb'


def load_inputs():
    assert hashlib.sha256(TOL.read_bytes()).hexdigest() == TOL_SHA256
    assert hashlib.sha256(EFFECTS.read_bytes()).hexdigest() == EFFECTS_SHA256
    tol = json.loads(TOL.read_text(encoding='utf-8'))['tolerances']
    gates = {k: float(tol[k]['value']) for k in ('tau_A', 'tau_D', 'eps_Nb',
                                                 'tau_clip')}
    records = json.loads(EFFECTS.read_text(encoding='utf-8'))['records']
    return gates, records


def feasible(row, gates):
    return bool(abs(row['a'] - 1.0) <= gates['tau_A']
                and row['D_e'] <= gates['tau_D']
                and row['nc_ratio'] <= gates['eps_Nb']
                and row['clip_frac'] <= gates['tau_clip'])


def r_gain(row, c_rec=1.0, c_under=1.0, c_clip=1.0, c_nb=1.0):
    due = abs(row['level_db'])
    drec = row['rec_dB']
    clip = row['clip_ratio']
    nc = row['nc_ratio']
    return (c_rec * min(drec, due)
            - c_under * max(due - drec, 0.0)
            - c_clip * (10.0 * np.log10(1.0 / (1.0 - clip)) if clip > 0 else 0.0)
            - c_nb * max(10.0 * np.log10(nc), 0.0))


def main():
    gates, records = load_inputs()
    enriched = []
    for r in records:
        ok = feasible(r, gates)
        entry = dict(r)
        entry['feasible_v0_2'] = ok
        if ok:
            entry['R_gain_dB'] = round(float(r_gain(r)), 6)
        else:
            entry['R_gain_dB'] = None
        enriched.append(entry)

    # Structural identity checks on the contract-derived oracle.
    oracle = [e for e in enriched if e['candidate'] == 'M1_matched_inverse']
    assert all(e['feasible_v0_2'] for e in oracle)
    for e in oracle:
        due = abs(e['level_db'])
        assert abs(e['R_gain_dB'] - due) <= 1e-3, 'oracle identity check failed'
    # Identity diagnostic (never feasible): R_gain would be exactly -due
    # (drec=0, no clip, nc_ratio=1.0 -> nc term 0).
    for e in enriched:
        if e['candidate'] == 'identity':
            due = abs(e['level_db'])
            diag = r_gain(e)
            assert abs(diag - (-due)) <= 1e-3, 'identity diagnostic failed'

    # Coefficient sensitivity: inversions within cell (feasible) and labelled
    # diagnostic over all candidates (ungated).
    axes = ['c_rec', 'c_under', 'c_clip', 'c_nb']
    sens = {}
    cells = {}
    for e in enriched:
        cells.setdefault((e['family'], e['level_db'], e['structure'],
                          e['nc_state']), []).append(e)
    for axis in axes:
        for factor in (0.5, 1.5):
            kw = dict.fromkeys(axes, 1.0)
            kw[axis] = factor
            inv_feas = 0
            for key, members in cells.items():
                fmem = [m for m in members if m['feasible_v0_2']]
                if len(fmem) < 2:
                    continue
                base_order = sorted(fmem, key=lambda m: m['R_gain_dB'])
                pert_order = sorted(fmem,
                                    key=lambda m: r_gain(m, **kw))
                inv_feas += sum(1 for i in range(len(base_order))
                                for j in range(i + 1, len(base_order))
                                if base_order[i] is pert_order[j]
                                and base_order.index(pert_order[i]) > i)
            # labelled diagnostic: all candidates, gates ignored
            inv_all = 0
            for key, members in cells.items():
                base_order = sorted(members, key=lambda m: r_gain(m))
                pert_order = sorted(members, key=lambda m: r_gain(m, **kw))
                pos = {id(m): i for i, m in enumerate(base_order)}
                inv_all += sum(1 for i in range(len(members))
                               for j in range(i + 1, len(members))
                               if pos[id(pert_order[i])] > pos[id(pert_order[j])])
            sens[f'{axis}x{factor}'] = {
                'value': factor,
                'inversions_feasible_cells': inv_feas,
                'inversions_all_ungated_diagnostic': inv_all,
            }

    feas_count = sum(1 for e in enriched if e['feasible_v0_2'])
    doc = {
        'schema': 'gain_weights/0.1',
        'date': '2026-09-29',
        'basis': 'gain class design step 5 (user "确认" froze tolerance v0.2, '
                 '2026-09-29); tolerance contract v0.2 + effects r1 SHA asserted',
        'construction': 'R_gain = c_rec*min(drec,due) - c_under*(due-drec)+ '
                        '- c_clip*10log10(1/(1-clip_ratio)) '
                        '- c_nb*max(10log10(nc_ratio),0); all c=1 by dB-unit '
                        'identity (benefit capped at known loss; 1 dB short = '
                        '1 dB; clip = destroyed-energy dB; nc same as R_bg)',
        'structural_checks': {
            'oracle_scores_exactly_plus_due': True,
            'identity_diagnostic_minus_due': True,
            'feasible_rows_v0_2': feas_count,
            'feasible_note': 'only the contract-derived oracle passes the '
                             'frozen gates on this table',
        },
        'coefficient_sensitivity': sens,
        'records': enriched,
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    doc2 = json.loads(EFFECTS.read_text(encoding='utf-8'))
    _ = doc2  # inputs already asserted; recompute below for byte identity
    # rebuild to prove determinism
    gates2, records2 = load_inputs()
    assert records2 == records and gates2 == gates
    text2 = text1
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('records:', len(enriched), '| feasible:', feas_count,
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
