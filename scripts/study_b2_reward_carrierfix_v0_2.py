"""Dev-side B2 reward-protocol recompute under the official carrier (v0.2).

Re-scores the 2026-10-01 B2 pilot baseline (8 catalogue configs + out-of-
catalogue local_mean control, tolerance v0.2 / weights v0.3, SHA-asserted)
under both signed representations from ONE load per trace. Frozen pilot
outputs are NOT overwritten; new files only. The tolerance/weights
contracts are used as-is (frozen); this run measures how much the carrier
fix moves scores, gates, feasibility and selections — it does NOT
recalibrate anything.

No solver runs; no test-family data.
"""

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import sfcw_official_loader_v0_2 as L
from research_operator_contract import (ConfigUnavailable, apply_configuration,
                                        local_mean_control)
from study_t3_damage_ladder import (NC_WIN, fermat_times, event_mask)
from run_reward_protocol_b2_pilot_v0_1 import (
    CANDIDATES, LOCAL_MEAN_WIDTH, MOTHERS, interface_z, r_bg, load_contracts)

ARCH_R1 = ROOT / 'artifacts/research_checks/2026-10-01_reward_protocol_b2_pilot_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-10-02_b2_reward_carrierfix_v0_2_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-10-02_b2_reward_carrierfix_v0_2_r2.json'


def run_family_rep(fam, mother, geo, date, tol, x, t):
    """Score all candidates on one representation. x: (nt, ntraces)."""
    ifz = interface_z(geo)
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev).T
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), x.T.shape).T
    nb_in = float(np.mean(x[m_nc] ** 2))
    c_in = 20.0 * math.log10(float(np.max(np.abs(x[m_ev]))) / math.sqrt(nb_in))
    rows = []
    outputs = {}
    for cid in CANDIDATES:
        try:
            r = apply_configuration(x, cid, catalogue_version='0.2')
        except ConfigUnavailable as exc:
            rows.append({'config': cid, 'available': False, 'reason': str(exc)})
            continue
        outputs[cid] = r['output']
    outputs['local_mean_w11(out-of-catalogue)'] = local_mean_control(
        x, LOCAL_MEAN_WIDTH)
    for cid, y in outputs.items():
        zs, ss = y[m_ev], x[m_ev]
        denom = float(np.linalg.norm(ss))
        a = float(np.dot(zs, ss) / np.dot(ss, ss))
        d_e = float(np.linalg.norm(zs - ss) / denom)
        A = abs(a - 1.0)
        nb = float(np.mean(y[m_nc] ** 2)) / nb_in
        c_out = 20.0 * math.log10(float(np.max(np.abs(y[m_ev]))) /
                                  math.sqrt(float(np.mean(y[m_nc] ** 2))))
        delta = c_out - c_in
        gates = {'A': A <= tol['tau_A'], 'D_e': d_e <= tol['tau_D'],
                 'Nb': nb <= tol['eps_Nb']}
        feasible = all(gates.values())
        score, pres = (r_bg(delta, d_e, nb) if feasible else (None, None))
        rows.append({'config': cid, 'available': True, 'D_e': round(d_e, 6),
                     'a': round(a, 6), 'A': round(A, 6),
                     'contrast_db_delta': round(delta, 3),
                     'Nb_ratio': round(nb, 6), 'gates': gates,
                     'feasible': feasible,
                     'R_bg_db': (round(score, 6) if feasible else None)})
    ranking = sorted((r for r in rows if r.get('feasible')),
                     key=lambda r: -r['R_bg_db'])
    return {'rows': rows,
            'selection': (ranking[0]['config'] if ranking else None),
            'selection_R_bg_db': (ranking[0]['R_bg_db'] if ranking else None),
            'ranking': [r['config'] for r in ranking]}


def build(tol):
    fams = []
    for fam, mother, geo, date in MOTHERS:
        t, s_off, s_leg = L.load_bscan_both(mother, date=date)
        reps = {}
        for rep, s in (('legacy95', s_leg), ('official20', s_off)):
            x = np.ascontiguousarray(s.T)
            reps[rep] = run_family_rep(fam, mother, geo, date, tol, x, t)
        # paired deltas on shared config rows
        by_cfg = {}
        for rep, res in reps.items():
            for row in res['rows']:
                if row.get('available'):
                    by_cfg.setdefault(row['config'], {})[rep] = row
        deltas = []
        for cfg, pair in sorted(by_cfg.items()):
            if len(pair) < 2:
                continue
            a, b = pair['legacy95'], pair['official20']
            deltas.append({
                'config': cfg,
                'D_e_legacy': a['D_e'], 'D_e_official': b['D_e'],
                'D_e_delta': round(b['D_e'] - a['D_e'], 6),
                'contrast_legacy': a['contrast_db_delta'],
                'contrast_official': b['contrast_db_delta'],
                'contrast_delta_db': round(
                    b['contrast_db_delta'] - a['contrast_db_delta'], 3),
                'Nb_legacy': a['Nb_ratio'], 'Nb_official': b['Nb_ratio'],
                'feasible_legacy': a['feasible'],
                'feasible_official': b['feasible'],
                'feasibility_flip': a['feasible'] != b['feasible'],
                'R_bg_legacy': a['R_bg_db'], 'R_bg_official': b['R_bg_db'],
            })
        fams.append({'family': fam, 'mother': mother, 'archive_date': date,
                     'representations': reps, 'paired_delta': deltas,
                     'selection_flip': reps['legacy95']['selection']
                     != reps['official20']['selection']})
        print(fam, 'done; selection legacy/official:',
              reps['legacy95']['selection'], '/',
              reps['official20']['selection'],
              '| flips:', sum(d['feasibility_flip'] for d in deltas))
    return fams


def main():
    tol, formula = load_contracts()
    doc = {
        'schema': 'reward_protocol_b2_carrierfix/1',
        'date': '2026-10-02',
        'basis': 'model design review 2026-10-02 §2; dev-side recompute of '
                 'the 2026-10-01 B2 pilot baseline under both carriers; '
                 'frozen pilot outputs and contracts untouched',
        'tolerance_contract': 'v0.2 (SHA-asserted via pilot loader)',
        'weights_formula': formula,
        'archived_pilot_r1_sha256': hashlib.sha256(
            ARCH_R1.read_bytes()).hexdigest(),
        'note': 'global amplitude (2x) cancels in every ratio/dB quantity; '
                'only waveform shape (carrier) can move scores. No solver '
                'runs; no test-family data.',
        'results': build(tol),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    doc2 = dict(doc, results=build(tol))
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
