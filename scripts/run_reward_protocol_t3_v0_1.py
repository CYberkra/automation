"""First executable run of the frozen reward protocol on the t3 development data.

Reward draft §4 calibration route is complete as of the weights freeze
(reward_weights_contract_v0.2.json, user "确认这版" 2026-09-29). This script
demonstrates the protocol as an executable pipeline, driven ONLY by the two
frozen contracts (SHA-asserted), independently of the study scripts:

  archived t3 BG CO B-scans (official SFCW chain, Hann window)
    -> apply background-class candidates (catalogue v0.2, G1/BG: B0,B2,B3,
       B4,B5,B7,B8,B9 + out-of-catalogue local_mean reference)
    -> v0.2 battery (D_e/a/A on Fermat window; contrast vs NC floor; N_b)
    -> hard gates (reward_tolerance_contract_v0.1)
    -> R_bg (reward_weights_contract_v0.2 formula)
    -> per-family admissible ranking and selection

This is a protocol-mechanics demonstration on archived development data: the
expected output reproduces the user-confirmed weights v0.2 ranking. No solver
runs; no test-family data; r1/r2 byte-identity asserted.

Run: artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe
     scripts/run_reward_protocol_t3_v0_1.py
"""

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_operator_contract import (ConfigUnavailable, apply_configuration,
                                        local_mean_control)
from study_t3_damage_ladder import (MOTHERS, NC_WIN, load_bscan, fermat_times,
                                    event_mask)

TOL_C = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
W_C = ROOT / 'configs/research/reward_weights_contract_v0.2.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_reward_protocol_t3_first_run_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_reward_protocol_t3_first_run_r2.json'

TOL_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'
W_SHA256 = '0427b6ed58328e5ae01d2e0bfbda90460ef1a68beea318808a3f2ae2ed69bfb0'

CANDIDATES = ['B0_G1_BG', 'B2_G1_BG', 'B3_G1_BG', 'B4_G1_BG', 'B5_G1_BG',
              'B7_G1_BG', 'B8_G1_BG', 'B9_G1_BG']
LOCAL_MEAN_WIDTH = 11


def load_contracts():
    assert hashlib.sha256(TOL_C.read_bytes()).hexdigest() == TOL_SHA256
    assert hashlib.sha256(W_C.read_bytes()).hexdigest() == W_SHA256
    tol = {k: v['value'] for k, v in
           json.loads(TOL_C.read_text(encoding='utf-8'))['tolerances'].items()
           if k in ('tau_A', 'tau_D', 'eps_Nb')}
    formula = json.loads(W_C.read_text(encoding='utf-8'))['reward_background_class']['formula']
    assert formula == ('R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) '
                       '- max(10*log10(Nb_ratio), 0)')
    return tol, formula


def r_bg(delta, d_e, nb):
    pres = -20.0 * math.log10(1.0 - d_e) if d_e < 1.0 else float('inf')
    return delta - pres - max(10.0 * math.log10(nb), 0.0), pres


def run_family(fam, mother, geo, tol):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev).T
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s.shape).T
    x = np.ascontiguousarray(s.T)
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
    outputs['local_mean_w11(out-of-catalogue)'] = local_mean_control(x, LOCAL_MEAN_WIDTH)
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
                     'feasible': feasible, 'R_bg_db': (round(score, 6) if feasible else None),
                     'preservation_db': (round(pres, 6) if feasible else None)})
    ranking = sorted((r for r in rows if r.get('feasible')),
                     key=lambda r: -r['R_bg_db'])
    excluded = [{'config': r['config'], 'reason': r.get('reason') or
                 'gate:' + ','.join(k for k, v in r['gates'].items() if not v)}
                for r in rows if not r.get('feasible')]
    return {'family': fam, 'mother': mother, 'rows': rows,
            'selection': (ranking[0]['config'] if ranking else None),
            'selection_R_bg_db': (ranking[0]['R_bg_db'] if ranking else None),
            'excluded': excluded, 'ranking': [r['config'] for r in ranking]}


def build():
    tol, _formula = load_contracts()
    return {'t3_families': [run_family(f, m, g, tol) for f, m, g in MOTHERS]}


def main():
    tol, formula = load_contracts()
    doc = {
        'schema': 'reward_protocol_first_run/0.1',
        'date': '2026-09-29',
        'basis': 'reward draft §4 complete: tolerances frozen (v0.1, user '
                 '"确认" 2026-09-28) + weights frozen (v0.2, user "确认这版" '
                 '2026-09-29); this run executes the protocol end-to-end from '
                 'archived t3 inputs through both frozen contracts',
        'contracts': {'tolerance': {'path': TOL_C.name, 'sha256': TOL_SHA256},
                      'weights': {'path': W_C.name, 'sha256': W_SHA256,
                                  'formula_applied': formula}},
        'candidates': CANDIDATES + ['local_mean_w11(out-of-catalogue)'],
        'note': 'protocol-mechanics demonstration on archived development '
                'data; expected to reproduce the confirmed weights v0.2 '
                'ranking; no solver runs; no test-family data',
        'results': build()['t3_families'],
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = dict(doc, results=build()['t3_families'])
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    for fam in doc['results']:
        print(f"{fam['family']:6s} selection: {fam['selection']} "
              f"R={fam['selection_R_bg_db']} | excluded: "
              + ', '.join(e['config'] for e in fam['excluded']))


if __name__ == '__main__':
    main()
