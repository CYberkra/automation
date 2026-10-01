"""Reward-protocol baseline scoring run on the B2 pilot batch (2026-10-01).

Mirrors scripts/run_reward_protocol_t3_v0_1.py (2026-09-29 first run) with
three deliberate changes:

  1. families: adds the new flat 1 m-cover family C1mX (BG archived
     2026-10-01 in the B2 pilot batch) alongside the three t3 calibration
     families S2X / S2TZX / C3mX (BG archived 2026-09-28), so the output is a
     full side-by-side baseline table;
  2. contracts: tolerance contract v0.2 (supersedes v0.1, background class
     unchanged) and weights contract v0.3, whose supersedes field declares the
     v0.2 background-class section unchanged - the R_bg formula is therefore
     taken from weights v0.2 and both SHAs are asserted;
  3. load_bscan takes an explicit archive-date parameter (the t3 loader
     hardcodes 2026-09-28).

Caveat recorded for the review: tolerances were calibrated on the t3
families (S2X/S2TZX/C3mX). C1mX is a new development-family geometry; its
scores here are a diagnostic baseline under t3-frozen tolerances, pending
its own calibration check (prerequisite ledger P1-2). The four B-layer probe
mothers (C1p5mS1X/C1p5mS3X/C2mS3X/C2mS3TZX) are NOT scored here: per the
task-definition contract scope limits, new scenarios need their own
reference windows and tolerance recalibration before the frozen reward
applies.

No solver runs; no test-family data; r1/r2 byte-identity asserted.

Run: artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe
     scripts/run_reward_protocol_b2_pilot_v0_1.py
"""

import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_operator_contract import (ConfigUnavailable, apply_configuration,
                                        local_mean_control)
from study_t3_damage_ladder import (FREQ, FC, NC_WIN, N_TRACES, fermat_times,
                                    event_mask)
from gprMax.toolboxes.SFCW.processing import (
    load_source, load_receiver, direct_frequency_response,
    reconstruct_time_response)

SIMS = ROOT / 'artifacts/simulations'
TOL_C = ROOT / 'configs/research/reward_tolerance_contract_v0.2.json'
W2_C = ROOT / 'configs/research/reward_weights_contract_v0.2.json'
W3_C = ROOT / 'configs/research/reward_weights_contract_v0.3.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-10-01_reward_protocol_b2_pilot_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-10-01_reward_protocol_b2_pilot_r2.json'

TOL_SHA256 = '2b2b6aff6a6bb3710f6d361a61ca710946a52e6d1cfeb4dc8828a64d4c9a5fee'
W2_SHA256 = '0427b6ed58328e5ae01d2e0bfbda90460ef1a68beea318808a3f2ae2ed69bfb0'
W3_SHA256 = '282d083b34477bafff25473f853e3e8d47c6c44807c82e3c5f9d41170399e0d0'

CANDIDATES = ['B0_G1_BG', 'B2_G1_BG', 'B3_G1_BG', 'B4_G1_BG', 'B5_G1_BG',
              'B7_G1_BG', 'B8_G1_BG', 'B9_G1_BG']
LOCAL_MEAN_WIDTH = 11

# (family label, mother model id, geometry, archive date)
MOTHERS = [
    ('C1mX', 'B2D-C1mX-BG', 'flat29', '2026-10-01'),
    ('S2X', 'B2D-C3mS2X-BG', 'slope', '2026-09-28'),
    ('S2TZX', 'B2D-C3mS2TZX-BG', 'slope', '2026-09-28'),
    ('C3mX', 'B2D-C3mX-BG', 'flat27', '2026-09-28'),
]


def load_bscan(mother, date):
    """Official SFCW chain, identical to study_t3_damage_ladder.load_bscan
    except the archive date is a parameter."""
    sigs, t = [], None
    for k in range(N_TRACES):
        rid = f'{mother}-CO33-t{k + 1:02d}'
        h5 = SIMS / f'{date}_{rid}' / f'{rid}.h5'
        with h5py.File(h5, 'r') as h:
            dt = float(h.attrs['dt'])
            items = list(h['rxs'].items())
            assert len(items) == 1
            rx_name = items[0][1].attrs['Name']
            raw = items[0][1]['Ex'][:]
        src = load_source(h5)
        rxt = load_receiver(h5, receiver_path='name:' + rx_name, component='Ex')
        n = min(len(raw), int(np.floor(1200e-9 / dt)) + 1)
        taper = (round(200.0e-9 / dt) - .25) / n
        rx = replace(rxt, samples=raw[:n])
        r = direct_frequency_response(src, rx, FREQ, tail_taper_fraction=taper)
        tr = reconstruct_time_response(r, zero_pad_factor=8, window='hann')
        tt = np.asarray(tr.time, dtype=float) * 1e9
        env = np.abs(np.asarray(tr.complex_envelope, dtype=np.complex128))
        sig = env * np.cos(2 * np.pi * FC * tr.time + np.angle(tr.complex_envelope))
        sigs.append(sig)
        t = tt
    return np.stack(sigs), t


def load_contracts():
    assert hashlib.sha256(TOL_C.read_bytes()).hexdigest() == TOL_SHA256
    assert hashlib.sha256(W2_C.read_bytes()).hexdigest() == W2_SHA256
    assert hashlib.sha256(W3_C.read_bytes()).hexdigest() == W3_SHA256
    tol = {k: v['value'] for k, v in
           json.loads(TOL_C.read_text(encoding='utf-8'))['tolerances'].items()
           if k in ('tau_A', 'tau_D', 'eps_Nb')}
    w2 = json.loads(W2_C.read_text(encoding='utf-8'))
    formula = w2['reward_background_class']['formula']
    assert formula == ('R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) '
                       '- max(10*log10(Nb_ratio), 0)')
    w3 = json.loads(W3_C.read_text(encoding='utf-8'))
    assert w3['supersedes']['sha256'] == W2_SHA256, \
        'v0.3 must declare v0.2 background-class section unchanged'
    return tol, formula


def r_bg(delta, d_e, nb):
    pres = -20.0 * math.log10(1.0 - d_e) if d_e < 1.0 else float('inf')
    return delta - pres - max(10.0 * math.log10(nb), 0.0), pres


def interface_z(geo):
    if geo == 'slope':
        return lambda y: np.minimum(0.2 * y + 22.625, 30.0)
    if geo == 'flat27':
        return lambda y: np.zeros_like(np.asarray(y, float)) + 27.0
    if geo == 'flat29':
        return lambda y: np.zeros_like(np.asarray(y, float)) + 29.0
    raise ValueError(geo)


def run_family(fam, mother, geo, date, tol):
    ifz = interface_z(geo)
    s, t = load_bscan(mother, date)
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
    return {'family': fam, 'mother': mother, 'archive_date': date, 'rows': rows,
            'selection': (ranking[0]['config'] if ranking else None),
            'selection_R_bg_db': (ranking[0]['R_bg_db'] if ranking else None),
            'excluded': excluded, 'ranking': [r['config'] for r in ranking]}


def build(tol):
    return [run_family(f, m, g, d, tol) for f, m, g, d in MOTHERS]


def main():
    tol, formula = load_contracts()
    doc = {
        'schema': 'reward_protocol_b2_pilot/0.1',
        'date': '2026-10-01',
        'basis': 'prerequisite ledger P4-1: baseline scoring of the operator '
                 'catalogue on the B2 pilot batch; mirrors the 2026-09-29 t3 '
                 'first run with tolerance v0.2 + weights v0.3 (background '
                 'formula inherited unchanged from v0.2, SHA-asserted)',
        'contracts': {'tolerance': {'path': TOL_C.name, 'sha256': TOL_SHA256},
                      'weights': {'path': W3_C.name, 'sha256': W3_SHA256,
                                  'background_formula_from': W2_C.name,
                                  'background_formula_sha256': W2_SHA256,
                                  'formula_applied': formula}},
        'candidates': CANDIDATES + ['local_mean_w11(out-of-catalogue)'],
        'note': 'C1mX is a new development-family geometry: its scores are a '
                'diagnostic baseline under t3-frozen tolerances, pending its '
                'own calibration check (ledger P1-2). B-layer probe mothers '
                'not scored (task contract scope limits). No solver runs; no '
                'test-family data.',
        'results': build(tol),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = dict(doc, results=build(tol))
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
