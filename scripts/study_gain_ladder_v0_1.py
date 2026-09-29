"""Gain-class damage ladder v0.1: known multiplicative attenuation controls.

Imposes the frozen gain_control_cases_v0.1 contract (SHA-asserted) on the t3
development families: five attenuation levels x three structures x two NC
states, all evaluated at IDENTITY gain (the attenuated radargram itself), with
the ORIGINAL unattenuated event as reference (constructed_reference; opposite
direction of the background-class ladder).

Identity-anchor self-consistency (contract §identity_anchor_expectations):
  event_only      -> measured a equals the linear factor exactly
  nc_untouched    -> NC energy ratio vs original equals 1.0
  nc_attenuated   -> NC energy ratio vs original equals k^2
  control (0 dB)  -> a=1, D_e=0
Any anchor failure aborts the run (the ladder scale is invalid otherwise).

Deterministic; r1/r2 byte-identity asserted. No solver; no test families; the
clip tolerance is NOT derived here (effects table, design step 3-4).
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from study_t3_damage_ladder import (MOTHERS, NC_WIN, load_bscan, fermat_times,
                                    event_mask)

CONTRACT = ROOT / 'configs/research/gain_control_cases_v0.1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_gain_ladder_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_gain_ladder_r2.json'
CONTRACT_SHA256 = 'c027798904a4506bf3bd74d12fef00a4dc3df94b3cf6a938af170960b6689a18'

BAND_EXT_NS = 20.0  # raised-cosine taper extent beyond the event window


def load_contract():
    assert hashlib.sha256(CONTRACT.read_bytes()).hexdigest() == CONTRACT_SHA256
    c = json.loads(CONTRACT.read_text(encoding='utf-8'))
    levels = [float(x) for x in c['attenuation_levels_db']]
    return levels, c['structures'].keys(), c['nc_states'].keys()


def factor_map(t, t_ev, m_ev, k, structure, nc_state):
    """Multiplicative factor map [trace, sample]; k = linear attenuation."""
    n_tr, n_s = m_ev.shape
    f = np.ones((n_tr, n_s))
    if structure == 'event_only':
        f[m_ev] = k
    elif structure == 'event_band':
        f[m_ev] = k
        for i in range(n_tr):
            lo = (t >= t_ev[i] - 10.0 - BAND_EXT_NS) & (t < t_ev[i] - 10.0)
            hi = (t > t_ev[i] + 10.0) & (t <= t_ev[i] + 10.0 + BAND_EXT_NS)
            d_lo = (t_ev[i] - 10.0) - t[lo]
            d_hi = t[hi] - (t_ev[i] + 10.0)
            f[i, lo] = 1.0 + (k - 1.0) * 0.5 * (1.0 + np.cos(np.pi * d_lo / BAND_EXT_NS))
            f[i, hi] = 1.0 + (k - 1.0) * 0.5 * (1.0 + np.cos(np.pi * d_hi / BAND_EXT_NS))
    elif structure == 'time_global':
        t_start = float(np.min(t_ev) - 10.0)
        f[:, t >= t_start] = k
    else:
        raise ValueError(structure)
    # NC dual state: NC must carry factor k (attenuated) or 1 (untouched).
    m_nc = (t >= NC_WIN[0]) & (t <= NC_WIN[1])
    if nc_state == 'nc_attenuated':
        f[:, m_nc] = k
    elif nc_state == 'nc_untouched':
        f[:, m_nc] = 1.0
    else:
        raise ValueError(nc_state)
    return f


def battery(z, s, m_ev):
    zs, ss = z[m_ev], s[m_ev]
    denom = float(np.linalg.norm(ss))
    a = float(np.dot(zs, ss) / np.dot(ss, ss))
    zn = float(np.linalg.norm(zs))
    rho = float(np.dot(zs, ss) / (zn * denom)) if zn > 0 else None
    return {'a': round(a, 9),
            'D_e': round(float(np.linalg.norm(zs - ss) / denom), 6),
            'A': round(abs(a - 1.0), 6),
            'H': round(float(np.linalg.norm(zs - a * ss) / denom), 6),
            'rho': (round(rho, 6) if rho is not None else None)}


def run_family(fam, mother, geo, levels, structures, nc_states):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev)
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s.shape)
    peak0 = float(np.max(np.abs(s)))
    rows = []
    # Control: no attenuation, identity.
    b = battery(s, s, m_ev)
    assert b['a'] == 1.0 and b['D_e'] == 0.0, 'control anchor failed'
    rows.append({'family': fam, 'level_db': 0.0, 'structure': 'control',
                 'nc_state': 'nc_untouched', **b,
                 'nc_ratio': 1.0, 'peak_ratio': 1.0, 'anchor_pass': True})
    for db in levels:
        k = 10.0 ** (db / 20.0)
        for st in structures:
            for nc in nc_states:
                f = factor_map(t, t_ev, m_ev, k, st, nc)
                z = f * s
                b = battery(z, s, m_ev)
                # Contract identity anchors: assert on UNROUNDED quantities
                # (battery() rounds to 9 decimals; rounding error ~5e-10
                # exceeds the 1e-12 relative anchor tolerance).
                zs, ss = z[m_ev], s[m_ev]
                a_raw = float(np.dot(zs, ss) / np.dot(ss, ss))
                nc_ratio = float(np.mean(z[m_nc] ** 2) / np.mean(s[m_nc] ** 2))
                ok = True
                if st == 'event_only':
                    ok &= abs(a_raw - k) <= 1e-12 * max(k, 1e-12)
                if nc == 'nc_untouched':
                    ok &= abs(nc_ratio - 1.0) <= 1e-12
                else:
                    ok &= abs(nc_ratio - k ** 2) <= 1e-12 * max(k ** 2, 1e-12)
                assert ok, f'anchor failed: {fam} {db} {st} {nc}'
                rows.append({'family': fam, 'level_db': db, 'structure': st,
                             'nc_state': nc, **b, 'k': round(k, 9),
                             'nc_ratio': round(nc_ratio, 9),
                             'peak_ratio': round(float(np.max(np.abs(z)) / peak0), 9),
                             'anchor_pass': ok})
    print(fam, 'done:', len(rows), 'ladder rows (all identity anchors pass)')
    return rows


def build_records():
    levels, structures, nc_states = load_contract()
    records = []
    for fam, mother, geo in MOTHERS:
        records.extend(run_family(fam, mother, geo, levels, structures, nc_states))
    return records


def main():
    doc = {
        'schema': 'gain_ladder/0.1',
        'date': '2026-09-29',
        'basis': 'gain class design (user "批准" 2026-09-29) step 1; contract '
                 'gain_control_cases_v0.1.json (SHA asserted) drives every '
                 'pre-registered value',
        'gain_applied': 'identity (the attenuated radargram itself); recovery '
                        'behaviour of gain operators is the effects-table step '
                        '(step 3), not this ladder',
        'reference': 'original unattenuated event (constructed_reference)',
        'records': build_records(),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    rec2 = build_records()
    doc2 = dict(doc, records=rec2)
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('records:', len(doc['records']),
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
