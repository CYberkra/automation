"""Clip-scale ladder + tolerance sensitivity for the gain class (design step 4).

Two tasks in one deterministic script:

1. METERED CLIP LADDER. The clip threshold may not be derived from candidate
   operator performance (Cawley-Talbot), so we build construction rungs in exact
   analogy to the damage ladder: on the FS-normalised original signature, amplify
   with a metered gain g(p) = 1/quantile(|x|, 1-p) so that EXACTLY the rung
   fraction p of samples exceeds full scale, then hard-clip at |y|<=1. Rung
   fractions are physical-language anchors (record share destroyed by clipping),
   not candidate properties. The measured battery (a/D_e on the event window)
   records how much metered clipping corrupts the event itself.

2. SENSITIVITY. Single-axis x0.5 / x1.5 perturbation of the reused frozen gates
   (tau_A, tau_D, eps_Nb) and of the PROPOSED clip threshold over the 864-row
   gain effects table (SHA asserted), counting feasibility flips vs baseline.
   No refitting: flips are reported, never used to retune.

Proposal (research note, NOT frozen): tau_clip on the full-record clip_frac
scale. Full-record is the recorded scale: every operator sees the whole record,
and the deep-time tail (0.57 FS content, see effects results R3) is exactly the
headroom a record-level curve can destroy. The tail-artifact caveat (catalogue
baselines clip only in the tail) is stated in the proposal, not hidden.

Deterministic; r1/r2 byte-identity asserted. No solver; no test families.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from study_t3_damage_ladder import MOTHERS, load_bscan, fermat_times, event_mask
from study_gain_ladder_v0_1 import battery

EFFECTS = ROOT / 'artifacts/research_checks/2026-09-29_gain_effects_r1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_clip_scale_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_clip_scale_r2.json'
EFFECTS_SHA256 = '057ede4c7223ae0e5082fa93b2956a67b8eeb6123de751da690310412d819abb'

TOL = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
TOL_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'

# Metered clip rungs (record share above FS), physical-language anchors:
# 1e-4 = 1 in 10,000 samples; 1e-3 = 0.1%; 5e-3 = 0.5%; 1e-2 = 1%;
# 5e-2 = 5% (record dominated by clipped content).
CLIP_RUNGS = [0.0, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2]

# Proposed clip threshold (gain class, new; everything else reused frozen).
# Rationale in the results doc: admits the incumbent catalogue baselines whose
# clipping is tail-only (documented in effects results), rejects the +40 dB and
# above hard-pull class. Placed between metered construction rungs, not fitted
# to any candidate ranking.
TAU_CLIP_PROPOSAL = 0.002


def load_inputs():
    assert hashlib.sha256(EFFECTS.read_bytes()).hexdigest() == EFFECTS_SHA256
    assert hashlib.sha256(TOL.read_bytes()).hexdigest() == TOL_SHA256
    tol = json.loads(TOL.read_text(encoding='utf-8'))['tolerances']
    eff = json.loads(EFFECTS.read_text(encoding='utf-8'))
    return (float(tol['tau_A']['value']), float(tol['tau_D']['value']),
            float(tol['eps_Nb']['value'])), eff['records']


def clip_ladder_family(fam, mother, geo):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev)
    fs = float(np.max(np.abs(s)))
    x = s / fs
    rows = []
    for p in CLIP_RUNGS:
        if p == 0.0:
            y = x.copy()
            g_db = 0.0
        else:
            amp = np.sort(np.abs(x).ravel())[::-1]   # descending
            k = max(1, int(np.ceil(p * amp.size)))
            amp_k = amp[k - 1]                        # k-th largest
            g = (1.0 + 1e-9) / amp_k               # strict excess with margin
            g_db = float(20.0 * np.log10(1.0 / amp_k))
            pre = g * x
            over_pre = np.abs(pre) > 1.0
            y = np.clip(pre, -1.0, 1.0)
        if p == 0.0:
            over = np.zeros_like(x, dtype=bool)
        else:
            over = over_pre
        clip_frac = float(np.mean(over))
        # Metering assertion: the construction must clip AT LEAST the rung
        # fraction (the top-k order statistic), and at most the rung plus the
        # boundary band (samples within 1e-6 relative of the k-th largest
        # amplitude ride the same construction gain).
        k = max(1, int(np.ceil(p * y.size))) if p > 0 else 0
        if p > 0:
            band = float(np.mean(amp >= amp[k - 1] * (1.0 - 1e-6)))
            assert clip_frac >= p - 1e-12 and clip_frac <= band + 2.0 / y.size, \
                f'clip metering failed: {fam} rung {p} measured {clip_frac}'
        else:
            assert clip_frac == 0.0
        clip_ratio = float(np.sum(y[over] ** 2) / np.sum(y ** 2)) \
            if over.any() else 0.0
        b = battery(y, x, m_ev)
        # Where does the clipping live? Earliest/latest clipped time.
        t_over = t[np.any(over, axis=0)] if over.any() else np.array([])
        rows.append({'family': fam, 'metered_clip_frac': p,
                     'gain_to_clip_db': round(g_db, 4),
                     'clip_frac': round(clip_frac, 9),
                     'clip_ratio': round(clip_ratio, 9),
                     'earliest_clip_ns': round(float(t_over.min()), 3) if t_over.size else None,
                     'latest_clip_ns': round(float(t_over.max()), 3) if t_over.size else None,
                     **b})
    print(fam, 'clip ladder done:', len(rows), 'rungs')
    return rows


def feasibility(row, tau_A, tau_D, eps_Nb, tau_clip):
    return bool(abs(row['a'] - 1.0) <= tau_A and row['D_e'] <= tau_D
                and row['nc_ratio'] <= eps_Nb
                and row['clip_frac'] <= tau_clip)


def sensitivity(records, tau_A, tau_D, eps_Nb):
    base = [feasibility(r, tau_A, tau_D, eps_Nb, TAU_CLIP_PROPOSAL)
            for r in records]
    axes = {'tau_A': tau_A, 'tau_D': tau_D, 'eps_Nb': eps_Nb,
            'tau_clip_proposal': TAU_CLIP_PROPOSAL}
    out = {}
    for name, val in axes.items():
        for factor in (0.5, 1.5):
            pert = {k: v for k, v in [('tau_A', tau_A), ('tau_D', tau_D),
                                      ('eps_Nb', eps_Nb),
                                      ('tau_clip', TAU_CLIP_PROPOSAL)]}
            if name == 'tau_clip_proposal':
                pert['tau_clip'] = val * factor
            else:
                pert[name] = val * factor
            flips = []
            for i, r in enumerate(records):
                f = feasibility(r, pert['tau_A'], pert['tau_D'],
                                pert['eps_Nb'], pert['tau_clip'])
                if f != base[i]:
                    flips.append({'family': r['family'],
                                  'level_db': r['level_db'],
                                  'structure': r['structure'],
                                  'nc_state': r['nc_state'],
                                  'candidate': r['candidate'],
                                  'baseline_feasible': base[i],
                                  'perturbed_feasible': f})
            out[f'{name}x{factor}'] = {
                'value': round(val * factor, 6),
                'flip_count': len(flips),
                'flips': flips[:40],  # cap detail; count is the statistic
                'detail_truncated': len(flips) > 40,
            }
    base_count = sum(base)
    return {'baseline_feasible_count': base_count,
            'baseline_note': 'only the M1 contract-derived oracle is feasible '
                             'at baseline' if base_count == 72 else '',
            'axes': out}


def clip_distribution(records):
    by_cand = {}
    for r in records:
        by_cand.setdefault(r['candidate'], []).append(r['clip_frac'])
    return {c: {'min': round(min(v), 9), 'median': round(float(np.median(v)), 9),
                'max': round(max(v), 9), 'rows': len(v)}
            for c, v in sorted(by_cand.items())}


def build_doc():
    (tau_A, tau_D, eps_Nb), records = load_inputs()
    ladder = []
    for fam, mother, geo in MOTHERS:
        ladder.extend(clip_ladder_family(fam, mother, geo))
    return {
        'schema': 'clip_scale/0.1',
        'date': '2026-09-29',
        'basis': 'gain class design step 4 (user "继续", 2026-09-29); effects '
                 'table 2026-09-29_gain_effects_r1.json (SHA asserted) + frozen '
                 'reward_tolerance_contract_v0.1.json (SHA asserted)',
        'tau_clip_proposal': TAU_CLIP_PROPOSAL,
        'tau_clip_proposal_status': 'PROPOSAL ONLY - not frozen; user '
                                    'confirmation required before entering any '
                                    'contract; no existing number reused',
        'metered_clip_rungs': CLIP_RUNGS,
        'reused_gates': {'tau_A': tau_A, 'tau_D': tau_D, 'eps_Nb': eps_Nb},
        'clip_ladder': ladder,
        'effects_clip_frac_distribution': clip_distribution(records),
        'sensitivity': sensitivity(records, tau_A, tau_D, eps_Nb),
    }


def main():
    doc = build_doc()
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    doc2 = build_doc()
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
