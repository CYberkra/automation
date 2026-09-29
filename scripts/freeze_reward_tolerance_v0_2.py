"""Freeze the reward tolerance contract v0.2 (gain class step 4: user-confirmed proposal).

User confirmed the tau_clip proposal ("确认", 2026-09-29) in docs/research/
2026-09-29_tolerance_v0_2_proposal.md. v0.2 = v0.1 verbatim (hash-locked) plus ONE
new tolerance, tau_clip = 0.002 (full-record clipped-sample fraction, gain class).

Every anchor claim is asserted against the committed evidence files before
writing; the document is built twice and byte-compared (freeze discipline).
Anchor claims:
  * v0.1 contract SHA-256 locked; reused tolerances copied verbatim;
  * metered clip ladder: zero rung exact, 0.1% rung confined to the direct
    wave with the event battery untouched, 0.5% rung record-wide;
  * effects-table distribution: identity/oracle exactly 0, catalogue baselines
    below the proposal, every mechanism-reference hard pull above it;
  * sensitivity: tau_clip +/-50% zero flips; eps x0.5 flips exactly the 72
    oracle rows (structural floor 1.0 reproduced in the gain class);
    tau_A x1.5 flips 36 (ladder boundary property, not refitted).
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/reward_tolerance_contract_v0.2.json'
V01 = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
CLIP = ROOT / 'artifacts/research_checks/2026-09-29_clip_scale_r1.json'
EFFECTS = ROOT / 'artifacts/research_checks/2026-09-29_gain_effects_r1.json'

V01_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'
CLIP_SHA256 = None    # filled after first run is committed; asserted if present
EFFECTS_SHA256 = '057ede4c7223ae0e5082fa93b2956a67b8eeb6123de751da690310412d819abb'

TAU_CLIP = 0.002


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_anchors(v01, clip):
    lad = clip['clip_ladder']
    dist = clip['effects_clip_frac_distribution']
    sens = clip['sensitivity']
    axes = sens['axes']

    # Metered ladder anchors (all three families). The uniform metered gain
    # rescales the event window (a = g, H = 0); clipping reaches the event
    # window itself only at the 5% rung (+69 dB), so H stays the sentinel.
    for fam in ('S2X', 'S2TZX', 'C3mX'):
        rows = {r['metered_clip_frac']: r for r in lad if r['family'] == fam}
        assert rows[0.0]['clip_frac'] == 0.0
        assert rows[1e-4]['clip_frac'] >= 1e-4 - 1e-12
        assert rows[1e-4]['a'] <= 1.0 + 1e-6 and rows[1e-4]['H'] == 0.0
        assert rows[1e-3]['H'] == 0.0                     # event window unclipped
        assert rows[1e-3]['latest_clip_ns'] <= 10.0       # direct-wave neighbourhood only
        assert rows[5e-3]['latest_clip_ns'] > 3000.0      # record-wide destruction
        assert rows[5e-3]['H'] == 0.0                     # event window still pure rescale
        assert rows[5e-2]['D_e'] > 10.0                   # event window itself clipped

    # Distribution anchors around the proposal value. The hard-pull class is
    # separated on the MEDIAN cell reading (per-cell minima vary structurally:
    # time_global cells attenuate the deep tail too, e.g. M2 min 0.000998).
    assert dist['identity']['max'] == 0.0
    assert dist['M1_matched_inverse']['max'] == 0.0
    assert dist['cat_G4']['max'] < TAU_CLIP
    mech = ('M2_exp_end_40dB', 'M3_exp_end_60dB', 'M4_sec_model', 'M5_sec_2a',
            'M6_power_t1', 'M7_rms_agc_50ns', 'M8_env_gain_50ns',
            'M9_smart_interp_16pt')
    assert TAU_CLIP < min(dist[c]['median'] for c in mech)

    # Sensitivity anchors.
    assert sens['baseline_feasible_count'] == 72
    assert axes['tau_clip_proposalx0.5']['flip_count'] == 0
    assert axes['tau_clip_proposalx1.5']['flip_count'] == 0
    assert axes['eps_Nbx0.5']['flip_count'] == 72
    assert axes['eps_Nbx1.5']['flip_count'] == 0
    assert axes['tau_Ax1.5']['flip_count'] == 36
    assert axes['tau_Dx0.5']['flip_count'] == 0
    assert axes['tau_Dx1.5']['flip_count'] == 0

    # v0.1 reuse: copied verbatim.
    t = v01['tolerances']
    for k in ('tau_A', 'tau_D', 'eps_Nb'):
        assert t[k]['value'] == {'tau_A': 0.20, 'tau_D': 0.95, 'eps_Nb': 1.5}[k]
    return t


def build_doc():
    assert sha256(V01) == V01_SHA256, 'v0.1 hash lock broken'
    assert sha256(EFFECTS) == EFFECTS_SHA256, 'effects table hash changed'
    if CLIP_SHA256 is not None:
        assert sha256(CLIP) == CLIP_SHA256
    v01 = json.loads(V01.read_text(encoding='utf-8'))
    clip = json.loads(CLIP.read_text(encoding='utf-8'))
    t = verify_anchors(v01, clip)
    doc = {
        'schema': 'reward-tolerance-contract/0.2',
        'status': 'frozen',
        'date': '2026-09-29',
        'frozen_by': 'user confirmation "确认" (2026-09-29) of docs/research/'
                     '2026-09-29_tolerance_v0_2_proposal.md (gain class step 4)',
        'supersedes': {'path': 'reward_tolerance_contract_v0.1.json',
                       'sha256': V01_SHA256,
                       'note': 'v0.1 remains valid for the background class; v0.2 '
                               'adds one gain-class tolerance, reuses the rest '
                               'verbatim'},
        'tolerances': {
            'tau_A': t['tau_A'],
            'tau_D': t['tau_D'],
            'eps_Nb': t['eps_Nb'],
            'arrival_drift': t['arrival_drift'],
            'tau_clip': {'value': TAU_CLIP,
                         'meaning': 'gain-class only: full-record clipped-sample '
                                    'fraction (|y|>1 vs simulated_full_scale) on '
                                    'the operator output',
                         'anchor': 'metered clip ladder (0.1% rung = direct-wave '
                                   'neighbourhood only, 0.5% rung = record-wide '
                                   'headroom destruction); incumbent catalogue '
                                   'baselines 0.0005-0.00125 below, all mechanism-'
                                   'reference hard pulls >=0.0022 above; +/-50% '
                                   'sensitivity zero ranking effect on the current '
                                   'effects table; guards the future high-recovery '
                                   '+over-drive attack surface'},
        },
        'judgement_rules': {
            'implementation': 'scripts/study_clip_scale_v0_1.py (feasibility) '
                              '+ scripts/study_gain_effects_v0_1.py (metrics)',
            'background_class': 'v0.1 rules unchanged (no FS convention there)',
            'gain_class_feasible': 'D_e<=tau_D and |a-1|<=tau_A and '
                                   'nc_ratio<=eps_Nb and clip_frac<=tau_clip',
            'hard_constraints_first': True,
            'scalar_weights': 'R_bg weights frozen in '
                              'reward_weights_contract_v0.2.json; R_gain '
                              'construction proposal at gain class step 5 '
                              '(dB-identity coefficients, user confirmation '
                              'pending at v0.2 freeze time)',
        },
        'evidence': {
            'tolerance_contract_v0_1': {'path': V01.name, 'sha256': V01_SHA256},
            'clip_scale_r1': {'path': CLIP.name, 'sha256': sha256(CLIP)},
            'gain_effects_r1': {'path': EFFECTS.name, 'sha256': EFFECTS_SHA256},
        },
        'scope_limits': v01['scope_limits'] + [
            'tau_clip is a property of the simulated_full_scale array convention '
            '(FS = max|original signature|, fixed per family); not a device ADC '
            'specification; conclusions do not extrapolate to hardware',
        ],
    }
    return doc


def main():
    text1 = json.dumps(build_doc(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    text2 = json.dumps(build_doc(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text1 == text2, 'double-build mismatch'
    OUT.write_text(text1, encoding='utf-8', newline='\n')
    print('anchors verified | wrote', OUT.relative_to(ROOT))
    print('sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest())


if __name__ == '__main__':
    main()
