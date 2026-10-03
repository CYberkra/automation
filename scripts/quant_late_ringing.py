# -*- coding: utf-8 -*-
"""Quantify late-time ringing residue: stock V6 (tau=73) vs revised (tau=58)."""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from augment_sim2real_v0_6 import calibrate_q, op1_ringing_v06, SEED  # noqa: E402
from augment_sim2real_v0_3 import load_line9, T_WINDOW_NS, N_SAMPLES  # noqa: E402
sys.path.insert(0, r'E:\automation_djh\automation_repo\scripts')
from render_sim_v6r_allgroups import op1_ringing_v07  # noqa: E402

SFCW = Path(r'E:\automation_djh\artifacts_check\sfcw_r1\sfcw_responses.npz')
QDIR = Path(r'E:\automation_djh\artifacts_check\ringing_field')

real = load_line9()
calib = calibrate_q(real, np.linspace(0, T_WINDOW_NS, N_SAMPLES))
q_all = np.concatenate([np.load(p) for p in sorted(QDIR.glob('*_q.npy'))])
q_lo, q_hi = float(np.percentile(q_all, 25)), float(np.percentile(q_all, 75))

d = np.load(SFCW)
for tag in ['3D_5cm', '2D_5cm', '2D_2p5cm']:
    mats, tm = [], None
    for k in range(1, 14):
        t = f't{k:02d}'
        mats.append(d[f'{tag}_{t}_tail_200ns_hann_real_bandpass'])
        tm = d[f'{tag}_{t}_tail_200ns_hann_envelope_time_s']
    B0 = np.stack(mats)
    tfull = tm * 1e9
    pk = np.max(np.abs(B0), axis=1)

    B1a, _ = op1_ringing_v06(B0, tfull, np.random.default_rng(SEED),
                             calib['p25'], calib['p75'])
    B1b, _ = op1_ringing_v07(B0, tfull, np.random.default_rng(SEED),
                             q_lo, q_hi, 58.0, 0.20)
    out = [tag]
    for lo, hi in [(45.0, 70.0), (160.0, 210.0)]:
        m = (tfull >= lo) & (tfull <= hi)
        ra = np.sqrt(((B1a - B0)[:, m] ** 2).mean(axis=1)) / pk * 100
        rb = np.sqrt(((B1b - B0)[:, m] ** 2).mean(axis=1)) / pk * 100
        out.append(f'{lo:.0f}-{hi:.0f}ns: v6(73) med={np.median(ra):.2f}% '
                   f'v6r(58) med={np.median(rb):.2f}%  ratio={np.median(rb/ra):.2f}')
    print(' | '.join(out))
