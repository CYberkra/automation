"""Ablation batch analysis: official SFCW responses + band metrics; numbers only, no figures.

Runs in the verified v401 env (needs gprMax.toolboxes.SFCW). Figures are rendered locally
from the NPZ/JSON this writes.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import hilbert
from gprMax.toolboxes.SFCW import processing as sf

from hs_capsule_identity import sha256 as sha

FREQ = 20e6 + np.arange(501) * 300e3
C0 = 299792458.
BANDS = {'direct_surface': [20e-9, 120e-9], 'multiple': [150e-9, 450e-9], 'deep': [450e-9, 900e-9]}
# Predicted C3 timing: subsurface leg scales by sqrt(18.4/11); 5 m air two-way fixed.
AIR_TWO_WAY_S = 2 * 5.05 / C0
C3_SCALE = float(np.sqrt(18.4 / 11.0))


def load_trace(path):
    with h5py.File(path) as h:
        x = h['rxs/rx1/Ez'][:]
        dt = float(h.attrs['dt'])
    t = np.arange(len(x)) * dt + 0.5 * dt
    return t, x


def band_peak(t, env, band):
    m = (t >= band[0]) & (t <= band[1])
    i = np.flatnonzero(m)[np.argmax(env[m])]
    return dict(t_ns=float(t[i] * 1e9), amp=float(env[i]))


def main(a):
    c = json.loads((a.execution / 'execution_contract.json').read_text('utf-8'))
    v = json.loads((a.execution / 'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256'] == sha(a.execution / 'execution_contract.json')
    receipts = {r['id']: r for r in v['groups']}
    out = dict(cases={}, bands=BANDS, air_two_way_ns=AIR_TWO_WAY_S * 1e9, c3_scale=C3_SCALE)
    traces = {}
    for g in c['groups']:
        cid = g['id']
        raw = Path(receipts[cid]['native_path'])
        assert sha(raw) == receipts[cid]['native_sha256']
        t, x = load_trace(raw)
        env = np.abs(hilbert(x))
        s = sf.load_source(raw)
        r = sf.load_receiver(raw, '/rxs/rx1', 'Ez')
        q = sf.direct_frequency_response(s, r, FREQ, tail_taper_fraction=0)
        assert q.source_valid.all()
        peaks = {name: band_peak(t, env, band) for name, band in BANDS.items()}
        out['cases'][cid] = dict(role=g['role'], variant=g['variant'],
                                 peaks=peaks,
                                 rms={n: float(np.sqrt(np.mean(x[(t >= b[0]) & (t <= b[1])] ** 2)))
                                      for n, b in BANDS.items()},
                                 response_abs_max=float(np.max(np.abs(q.response))))
        traces[cid] = (t, x, env, np.asarray(q.response))
    # H1-H0 differential per variant.
    out['differential'] = {}
    for g in c['groups']:
        if g['role'] != 'H1':
            continue
        v1 = g['variant']
        t1, x1, e1, _ = traces[f'abl_{v1}_x190_H1']
        t0, x0, e0, _ = traces[f'abl_{v1}_x190_H0']
        assert np.array_equal(t1, t0)
        d = x1 - x0
        out['differential'][v1] = {n: float(np.sqrt(np.mean(d[(t1 >= b[0]) & (t1 <= b[1])] ** 2)))
                                   for n, b in BANDS.items()}
    # C3 shift test on the multiple band: measured envelope peak vs scaled prediction.
    base = out['cases']['abl_base_x190_H1']['peaks']['multiple']['t_ns']
    c3 = out['cases']['abl_c3_cover_er184_x190_H1']['peaks']['multiple']['t_ns']
    pred = AIR_TWO_WAY_S * 1e9 + C3_SCALE * (base - AIR_TWO_WAY_S * 1e9)
    out['c3_shift_test'] = dict(base_multiple_peak_ns=base, c3_multiple_peak_ns=c3,
                                predicted_c3_ns=float(pred), residual_ns=float(c3 - pred))
    a.out.mkdir(parents=True)
    (a.out / 'ablation_analysis.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(a.out / 'ablation_traces.npz',
                        **{cid: np.vstack(traces[cid]) for cid in traces})
    print(json.dumps(dict(cases=len(out['cases']), c3_test=out['c3_shift_test'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    main(a)
