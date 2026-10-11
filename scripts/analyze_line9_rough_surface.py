"""Rough-vs-smooth sensitivity analysis at x190 AGL5: rough pair vs ablation base pair (numbers only).

Runs in the verified v401 env (needs gprMax.toolboxes.SFCW). Smooth reference:
line9_material_ablation_20261011_r1 combined view (abl_base_x190_H0/H1), not re-run.
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
BANDS = {'direct_surface': [20e-9, 120e-9], 'multiple': [150e-9, 450e-9], 'deep': [450e-9, 900e-9]}


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


def band_rms(t, x, band):
    m = (t >= band[0]) & (t <= band[1])
    return float(np.sqrt(np.mean(x[m] ** 2)))


def receipts(execution):
    c = json.loads((execution / 'execution_contract.json').read_text('utf-8'))
    v = json.loads((execution / 'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256'] == sha(execution / 'execution_contract.json')
    return {r['id']: r for r in v['groups']}


def main(a):
    rough = receipts(a.rough)
    smooth = receipts(a.smooth)
    pairs = {'H0': ('rough_x190_H0', 'abl_base_x190_H0'), 'H1': ('rough_x190_H1', 'abl_base_x190_H1')}
    out = dict(cases={}, sensitivity={}, bands=BANDS,
               smooth_reference='line9_material_ablation_20261011_r1 abl_base_x190 pair')
    traces = {}
    for role, (rid, sid) in pairs.items():
        rr, sr = rough[rid], smooth[sid]
        assert sha(rr['native_path']) == rr['native_sha256'] and sha(sr['native_path']) == sr['native_sha256']
        t_r, x_r = load_trace(rr['native_path'])
        t_s, x_s = load_trace(sr['native_path'])
        assert np.array_equal(t_r, t_s)
        e_r, e_s = np.abs(hilbert(x_r)), np.abs(hilbert(x_s))
        q_r = sf.direct_frequency_response(sf.load_source(rr['native_path']), sf.load_receiver(rr['native_path'], '/rxs/rx1', 'Ez'), FREQ, tail_taper_fraction=0)
        q_s = sf.direct_frequency_response(sf.load_source(sr['native_path']), sf.load_receiver(sr['native_path'], '/rxs/rx1', 'Ez'), FREQ, tail_taper_fraction=0)
        assert q_r.source_valid.all() and q_s.source_valid.all()
        out['cases'][role] = dict(
            rough={n: band_peak(t_r, e_r, b) for n, b in BANDS.items()},
            smooth={n: band_peak(t_s, e_s, b) for n, b in BANDS.items()})
        out['sensitivity'][role] = {n: dict(
            delta_rms=band_rms(t_r, x_r - x_s, b),
            rel_delta_rms=band_rms(t_r, x_r - x_s, b) / band_rms(t_s, x_s, b),
            delta_env_rms=band_rms(t_r, e_r - e_s, b),
            peak_time_shift_ns=out['cases'][role]['rough'][n]['t_ns'] - out['cases'][role]['smooth'][n]['t_ns'])
            for n, b in BANDS.items()}
        traces[f'{role}_rough'] = (t_r, x_r, e_r, np.asarray(q_r.response))
        traces[f'{role}_smooth'] = (t_s, x_s, e_s, np.asarray(q_s.response))
    # H1-H0 differential under each surface.
    diff = {}
    for surf in ['rough', 'smooth']:
        t1, x1 = traces[f'H1_{surf}'][0], traces[f'H1_{surf}'][1]
        x0 = traces[f'H0_{surf}'][1]
        diff[surf] = {n: band_rms(t1, x1 - x0, b) for n, b in BANDS.items()}
    out['differential_rms'] = diff
    out['differential_change'] = {n: dict(
        rel=(diff['rough'][n] - diff['smooth'][n]) / diff['smooth'][n]
        if diff['smooth'][n] != 0 else None,
        note=None if diff['smooth'][n] != 0 else 'smooth differential RMS is exactly zero in this band')
        for n in BANDS}
    a.out.mkdir(parents=True)
    (a.out / 'rough_analysis.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(a.out / 'rough_traces.npz',
                        **{f'{cid}__{name}': arr for cid, arrays in traces.items()
                           for name, arr in zip(['t', 'x', 'env', 'response'], arrays)})
    print(json.dumps(dict(cases=len(out['cases']), differential_change=out['differential_change'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rough', type=Path, required=True)
    p.add_argument('--smooth', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
