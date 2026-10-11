"""Staircase ablation analysis: cross-grid SFCW response errors + band metrics; numbers only.

Pre-registered metrics (batch manifest analysis_prereg):
- E_s(d) = RMS_tones|S_{s,d}(f) - S_{s,d625}(f)| / RMS_tones|S_{s,d625}(f)|, d in {d25,d125},
  computed on raw responses and on dl-normalised responses (source discretisation scaling check).
- Convergence ratio E(d25)/E(d125) per surface.
- Staircase attribution: excess of E_slope / E_rough over E_flat (flat carries only
  dispersion/source/PML grid effects).
- Secondary: Hann-envelope band rel_delta_rms on [20,120] and [150,450] ns vs d625,
  envelope peak tracking per band.

Runs in the verified v401 env (needs gprMax.toolboxes.SFCW). Figures render locally from NPZ.
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
BANDS = {'direct_surface': [20e-9, 120e-9], 'multiple': [150e-9, 450e-9]}
REF = 'd625'
GRIDS = {'d25': 0.025, 'd125': 0.0125, 'd625': 0.00625}


def load_trace(path):
    with h5py.File(path) as h:
        x = h['rxs/rx1/Ez'][:]
        dt = float(h.attrs['dt'])
    t = np.arange(len(x)) * dt + 0.5 * dt
    return t, x


def rms(v):
    return float(np.sqrt(np.mean(v ** 2)))


def main(a):
    c = json.loads((a.execution / 'execution_contract.json').read_text('utf-8'))
    v = json.loads((a.execution / 'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256'] == sha(a.execution / 'execution_contract.json')
    assert len(v['groups']) == 9
    receipts = {r['id']: r for r in v['groups']}
    traces, spectra = {}, {}
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
        traces[cid] = (t, x, env)
        spectra[cid] = np.asarray(q.response)
    out = dict(bands={k: [b[0] * 1e9, b[1] * 1e9] for k, b in BANDS.items()},
               reference_grid=REF, surfaces={})
    for surface in ['flat', 'slope', 'rough']:
        entry = dict(frequency={}, time_bands={}, peaks={})
        t_ref, x_ref, env_ref = traces[f'stair_{surface}_{REF}']
        s_ref = spectra[f'stair_{surface}_{REF}']
        den = float(np.sqrt(np.mean(np.abs(s_ref) ** 2)))
        den_norm = float(np.sqrt(np.mean(np.abs(s_ref / GRIDS[REF]) ** 2)))
        for label, dl in GRIDS.items():
            t, x, env = traces[f'stair_{surface}_{label}']
            s = spectra[f'stair_{surface}_{label}']
            e_raw = float(np.sqrt(np.mean(np.abs(s - s_ref) ** 2)) / den)
            e_norm = float(np.sqrt(np.mean(np.abs(s / dl - s_ref / GRIDS[REF]) ** 2)) / den_norm)
            entry['frequency'][label] = dict(E_raw=e_raw, E_dl_normalised=e_norm)
            env_i = np.interp(t_ref, t, env)
            bands = {}
            for name, b in BANDS.items():
                m = (t_ref >= b[0]) & (t_ref <= b[1])
                d = rms(env_i[m] - env_ref[m])
                dn = rms(env_ref[m])
                bands[name] = dict(rel_delta_rms=(d / dn if dn else None),
                                   ref_rms=dn, delta_rms=d)
                i = np.flatnonzero(m)[np.argmax(env_i[m])]
                ir = np.flatnonzero(m)[np.argmax(env_ref[m])]
                bands[name]['peak_t_ns'] = float(t_ref[i] * 1e9)
                bands[name]['peak_shift_ns'] = float((t_ref[i] - t_ref[ir]) * 1e9)
            entry['time_bands'][label] = bands
        entry['convergence_ratio_d25_over_d125'] = (
            entry['frequency']['d25']['E_raw'] / entry['frequency']['d125']['E_raw']
            if entry['frequency']['d125']['E_raw'] else None)
        out['surfaces'][surface] = entry
    out['staircase_attribution'] = {
        s: dict(excess_slope_over_flat=None) for s in ['d25', 'd125']}
    for label in ['d25', 'd125']:
        ef = out['surfaces']['flat']['frequency'][label]['E_raw']
        out['staircase_attribution'][label] = dict(
            E_flat=ef,
            E_slope=out['surfaces']['slope']['frequency'][label]['E_raw'],
            E_rough=out['surfaces']['rough']['frequency'][label]['E_raw'],
            slope_excess=out['surfaces']['slope']['frequency'][label]['E_raw'] - ef,
            rough_excess=out['surfaces']['rough']['frequency'][label]['E_raw'] - ef)
    a.out.mkdir(parents=True)
    (a.out / 'staircase_analysis.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(
        a.out / 'staircase_traces.npz',
        frequency_Hz=FREQ,
        **{f'{cid}__{name}': arr
           for cid, arrays in traces.items() for name, arr in zip(['t', 'x', 'env'], arrays)},
        **{f'{cid}__response': spectra[cid] for cid in spectra})
    print(json.dumps(out['staircase_attribution']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    main(a)
