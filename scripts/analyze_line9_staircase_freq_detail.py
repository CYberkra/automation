"""Staircase ablation frequency-resolved follow-up (local CPU, read-only on retrieved NPZ).

Motivation: the batch-level E metric mixes magnitude and wrapped phase; flat/rough showed
non-monotone E(d25) vs E(d125), consistent with phase wrapping once the accumulated
dispersion phase exceeds pi at high tones. This addendum separates:
- E restricted to sub-bands [20,100] and [100,170] MHz (low band: phase error < pi expected),
- magnitude-only relative error (phase-insensitive),
- per-tone phase difference statistics (wrapping check),
- band-limited impulse responses reconstructed from the 501-tone complex responses on a
  common 0.5 ns time base (no cross-grid interpolation), with envelope band rel_delta_rms
  on [20,120] and [150,450] ns.
"""
import argparse
import json
from pathlib import Path

import numpy as np

GRIDS = {'d25': 0.025, 'd125': 0.0125, 'd625': 0.00625}
REF = 'd625'
BANDS_NS = {'direct_surface': [20., 120.], 'multiple': [150., 450.]}
SUBBANDS_HZ = {'low_20_100MHz': (20e6, 100e6), 'high_100_170MHz': (100e6, 170.3e6)}
PROBE_TONES_MHZ = [30, 60, 90, 120, 150, 170]


def band_limited(response, freq, t):
    """x(t) = sum_k S_k exp(2 pi i f_k t); analytic band-limited signal, common time base."""
    return np.exp(2j * np.pi * np.outer(t, freq)) @ response


def rms(v):
    return float(np.sqrt(np.mean(np.abs(v) ** 2)))


def main(a):
    z = np.load(a.npz)
    freq = z['frequency_Hz']
    t = np.arange(0., 900e-9, 0.5e-9)
    out = dict(reference_grid=REF, subbands_Hz=SUBBANDS_HZ, surfaces={})
    for surface in ['flat', 'slope', 'rough']:
        resp = {g: z[f'stair_{surface}_{g}__response'] for g in GRIDS}
        ref = resp[REF]
        signals = {g: band_limited(resp[g], freq, t) for g in GRIDS}
        env_ref = np.abs(signals[REF])
        entry = dict(subband_E={}, magnitude_only={}, phase_deg={}, time_bands={})
        for g in ['d25', 'd125']:
            d = resp[g] - ref
            entry['subband_E'][g] = {
                'full_20_170MHz': rms(d) / rms(ref),
                **{name: rms(d[(freq >= lo) & (freq <= hi)]) / rms(ref[(freq >= lo) & (freq <= hi)])
                   for name, (lo, hi) in SUBBANDS_HZ.items()}}
            entry['magnitude_only'][g] = (
                rms(np.abs(resp[g]) - np.abs(ref)) / rms(np.abs(ref)))
            phase = np.angle(resp[g] * np.conj(ref))
            entry['phase_deg'][g] = {
                'rms_all_tones': float(np.degrees(rms(phase))),
                **{f'{f0}MHz': float(np.degrees(phase[np.argmin(np.abs(freq - f0 * 1e6))]))
                   for f0 in PROBE_TONES_MHZ}}
            env = np.abs(signals[g])
            bands = {}
            for name, (b0, b1) in BANDS_NS.items():
                m = (t >= b0 * 1e-9) & (t <= b1 * 1e-9)
                dn = rms(env_ref[m])
                bands[name] = dict(rel_delta_rms=(rms(env[m] - env_ref[m]) / dn if dn else None),
                                   ref_rms=dn,
                                   peak_shift_ns=float(
                                       (t[m][np.argmax(env[m])] - t[m][np.argmax(env_ref[m])]) * 1e9))
            entry['time_bands'][g] = bands
        out['surfaces'][surface] = entry
    a.out.mkdir(parents=True, exist_ok=True)
    p = a.out / 'staircase_freq_detail.json'
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(a.out / 'staircase_bandlimited.npz',
                        time_s=t, frequency_Hz=freq,
                        **{f'{s}_{g}__bl': v for s in ['flat', 'slope', 'rough']
                           for g, v in
                           [(g, band_limited(z[f'stair_{s}_{g}__response'], freq, t)) for g in GRIDS]})
    print(json.dumps({s: {g: out['surfaces'][s]['subband_E'][g] for g in ['d25', 'd125']}
                      for s in out['surfaces']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--npz', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    main(a)
