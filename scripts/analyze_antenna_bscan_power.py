"""Antenna-B-scan analysis + 5 W / 20 W power post-processing (P2 batch).

User (2026-10-05): 后处理给我看看 5W、20W 功率的 B-scan.
Reads the frozen P2 batch (13 stations x rough/halfspace, x4-scaled GSSI-400
antenna at 15 m altitude), runs the official SFCW actual-source chain per
trace, then applies the SAME declared noise model as
analyze_hs4_power_budget.py (per-frequency-bin iid complex Gaussian,
sigma = max|H_total_centre| / 10^(D_eff/20), contrast sigma x sqrt(2),
D_eff = D_REF + (P - 36 dBm), +1 dB power = +1 dB effective dynamic range).
37 dBm = 5 W, 43 dBm = 20 W. Diagnostic only; the noise floor is a declared
assumption, not hardware-measured; not a physical acceptance run.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from hs_capsule_identity import sha256
from analyze_hs4_power_budget import (relief_twoway, noise_spectrum, P_REF_DBM,
                                      PANEL_D, SEED, COVER_EPS_R)
from plot_hs4_permittivity_bscans import C_AIR, interface_relief
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime
from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
                                              direct_frequency_response,
                                              reconstruct_time_response)

P2 = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p2'
N_STATIONS = 13
STEP = 0.5
POWERS_DBM = {37.0: '5W', 43.0: '20W'}
TAIL_TAPER = 0.2
ANT_Z_ABS = 27.0
VIEW_NS = (80.0, 240.0)
# Declared range gate on the CONTRAST product only: removes direct/surface waves
# and the cross-scene numerical floor observed at their times (see
# docs/research/2026-10-05_cross_scene_subtraction_floor.md); causal floor for
# any bedrock-relief influence at the receiver is ~118 ns.
GATE_START_NS, GATE_FULL_NS = 120.0, 140.0


def gate(tr):
    """Zero the contrast before GATE_START_NS with a cos ramp to GATE_FULL_NS."""
    t = tr.time * 1e9
    w = np.ones_like(t)
    w[t < GATE_START_NS] = 0.0
    m = (t >= GATE_START_NS) & (t < GATE_FULL_NS)
    w[m] = 0.5 * (1 - np.cos(np.pi * (t[m] - GATE_START_NS)
                             / (GATE_FULL_NS - GATE_START_NS)))
    return replace(tr, complex_envelope=tr.complex_envelope * w,
                   complex_bandpass=tr.complex_bandpass * w,
                   real_bandpass=tr.real_bandpass * w)


def station_x(k):
    return 15.25 + 0.5 * k


def load_traces():
    index = json.loads((P2 / 'batch_index.json').read_text('utf-8'))
    if len(index) != 2 * N_STATIONS or any(e['status'] != 'PASS' for e in index):
        raise ValueError('P2 batch incomplete or not PASS')
    spectra, midx = {}, {}
    for e in index:
        h5 = P2 / e['tag'] / 'p2.h5'
        if sha256(h5) != e['h5_sha256']:
            raise ValueError(f'raw identity differs: {e["tag"]}')
        src = load_source(h5)
        rx = load_receiver(h5, receiver_path='/rxs/rx1', component='Ey')
        r = direct_frequency_response(src, rx, FREQ, tail_taper_fraction=TAIL_TAPER)
        if not np.all(r.source_valid):
            raise ValueError(f'source invalid: {e["tag"]}')
        spectra[e['tag']] = r
        midx[e['tag']] = 0.5 * (e['src_Position'][0] + e['rx_Position'][0]) + 11.6
    return spectra, midx


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    spectra, midx = load_traces()
    template = spectra['rough_k06']
    nfreq = len(template.response)

    xmids = np.array([midx[f'rough_k{k:02d}'] for k in range(N_STATIONS)])
    fx, fz = interface_relief()
    relief = relief_twoway(ANT_Z_ABS, xmids, fx, fz, C_AIR / np.sqrt(COVER_EPS_R))

    def tag(role, k):
        return f'{role}_k{k:02d}'

    def product(spec):
        return reconstruct_time_response(replace(template, response=spec),
                                         window='hann', zero_pad_factor=8)
    total, contrast, contrast_g = {}, {}, {}
    for k in range(N_STATIONS):
        total[k] = product(spectra[tag('rough', k)].response)
        contrast[k] = product(spectra[tag('rough', k)].response
                              - spectra[tag('halfspace', k)].response)
        contrast_g[k] = gate(contrast[k])
    time = total[0].time * 1e9

    # noise-free levels (gated contrast: echo window is well past the gate)
    levels = []
    for i, k in enumerate(range(N_STATIONS)):
        m = np.abs(time - relief[i]) <= 25.0
        a_echo = float(np.max(np.abs(contrast_g[k].complex_envelope[m])))
        a_direct = float(np.max(np.abs(total[k].complex_envelope)))
        levels.append({'station_k': k, 'echo_peak': a_echo, 'direct_peak': a_direct,
                       'echo_over_direct_dB': 20 * np.log10(a_echo / a_direct)})

    sigma0 = float(np.max(np.abs(template.response)))
    margins = {}
    i6 = 6
    m6 = np.abs(time - relief[i6]) <= 25.0
    a_echo6 = levels[i6]['echo_peak']
    for P in POWERS_DBM:
        d_eff = PANEL_D + (P - P_REF_DBM)
        sigma = sigma0 / 10 ** (d_eff / 20.0) * np.sqrt(2)
        rms = []
        for real in range(8):
            rng = np.random.default_rng((SEED, i6, int(P * 10), real))
            ns = noise_spectrum(rng, nfreq, sigma)
            prod = gate(product(ns))
            rms.append(float(np.sqrt(np.mean(np.abs(prod.complex_envelope[m6]) ** 2))))
        margins[f'P{P:g}'] = {'d_eff': d_eff,
                              'noise_rms_median': float(np.median(rms)),
                              'margin_dB': 20 * np.log10(a_echo6 / np.median(rms))}

    # figure: rows = total clean, contrast clean raw, contrast clean gated,
    #               gated contrast + noise per power
    vmax_t = max(np.max(np.abs(total[k].real_bandpass)) for k in range(N_STATIONS))
    vmax_c = max(np.max(np.abs(contrast[k].real_bandpass)) for k in range(N_STATIONS))
    vmax_g = max(np.max(np.abs(contrast_g[k].real_bandpass)) for k in range(N_STATIONS))
    lo, hi = VIEW_NS
    mt = (time >= lo) & (time <= hi)
    t_edge = np.concatenate(([time[mt][0] - (time[1] - time[0]) / 2],
                             (time[mt][:-1] + time[mt][1:]) / 2,
                             [time[mt][-1] + (time[1] - time[0]) / 2]))
    rows = [('total noise-free', 'total', None, vmax_t),
            ('contrast noise-free (raw)', 'raw', None, vmax_c),
            ('contrast noise-free (gated 120-140ns)', 'gated', None, vmax_g)]
    rows += [(f'contrast gated+noise {POWERS_DBM[P]} (P={P:g}dBm)', 'noisy', P, vmax_g)
             for P in POWERS_DBM]
    fig, axes = plt.subplots(len(rows), 1, figsize=(13, 14), sharex=True,
                             layout='constrained')
    for ax, (label, kind, P, vmax) in zip(axes, rows):
        for k in range(N_STATIONS):
            if kind == 'total':
                band = total[k].real_bandpass
            elif kind == 'raw':
                band = contrast[k].real_bandpass
            elif kind == 'gated':
                band = contrast_g[k].real_bandpass
            else:
                d_eff = PANEL_D + (P - P_REF_DBM)
                sigma = sigma0 / 10 ** (d_eff / 20.0) * np.sqrt(2)
                rng = np.random.default_rng((SEED, k, int(P * 10), 0))
                ns = noise_spectrum(rng, nfreq, sigma)
                band = gate(product(spectra[tag('rough', k)].response
                                    - spectra[tag('halfspace', k)].response
                                    + ns)).real_bandpass
            ax.pcolormesh([station_x(k) - STEP / 2, station_x(k) + STEP / 2], t_edge,
                          band[mt][:, None], cmap='RdBu_r', shading='flat',
                          norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * vmax,
                                                            vmin=-vmax, vmax=vmax))
        ax.plot(xmids, relief, 'k-', lw=1.0)
        ax.set_xlim(station_x(0) - STEP / 2, station_x(N_STATIONS - 1) + STEP / 2)
        ax.set_ylim(hi, lo)
        ax.set_ylabel('ns')
        ax.set_title(label, fontsize=10)
    axes[-1].set_xlabel('source x (m), shifted-back absolute frame')
    fig.suptitle('3D x4-scaled GSSI-400-like antenna, 15 m altitude; declared noise '
                 'model (D=100 dB at 36 dBm, +1 dB/dB); diagnostic only')

    a.out.mkdir(parents=True)
    fig.savefig(a.out / 'antenna_bscan_power_panels.png', dpi=150)
    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'code_sha256': sha256(__file__),
               'p2_batch': str(P2.relative_to(ROOT)),
               'sfcw_chain': {'freqs': '20-170 MHz, 501 bins', 'tail_taper': TAIL_TAPER,
                              'window': 'hann', 'zero_pad_factor': 8},
               'noise_model': {'same_as': 'analyze_hs4_power_budget.py',
                               'D_ref_dB_at_36dBm': PANEL_D, 'seed': SEED},
               'range_gate_ns': [GATE_START_NS, GATE_FULL_NS],
               'range_gate_reason': ('cross-scene subtraction floor at direct/surface-wave '
                                     'times (see 2026-10-05_cross_scene_subtraction_floor.md); '
                                     'causal floor for relief influence ~118 ns'),
               'powers': {f'{P:g}dBm': POWERS_DBM[P] for P in POWERS_DBM},
               'echo_over_direct_dB_per_station': [l['echo_over_direct_dB'] for l in levels],
               'margins_centre_station': margins,
               'relief_two_way_ns': list(map(float, relief)),
               'physical_attribution_certified': False,
               'scope': ('Antenna is an ASSUMED x4-scaled design on a coarse 4 cm grid '
                         '(thin details quantized); noise floor declared, not measured; '
                         'no field claims.')}
    (a.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                                   allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'echo_over_direct_dB_k06': round(levels[6]['echo_over_direct_dB'], 2),
                      'margins': {k: round(v['margin_dB'], 2) for k, v in margins.items()}},
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
