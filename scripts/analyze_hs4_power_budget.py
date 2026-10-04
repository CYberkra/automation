"""Power-budget post-processing study on archived noise-free traces (no new solver runs).

Physics: in linear-media FDTD without a noise model, transmit power is a scalar
multiplier and changes nothing relative. Power only matters against a receiver
noise floor. This script therefore takes the archived noise-free rough/halfspace
pairs at four antenna heights (ground / 2 m / 8 m / 15 m) and adds a DECLARED
noise model in post-processing, then sweeps transmit power.

Declared noise model (assumption, not measured):
- per-frequency-bin i.i.d. complex Gaussian noise, sigma = max_f|H_total(f)| / 10^(D/20)
- D = system dynamic range at P_ref = 36 dBm (Wang 2026 best field power)
- receiver-dominated floor: +1 dB transmit power -> +1 dB effective D
- contrast noise sigma_contrast = sqrt(2)*sigma (two independent noisy acquisitions,
  perfect coherent subtraction); total noise sigma
Detection metric: margin_dB = 20*log10(A_echo / RMS_noise_chain_output) inside a
+/-25 ns window around the true-relief two-way time; 6 dB margin = detectable.
Diagnostic only; not a physical acceptance criterion.
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
from hs4_large_domain_dense import N_STATIONS, STEP, RX_OFFSET, station_x
from hs4_ground_gpr_dense import GROUND_Z
from analyze_hs4_height_wavefield import response
from check_hs4_v4_factor_evidence import direct_response
from plot_hs4_permittivity_bscans import C_AIR, SURFACE_Z, interface_relief
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

CHECKS = ROOT / 'artifacts/research_checks'
TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
P_REF_DBM = 36.0
D_LEVELS = [60.0, 80.0, 100.0, 120.0]
P_LEVELS = [10.0, 20.0, 30.0, 36.0, 46.0]
PANEL_POWERS = [20.0, 30.0, 36.0, 46.0]
PANEL_D = 100.0
MARGIN_DB = 6.0
K_NOISE = 8
SEED = 20261005
COVER_EPS_R = 18.017

# height label -> (capsule dir, verification file, group id pattern, antenna z, dense 13 stations?)
HEIGHTS = {
    'ground': ('2026-10-04_hs4_ground_gpr_r1', 'completed_verification.json', 'g{k:02d}_{role}', GROUND_Z, True),
    '2m': ('2026-10-04_hs4_height_wavefield_continuation', 'independent_verification.json', 'low_{role}', 14.0, False),
    '8m': ('2026-10-04_hs4_height8m_wavefield_c', 'completed_verification.json', 'mid8_{role}', 20.0, False),
    '15m': ('2026-10-04_hs4_large_domain_dense_r2', 'completed_verification.json', 's{k:02d}_{role}', 27.0, True),
}
STATIONS_SINGLE = [6]  # sparse capsules contain only the centre station


def relief_twoway(ant_z, xmids, fx, fz, v_cover):
    out = []
    for xmid in xmids:
        z = float(np.interp(xmid, fx, fz))
        out.append(2 * np.hypot(0.65, ant_z - SURFACE_Z) / C_AIR
                   + 2 * (SURFACE_Z - z) / v_cover)
    return np.array(out)


def noise_spectrum(rng, shape, sigma):
    return sigma * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    fx, fz = interface_relief()
    v_cover = C_AIR / np.sqrt(COVER_EPS_R)

    spectra, meta, dft = {}, {}, {}
    for label, (cap, vfile, pat, ant_z, dense) in HEIGHTS.items():
        v = json.loads((CHECKS / cap / vfile).read_text('utf-8'))
        if v['status'] != 'PASS':
            raise ValueError(f'{cap} verification not PASS')
        rows = {g['id']: g for g in v['groups']}
        ks = range(N_STATIONS) if dense else STATIONS_SINGLE
        for k in ks:
            for role in ('rough', 'halfspace'):
                gid = pat.format(k=k, role=role)
                raw = CHECKS / cap / gid / 'profile.h5'
                if sha256(raw) != rows[gid]['raw_sha256']:
                    raise ValueError(f'raw identity differs: {gid}')
                r = response(raw)
                spectra[label, k, role] = r
                exact = direct_response(raw, r.frequency[TAKE])
                err = float(np.linalg.norm(exact - r.response[TAKE]) / np.linalg.norm(exact))
                if err > 1e-9:
                    raise ValueError(f'independent DFT failed: {gid}')
                dft[f'{label}_{gid}'] = err
        meta[label] = {'capsule': cap, 'antenna_z_m': ant_z, 'dense': dense,
                       'source_waveform': 'impulse' if 'ground' in cap or 'large_domain' in cap else 'ricker95'}

    template = spectra['ground', 6, 'rough']
    for key, r in spectra.items():
        if not np.array_equal(r.frequency, template.frequency):
            raise ValueError('frequency grids differ')

    products, relief, xmids_all = {}, {}, {}
    for label, (_, _, _, ant_z, dense) in HEIGHTS.items():
        ks = range(N_STATIONS) if dense else STATIONS_SINGLE
        xmids = np.array([station_x(k) + RX_OFFSET / 2 for k in ks])
        xmids_all[label] = xmids
        relief[label] = relief_twoway(ant_z, xmids, fx, fz, v_cover)
        for k in ks:
            rough = spectra[label, k, 'rough']
            contrast_spec = rough.response - spectra[label, k, 'halfspace'].response
            total = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            contrast = reconstruct_time_response(replace(template, response=contrast_spec),
                                                 window='hann', zero_pad_factor=8)
            products[label, k, 'total'] = total
            products[label, k, 'contrast'] = contrast
    time = products['ground', 6, 'total'].time * 1e9
    nfreq = len(template.frequency)

    # echo/direct levels (noise-free)
    levels = {}
    for label in HEIGHTS:
        ks = range(N_STATIONS) if HEIGHTS[label][4] else STATIONS_SINGLE
        per_k = []
        for i, k in enumerate(ks):
            m = np.abs(time - relief[label][i]) <= 25.0
            a_echo = float(np.max(np.abs(products[label, k, 'contrast'].complex_envelope[m])))
            a_direct = float(np.max(np.abs(products[label, k, 'total'].complex_envelope)))
            per_k.append({'echo_peak': a_echo, 'direct_peak': a_direct,
                          'echo_over_direct_dB': 20 * np.log10(a_echo / a_direct)})
        levels[label] = per_k

    # noise RMS after the chain per (label, D_eff): analytic-free Monte Carlo, K realizations
    noise_rms_cache = {}
    sigma0 = {}  # per label: reference sigma at D=0 dB (= max |H_total| at centre)
    for label in HEIGHTS:
        sigma0[label] = float(np.max(np.abs(spectra[label, 6, 'rough'].response)))
    def chain_noise_rms(label, i, window_mask, d_eff, kind):
        key = (label, i, d_eff, kind)
        if key in noise_rms_cache:
            return noise_rms_cache[key]
        sigma = sigma0[label] / 10 ** (d_eff / 20.0)
        if kind == 'contrast':
            sigma *= np.sqrt(2)  # two independent noisy acquisitions
        vals = []
        for j in range(K_NOISE):
            rng = np.random.default_rng((SEED, hash(key) % (2**31), j))
            ns = noise_spectrum(rng, nfreq, sigma)
            prod = reconstruct_time_response(replace(template, response=ns),
                                             window='hann', zero_pad_factor=8)
            vals.append(np.sqrt(np.mean(np.abs(prod.complex_envelope[window_mask]) ** 2)))
        out = float(np.mean(vals))
        noise_rms_cache[key] = out
        return out

    margins = {}
    for label in HEIGHTS:
        ks = range(N_STATIONS) if HEIGHTS[label][4] else STATIONS_SINGLE
        margins[label] = {}
        i = ks.index(6) if not isinstance(ks, list) else ks.index(6)
        m = np.abs(time - relief[label][i]) <= 25.0
        a_echo = levels[label][i]['echo_peak']
        for D in D_LEVELS:
            for P in P_LEVELS:
                d_eff = D + (P - P_REF_DBM)
                rms = chain_noise_rms(label, i, m, d_eff, 'contrast')
                margins[label][f'D{D:g}_P{P:g}'] = 20 * np.log10(a_echo / rms)

    # minimum requirements
    requirements = {}
    for label in HEIGHTS:
        need_D = {}
        for P in P_LEVELS:
            ok = [D for D in D_LEVELS if margins[label][f'D{D:g}_P{P:g}'] >= MARGIN_DB]
            need_D[f'P{P:g}'] = min(ok) if ok else None
        need_P = {}
        for D in D_LEVELS:
            ok = [P for P in P_LEVELS if margins[label][f'D{D:g}_P{P:g}'] >= MARGIN_DB]
            need_P[f'D{D:g}'] = min(ok) if ok else None
        requirements[label] = {'min_dynamic_range_dB_for_6dB_margin': need_D,
                               'min_power_dBm_for_6dB_margin': need_P}

    # figure: noisy contrast B-scans (ground, 15m) x powers + total column at best case
    vmax_c = max(np.max(np.abs(products[l, k, 'contrast'].real_bandpass))
                 for l in ('ground', '15m') for k in range(N_STATIONS))
    vmax_t = max(np.max(np.abs(products[l, k, 'total'].real_bandpass))
                 for l in ('ground', '15m') for k in range(N_STATIONS))
    view = {'ground': (10.0, 160.0), '15m': (80.0, 240.0)}
    fig, axes = plt.subplots(2, len(PANEL_POWERS) + 1, figsize=(19, 7), sharey=False,
                             layout='constrained')
    for row, label in enumerate(('ground', '15m')):
        lo, hi = view[label]
        mt = (time >= lo) & (time <= hi)
        t_edge = np.concatenate(([time[mt][0] - (time[1] - time[0]) / 2],
                                 (time[mt][:-1] + time[mt][1:]) / 2,
                                 [time[mt][-1] + (time[1] - time[0]) / 2]))
        for col, P in enumerate(PANEL_POWERS):
            d_eff = PANEL_D + (P - P_REF_DBM)
            sigma = sigma0[label] / 10 ** (d_eff / 20.0) * np.sqrt(2)
            ax = axes[row, col]
            for k in range(N_STATIONS):
                rng = np.random.default_rng((SEED, k, int(P * 10), row))
                ns = noise_spectrum(rng, nfreq, sigma)
                spec = (spectra[label, k, 'rough'].response
                        - spectra[label, k, 'halfspace'].response) + ns
                prod = reconstruct_time_response(replace(template, response=spec),
                                                 window='hann', zero_pad_factor=8)
                ax.pcolormesh([station_x(k) - STEP / 2, station_x(k) + STEP / 2], t_edge,
                              prod.real_bandpass[mt][:, None], cmap='RdBu_r', shading='flat',
                              norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * vmax_c,
                                                                vmin=-vmax_c, vmax=vmax_c))
            ax.plot(xmids_all[label], relief[label], 'k-', lw=1.0)
            ax.set_xlim(station_x(0) - STEP / 2, station_x(N_STATIONS - 1) + RX_OFFSET + STEP / 2)
            ax.set_ylim(hi, lo)
            ax.set_title(f'{label} contrast+noise  P={P:g}dBm (D={PANEL_D:g}dB)', fontsize=9)
            ax.set_xlabel('source x (m)')
            if col == 0:
                ax.set_ylabel('actual receiver time (ns)')
        ax = axes[row, len(PANEL_POWERS)]
        for k in range(N_STATIONS):
            ax.pcolormesh([station_x(k) - STEP / 2, station_x(k) + STEP / 2], t_edge,
                          products[label, k, 'total'].real_bandpass[mt][:, None],
                          cmap='RdBu_r', shading='flat',
                          norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * vmax_t,
                                                            vmin=-vmax_t, vmax=vmax_t))
        ax.plot(xmids_all[label], relief[label], 'k-', lw=1.0)
        ax.set_xlim(station_x(0) - STEP / 2, station_x(N_STATIONS - 1) + RX_OFFSET + STEP / 2)
        ax.set_ylim(hi, lo)
        ax.set_title(f'{label} total noise-free (coherent limit)', fontsize=9)
        ax.set_xlabel('source x (m)')
    fig.suptitle('Power only fights the noise floor (declared model, D at 36 dBm); '
                 'contrast rows share one symlog scale; total column shows the coherent '
                 'direct/surface-wave limit that no power level fixes')

    a.out.mkdir(parents=True)
    fig.savefig(a.out / 'power_sweep_panels.png', dpi=150)

    # margin curves
    fig2, ax2 = plt.subplots(figsize=(9, 5.5), layout='constrained')
    for label in HEIGHTS:
        for D in (80.0, 100.0):
            ms = [margins[label][f'D{D:g}_P{P:g}'] for P in P_LEVELS]
            ax2.plot(P_LEVELS, ms, marker='o', label=f'{label}, D={D:g} dB')
    ax2.axhline(MARGIN_DB, color='k', ls='--', lw=1, label=f'{MARGIN_DB:g} dB detection margin')
    ax2.set_xlabel('transmit power (dBm)')
    ax2.set_ylabel('interface-echo margin over noise floor (dB)')
    ax2.legend(fontsize=8, ncol=2)
    ax2.set_title('Declared noise model; margins include ~27 dB SFCW chain processing gain')
    fig2.savefig(a.out / 'power_margin_curves.png', dpi=150)

    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'code_sha256': sha256(__file__),
               'noise_model': {'type': 'per-frequency-bin iid complex Gaussian, flat spectrum',
                               'sigma_definition': 'max_f|H_total_centre(f)| / 10^(D/20); contrast sigma x sqrt(2)',
                               'D_reference_power_dBm': P_REF_DBM,
                               'power_scaling': '+1 dB transmit power -> +1 dB effective dynamic range (receiver-dominated floor)',
                               'seed': SEED, 'K_realizations': K_NOISE},
               'heights': meta,
               'echo_over_direct_dB_per_station': {l: [p['echo_over_direct_dB'] for p in v]
                                                   for l, v in levels.items()},
               'margin_dB': margins,
               'requirements': requirements,
               'relief_two_way_ns': {l: list(map(float, r)) for l, r in relief.items()},
               'independent_direct_DFT_relative_L2': dft,
               'coherent_limit_statement': 'In the raw total B-scan the interface echo sits under coherent direct/surface waves (48.8 dB down at 15 m); transmit power cannot change coherent ratios.',
               'physical_attribution_certified': False,
               'scope': 'Post-processing of archived noise-free traces with a DECLARED noise model; no new solver runs; noise floor not measured from hardware; 2D line source; no field claims.'}
    (a.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                                   allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'echo_over_direct_dB': {l: round(levels[l][0 if not HEIGHTS[l][4] else 6]['echo_over_direct_dB'], 2) for l in HEIGHTS},
                      'requirements': requirements}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
