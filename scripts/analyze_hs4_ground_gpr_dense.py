"""Ground-coupled vs 15 m airborne dense B-scan comparison (13 stations).

Reads the completed ground capsule (antenna z=12.025 m, one grid step above the
flat surface) and pairs every trace with its archived 15 m raw declared in the
contract (identity-checked). Both run through the identical source-normalized
SFCW chain. Figure: 2 rows (total / contrast) x 2 cols (ground / 15 m), one
symlog norm per row shared by both columns so the dynamic-range change is
visible, with the true-relief two-way-time curve overlaid per column
(dispersion-ignored straight-ray estimate, per-column antenna height).
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
from plot_hs4_permittivity_bscans import C_AIR, SURFACE_Z, ANTENNA_Z, interface_relief
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
COVER_EPS_R = 18.017
VIEW = {'ground': (10.0, 160.0), '15m': (80.0, 240.0)}        # figure y limits (ns)
PEAK_WIN = {'ground': (40.0, 160.0), '15m': (160.0, 220.0)}   # contrast envelope peak search (ns)
ANT_Z = {'ground': GROUND_Z, '15m': float(ANTENNA_Z)}


def relief_two_way(label, v_cover, fx, fz):
    out = []
    for k in range(N_STATIONS):
        tx, rx = station_x(k), station_x(k) + RX_OFFSET
        xmid = 0.5 * (tx + rx)
        z = float(np.interp(xmid, fx, fz))
        out.append((np.hypot(tx - xmid, ANT_Z[label] - SURFACE_Z)
                    + np.hypot(rx - xmid, ANT_Z[label] - SURFACE_Z)) / C_AIR
                   + 2 * (SURFACE_Z - z) / v_cover)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    cp = a.capsule / 'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    v = json.loads((a.capsule / 'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (a.capsule / 'execution.jsonl').read_text('utf-8').splitlines()]
    if (v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp)
            or events[-1]['status'] != 'COMPLETED' or events[-1]['traces'] != 2 * N_STATIONS
            or events[-1]['verification_sha256'] != sha256(a.capsule / 'completed_verification.json')):
        raise ValueError('completed 26-trace audited ground capsule required')
    verify_official_runtime()
    rows = {g['id']: g for g in v['groups']}
    responses, dft, arrays = {}, {}, {}
    for g in c['groups']:
        new = Path(g['input']).with_suffix('.h5')
        old = Path(g['reference_raw'])
        if sha256(new) != rows[g['id']]['raw_sha256'] or sha256(old) != g['reference_raw_sha256']:
            raise ValueError('new/archived raw identity differs')
        for label, raw in [('ground', new), ('15m', old)]:
            r = response(raw)
            responses[label, g['station_index'], g['role']] = r
            exact = direct_response(raw, r.frequency[TAKE])
            error = float(np.linalg.norm(exact - r.response[TAKE]) / np.linalg.norm(exact))
            if error > 1e-9:
                raise ValueError('independent actual-source direct DFT failed')
            dft[label + '_' + g['id']] = error
    template = responses['ground', 6, 'rough']
    products = {}
    for k in range(N_STATIONS):
        for label in ('ground', '15m'):
            rough = responses[label, k, 'rough']
            half = responses[label, k, 'halfspace']
            if (not np.array_equal(rough.frequency, template.frequency)
                    or rough.source.quantity != template.source.quantity
                    or rough.source.spatial_scale != template.source.spatial_scale):
                raise ValueError('source/frequency conventions differ')
            total = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            contrast = reconstruct_time_response(replace(template, response=rough.response - half.response),
                                                 window='hann', zero_pad_factor=8)
            products[label, k, 'total'] = total
            products[label, k, 'contrast'] = contrast
            arrays[f'{label}_s{k:02d}_total_signed'] = total.real_bandpass
            arrays[f'{label}_s{k:02d}_contrast_signed'] = contrast.real_bandpass
            arrays[f'{label}_s{k:02d}_contrast_envelope'] = contrast.complex_envelope
    time = products['ground', 6, 'total'].time * 1e9
    for k in range(N_STATIONS):
        for label in ('ground', '15m'):
            if not np.array_equal(products[label, k, 'total'].time * 1e9, time):
                raise ValueError('time grids differ across cases')
    arrays['time_ns'] = time
    arrays['frequency_hz'] = template.frequency

    fx, fz = interface_relief()
    v_cover = C_AIR / np.sqrt(COVER_EPS_R)
    relief = {label: relief_two_way(label, v_cover, fx, fz) for label in ('ground', '15m')}
    xmids = np.array([station_x(k) + RX_OFFSET / 2 for k in range(N_STATIONS)])

    masks = {label: (time >= lo) & (time <= hi) for label, (lo, hi) in VIEW.items()}
    vmax = {kind: max(np.max(np.abs(products[l, k, kind].real_bandpass[masks[l]]))
                      for l in ('ground', '15m') for k in range(N_STATIONS))
            for kind in ('total', 'contrast')}

    peaks, peak_corr, span = {}, {}, {}
    for label in ('ground', '15m'):
        lo, hi = PEAK_WIN[label]
        m = (time >= lo) & (time <= hi)
        peaks[label] = [float(time[m][np.argmax(np.abs(products[label, k, 'contrast'].complex_envelope[m]))])
                        for k in range(N_STATIONS)]
        peak_corr[label] = float(np.corrcoef(peaks[label], relief[label])[0, 1])
        span[label] = max(peaks[label]) - min(peaks[label])
    dyn_db = {label: float(20 * np.log10(
        max(np.max(np.abs(products[label, k, 'total'].real_bandpass[masks[label]])) for k in range(N_STATIONS))
        / max(np.max(np.abs(products[label, k, 'contrast'].real_bandpass[masks[label]])) for k in range(N_STATIONS))))
        for label in ('ground', '15m')}
    contrast_vmax = {label: max(np.max(np.abs(products[label, k, 'contrast'].real_bandpass[masks[label]]))
                                for k in range(N_STATIONS)) for label in ('ground', '15m')}
    gain_db = float(20 * np.log10(contrast_vmax['ground'] / contrast_vmax['15m']))

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharey=False, layout='constrained')
    for row, kind in enumerate(('total', 'contrast')):
        for col, label in enumerate(('ground', '15m')):
            ax = axes[row, col]
            lo, hi = VIEW[label]
            m = masks[label]
            t_edge = np.concatenate(([time[m][0] - (time[1] - time[0]) / 2],
                                     (time[m][:-1] + time[m][1:]) / 2,
                                     [time[m][-1] + (time[1] - time[0]) / 2]))
            for k in range(N_STATIONS):
                x = station_x(k)
                strip = products[label, k, kind].real_bandpass[m]
                ax.pcolormesh([x - STEP / 2, x + STEP / 2], t_edge,
                              strip[:, None], cmap='RdBu_r', shading='flat',
                              norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * vmax[kind],
                                                                vmin=-vmax[kind], vmax=vmax[kind]))
            ax.plot(xmids, relief[label], 'k-', lw=1.0,
                    label='true relief two-way time (no dispersion)')
            ax.set_xlim(station_x(0) - STEP / 2, station_x(N_STATIONS - 1) + RX_OFFSET + STEP / 2)
            ax.set_ylim(hi, lo)
            ax.set_title(f'{label}  {kind} (row symlog |max|={vmax[kind]:.2e})', fontsize=10)
            ax.set_xlabel('source x (m)')
            if col == 0:
                ax.set_ylabel('actual receiver time (ns)')
            ax.legend(fontsize=7, loc='lower right')
    fig.suptitle('Ground-coupled (z=12.025 m) vs 15 m airborne, same 36 m model/grid/materials; '
                 'V4.0.0 double, source-normalized 20-170 MHz; 13 adjacent 0.5 m columns, no interpolation')
    a.out.mkdir(parents=True)
    fig.savefig(a.out / 'ground_vs_15m_bscan.png', dpi=150)
    np.savez_compressed(a.out / 'comparison_arrays.npz', **arrays)
    completed = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'contract_sha256': sha256(cp),
               'code_sha256': sha256(__file__),
               'verification_sha256': sha256(a.capsule / 'completed_verification.json'),
               'completed_cases': len(c['groups']),
               'factor': 'antenna height only: z 27 m -> 12.025 m; archived 36 m model otherwise byte-identical inputs',
               'cover_eps_r_assumed': COVER_EPS_R,
               'relief_two_way_ns': {l: list(map(float, relief[l])) for l in relief},
               'contrast_envelope_peak_ns_per_station': peaks,
               'peak_time_vs_relief_pearson': peak_corr,
               'peak_span_ns': span,
               'contrast_vmax': contrast_vmax,
               'contrast_amplitude_gain_ground_over_15m_dB': gain_db,
               'dynamic_range_total_over_contrast_dB': dyn_db,
               'figure_vmax': vmax,
               'view_windows_ns': VIEW, 'peak_windows_ns': PEAK_WIN,
               'independent_direct_DFT_relative_L2': dft,
               'sum_solver_elapsed_s': sum(e['elapsed_s'] for e in completed),
               'maximum_owned_RSS_GiB': max(e['peak_owned_RSS_bytes'] for e in completed) / 2**30,
               'physical_attribution_certified': False, 'grid_convergence_certified': False,
               'scope': 'Antenna-height factor control; 13-station segment; 2D line source; not a real ground-antenna coupling model; no field claims.'}
    (a.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                                   allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'peak_corr': peak_corr, 'span_ns': span, 'gain_dB': gain_db,
                      'dyn_range_dB': dyn_db,
                      'peaks_ground': peaks['ground'], 'peaks_15m': peaks['15m']}, indent=1))


if __name__ == '__main__':
    main()
