"""Per-station 36m/108m comparison and dense B-scan figure for the dense capsule.

All 13 stations pair with their archived 36 m raws (identity checked against
the capsule contract and the historical audit.json rows). Fixed-window metrics
reuse the same source-normalized SFCW chain as the six-case large-domain
analysis. The figure draws all 13 adjacent 0.5 m columns -- a true dense
B-scan over the instrumented segment, no masking needed.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hs_capsule_identity import sha256
from hs4_large_domain_dense import N_STATIONS, SRC_X0, STEP, SHIFT, station_x
from analyze_hs4_height_wavefield import response
from analyze_hs4_v4_factor_controls import compare
from check_hs4_v4_factor_evidence import direct_response
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
T_LO, T_HI = 80.0, 240.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--failed-attempt', type=Path,
                    help='optional preserved failed capsule whose group raws must be bitwise identical')
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
        raise ValueError('completed 26-trace audited capsule required')
    verify_official_runtime()
    rows = {g['id']: g for g in v['groups']}
    responses, dft, arrays = {}, {}, {}
    for g in c['groups']:
        p = Path(g['input']).with_suffix('.h5')
        old = Path(g['reference_raw'])
        if sha256(p) != rows[g['id']]['raw_sha256'] or sha256(old) != g['reference_raw_sha256']:
            raise ValueError('new/historical raw identity differs')
        for label, raw in [('36m', old), ('108m', p)]:
            r = response(raw)
            responses[label, g['station_index'], g['role']] = r
            exact = direct_response(raw, r.frequency[TAKE])
            error = float(np.linalg.norm(exact - r.response[TAKE]) / np.linalg.norm(exact))
            if error > 1e-9:
                raise ValueError('independent actual-source direct DFT failed')
            dft[label + '_' + g['id']] = error
    model = responses['36m', 6, 'rough']
    metrics, peaks = {}, {'36m': [], '108m': []}
    products = {}
    for k in range(N_STATIONS):
        for label in ('36m', '108m'):
            rough = responses[label, k, 'rough']
            half = responses[label, k, 'halfspace']
            if (not np.array_equal(rough.frequency, model.frequency)
                    or rough.source.quantity != model.source.quantity
                    or rough.source.spatial_scale != model.source.spatial_scale):
                raise ValueError('source/frequency conventions differ')
            total = reconstruct_time_response(rough, window='hann', zero_pad_factor=8)
            contrast = reconstruct_time_response(replace(model, response=rough.response - half.response),
                                                 window='hann', zero_pad_factor=8)
            products[label, k, 'total'] = total
            products[label, k, 'contrast'] = contrast
            arrays[f'{label}_s{k:02d}_total_signed'] = total.real_bandpass
            arrays[f'{label}_s{k:02d}_contrast_signed'] = contrast.real_bandpass
            arrays[f'{label}_s{k:02d}_contrast_envelope'] = contrast.complex_envelope
        time = products['36m', k, 'contrast'].time * 1e9
        metrics[f's{k:02d}'] = {
            'contrast': compare(products['36m', k, 'contrast'], products['108m', k, 'contrast'], time),
            'total': compare(products['36m', k, 'total'], products['108m', k, 'total'], time)}
        for label in ('36m', '108m'):
            m = (time >= 160) & (time <= 220)
            peaks[label].append(float(time[m][np.argmax(np.abs(products[label, k, 'contrast'].complex_envelope[m]))]))
    time = products['36m', 6, 'total'].time * 1e9
    arrays['time_ns'] = time
    arrays['frequency_hz'] = model.frequency
    mask = (time >= T_LO) & (time <= T_HI)
    vmax = {kind: max(np.max(np.abs(products[l, k, kind].real_bandpass[mask]))
                      for l in ('36m', '108m') for k in range(N_STATIONS))
            for kind in ('total', 'contrast')}
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharey=True, layout='constrained')
    t_edge = np.concatenate(([time[mask][0] - (time[1] - time[0]) / 2],
                             (time[mask][:-1] + time[mask][1:]) / 2,
                             [time[mask][-1] + (time[1] - time[0]) / 2]))
    for row, kind in enumerate(('total', 'contrast')):
        for col, label in enumerate(('36m', '108m')):
            ax = axes[row, col]
            shift = 0.0 if label == '36m' else SHIFT
            x0, x1 = (10.0, 26.0) if label == '36m' else (10.0 + SHIFT, 26.0 + SHIFT)
            ax.axvspan(x0, x1, color='0.92', zorder=0)
            for k in range(N_STATIONS):
                x = station_x(k) + shift
                strip = products[label, k, kind].real_bandpass[mask]
                ax.pcolormesh([x - STEP / 2, x + STEP / 2], t_edge, strip[:, None],
                              cmap='RdBu_r', shading='flat',
                              norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * vmax[kind],
                                                                vmin=-vmax[kind], vmax=vmax[kind]))
            ax.set_xlim(x0, x1)
            ax.set_ylim(T_HI, T_LO)
            ax.set_title(f'{label}  {kind} (fixed symlog |max|={vmax[kind]:.2e})', fontsize=10)
            ax.set_xlabel('source x (m)')
            if col == 0:
                ax.set_ylabel('actual receiver time (ns)')
    fig.suptitle('Dense 13-station B-scan, V4.0.0 double, source-normalized 20-170 MHz; '
                 'all 13 adjacent 0.5 m columns computed, no interpolation')
    a.out.mkdir(parents=True)
    fig.savefig(a.out / 'dense_bscan.png', dpi=150)
    np.savez_compressed(a.out / 'comparison_arrays.npz', **arrays)
    completed = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    bitwise = None
    if a.failed_attempt is not None:
        prior = {}
        for line in (a.failed_attempt / 'execution.jsonl').read_text('utf-8').splitlines():
            e = json.loads(line)
            if e['status'] == 'COMPLETED' and 'group' in e:
                prior[e['group']] = e['raw_sha256']
        bitwise = {}
        for g in c['groups']:
            current = rows[g['id']]['raw_sha256']
            if g['id'] not in prior:
                raise ValueError('failed attempt lacks group ' + g['id'])
            bitwise[g['id']] = bool(prior[g['id']] == current)
        if not all(bitwise.values()):
            raise ValueError('failed-attempt raws differ bitwise: '
                             + ', '.join(k for k, ok in bitwise.items() if not ok))
    span = {label: max(p) - min(p) for label, p in peaks.items()}
    worst = max(max(w['signed_waveform_relative_L2'] for w in m[kind].values())
                for m in metrics.values() for kind in ('contrast', 'total'))
    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'contract_sha256': sha256(cp),
               'code_sha256': sha256(__file__),
               'verification_sha256': sha256(a.capsule / 'completed_verification.json'),
               'completed_cases': len(c['groups']),
               'per_station_metrics': metrics,
               'worst_fixed_window_signed_relative_L2_all_stations': worst,
               'contrast_envelope_peak_ns_per_station': peaks,
               'peak_span_ns': span,
               'independent_direct_DFT_relative_L2': dft,
               'sum_solver_elapsed_s': sum(e['elapsed_s'] for e in completed),
               'maximum_owned_RSS_GiB': max(e['peak_owned_RSS_bytes'] for e in completed) / 2**30,
               'failed_attempt_bitwise_equal': bitwise,
               'figure_vmax': vmax,
               'physical_attribution_certified': False, 'grid_convergence_certified': False,
               'scope': 'Dense 13-station segment in enlarged X domain; 2D line source; unchanged top/bottom; no full-line/3D/field claims.'}
    (a.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                                   allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'worst': worst, 'span_ns': span,
                      'peaks_36m': peaks['36m'], 'peaks_108m': peaks['108m']}, indent=1))


if __name__ == '__main__':
    main()
