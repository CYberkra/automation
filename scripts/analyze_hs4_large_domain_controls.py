"""Same-chain comparison and independent DFT checks for 108m domain controls."""
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
from hs4_large_domain_controls import STATIONS
from analyze_hs4_height_wavefield import response
from analyze_hs4_v4_factor_controls import compare
from check_hs4_v4_factor_evidence import direct_response
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    cp = a.capsule/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    v = json.loads((a.capsule/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (a.capsule/'execution.jsonl').read_text('utf-8').splitlines()]
    if (v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp)
            or events[-1]['status'] != 'COMPLETED' or events[-1]['traces'] != 6
            or events[-1]['verification_sha256'] != sha256(a.capsule/'completed_verification.json')):
        raise ValueError('completed six-case audited capsule required')
    verify_official_runtime()
    rows = {g['id']: g for g in v['groups']}
    responses, native, dft, products, metrics, totals, independent = {}, {}, {}, {}, {}, {}, {}
    take = np.array([0, 83, 167, 250, 333, 417, 500])
    arrays = {}
    for g in c['groups']:
        p = a.capsule/g['id']/'profile.h5'; old = Path(g['reference_raw'])
        if sha256(p) != rows[g['id']]['raw_sha256'] or sha256(old) != g['reference_raw_sha256']:
            raise ValueError('new/historical raw identity differs')
        old_row = next(r for r in json.loads((old.parent/'audit.json').read_text('utf-8')) if r['file'] == old.name)
        if old_row['sha256'] != sha256(old):
            raise ValueError('historical original audit differs')
        for label, raw in [('36m', old), ('108m', p)]:
            key = label+'_'+g['id']; r = response(raw); responses[label, g['station'], g['role']] = r
            exact = direct_response(raw, r.frequency[take])
            error = float(np.linalg.norm(exact-r.response[take])/np.linalg.norm(exact))
            if error > 1e-9:
                raise ValueError('independent actual-source direct DFT failed')
            dft[key] = error
            arrays[key+'_frequency_response'] = r.response
            products[key] = reconstruct_time_response(r, window='hann', zero_pad_factor=8)
            with h5py.File(raw) as h:
                y = h['rxs/rx1/Ey'][:]; dt = float(h.attrs['dt'])
            native[label, g['station'], g['role']] = y
        total_a, total_b = products['36m_'+g['id']], products['108m_'+g['id']]
        totals[g['id']] = compare(total_a, total_b, total_a.time*1e9)
    model = responses['36m', 'centre', 'rough']
    for station in STATIONS:
        for label in ('36m', '108m'):
            rough, half = responses[label, station, 'rough'], responses[label, station, 'halfspace']
            if (not np.array_equal(rough.frequency, model.frequency)
                    or rough.source.quantity != model.source.quantity
                    or rough.source.spatial_scale != model.source.spatial_scale):
                raise ValueError('source/frequency conventions differ')
            diff = rough.response-half.response
            key = label+'_'+station
            products[key] = reconstruct_time_response(replace(model, response=diff), window='hann', zero_pad_factor=8)
            arrays[key+'_frequency_response'] = diff
            arrays[key+'_signed_waveform'] = products[key].real_bandpass
            arrays[key+'_complex_envelope'] = products[key].complex_envelope
        time = products['36m_'+station].time*1e9
        metrics[station] = compare(products['36m_'+station], products['108m_'+station], time)
        errors = []
        for window, (lo, hi) in {'underground': (160, 220), 'early': (160, 180), 'late': (180, 220)}.items():
            mask = (time >= lo)&(time <= hi)
            x, y = arrays['36m_'+station+'_signed_waveform'][mask], arrays['108m_'+station+'_signed_waveform'][mask]
            value = float(np.sqrt(sum(float(z)*float(z) for z in y-x)/sum(float(z)*float(z) for z in x)))
            errors.append(abs(value-metrics[station][window]['signed_waveform_relative_L2']))
        independent[station] = max(errors)
    arrays['time_ns'] = time; arrays['frequency_hz'] = model.frequency
    raw_metrics = {}
    for station in STATIONS:
        base = native['36m', station, 'rough']-native['36m', station, 'halfspace']
        wide = native['108m', station, 'rough']-native['108m', station, 'halfspace']
        raw_time = np.arange(len(base))*dt*1e9
        raw_metrics[station] = {}
        for name, (lo, hi) in {'underground': (160, 220), 'late': (325, 400), 'tail': (400, 600)}.items():
            mask = (raw_time >= lo)&(raw_time <= hi)
            den = float(np.linalg.norm(base[mask]))
            raw_metrics[station][name] = float(np.linalg.norm(wide[mask]-base[mask])/den) if den else None
    shape = {}
    for label in ('36m', '108m'):
        peaks = [metrics[s]['underground'][('reference' if label == '36m' else 'candidate')+'_envelope_peak_ns'] for s in STATIONS]
        shape[label] = {'peak_ns_left_centre_right': peaks, 'span_ns': max(peaks)-min(peaks)}
    completed = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    result = {'status': 'COMPLETED_SIX_108M_DOMAIN_CONTROLS', 'contract_sha256': sha256(cp),
              'code_sha256': sha256(__file__), 'verification_sha256': sha256(a.capsule/'completed_verification.json'),
              'completed_cases': 6, 'fixed_window_interface_comparisons': metrics,
              'fixed_window_total_response_comparisons': totals, 'raw_time_interface_relative_L2': raw_metrics,
              'three_station_shapes': shape, 'independent_direct_DFT_relative_L2': dft,
              'independent_scalar_metric_max_abs_difference': independent,
              'sum_solver_elapsed_s': sum(e['elapsed_s'] for e in completed),
              'maximum_owned_RSS_GiB': max(e['peak_owned_RSS_bytes'] for e in completed)/2**30,
              'physical_attribution_certified': False, 'grid_convergence_certified': False,
              'scope': 'Only X enlarged, with old ROI/edge continuation and actual source equal; 3 anchors, not full B-scan; unchanged top/bottom; no3D/field claims.'}
    if max(independent.values()) > 1e-12:
        raise ValueError('independent scalar metrics differ')
    a.out.mkdir(parents=True)
    np.savez_compressed(a.out/'comparison_arrays.npz', **arrays)
    fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True, layout='constrained')
    for col, station in enumerate(STATIONS):
        for label in ('36m', '108m'):
            p = products[label+'_'+station]; mask = (time >= 140)&(time <= 240)
            axes[0, col].plot(time[mask], p.real_bandpass[mask], label=label,
                              linestyle='-' if label == '36m' else '--')
            axes[1, col].plot(time[mask], 2*abs(p.complex_envelope[mask]), label=label,
                              linestyle='-' if label == '36m' else '--')
        axes[0, col].set_title(station)
        axes[1, col].set_xlabel('actual receiver time (ns)')
        for ax in axes[:, col]:
            ax.axvspan(160, 220, alpha=.08, color='grey'); ax.grid(alpha=.2)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    for row in axes:
        lo = min(ax.get_ylim()[0] for ax in row); hi = max(ax.get_ylim()[1] for ax in row)
        for ax in row:
            ax.set_ylim(lo, hi)
    axes[0, 0].set_ylabel('signed interface contrast ((V/m)/A)')
    axes[1, 0].set_ylabel('2|complex envelope| ((V/m)/A)')
    axes[0, 1].legend()
    fig.suptitle('V4.0.0 native double; 15m altitude; X36 vs108m; same grid/material/PML/local geometry')
    fig.savefig(a.out/'large_domain_ascans.png', dpi=140); plt.close(fig)
    (a.out/'summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
