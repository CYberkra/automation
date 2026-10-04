"""Source-normalized, phase-preserving comparisons of local V4 factor controls."""
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
from hs4_v4_factor_controls import STATIONS, reference_raw, ROOT, PROFILE
from analyze_hs4_height_wavefield import response
from hs4t2d_path_diagnostic import cover_index
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

WINDOWS = {'underground': (160, 220), 'early': (160, 180), 'late': (180, 220)}


def relative(a, b):
    denominator = np.linalg.norm(a)
    return None if denominator == 0 else float(np.linalg.norm(b-a)/denominator)


def archived_fine_raw(station, role):
    tag = 'base' if role == 'rough' else 'halfspace'
    study = ROOT/'artifacts/research_checks'/('2026-10-04_hs4_patch_finest' if station == 'centre' else '2026-10-04_hs4_station_grid_endpoints')
    name = tag if station == 'centre' else station+'_'+tag
    p = study/name/'profile.h5'
    verification = json.loads((study/'completed_verification.json').read_text('utf-8'))
    cp = study/'execution_contract.json'
    if verification['status'] != 'PASS' or not verification['completed'] or verification['contract_sha256'] != sha256(cp):
        raise ValueError('archived finer grid must have matching completed audit')
    group = next(g for g in verification['groups'] if g.get('id', g.get('group')) == name)
    digest = group.get('raw_sha256') if 'raw_sha256' in group else next(r['sha256'] for r in group['outputs'] if r['file'] == p.name)
    if sha256(p) != digest:
        raise ValueError('archived finer-grid native identity differs')
    return p


def compare(reference, candidate, time):
    result = {}
    for name, (lo, hi) in WINDOWS.items():
        mask = (time >= lo) & (time <= hi)
        a, b = reference.real_bandpass[mask], candidate.real_bandpass[mask]
        za, zb = reference.complex_envelope[mask], candidate.complex_envelope[mask]
        ca = float(np.linalg.norm(za)); cb = float(np.linalg.norm(zb))
        if ca == 0 or cb == 0:
            raise ValueError('zero complex event denominator')
        coherence = np.vdot(za, zb)/(ca*cb)
        peak_a = float(time[mask][np.argmax(abs(za))])
        peak_b = float(time[mask][np.argmax(abs(zb))])
        result[name] = {'signed_waveform_relative_L2': relative(a, b),
                        'complex_envelope_relative_L2': relative(za, zb),
                        'envelope_magnitude_relative_L2': relative(abs(za), abs(zb)),
                        'amplitude_norm_ratio': cb/ca, 'amplitude_norm_change_dB': float(20*np.log10(cb/ca)),
                        'complex_coherence_real': float(coherence.real), 'complex_coherence_imag': float(coherence.imag),
                        'reference_envelope_peak_ns': peak_a, 'candidate_envelope_peak_ns': peak_b,
                        'peak_shift_ns': peak_b-peak_a}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new output directory required')
    cp = args.capsule/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    v = json.loads((args.capsule/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (args.capsule/'execution.jsonl').read_text('utf-8').splitlines()]
    if (v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp)
            or events[-1]['status'] != 'COMPLETED' or events[-1]['traces'] != 11
            or events[-1]['verification_sha256'] != sha256(args.capsule/'completed_verification.json')):
        raise ValueError('complete identity-verified eleven-case capsule required')
    verify_official_runtime()
    verified = {g['id']: g for g in v['groups']}
    spectra, products, refs, audits = {}, {}, {}, {}
    for station in STATIONS:
        for role in ('rough', 'halfspace'):
            path = reference_raw(station, role)
            row = next(x for x in json.loads((path.parent/'audit.json').read_text('utf-8')) if x['file'] == path.name)
            if sha256(path) != row['sha256']:
                raise ValueError('historical receiver differs from original audit')
            refs[station, role] = response(path)
            audits[f'reference_{station}_{role}'] = {'raw_sha256': sha256(path), 'source_all501_valid': True}
    for group in c['groups']:
        name = group['id']; path = args.capsule/name/'profile.h5'
        if sha256(path) != verified[name]['raw_sha256']:
            raise ValueError('new raw receiver identity differs')
        r = response(path)
        old = refs[group['station'], group['role']]
        if (not np.array_equal(r.frequency, old.frequency) or r.source.spatial_scale != old.source.spatial_scale
                or r.source.quantity != old.source.quantity or r.source.units != old.source.units):
            raise ValueError('source-normalization or frequency semantics changed')
        spectra[name] = r
        audits[name] = {'raw_sha256': sha256(path), 'source_all501_valid': True,
                        'source_quantity': r.source.quantity, 'source_units': r.source.units,
                        'source_spatial_scale': r.source.spatial_scale, 'native_dt_s': r.source.dt}
    model = refs['centre', 'rough']
    visibility = {}
    for station in STATIONS:
        base = refs[station, 'rough'].response-refs[station, 'halfspace'].response
        avg = spectra[f'averaging_y_{station}_rough'].response-spectra[f'averaging_y_{station}_halfspace'].response
        flat = spectra[f'flat_{station}_rough'].response-refs[station, 'halfspace'].response
        for tag, array in (('baseline', base), ('averaging_y', avg), ('flat', flat)):
            key = f'{tag}_{station}'
            products[key] = reconstruct_time_response(replace(model, response=array), window='hann', zero_pad_factor=8)
        fine = [response(archived_fine_raw(station, role)) for role in ('rough', 'halfspace')]
        for r in fine:
            if r.source.quantity != model.source.quantity or r.source.spatial_scale != model.source.spatial_scale:
                raise ValueError('archived grid source normalization differs')
        products[f'spatial_fine_{station}'] = reconstruct_time_response(replace(model, response=fine[0].response-fine[1].response), window='hann', zero_pad_factor=8)
        for tag, raw_product in (('baseline', refs[station, 'rough']), ('averaging_y', spectra[f'averaging_y_{station}_rough'])):
            total = reconstruct_time_response(raw_product, window='hann', zero_pad_factor=8)
            p = products[f'{tag}_{station}']
            mask = (p.time*1e9 >= 160)&(p.time*1e9 <= 220)
            ratio = float(abs(p.complex_envelope[mask]).max()/abs(total.complex_envelope).max())
            visibility[f'{tag}_{station}'] = {'interface_event_over_full_total_peak': ratio,
                                             'interface_event_over_full_total_peak_dB': float(20*np.log10(ratio))}
    half_dt = spectra['dt_half_centre_rough'].response-spectra['dt_half_centre_halfspace'].response
    products['dt_half_centre'] = reconstruct_time_response(replace(model, response=half_dt), window='hann', zero_pad_factor=8)
    time = products['baseline_centre'].time*1e9
    metrics = {}
    for key, product in products.items():
        tag, station = key.rsplit('_', 1)
        if tag != 'baseline':
            metrics[key] = compare(products[f'baseline_{station}'], product, time)
    profile = np.genfromtxt(PROFILE, delimiter=',', names=True)
    shape = {}
    for tag in ('baseline', 'averaging_y', 'flat'):
        peaks = [float(time[(time >= 160)&(time <= 220)][np.argmax(abs(products[f'{tag}_{s}'].complex_envelope[(time >= 160)&(time <= 220)]))]) for s in STATIONS]
        shape[tag] = {'envelope_peak_ns_by_left_centre_right': peaks, 'three_station_peak_span_ns': max(peaks)-min(peaks)}
    frequencies = np.linspace(20e6, 170e6, 501)
    n = cover_index(frequencies)
    wavelengths = 299792458/(frequencies*n.real)
    mesh = {'minimum_phase_wavelength_m_in_20_170MHz': float(wavelengths.min()),
            'minimum_cells_per_phase_wavelength_at_2p5cm': float((wavelengths/.025).min()),
            'tau_over_actual_dt': 6.4567e-9/c['base_dt_s'],
            'note': 'Band-limited wavelength guideline, not convergence. Impulse has no finite maximum source frequency.'}
    result = {'status': 'COMPLETED_CONTROLLED_NUMERICAL_AND_GEOMETRY_DIAGNOSTIC',
              'code_sha256': sha256(__file__), 'contract_sha256': sha256(cp),
              'completed_cases': 11, 'processing_audit': audits, 'fixed_window_comparisons': metrics,
              'three_station_shapes': shape, 'mesh_guideline': mesh,
              'visibility_diagnostic': visibility,
              'physical_attribution_certified': False,
              'limitations': ['Sensitivity is not a bound on true error; averaging changes component assignment, not bulk cell map.',
                              'Three stations are not a complete B-scan.', 'No new 3D/finite-antenna or field validation.']}
    args.out.mkdir(parents=True)
    arrays = {'time_ns': time, 'frequency_hz': model.frequency}
    for key, p in products.items():
        arrays[key+'_complex_envelope'] = p.complex_envelope
        arrays[key+'_signed_waveform'] = p.real_bandpass
    for key, r in spectra.items():
        arrays[key+'_frequency_response'] = r.response
    np.savez_compressed(args.out/'comparison_arrays.npz', **arrays)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, layout='constrained')
    for column, station in enumerate(STATIONS):
        for tag in ('baseline', 'averaging_y', 'flat'):
            p = products[f'{tag}_{station}']
            mask = (time >= 140)&(time <= 240)
            axes[0, column].plot(time[mask], p.real_bandpass[mask], label=tag)
            axes[1, column].plot(time[mask], 2*abs(p.complex_envelope[mask]), label=tag)
        if station == 'centre':
            p = products['dt_half_centre']
            axes[0, column].plot(time[mask], p.real_bandpass[mask], '--', label='dt/2')
            axes[1, column].plot(time[mask], 2*abs(p.complex_envelope[mask]), '--', label='dt/2')
        axes[0, column].set_title(f'{station}; Tx x={STATIONS[station]:g} m')
        axes[1, column].set_xlabel('actual receiver time (ns)')
        for ax in axes[:, column]:
            ax.axvspan(160, 220, color='grey', alpha=.08); ax.grid(alpha=.2)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # Shared physical amplitude scales; never normalize individual traces.
    for row in axes:
        lo = min(ax.get_ylim()[0] for ax in row); hi = max(ax.get_ylim()[1] for ax in row)
        for ax in row:
            ax.set_ylim(lo, hi)
    axes[0, 0].set_ylabel('signed bandpass response (transfer units)')
    axes[1, 0].set_ylabel('2|complex envelope| (transfer units)')
    axes[0, 1].legend(fontsize=8)
    fig.suptitle('Local V4/double, 15 m; same spatial grid / Debye bulk; paired interface contrast')
    fig.savefig(args.out/'factor_control_ascans.png', dpi=130); plt.close(fig)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'mesh': mesh, 'metrics': metrics, 'shape': shape}, indent=2))


if __name__ == '__main__':
    main()
