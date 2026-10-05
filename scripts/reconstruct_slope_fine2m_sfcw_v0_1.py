"""Reconstruct the archived 26-trace slope experiment using official SFCW.

Read-only native inputs; no solver or field fitting. Source-normalised Ey is
not antenna-port S21. Rectangular and Hann branches share plotting ranges.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

import hs4_slope_bscan_fine2m_v0_1 as sim
from research_operator_contract import apply_configuration

PROCESSING_SHA = 'adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b'
FREQUENCY = np.linspace(20e6, 170e6, 501)
WINDOWS = ('rectangular', 'hann')
C = 299792458.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(a, b):
    den = np.linalg.norm(b)
    if den == 0:
        raise ValueError('zero reference norm')
    return float(np.linalg.norm(a-b)/den)


def main(study, out):
    if out.exists():
        raise ValueError('fresh output required; archived analysis cannot be overwritten')
    if sha(sf.__file__) != PROCESSING_SHA:
        raise ValueError('reviewed official V4 SFCW processing required')
    contract = study/'execution_contract.json'
    c = json.loads(contract.read_text('utf-8'))
    audit = json.loads((study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha(contract):
        raise ValueError('completed native identity required')
    identities = {g['id']: g['raw_sha256'] for g in audit['groups']}
    stations = np.array(c['geometry']['stations_tx_x_m'])
    count = len(stations)
    spectra, tapered, raw_mats, records, sources, seen = {}, {}, {}, [], [], set()
    max_dft_error = 0.
    for g in c['groups']:
        raw = study/g['id']/'profile.h5'
        if sha(raw) != identities[g['id']] or sha(raw.with_suffix('.in')) != g['input_sha256']:
            raise ValueError('archived native/input identity mismatch')
        key = (g['case'], g['station_index'])
        if key in seen:
            raise ValueError('duplicate case/station')
        seen.add(key)
        source = sf.load_source(raw)
        receiver = sf.load_receiver(raw, '/rxs/rx1', 'Ey')
        with h5py.File(raw, 'r') as h:
            spacing = np.array(h.attrs['dx_dy_dz'])
            if str(h.attrs['gprMax']) != '4.0.0' or h['rxs/rx1/Ey'].dtype != np.float64:
                raise ValueError('native V4.0.0 float64 required')
            for path, position in [('srcs/src1', g['tx_m']), ('rxs/rx1', g['receivers_m'][0])]:
                np.testing.assert_array_equal(h[path].attrs['GridPosition'], np.rint(np.array(position)/spacing).astype(int))
        if source.spatial_scale != spacing[1] or source.spatial_scale != .05:
            raise ValueError('native y-current element length differs')
        if not np.isfinite(source.samples).all() or not np.isfinite(receiver.samples).all():
            raise ValueError('nonfinite native values')
        if records:
            np.testing.assert_array_equal(source.samples, sources[0].samples)
            np.testing.assert_array_equal(receiver.times, native_time)
            if source.time_offset != sources[0].time_offset or source.dt != sources[0].dt:
                raise ValueError('source clocks differ')
        native_time = receiver.times
        product = sf.direct_frequency_response(source, receiver, FREQUENCY, tail_taper_fraction=0)
        fraction = (round(20e-9/receiver.dt)-.25)/len(receiver.samples)
        sensitivity = sf.direct_frequency_response(source, receiver, FREQUENCY, tail_taper_fraction=fraction)
        if not product.source_valid.all() or not sensitivity.source_valid.all():
            raise ValueError('unsupported source frequency')
        product = replace(product, response=product.response/source.spatial_scale)
        sensitivity = replace(sensitivity, response=sensitivity.response/source.spatial_scale)
        take = np.array([0, 83, 167, 250, 333, 417, 500])
        xf = source.dt*np.exp(-2j*np.pi*FREQUENCY[take, None]*source.times)@source.samples
        yf = receiver.dt*np.exp(-2j*np.pi*FREQUENCY[take, None]*receiver.times)@receiver.samples
        error = relative(product.response[take], yf/xf/source.spatial_scale)
        max_dft_error = max(max_dft_error, error)
        spectra.setdefault(g['case'], np.zeros((501, count), complex))[:, g['station_index']] = product.response
        tapered.setdefault(g['case'], np.zeros((501, count), complex))[:, g['station_index']] = sensitivity.response
        raw_mats.setdefault(g['case'], np.zeros((len(receiver.samples), count)))[:, g['station_index']] = receiver.samples
        records.append(dict(id=g['id'], raw_sha256=sha(raw), input_sha256=sha(raw.with_suffix('.in')),
            independent_DFT_relative_L2=error, source_length_m=source.spatial_scale,
            source_peak_ns=float(source.times[np.argmax(abs(source.samples))]*1e9),
            source_peak_current_moment_Am=float(np.max(abs(source.samples))*source.spatial_scale),
            historical_frequency_normalisation_over_correct=source.spatial_scale/g['spacing_m']))
        sources.append(source)
    if seen != {(case, j) for case in ('slope_rough', 'slope_fullcover') for j in range(count)}:
        raise ValueError('incomplete matched case/station set')
    if max_dft_error > 1e-9:
        raise ValueError('independent source/receiver transform mismatch')
    source = sources[0]
    iso_spectrum = spectra['slope_rough']-spectra['slope_fullcover']
    iso_tapered = tapered['slope_rough']-tapered['slope_fullcover']
    raw_iso = raw_mats['slope_rough']-raw_mats['slope_fullcover']
    ntail = round(20e-9/source.dt)
    tail_rms = np.sqrt(np.mean(raw_iso[-ntail:]**2, axis=0))/np.max(abs(raw_iso), axis=0)
    tail_last = abs(raw_iso[-1])/np.max(abs(raw_iso), axis=0)
    tail_spectral = [relative(iso_tapered[:, j], iso_spectrum[:, j]) for j in range(count)]
    # Baseline deliberately uses no tail taper. Tail diagnostics are not a
    # proof that every weak path has completely decayed within the 200 ns run.
    zsurf = np.rint(sim.surface_z(stations)/sim.DX)*sim.DX
    ziface = np.rint(sim.interface_z(stations)/sim.DX)*sim.DX
    w = 2*np.pi*95e6
    eps = 11+.5/(1+1j*w*6.4567e-9)-1j*.001/(w*8.854187817e-12)
    approximate_interface_ns = (2*sim.AGL+2*(zsurf-ziface)*np.sqrt(eps.real))/C*1e9
    profiles, arrays, metrics = {}, {}, {}
    max_inverse_error = 0.
    max_linear_error = 0.
    for window in WINDOWS:
        full = {}
        for case in ('slope_rough', 'slope_fullcover'):
            matrix_product = replace(product, response=spectra[case])
            tr = sf.reconstruct_time_response(matrix_product, window=window, zero_pad_factor=8, time_shift=0)
            full[case] = tr
            arrays[f'{window}_{case}_real'] = tr.real_bandpass
            arrays[f'{window}_{case}_complex_envelope'] = tr.complex_envelope
            arrays[f'{window}_weights'] = tr.weights
        target_tr = sf.reconstruct_time_response(replace(product, response=iso_spectrum), window=window, zero_pad_factor=8)
        # Independently sum physical-frequency tones at selected delays.
        indices = np.array([0, 20, 50, 90, 140, 220, 239])
        direct = (2/len(FREQUENCY))*np.real(np.exp(2j*np.pi*tr.time[indices, None]*FREQUENCY)@(iso_spectrum[:, count//2]*tr.weights))
        inverse_error = relative(target_tr.real_bandpass[indices, count//2], direct)
        max_inverse_error = max(max_inverse_error, inverse_error)
        difference = full['slope_rough'].real_bandpass-full['slope_fullcover'].real_bandpass
        linear_error = relative(difference, target_tr.real_bandpass)
        max_linear_error = max(max_linear_error, linear_error)
        arrays[f'{window}_ideal_real'] = target_tr.real_bandpass
        arrays[f'{window}_ideal_complex_envelope'] = target_tr.complex_envelope
        time_ns = tr.time*1e9
        mask = time_ns <= 200
        view_time = time_ns[mask]
        total = full['slope_rough'].real_bandpass[mask]
        ideal = target_tr.real_bandpass[mask]
        ops, diagnostics = {}, {}
        for cid, name in [('B3_G1_BG', 'mean_full'), ('B5_G1_BG', 'svd_rank2')]:
            try:
                result = apply_configuration(total, cid)
                ops[name] = result['output']
                diagnostics[name] = {'status': 'available', 'steps': [s['diagnostics'] for s in result['steps']]}
            except ValueError as exc:
                ops[name] = None
                diagnostics[name] = {'status': 'unavailable', 'reason': str(exc)}
        win = (view_time[:, None] >= approximate_interface_ns[None, :]-8) & (view_time[:, None] <= approximate_interface_ns[None, :]+12)
        def wrms(v):
            return np.sqrt((v*v*win).sum(axis=0)/win.sum(axis=0))
        row = {'independent_inverse_relative_L2': inverse_error, 'complex_subtraction_linearity_relative_L2': linear_error,
               'total_window_rms_over_ideal_median': float(np.median(wrms(total)/wrms(ideal))), 'operator_diagnostics': diagnostics}
        for name, value in ops.items():
            if value is None:
                row[name] = {'status': 'unavailable'}
                continue
            corr = (value*ideal*win).sum(axis=0)/np.sqrt((value*value*win).sum(axis=0)*(ideal*ideal*win).sum(axis=0))
            row[name] = {'window_rms_over_ideal_median': float(np.median(wrms(value)/wrms(ideal))),
                         'window_corr_with_ideal_median': float(np.median(corr))}
            arrays[f'{window}_{name}_view_real'] = value
        profiles[window] = {'total': total, 'ideal': ideal, **ops}
        metrics[window] = row
    if max_inverse_error > 1e-9 or max_linear_error > 1e-9:
        raise ValueError('inverse transform or matched complex subtraction mismatch')
    out.mkdir(parents=True)
    arrays.update(frequency_Hz=FREQUENCY, sfcw_time_ns=time_ns, view_time_ns=view_time,
        native_time_ns=native_time*1e9, stations_tx_x_m=stations,
        approximate_interface_ns=approximate_interface_ns, source_excitation=source.samples,
        source_time_ns=source.times*1e9, source_spectrum=product.source_spectrum,
        spectrum_slope_rough=spectra['slope_rough'], spectrum_slope_fullcover=spectra['slope_fullcover'],
        spectrum_ideal=iso_spectrum, raw_slope_rough=raw_mats['slope_rough'], raw_slope_fullcover=raw_mats['slope_fullcover'])
    np.savez_compressed(out/'sfcw_profiles.npz', **arrays)
    scale_total = max(np.max(abs(v['total'])) for v in profiles.values())*.005
    scale_target = max(np.max(abs(v['ideal'])) for v in profiles.values())
    make_plots(out, stations, view_time, approximate_interface_ns, profiles, scale_total, scale_target)
    summary = dict(status='PASS_SFCW_PROCESSING_NOT_FIELD_VALIDATION', calls_solver=False, calls_training=False,
        field_fitted=False, script_sha256=sha(__file__), official_processing_sha256=sha(sf.__file__),
        source_contract_sha256=sha(contract), native_records=records,
        frequency=dict(start_Hz=20e6, stop_Hz=170e6, step_Hz=.3e6, count=501),
        reconstruction=dict(method='official direct Y/X/native source length; retained complex phase',
            source_length_m=.05, response_units='(V/m)/(A*m), not antenna-port S21',
            windows=list(WINDOWS), window_normalisation='unit mean', zero_pad_factor=8,
            source_floor_db=-100, source_band_min_relative_db=float(np.min(20*np.log10(abs(product.source_spectrum)/np.max(abs(product.source_spectrum))))),
            tail_taper_ns=0, time_shift_ns=0, source_peak_delay_removed_by_complex_division=True,
            signed_product='official real_bandpass=2*Re(complex_envelope*exp(2j*pi*20MHz*t))',
            gain_applied=False, operator_input_domain='signed SFCW time_real, 0-200ns crop',
            time_step_ns=float(np.diff(time_ns)[0]), periodic_delay_ns=1/.3e6*1e9,
            native_record_end_ns=float(native_time[-1]*1e9), historical_instrument_window_unknown=True),
        checks=dict(independent_DFT_relative_L2_max=max_dft_error,
            independent_inverse_relative_L2_max=max_inverse_error,
            complex_subtraction_linearity_relative_L2_max=max_linear_error),
        tail=dict(isolated_last20ns_rms_over_own_peak_per_station=tail_rms.tolist(),
            isolated_last_sample_over_own_peak_per_station=tail_last.tolist(),
            isolated_20ns_taper_complex_spectrum_relative_L2_per_station=tail_spectral,
            completeness='Measured residual and taper sensitivity only; no missing-path recovery or full-tail certificate'),
        approximate_interface_ns_per_station=approximate_interface_ns.tolist(),
        event_window='95MHz vertical two-way approximation [-8,+12]ns, no source delay; bistatic/refraction peak truth not certified',
        metrics=metrics, metrics_role='processing diagnostics, not physical acceptance or training labels',
        plot=dict(cmap='gray', negative='black', zero='midgray', positive='white',
            total_symmetric_limit=scale_total, target_and_operators_shared_symmetric_limit=scale_target,
            interpolation='nearest; 13 measured simulation stations, no added traces',
            clipping='total intentionally clipped at 0.5% of common total peak; ideal/mean/SVD share one range across both windows'))
    summary['products'] = {p.name: sha(p) for p in out.iterdir() if p.is_file()}
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'out': str(out), 'status': summary['status'], 'checks': summary['checks'],
        'metrics': {w: {k: v for k, v in m.items() if k != 'operator_diagnostics'} for w, m in metrics.items()}}, indent=2))


def make_plots(out, stations, times, predicted, profiles, scale_total, scale_target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    labels = {'total': '总场（0.5%截幅）', 'ideal': '理想参考（复谱匹配背景差分）',
              'mean_full': '均值道去背景', 'svd_rank2': 'SVD rank-2去背景'}
    def panel(fig, ax, values, name, window):
        if values is None:
            ax.text(.5, .5, '算子不可用', ha='center'); return
        limit = scale_total if name == 'total' else scale_target
        step = float(np.diff(times)[0])
        im = ax.imshow(values, cmap='gray', aspect='auto', interpolation='nearest', vmin=-limit, vmax=limit,
            extent=[stations[0]-.25, stations[-1]+.25, times[-1]+step/2, times[0]-step/2])
        ax.plot(stations, predicted, 'k--', lw=.8, label='95MHz垂直路径近似')
        ax.set(title=f'{window}: {labels[name]}', xlabel='发射站位 x(m)', ylabel='源归一化延迟(ns)', ylim=(200, 0))
        ax.legend(fontsize=6, loc='lower right')
        fig.colorbar(im, ax=ax, shrink=.8, label='源电流矩归一化场响应')
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), layout='constrained')
    for row, window in enumerate(WINDOWS):
        for col, name in enumerate(labels):
            panel(fig, axes[row, col], profiles[window][name], name, window)
    fig.suptitle('真实SFCW重建：20–170MHz / 0.3MHz / 501频点；2m航高、1.25cm网格、13站位\n'
        '无增益；两种窗共用幅值范围（总场另用截幅范围）；源延迟已由Y/X去除；窗口为机制分支，非实测设置认证', fontsize=12)
    fig.savefig(out/'sfcw_windows_gray.png', dpi=150); plt.close(fig)
    for window in WINDOWS:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
        for ax, name in zip(axes.flat, labels):
            panel(fig, ax, profiles[window][name], name, window)
        fig.suptitle(f'SFCW {window} 灰度图：20–170MHz，2m航高；无增益\n理想参考/均值道/SVD共用幅值范围；虚线是近似路径预测', fontsize=12)
        fig.savefig(out/f'sfcw_{window}_gray.png', dpi=145); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    main(args.study, args.out)
