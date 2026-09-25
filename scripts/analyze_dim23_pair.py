"""Paired-difference analysis for the dim23_pair_v1 batch (CPU/NumPy only).

Per dimension (2D TMx, 3D), compute TGT-BG paired differences of the official
SFCW direct frequency response and of the raw time traces, then compare the
two dimensions by arrival time, phase slope and NORMALIZED amplitude trend.
Cross-dimension absolute amplitudes are never compared (different source
normalization and geometric spreading); the optional 3D->2D spreading-factor
transform is exploratory with constants omitted. No solver is invoked.
"""
import argparse
import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import hilbert
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response

FREQ = np.linspace(20e6, 170e6, 501)
EXPECTED = {'2D': dict(nx_ny_nz=[1, 200, 780]), '3D': dict(nx_ny_nz=[200, 200, 780])}


def load_group(root, dim):
    paths = {k: Path(root) / f'2026-09-25_DIM23_{dim}_{k}' / f'DIM23_{dim}_{k}.h5' for k in ('BG', 'TGT')}
    for k, v in paths.items():
        assert v.is_file(), v
    src = {k: load_source(v) for k, v in paths.items()}
    rx = {k: load_receiver(v, receiver_path='name:measurement', component='Ex') for k, v in paths.items()}
    for k, v in paths.items():
        with h5py.File(v, 'r') as h:
            assert np.array_equal(h.attrs['nx_ny_nz'], EXPECTED[dim]['nx_ny_nz']), (dim, h.attrs['nx_ny_nz'])
            assert np.allclose(h.attrs['dx_dy_dz'], [.05, .05, .05], rtol=0, atol=1e-12)
            assert np.allclose(h[rx[k].path].parent.attrs['Position'][1:], [5.65, 36], rtol=0, atol=1e-12)
            assert h[rx[k].path].dtype == np.float64
    assert src['BG'].dt == rx['BG'].dt == src['TGT'].dt == rx['TGT'].dt
    assert np.array_equal(src['BG'].samples, src['TGT'].samples), f'{dim} source mismatch'
    assert all(np.all(np.isfinite(r.samples)) for r in rx.values())
    return paths, src, rx


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--runs-root', type=Path, default=Path('artifacts/simulations'))
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)

    data = {}
    for dim in ('2D', '3D'):
        paths, src, rx = load_group(a.runs_root, dim)
        dt = rx['BG'].dt
        n = min(len(rx['BG'].samples), int(np.floor(800e-9 / dt)) + 1)
        taper = (round(200e-9 / dt) - .25) / n
        spec = {}
        for k in ('BG', 'TGT'):
            r = direct_frequency_response(src[k], replace(rx[k], samples=rx[k].samples[:n]), FREQ,
                                          tail_taper_fraction=taper)
            assert np.all(r.source_valid) and np.all(np.isfinite(r.response))
            spec[k] = r.response
        t = rx['BG'].times[:n]
        diff_t = rx['TGT'].samples[:n] - rx['BG'].samples[:n]
        diff_f = spec['TGT'] - spec['BG']
        env = np.abs(hilbert(diff_t))
        data[dim] = dict(paths=paths, dt=dt, n=n, t=t, diff_t=diff_t, diff_f=diff_f, env=env)

    rows = {}
    for dim, d in data.items():
        f = FREQ
        mag = np.abs(d['diff_f'])
        phase = np.unwrap(np.angle(d['diff_f']))
        i = int(np.argmax(mag))
        # Normalized amplitude trend: shape only, per-dimension max == 1.
        norm_db = 20 * np.log10(mag / mag.max())
        sel = np.array([20, 50, 80, 110, 140, 170]) * 1e6
        idx = [int(np.argmin(abs(f - s))) for s in sel]
        rows[dim] = dict(
            dt_s=d['dt'], transformed_samples=d['n'], taper_fraction=d['n'] and (round(200e-9 / d['dt']) - .25) / d['n'],
            time_diff_peak_V_m=float(np.max(np.abs(d['diff_t']))),
            envelope_peak_time_ns=float(d['t'][int(np.argmax(d['env']))] * 1e9),
            spectral_peak_frequency_MHz=float(f[i] / 1e6),
            spectral_peak_phase_deg=float(np.degrees(np.angle(d['diff_f'][i]))),
            normalized_mag_dB_at_selected_MHz={f'{s / 1e6:g}': float(norm_db[j]) for s, j in zip(sel, idx)})

    # Cross-dimension shape comparisons (arrival, phase slope, normalized trend only).
    d2, d3 = data['2D'], data['3D']
    t_grid = d2['t']  # 2D is the coarser time base (larger dt)
    diff3_on_2d = np.interp(t_grid, d3['t'], d3['diff_t'])
    a2 = d2['diff_t'] / np.max(np.abs(d2['diff_t']))
    a3 = diff3_on_2d / np.max(np.abs(diff3_on_2d))
    xcorr = np.correlate(a2, a3, mode='full')
    lag_samples = int(np.argmax(xcorr) - (len(a3) - 1))
    norm = float(np.sqrt(np.sum(a2**2) * np.sum(a3**2)))
    cross = dict(
        envelope_peak_time_difference_ns=float(rows['2D']['envelope_peak_time_ns'] - rows['3D']['envelope_peak_time_ns']),
        waveform_max_normalized_correlation=float(np.max(xcorr) / norm),
        waveform_lag_ns=float(lag_samples * d2['dt'] * 1e9),
        spectral_peak_frequency_difference_MHz=float(rows['2D']['spectral_peak_frequency_MHz'] - rows['3D']['spectral_peak_frequency_MHz']))
    ph2 = np.unwrap(np.angle(d2['diff_f']))
    ph3 = np.unwrap(np.angle(d3['diff_f']))
    dphi = ph2 - ph3
    slope, intercept = np.polyfit(FREQ, dphi, 1)
    cross['phase_difference_linear_fit'] = dict(
        slope_rad_per_Hz=float(slope), equivalent_time_offset_ns=float(-slope / (2 * np.pi) * 1e9),
        residual_std_deg=float(np.degrees(np.std(dphi - (slope * FREQ + intercept)))))
    sel = np.array([20, 50, 80, 110, 140, 170]) * 1e6
    cross['normalized_mag_trend_difference_dB'] = {
        f'{s / 1e6:g}': float(rows['2D']['normalized_mag_dB_at_selected_MHz'][f'{s / 1e6:g}']
                              - rows['3D']['normalized_mag_dB_at_selected_MHz'][f'{s / 1e6:g}']) for s in sel}
    n2 = np.abs(d2['diff_f']) / np.abs(d2['diff_f']).max()
    n3 = np.abs(d3['diff_f']) / np.abs(d3['diff_f']).max()
    cross['normalized_spectrum_shape_correlation'] = float(np.corrcoef(n2, n3)[0, 1])

    # Exploratory: far-field 3D->2D frequency factor 1/sqrt(i*2*pi*f), constants omitted
    # (Bleistein/Forbriger form); near-field caveats recorded in the report.
    g = 1 / np.sqrt(1j * 2 * np.pi * FREQ)
    t3 = np.abs(d3['diff_f'] * g)
    t3 /= t3.max()
    exploratory = dict(form='D_2D_equiv(f) ∝ D_3D(f) / sqrt(i*2*pi*f); constants omitted, shape only',
                       transformed_vs_2D_shape_correlation=float(np.corrcoef(n2, t3)[0, 1]),
                       transformed_normalized_mag_dB_at_selected_MHz={
                           f'{s / 1e6:g}': float(20 * np.log10(t3[int(np.argmin(abs(FREQ - s)))])) for s in sel})

    result = dict(batch='dim23_pair_v1', per_dimension=rows, cross_dimension_shape_only=cross,
                  exploratory_3d_to_2d_transform=exploratory,
                  inputs={str(v): hashlib.sha256(v.read_bytes()).hexdigest()
                          for d in data.values() for v in d['paths'].values()},
                  frequency_Hz='linspace(20e6,170e6,501)', taper_ns=200,
                  cross_dimension_absolute_amplitude_compared=False,
                  grid_convergence_certified=False, physical_accuracy_threshold=None,
                  solver_invoked=False)
    (a.output / 'results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(a.output / 'arrays.npz', frequency_Hz=FREQ,
                        diff_2D=d2['diff_f'], diff_3D=d3['diff_f'],
                        time_2D_s=d2['t'], diff_t_2D=d2['diff_t'],
                        time_3D_s=d3['t'], diff_t_3D=d3['diff_t'],
                        transformed_3d_to_2d_normalized=t3)
    with (a.output / 'difference_spectra.csv').open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['frequency_Hz', 'diff_2D_re', 'diff_2D_im', 'diff_3D_re', 'diff_3D_im'])
        w.writerows(zip(FREQ, d2['diff_f'].real, d2['diff_f'].imag, d3['diff_f'].real, d3['diff_f'].imag))
    print(json.dumps(result))


if __name__ == '__main__':
    main()
