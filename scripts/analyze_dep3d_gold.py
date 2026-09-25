"""Paired-difference analysis for the dep3d_gold_v1 batch (CPU/NumPy only, no solver).

Adapted from scripts/analyze_dim23_pair.py per docs/research/2026-09-26_dep3d_gold_design.md §9:
four TGT-BG paired differences (B2D5CM, B2DANISO, DEP3D_5CM, DEP3D_ANISO), each with the
official SFCW direct frequency response on 501 points 20-170 MHz (200 ns tail window as the
main result, 400 ns tail window as a robustness check, 440-600 ns real-waveform diagnostic
window). Comparisons: 3D internal dz-refinement convergence (DEP3D_ANISO vs DEP3D_5CM),
same-grid-spacing cross-dimension shape comparisons (B2D5CM vs DEP3D_5CM, B2DANISO vs
DEP3D_ANISO), an exploratory far-field 3D->2D transform, and a shape-only comparison to the
DEP 2D ZFINE2 chain (ZFINE2_200 in 2026-09-25_deep_zfine2_analysis_r1/arrays.npz).

Hard limits enforced here: cross-dimension absolute amplitudes are never compared (different
source normalization and geometric spreading); phase is reported only as a linear-fit formal
equivalent time shift with residual std, never as a transferable physical delay; no physical
thresholds, training labels or clean truths are produced. Running twice yields byte-identical
outputs (no randomness, no timestamps; inputs are identified by SHA-256 only).
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
SEL_MHZ = (20, 50, 80, 110, 140, 170)
TIME_WINDOW_S = 1200e-9
TAIL_NS = {'main_200': 200.0, 'robust_400': 400.0}
DIAG_WINDOW_NS = (440.0, 600.0)
PAIRS = ('B2D5CM', 'B2DANISO', 'DEP3D_5CM', 'DEP3D_ANISO')
CROSS_PAIRS = (('B2D5CM', 'DEP3D_5CM'), ('B2DANISO', 'DEP3D_ANISO'))
# Grid identity per design §3/§8; dt per design §7.1 (round_float convention, solver-verified).
CASES = {
    'B2D5CM': dict(dim='2D', shape=[1, 640, 1000], spacing=[.05, .05, .05],
                   rx_yz=[16.65, 45.], dt_s=117.9327e-12),
    'B2DANISO': dict(dim='2D', shape=[1, 640, 2000], spacing=[.05, .05, .025],
                     rx_yz=[16.65, 45.], dt_s=74.5872e-12),
    'DEP3D_5CM': dict(dim='3D', shape=[160, 200, 1000], spacing=[.05, .05, .05],
                      rx_yz=[5.65, 45.], dt_s=96.2917e-12),
    'DEP3D_ANISO': dict(dim='3D', shape=[160, 200, 2000], spacing=[.05, .05, .025],
                        rx_yz=[5.65, 45.], dt_s=68.0885e-12),
}
ZFINE2_DEFAULT = Path('artifacts/research_checks/2026-09-25_deep_zfine2_analysis_r1/arrays.npz')


def load_pair(runs_root, pair):
    """Load one BG/TGT pair, assert grid identity, receiver position, dtype, finiteness,
    identical BG/TGT source samples (design §9.3 pre-analysis checks)."""
    case = CASES[pair]
    paths = {k: Path(runs_root) / f'2026-09-26_{pair}_{k}' / f'{pair}_{k}.h5' for k in ('BG', 'TGT')}
    for v in paths.values():
        assert v.is_file(), v
    src, rx = {}, {}
    for k, v in paths.items():
        src[k] = load_source(v)
        rx[k] = load_receiver(v, receiver_path='name:measurement', component='Ex')
        with h5py.File(v, 'r') as h:
            assert np.array_equal(h.attrs['nx_ny_nz'], case['shape']), (pair, h.attrs['nx_ny_nz'])
            assert np.allclose(h.attrs['dx_dy_dz'], case['spacing'], rtol=0, atol=1e-12), pair
            assert np.allclose(h[rx[k].path].parent.attrs['Position'][1:], case['rx_yz'],
                               rtol=0, atol=1e-12), pair
            assert h[rx[k].path].dtype == np.float64, pair
        # dt sanity: solver-reported value (design §7.1), float32-rounded in storage.
        assert np.isclose(src[k].dt, rx[k].dt, rtol=0, atol=0), pair
        assert np.isclose(rx[k].dt, case['dt_s'], rtol=1e-4, atol=0), (pair, rx[k].dt)
        assert src[k].time_offset == src[k].dt / 2 and rx[k].time_offset == 0, pair
        assert np.all(np.isfinite(rx[k].samples)), pair
    assert np.array_equal(src['BG'].samples, src['TGT'].samples), f'{pair} source mismatch'
    dt = rx['BG'].dt
    n = min(len(rx['BG'].samples), int(np.floor(TIME_WINDOW_S / dt)) + 1)
    tapers = {tag: (round(ns * 1e-9 / dt) - .25) / n for tag, ns in TAIL_NS.items()}
    t = rx['BG'].times[:n]
    diff_t = rx['TGT'].samples[:n] - rx['BG'].samples[:n]
    spec = {}
    for tag in TAIL_NS:
        resp = {}
        for k in ('BG', 'TGT'):
            r = direct_frequency_response(src[k], replace(rx[k], samples=rx[k].samples[:n]), FREQ,
                                          tail_taper_fraction=tapers[tag])
            assert np.all(r.source_valid) and np.all(np.isfinite(r.response)), (pair, tag)
            resp[k] = r.response
        spec[tag] = resp['TGT'] - resp['BG']
    return dict(paths=paths, dim=case['dim'], shape=case['shape'], spacing=case['spacing'],
                dt=dt, n=n, t=t, tapers=tapers, diff_t=diff_t, spec=spec)


def selected_norm_db(mag):
    norm_db = 20 * np.log10(mag / mag.max())
    return {f'{s:g}': float(norm_db[int(np.argmin(abs(FREQ - s * 1e6)))]) for s in SEL_MHZ}


def spectrum_shape_corr(mag_a, mag_b):
    return float(np.corrcoef(mag_a / mag_a.max(), mag_b / mag_b.max())[0, 1])


def waveform_corr(t_coarse, x_coarse, t_fine, x_fine):
    """Max normalized cross-correlation after resampling the finer-dt trace onto the coarser
    time base, plus the optimal lag (dim23 convention)."""
    x_fine_on = np.interp(t_coarse, t_fine, x_fine)
    a = x_coarse / np.max(np.abs(x_coarse))
    b = x_fine_on / np.max(np.abs(x_fine_on))
    xcorr = np.correlate(a, b, mode='full')
    lag = int(np.argmax(xcorr) - (len(b) - 1))
    denom = float(np.sqrt(np.sum(a ** 2) * np.sum(b ** 2)))
    return dict(waveform_max_normalized_correlation=float(np.max(xcorr) / denom),
                waveform_optimal_lag_ns=float(lag * (t_coarse[1] - t_coarse[0]) * 1e9))


def phase_linear_fit(diff_a, diff_b):
    """Phase difference vs frequency, linear fit only. The equivalent time offset is a formal
    quantity; residual std quantifies non-linearity (dim23 measured 708 deg). Never a
    transferable physical delay."""
    dphi = np.unwrap(np.angle(diff_a)) - np.unwrap(np.angle(diff_b))
    slope, intercept = np.polyfit(FREQ, dphi, 1)
    return dict(slope_rad_per_Hz=float(slope),
                equivalent_time_offset_ns=float(-slope / (2 * np.pi) * 1e9),
                residual_std_deg=float(np.degrees(np.std(dphi - (slope * FREQ + intercept)))))


def far_field_transform(diff_3d, mag_2d):
    """Exploratory far-field 3D->2D factor D_2D_equiv(f) ∝ D_3D(f)/sqrt(i*2*pi*f)
    (Bleistein/Forbriger form); constants omitted, shape only, near-field caveats apply.
    Exploratory evidence only: never a conclusion basis, never a physical/training label."""
    g = 1 / np.sqrt(1j * 2 * np.pi * FREQ)
    transformed = np.abs(diff_3d * g)
    transformed /= transformed.max()
    return dict(form='D_2D_equiv(f) ∝ D_3D(f) / sqrt(i*2*pi*f); constants omitted, shape only',
                transformed_vs_2D_shape_correlation=spectrum_shape_corr(transformed, mag_2d),
                transformed_normalized_mag_dB_at_selected_MHz=selected_norm_db(transformed))


def diagnostic_window(t, env):
    lo, hi = DIAG_WINDOW_NS
    mask = (t >= lo * 1e-9) & (t <= hi * 1e-9)
    i = int(np.argmax(np.where(mask, env, -np.inf)))
    return dict(window_envelope_peak_time_ns=float(t[i] * 1e9),
                window_envelope_peak_to_global_peak_ratio=float(env[i] / env.max()),
                predicted_two_way_arrival_ns_note='air 100.1 + cover 80.0 + rock 340.2 ≈ 520 ns (design §9.3)')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--runs-root', type=Path, default=Path('artifacts/research_checks'))
    p.add_argument('--zfine2-arrays', type=Path, default=ZFINE2_DEFAULT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    assert a.zfine2_arrays.is_file(), a.zfine2_arrays

    data = {pair: load_pair(a.runs_root, pair) for pair in PAIRS}
    inputs = {str(v): hashlib.sha256(v.read_bytes()).hexdigest()
              for d in data.values() for v in d['paths'].values()}
    inputs[str(a.zfine2_arrays)] = hashlib.sha256(a.zfine2_arrays.read_bytes()).hexdigest()

    # Per-pair metrics: 200 ns main, 400 ns robustness, 440-600 ns diagnostic window.
    per_pair = {}
    for pair, d in data.items():
        mag = np.abs(d['spec']['main_200'])
        i = int(np.argmax(mag))
        env = np.abs(hilbert(d['diff_t']))
        robust_rel_l2 = float(np.linalg.norm(d['spec']['main_200'] - d['spec']['robust_400'])
                              / np.linalg.norm(d['spec']['main_200']))
        mag_rob = np.abs(d['spec']['robust_400'])
        per_pair[pair] = dict(
            dimension=d['dim'], grid=dict(nx_ny_nz=d['shape'], dx_dy_dz=d['spacing']),
            dt_s=float(d['dt']), transformed_samples=int(d['n']),
            tail_taper_fraction={tag: float(f) for tag, f in d['tapers'].items()},
            main_200ns=dict(
                time_diff_peak_V_m=float(np.max(np.abs(d['diff_t']))),
                envelope_peak_time_ns=float(d['t'][int(np.argmax(env))] * 1e9),
                spectral_peak_frequency_MHz=float(FREQ[i] / 1e6),
                spectral_peak_phase_deg=float(np.degrees(np.angle(d['spec']['main_200'][i]))),
                normalized_mag_dB_at_selected_MHz=selected_norm_db(mag)),
            robustness_400ns=dict(spectral_peak_frequency_MHz=float(FREQ[int(np.argmax(mag_rob))] / 1e6),
                                  relative_L2_vs_main=robust_rel_l2),
            diagnostic_window_440_600ns=diagnostic_window(d['t'], env),
            note='time_diff_peak_V_m carries the pair-internal source normalization; not comparable across dimensions')

    # 3D internal dz-refinement convergence (design §9.2): pure dz contrast, identical domains.
    d5, da = data['DEP3D_5CM'], data['DEP3D_ANISO']
    m5, ma = np.abs(d5['spec']['main_200']), np.abs(da['spec']['main_200'])
    e5, ea = np.abs(hilbert(d5['diff_t'])), np.abs(hilbert(da['diff_t']))
    convergence_3d = dict(
        comparison='DEP3D_ANISO (dz 2.5 cm) vs DEP3D_5CM (dz 5 cm), TGT-BG differences, 200 ns tail',
        envelope_peak_time_difference_ns=float(d['t'][int(np.argmax(ea))] * 1e9 - d['t'][int(np.argmax(e5))] * 1e9),
        spectral_peak_frequency_difference_MHz=float(FREQ[int(np.argmax(ma))] / 1e6 - FREQ[int(np.argmax(m5))] / 1e6),
        normalized_spectrum_shape_correlation=spectrum_shape_corr(ma, m5))
    convergence_3d.update(waveform_corr(d5['t'], d5['diff_t'], da['t'], da['diff_t']))
    convergence_3d['threshold_commitment'] = ('none: differences and their direction are reported; '
                                              'no pre-declared convergence criterion (design §9.2)')

    # Same-grid-spacing cross-dimension shape comparisons (design §9.3).
    cross = {}
    for p2d, p3d in CROSS_PAIRS:
        d2, d3 = data[p2d], data[p3d]
        m2, m3 = np.abs(d2['spec']['main_200']), np.abs(d3['spec']['main_200'])
        e2, e3 = np.abs(hilbert(d2['diff_t'])), np.abs(hilbert(d3['diff_t']))
        db2, db3 = selected_norm_db(m2), selected_norm_db(m3)
        row = dict(
            pair=f'{p2d}_vs_{p3d}',
            envelope_peak_time_difference_ns=float(d['t'][int(np.argmax(e2))] * 1e9 - d['t'][int(np.argmax(e3))] * 1e9),
            spectral_peak_frequency_difference_MHz=float(FREQ[int(np.argmax(m2))] / 1e6 - FREQ[int(np.argmax(m3))] / 1e6),
            normalized_spectrum_shape_correlation=spectrum_shape_corr(m2, m3),
            normalized_mag_trend_difference_dB={k: float(db2[k] - db3[k]) for k in db2})
        row.update(waveform_corr(d2['t'], d2['diff_t'], d3['t'], d3['diff_t']))
        row['phase_difference_linear_fit'] = phase_linear_fit(d2['spec']['main_200'], d3['spec']['main_200'])
        row['far_field_3d_to_2d_exploratory'] = far_field_transform(d3['spec']['main_200'], m2)
        cross[row['pair']] = row
    cross['confounds'] = (
        'even with aligned grid spacing, each cross-dimension pair mixes (i) the measured dimension '
        'effect (2D x-infinite target/line source/cylindrical spreading vs 3D 2 m finite target/point '
        'source/spherical spreading), (ii) lateral domain width (2D 32 m vs 3D 10 m) and (iii) target y '
        'position (2D y 14-18 vs 3D y 3-7; same relative position, different absolute distance to the '
        'lateral PML). These are not clean dimension factors; the two pair differences must not be '
        'subtracted from each other as if additive (design §9.3).')

    # Comparison to the DEP 2D ZFINE2 fine chain (shape + relative L2 on the shared 501-point grid;
    # grid spacing NOT aligned: dy 12.5 mm / dz 3.125 mm — reference chain, not a convergence pair).
    with np.load(a.zfine2_arrays) as z:
        assert np.array_equal(z['frequency_Hz'], FREQ), 'ZFINE2 frequency grid mismatch'
        zref = z['ZFINE2_200'].copy()
    vs_zfine2 = dict(
        reference='ZFINE2_200 (DEP 2D chain TGT-BG, 200 ns tail) from ' + str(a.zfine2_arrays),
        grid_alignment='none (dy 12.5 mm / dz 3.125 mm vs 5 cm-class grids); spectra directly comparable on the 501-point grid',
        pairs={})
    for pair in ('B2D5CM', 'B2DANISO'):
        diff = data[pair]['spec']['main_200']
        vs_zfine2['pairs'][pair] = dict(
            relative_L2_difference=float(np.linalg.norm(diff - zref) / np.linalg.norm(zref)),
            normalized_spectrum_shape_correlation=spectrum_shape_corr(np.abs(diff), np.abs(zref)))

    result = dict(
        batch='dep3d_gold_v1',
        frequency_Hz='linspace(20e6,170e6,501)',
        time_window_ns=1200, tail_windows_ns=dict(main=200, robustness=400),
        diagnostic_window_ns=list(DIAG_WINDOW_NS),
        per_pair=per_pair,
        convergence_3d_dz_refinement=convergence_3d,
        cross_dimension_same_spacing=cross,
        vs_dep2d_zfine2_chain=vs_zfine2,
        hard_limits=dict(
            cross_dimension_absolute_amplitude_compared=False,
            phase_not_transferable_across_dimensions=True,
            phase_reported_only_as='linear fit of the phase difference vs frequency; the equivalent time offset is a formal quantity, not a physical delay',
            physical_acceptance_threshold=None,
            training_labels_generated=False,
            clean_truth_generated=False,
            grid_convergence_certified=False,
            far_field_transform_is_conclusion_basis=False,
            conclusions_scope='limited to the dep3d_gold_v1 scenario family and the grids analysed; not extrapolated to field data',
            solver_invoked=False),
        inputs=inputs)
    (a.output / 'results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')

    np.savez_compressed(a.output / 'arrays.npz', frequency_Hz=FREQ,
                        **{f'diff_{pair}_{tag}': data[pair]['spec'][tag]
                           for pair in PAIRS for tag in TAIL_NS},
                        **{f'time_{pair}_s': data[pair]['t'] for pair in PAIRS},
                        **{f'diff_t_{pair}': data[pair]['diff_t'] for pair in PAIRS},
                        zfine2_200=zref)
    with (a.output / 'difference_spectra.csv').open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['frequency_Hz'] + [f'{pair}_200_{part}' for pair in PAIRS for part in ('re', 'im')]
                   + ['zfine2_200_re', 'zfine2_200_im'])
        w.writerows(zip(FREQ, *[c for pair in PAIRS for c in (data[pair]['spec']['main_200'].real,
                                                              data[pair]['spec']['main_200'].imag)],
                        zref.real, zref.imag))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    labels = {'B2D5CM': '2D TMx dy=dz=5cm', 'B2DANISO': '2D TMx dy=5cm dz=2.5cm',
              'DEP3D_5CM': '3D 5cm isotropic', 'DEP3D_ANISO': '3D dz=2.5cm'}
    for pair in PAIRS:
        ax[0].plot(FREQ / 1e6, np.abs(data[pair]['spec']['main_200']), label=labels[pair])
    ax[0].plot(FREQ / 1e6, np.abs(zref), 'k--', label='DEP 2D ZFINE2 chain (200 ns)')
    ax[0].set(xlabel='Frequency [MHz]', ylabel='|TGT-BG difference| [(V/m)/A]', title='Difference spectra, 200 ns tail')
    for pair in PAIRS:
        env = np.abs(hilbert(data[pair]['diff_t']))
        ax[1].plot(data[pair]['t'] * 1e9, env / env.max(), label=labels[pair])
    ax[1].axvspan(*DIAG_WINDOW_NS, color='0.85', zorder=0)
    ax[1].annotate('440-600 ns diagnostic window\n(predicted two-way ≈ 520 ns,\ntarget centre depth 20 m)',
                   xy=(520, .96), xycoords='data', fontsize=7, ha='center', va='top')
    ax[1].set(xlabel='Time [ns]', ylabel='Normalized envelope of TGT-BG', xlim=(0, 700),
              title='Time-domain differences (native time axes)')
    for axis in ax:
        axis.legend(fontsize=7)
        axis.grid(alpha=.2)
    fig.suptitle('dep3d_gold_v1: TGT-BG paired differences (antennas z=45 m, 15 m above surface z=30 m; '
                 'grid spacing aligned within each cross pair)')
    fig.savefig(a.output / 'comparison.png', dpi=160)
    plt.close(fig)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
