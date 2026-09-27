"""Paired-difference analysis for the a0_3d_v1 3D cross-check batch (CPU/NumPy only; no solver).

Adapted from scripts/analyze_dep3d_gold.py (unchanged SFCW machinery, tail-taper convention and
metric definitions). The comparison design follows
docs/research/2026-09-26_a0_3d_validation_contract_proposal.md §1-§3 and §8 (its SHA-256 is recorded
as a self-check but never asserted: the contract text may be revised, while the structural assertions
on the h5 files are always asserted).

Structure (proposal §2.1):
  1. one 3D-internal paired difference, TGT - BG (same domain, same grid, same source samples);
  2. cross-dimension **shape-only** comparisons of that 3D difference against the two already
     archived 2D A0 anchor differences, kept as two separate tier controls (BASE dy=dz=25 mm, FINE2
     dy 12.5 / dz 3.125 mm). The two tiers are never merged, never subtracted from each other, and
     the 3D and 2D arrays are never subtracted from one another - each side is differenced inside
     itself first and only the resulting shapes are compared (proposal §2.2 items 2-3);
  3. official gprMax SFCW post-processing (gprMax.toolboxes.SFCW.processing
     load_source/load_receiver/direct_frequency_response) on 501 frequency points 20-170 MHz, tail
     window 200 ns (main) and 400 ns (robustness), plus a 240-400 ns real-waveform diagnostic window
     built with the same +/-80 ns rule around the D10m self-computed two-way time 320.2 ns
     (proposal §1.3 difference 4);
  4. the seven comparison quantities of proposal §2.1: (1) normalised spectrum shape correlation
     (the main quantity), (2) envelope-peak arrival direction and its difference in ns, (3) spectral
     peak frequency difference, (4) max normalised cross-correlation of the time-domain difference
     plus the optimal lag, (5) per-frequency (20/50/80/110/140/170 MHz) normalised magnitude trend
     difference in dB, (6) linear fit of the phase difference vs frequency: a formal equivalent time
     shift plus the residual std (never a transferable physical delay), (7) an exploratory far-field
     3D->2D transform reported separately with far_field_transform_is_conclusion_basis False.

Hard limits are the six "explicitly not done" items of proposal §2.2 copied into the result file
verbatim as fields (no cross-dimension absolute amplitude comparison, no cross-grid subtraction, no
cross-tier merging/subtraction, no cross-family differencing, no convergence certification and no
pass/fail threshold, no labels of any kind), plus solver_invoked False and the statement that the
dep3d_gold_v1 D20m-class numbers are a mutually quotable magnitude reference, not an expected value
for this D10m anchor. Two runs on identical inputs give byte-identical results.json (no randomness,
no timestamps; inputs are identified by SHA-256 only).
"""
import argparse
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
TAIL_NS = {'200ns': 200.0, '400ns': 400.0}
MAIN = '200ns'
ROBUST = '400ns'
DIAG_WINDOW_NS = (240.0, 400.0)
RX_NAME = 'name:measurement'
COMPONENT = 'Ex'

BG_3D = 'B2D-C3m-BG-3D5CM'
TGT_3D = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM'
DIR_PREFIX_3D = '2026-09-27_'
# Structural expectations asserted before any analysis (mirrors
# configs/research/a0_3d_v1/static_check.json; iterations follow the ceil(window/dt)+1 convention).
EXPECTED_3D = dict(nx_ny_nz=[160, 200, 1000], dx_dy_dz=[0.05, 0.05, 0.05], iterations=12464,
                   dt_s=9.629166007732351e-11, rx_yz=[5.65, 45.0])
# 2D side: the same A0 mother model on two grid tiers; dt is never hardcoded for the 2D side, it is
# read from the h5 metadata (proposal §2.2 item 3: tiers stay stratified).
TIERS_2D = {
    'BASE': dict(date_prefix='2026-09-26_',
                 bg='B2D-C3m-BG', tgt='B2D-C3m-D10m-W4m-T0.5m-E20-S0.02',
                 nx_ny_nz=[1, 1280, 2000], dx_dy_dz=[0.025, 0.025, 0.025],
                 roots_key=('archive', 'runs')),
    'FINE2': dict(date_prefix='2026-09-26_',
                  bg='B2D-C3m-BG-F2', tgt='B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-F2',
                  nx_ny_nz=[1, 2560, 16000], dx_dy_dz=[0.0125, 0.0125, 0.003125],
                  roots_key=('runs', 'archive')),
}
SELFTEST_PSEUDO = {'BG': Path('artifacts/simulations/2026-09-26_DEP3D_5CM_BG/DEP3D_5CM_BG.h5'),
                   'TGT': Path('artifacts/simulations/2026-09-26_DEP3D_5CM_TGT/DEP3D_5CM_TGT.h5')}
CONTRACT_FILES = ('docs/research/2026-09-26_a0_3d_validation_contract_proposal.md',
                  'configs/research/a0_3d_v1/cases.json',
                  'configs/research/a0_3d_v1/groups.json',
                  'configs/research/a0_3d_v1/static_check.json')
DIAG_NOTE = ('reference only: the same +/-80 ns rule applied to the D10m self-computed two-way time '
             '320.2 ns (air 30 m/c + cover 6 m/(c/4) + rock 14 m/(c/3)); not a physical timing '
             'calibration (proposal §1.3 difference 4)')
CONFOUNDS = (
    'Cross-dimension comparisons here mix at least four factors that are NOT separated by this '
    'analysis (proposal §8.4 and §8.5): (i) the dimension effect under study - 2D x-infinite target / '
    'line source / cylindrical spreading vs 3D x-limited 2 m target / point source / spherical '
    'spreading; (ii) lateral domain width - 2D 32 m vs 3D 10 m; (iii) target y position - 2D y 14-18 m '
    'vs 3D y 3-7 m, the same relative position but a different absolute distance to the lateral PML; '
    '(iv) grid spacing, which is not aligned - 3D 5 cm isotropic vs 2D BASE dy=dz=25 mm and 2D FINE2 '
    'dy 12.5 / dz 3.125 mm. None of these is claimed to be removed. Consequently the BASE and FINE2 '
    'rows must not be subtracted from one another, and no row is a grid-error or dimension-error '
    'decomposition.')
REFERENCE_CONTEXT = dict(
    note='read-only magnitude cross-reference recorded for interpretation; it never enters any '
         'judgement, threshold or expectation in this analysis (proposal §2.3)',
    dep3d_d20m_class_measured=dict(
        normalized_spectrum_shape_correlation={'B2D5CM_vs_DEP3D_5CM': 0.935,
                                               'B2DANISO_vs_DEP3D_ANISO': 0.958},
        far_field_transformed_shape_correlation={'B2D5CM_vs_DEP3D_5CM': 0.967,
                                                 'B2DANISO_vs_DEP3D_ANISO': 0.979},
        envelope_peak_time_difference_ns={'B2D5CM_vs_DEP3D_5CM': -71.02,
                                          'B2DANISO_vs_DEP3D_ANISO': -53.72},
        spectral_peak_frequency_difference_MHz={'B2D5CM_vs_DEP3D_5CM': -5.70,
                                                'B2DANISO_vs_DEP3D_ANISO': -0.30},
        max_normalized_waveform_correlation={'B2D5CM_vs_DEP3D_5CM': 0.738,
                                             'B2DANISO_vs_DEP3D_ANISO': 0.732},
        max_selected_freq_trend_difference_dB=17.43,
        status='dep3d_gold_v1 measured values on the D20m class; mutually quotable as a magnitude '
               'reference only, NOT an expected value for this D10m anchor (proposal §2.3/§8)'),
    two_d_a0_anchor=dict(
        envelope_peak_time_ns={'BASE': 331.8626651806806, 'FINE2': 316.1921349940624},
        in_band_energy_ratio_DIFF_over_BG={'BASE': 2.8635252336237045e-07,
                                           'FINE2': 3.092143334046171e-07},
        spectral_peak_frequency_MHz={'BASE': 169.7, 'FINE2': 170.0},
        coarse_to_fine_shape_correlation=0.9994178181733713,
        coarse_to_fine_sign_agreement=0.9222,
        coarse_to_fine_envelope_peak_FINE2_minus_BASE_ns=-15.67053018661818,
        source='proposal §2.3 table (batch2d_v1 BASE and fine2 analyses)'))


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve_h5(run_id, roots, date_prefix):
    """Locate <root>/<date_prefix><run_id>/<run_id>.h5; fall back to a single h5 inside the run
    directory. Returns (path, how, tried)."""
    tried = []
    for root in roots:
        d = Path(root) / (date_prefix + run_id)
        p = d / (run_id + '.h5')
        tried.append(p.as_posix())
        if p.is_file():
            return p, 'exact run_id stem', tried
        if d.is_dir():
            cands = sorted(d.glob('*.h5'))
            if len(cands) == 1:
                return cands[0], 'single h5 found in the run directory', tried
            if len(cands) > 1:
                raise AssertionError('ambiguous h5 files in %s: %s'
                                     % (d.as_posix(), [c.as_posix() for c in cands]))
    raise AssertionError('h5 for run %s not found; tried: %s' % (run_id, tried))


def load_case(path, expected, label):
    """Load one h5: source, Ex receiver, structural assertions, dtype, finiteness, then the official
    SFCW direct frequency response for both tail windows."""
    src = load_source(path)
    rx = load_receiver(path, receiver_path=RX_NAME, component=COMPONENT)
    with h5py.File(path, 'r') as h:
        nx_ny_nz = [int(v) for v in np.asarray(h.attrs['nx_ny_nz']).ravel()]
        dx_dy_dz = [float(v) for v in np.asarray(h.attrs['dx_dy_dz']).ravel()]
        iterations = int(h.attrs['Iterations'])
        rx_pos = [float(v) for v in np.asarray(h[rx.path].parent.attrs['Position']).ravel()]
        dtype = h[rx.path].dtype
    assert dtype == np.float64, (label, dtype)
    assert np.isclose(src.dt, rx.dt, rtol=0, atol=0), label
    assert src.time_offset == src.dt / 2 and rx.time_offset == 0, label
    assert np.all(np.isfinite(rx.samples)) and np.all(np.isfinite(src.samples)), label
    assert len(rx.samples) == iterations, (label, len(rx.samples), iterations)
    dt = float(rx.dt)
    if expected:
        if 'nx_ny_nz' in expected:
            assert nx_ny_nz == list(expected['nx_ny_nz']), (label, nx_ny_nz)
        if 'dx_dy_dz' in expected:
            assert np.allclose(dx_dy_dz, expected['dx_dy_dz'], rtol=0, atol=1e-12), (label, dx_dy_dz)
        if 'iterations' in expected:
            assert iterations == int(expected['iterations']), (label, iterations)
        if 'dt_s' in expected:
            assert np.isclose(dt, expected['dt_s'], rtol=1e-6, atol=0), (label, dt)
        if 'rx_yz' in expected:
            assert np.allclose(rx_pos[1:], expected['rx_yz'], rtol=0, atol=1e-12), (label, rx_pos)
    n = min(len(rx.samples), int(np.floor(TIME_WINDOW_S / dt)) + 1)
    spec, tapers = {}, {}
    for tag, ns in TAIL_NS.items():
        taper = (round(ns * 1e-9 / dt) - 0.25) / n
        r = direct_frequency_response(src, replace(rx, samples=rx.samples[:n]), FREQ,
                                      tail_taper_fraction=taper)
        assert bool(np.all(r.source_valid)) and bool(np.all(np.isfinite(r.response))), (label, tag)
        spec[tag] = r.response
        tapers[tag] = float(taper)
    return dict(path=Path(path), label=label, run_id=Path(path).stem,
                grid=dict(nx_ny_nz=nx_ny_nz, dx_dy_dz=dx_dy_dz), iterations=iterations, dt=dt,
                rx_position_m=rx_pos, dtype='float64', all_finite=True, n=n, t=rx.times[:n],
                samples=rx.samples[:n], src_samples=src.samples, spec=spec, tapers=tapers)


def case_entry(obj):
    return dict(path=obj['path'].as_posix(), sha256=sha256_file(obj['path']), run_id=obj['run_id'],
                nx_ny_nz=list(obj['grid']['nx_ny_nz']), dx_dy_dz=list(obj['grid']['dx_dy_dz']),
                dt_s=float(obj['dt']), iterations=obj['iterations'],
                transformed_samples=int(obj['n']), receiver_component=COMPONENT,
                receiver_position_m=[float(v) for v in obj['rx_position_m']],
                dtype=obj['dtype'], all_finite=obj['all_finite'],
                tail_taper_fraction={k: float(v) for k, v in obj['tapers'].items()})


def paired_difference(tgt, bg, label):
    """Same-side TGT - BG difference (identical domain, grid, source samples and time base)."""
    assert np.array_equal(tgt['src_samples'], bg['src_samples']), label
    assert np.isclose(tgt['dt'], bg['dt'], rtol=0, atol=0) and tgt['n'] == bg['n'], label
    return dict(diff={tag: tgt['spec'][tag] - bg['spec'][tag] for tag in TAIL_NS},
                diff_t=tgt['samples'] - bg['samples'], t=tgt['t'], dt=tgt['dt'], n=tgt['n'],
                source_samples_identical_to_bg=True)


def selected_norm_db(mag):
    norm_db = 20 * np.log10(mag / mag.max())
    return {f'{s:g}': float(norm_db[int(np.argmin(abs(FREQ - s * 1e6)))]) for s in SEL_MHZ}


def spectrum_shape_corr(mag_a, mag_b):
    return float(np.corrcoef(mag_a / mag_a.max(), mag_b / mag_b.max())[0, 1])


def waveform_corr(t_coarse, x_coarse, t_fine, x_fine):
    """Max normalised cross-correlation after resampling the finer-dt trace onto the coarser time
    base, plus the optimal lag (dim23/dep3d convention). Positive lag: the coarser-dt side is
    delayed by that amount relative to the finer-dt side."""
    x_fine_on = np.interp(t_coarse, t_fine, x_fine)
    a = x_coarse / np.max(np.abs(x_coarse))
    b = x_fine_on / np.max(np.abs(x_fine_on))
    xcorr = np.correlate(a, b, mode='full')
    lag = int(np.argmax(xcorr) - (len(b) - 1))
    denom = float(np.sqrt(np.sum(a ** 2) * np.sum(b ** 2)))
    return dict(waveform_max_normalized_correlation=float(np.max(xcorr) / denom),
                waveform_optimal_lag_ns=float(lag * (t_coarse[1] - t_coarse[0]) * 1e9))


def waveform_corr_pair(a, b):
    """dt-ordered wrapper: the larger-dt side is the correlation base, the finer one is resampled."""
    coarse, fine = (a, b) if a['dt'] >= b['dt'] else (b, a)
    out = waveform_corr(coarse['t'], coarse['diff_t'], fine['t'], fine['diff_t'])
    out['base_side_larger_dt'] = coarse['label']
    out['resampled_side_smaller_dt'] = fine['label']
    out['native_dt_s'] = {coarse['label']: float(coarse['dt']), fine['label']: float(fine['dt'])}
    return out


def phase_linear_fit(diff_a, diff_b):
    """Phase difference vs frequency, linear fit only. The equivalent time offset is a formal
    quantity; the residual std quantifies non-linearity. Never a transferable physical delay."""
    dphi = np.unwrap(np.angle(diff_a)) - np.unwrap(np.angle(diff_b))
    slope, intercept = np.polyfit(FREQ, dphi, 1)
    return dict(difference_is='angle(A) - angle(B) with A the 2D side and B the 3D side',
                slope_rad_per_Hz=float(slope),
                equivalent_time_offset_ns=float(-slope / (2 * np.pi) * 1e9),
                residual_std_deg=float(np.degrees(np.std(dphi - (slope * FREQ + intercept)))),
                is_physical_delay=False)


def far_field_transform(diff_3d, mag_2d):
    """Exploratory far-field 3D->2D factor D_2D_equiv(f) proportional to D_3D(f)/sqrt(i*2*pi*f)
    (Bleistein/Forbriger form); constants omitted, shape only, near-field caveats apply.
    Exploratory evidence only: never a conclusion basis, never a physical or training label."""
    g = 1 / np.sqrt(1j * 2 * np.pi * FREQ)
    transformed = np.abs(diff_3d * g)
    transformed /= transformed.max()
    return dict(form='D_2D_equiv(f) proportional to D_3D(f) / sqrt(i*2*pi*f); constants omitted, '
                     'shape only',
                transformed_vs_2D_shape_correlation=spectrum_shape_corr(transformed, mag_2d),
                transformed_normalized_mag_dB_at_selected_MHz=selected_norm_db(transformed),
                is_conclusion_basis=False)


def diagnostic_window(t, env):
    lo, hi = DIAG_WINDOW_NS
    mask = (t >= lo * 1e-9) & (t <= hi * 1e-9)
    i = int(np.argmax(np.where(mask, env, -np.inf)))
    return dict(window_envelope_peak_time_ns=float(t[i] * 1e9),
                window_envelope_peak_to_global_peak_ratio=float(env[i] / env.max()),
                note=DIAG_NOTE)


def difference_block(d, bg_run_id, tgt_run_id, tier):
    """Per-side metrics of one internal TGT-BG paired difference (spectral + real-waveform)."""
    mag = np.abs(d['diff'][MAIN])
    i = int(np.argmax(mag))
    env = np.abs(hilbert(d['diff_t']))
    mag_rob = np.abs(d['diff'][ROBUST])
    return dict(
        bg_run_id=bg_run_id, tgt_run_id=tgt_run_id, grid_tier=tier,
        source_samples_identical_to_bg=d['source_samples_identical_to_bg'],
        dt_s=float(d['dt']), transformed_samples=int(d['n']),
        time_diff_peak_V_m=float(np.max(np.abs(d['diff_t']))),
        envelope_peak_time_ns=float(d['t'][int(np.argmax(env))] * 1e9),
        spectral_peak_frequency_MHz=float(FREQ[i] / 1e6),
        spectral_peak_phase_deg=float(np.degrees(np.angle(d['diff'][MAIN][i]))),
        normalized_mag_dB_at_selected_MHz=selected_norm_db(mag),
        robustness_400ns=dict(spectral_peak_frequency_MHz=float(FREQ[int(np.argmax(mag_rob))] / 1e6),
                              relative_L2_vs_main=float(np.linalg.norm(d['diff'][MAIN] - d['diff'][ROBUST])
                                                        / np.linalg.norm(d['diff'][MAIN]))),
        diagnostic_window_240_400ns=diagnostic_window(d['t'], env),
        note='time_diff_peak_V_m carries the pair-internal excitation scale; it is never compared '
             'across dimensions or across grid tiers (proposal §2.2 item 1)')


def time_direction(delta_ns):
    if delta_ns > 0:
        return '3D envelope peak later than 2D'
    if delta_ns < 0:
        return '3D envelope peak earlier than 2D'
    return 'identical envelope peak index time'


def cross_row(tier, two_d, three_d, cfg, three_d_grid):
    """One tier-stratified cross-dimension shape row: the seven quantities of proposal §2.1.
    Shape-only; nothing here is a pass/fail decision and the two tiers are never combined."""
    m2 = np.abs(two_d['diff'][MAIN])
    m3 = np.abs(three_d['diff'][MAIN])
    e2 = np.abs(hilbert(two_d['diff_t']))
    e3 = np.abs(hilbert(three_d['diff_t']))
    db2, db3 = selected_norm_db(m2), selected_norm_db(m3)
    t2 = float(two_d['t'][int(np.argmax(e2))] * 1e9)
    t3 = float(three_d['t'][int(np.argmax(e3))] * 1e9)
    row = dict(
        tier=tier, computed=True,
        two_d=dict(bg_run_id=two_d['bg_run_id'], tgt_run_id=two_d['tgt_run_id'],
                   nx_ny_nz=list(cfg['nx_ny_nz']), dx_dy_dz=list(cfg['dx_dy_dz']),
                   dt_s=float(two_d['dt']), native_time_base='own dt, own sample count'),
        three_d=dict(bg_run_id=three_d['bg_run_id'], tgt_run_id=three_d['tgt_run_id'],
                     nx_ny_nz=list(three_d_grid['nx_ny_nz']),
                     dx_dy_dz=list(three_d_grid['dx_dy_dz']),
                     dt_s=float(three_d['dt']), native_time_base='own dt, own sample count'),
        grids_aligned=False,
        grids_alignment_note='3D is 5 cm isotropic; this 2D tier has its own spacing; the two sides '
                             'are only ever compared as shapes after each side is differenced inside '
                             'itself (proposal §2.2 items 2-3)',
        main_tail_ns=200, robustness_tail_ns=400,
        metric_1_normalized_spectrum_shape_correlation=dict(
            main_200ns=spectrum_shape_corr(m2, m3),
            robustness_400ns=spectrum_shape_corr(np.abs(two_d['diff'][ROBUST]),
                                                 np.abs(three_d['diff'][ROBUST])),
            is_main_quantity=True),
        metric_2_envelope_peak_arrival_direction=dict(
            two_d_ns=t2, three_d_ns=t3, three_d_minus_two_d_ns=float(t3 - t2),
            direction=time_direction(t3 - t2),
            note='each side is read on its own native time axis; the dt differ, so this compares the '
                 'peak index time. Reported as a direction and a difference only, never as an '
                 'offset error (proposal §2.1 item 2)'),
        metric_3_spectral_peak_frequency_difference_MHz=dict(
            two_d=float(FREQ[int(np.argmax(m2))] / 1e6), three_d=float(FREQ[int(np.argmax(m3))] / 1e6),
            three_d_minus_two_d=float(FREQ[int(np.argmax(m3))] / 1e6 - FREQ[int(np.argmax(m2))] / 1e6),
            threshold='none: reported as a number, no pass reading (proposal §2.1 item 3)'),
        metric_4_time_difference_waveform_correlation=waveform_corr_pair(two_d, three_d),
        metric_5_selected_freq_normalized_mag_trend_difference_dB=dict(
            two_d_minus_three_d_dB={k: float(db2[k] - db3[k]) for k in db2},
            selected_MHz=list(SEL_MHZ),
            note='each side normalised by its own maximum first; absolute amplitudes are never '
                 'compared (proposal §2.2 item 1)'),
        metric_6_phase_difference_linear_fit=phase_linear_fit(two_d['diff'][MAIN], three_d['diff'][MAIN]),
        metric_7_far_field_3d_to_2d_exploratory=far_field_transform(three_d['diff'][MAIN], m2),
        confounds=CONFOUNDS,
        thresholds=dict(direction_consistency_threshold=None, physical_acceptance_threshold=None,
                        shape_correlation_threshold=None,
                        note='no threshold of any kind is declared (proposal §2.2 item 5)'),
        interpretation='direction and magnitude evidence on paired-difference shapes; not a '
                       'convergence result, not an absolute-accuracy result, not a physical decision')
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=None,
                   help='default: artifacts/research_checks/2026-09-27_a0_3d_analysis_r1 '
                        '(artifacts/local_checks/a0_3d_analysis_selftest with --selftest)')
    p.add_argument('--runs-root', type=Path, default=Path('artifacts/simulations'))
    p.add_argument('--archive-root', type=Path, default=Path('artifacts/research_checks'))
    p.add_argument('--selftest', action='store_true',
                   help='run the whole pipeline on the archived dep3d 5 cm pair as a pseudo input; '
                        'structural expectations are read from the h5 instead of being asserted')
    p.add_argument('--allow-missing-2d', action='store_true',
                   help='record a missing 2D tier as unavailable instead of aborting')
    a = p.parse_args()
    output = a.output or (Path('artifacts/local_checks/a0_3d_analysis_selftest') if a.selftest
                          else Path('artifacts/research_checks/2026-09-27_a0_3d_analysis_r1'))
    output.mkdir(parents=True, exist_ok=a.selftest)

    cases, provenance, unavailable_2d = {}, {}, []

    # --- 3D side (mandatory): one BG/TGT pair, differenced inside the 3D domain.
    if a.selftest:
        for k, v in SELFTEST_PSEUDO.items():
            assert v.is_file(), 'selftest pseudo input missing: %s' % v.as_posix()
        three_d = {k: load_case(v, None, '3D_%s' % k) for k, v in SELFTEST_PSEUDO.items()}
        provenance['3D'] = dict(selftest_pseudo_input=True,
                                BG=SELFTEST_PSEUDO['BG'].as_posix(),
                                TGT=SELFTEST_PSEUDO['TGT'].as_posix(),
                                note='dep3d_gold_v1 5 cm pair stands in for the a0_3d_v1 pair; it is '
                                     'the D20m class, so every number here exercises the pipeline '
                                     'only and is NOT an A0 result',
                                structural_expectations='read from the h5, not asserted')
        bg_id, tgt_id = SELFTEST_PSEUDO['BG'].stem, SELFTEST_PSEUDO['TGT'].stem
    else:
        three_d, bg_id, tgt_id = {}, None, None
        for k, rid in (('BG', BG_3D), ('TGT', TGT_3D)):
            path, how, tried = resolve_h5(rid, (a.runs_root, a.archive_root), DIR_PREFIX_3D)
            three_d[k] = load_case(path, EXPECTED_3D, '3D_%s' % k)
            provenance['3D_%s' % k] = dict(run_id=rid, path=path.as_posix(), resolution=how,
                                           tried=tried)
            if k == 'BG':
                bg_id = rid
            else:
                tgt_id = rid
        provenance['3D_structural_expectations'] = dict(
            source='configs/research/a0_3d_v1/static_check.json (mirrored in this script)',
            asserted=dict(nx_ny_nz=EXPECTED_3D['nx_ny_nz'], dx_dy_dz=EXPECTED_3D['dx_dy_dz'],
                          iterations=EXPECTED_3D['iterations'], dt_s=EXPECTED_3D['dt_s'],
                          rx_yz=EXPECTED_3D['rx_yz'], dtype='float64', all_finite=True,
                          bg_tgt_source_samples_identical=True))
    for k, v in three_d.items():
        cases['3D_%s' % k] = case_entry(v)
        cases['3D_%s' % k]['side'] = '3D'
        cases['3D_%s' % k]['grid_tier'] = '3D5CM' if not a.selftest else 'selftest-pseudo-3D5CM'
    d3 = paired_difference(three_d['TGT'], three_d['BG'], '3D')
    d3['label'] = '3D'
    d3['bg_run_id'], d3['tgt_run_id'] = bg_id, tgt_id

    # --- 2D side: two separate tiers, each differenced inside itself.
    two_d, two_d_cases = {}, {}
    roots_by_key = {'runs': a.runs_root, 'archive': a.archive_root}
    for tier, cfg in TIERS_2D.items():
        expected = dict(nx_ny_nz=cfg['nx_ny_nz'], dx_dy_dz=cfg['dx_dy_dz'])
        tier_roots = tuple(roots_by_key[k] for k in cfg['roots_key'])
        loaded, why = {}, None
        for role, rid in (('BG', cfg['bg']), ('TGT', cfg['tgt'])):
            try:
                path, how, tried = resolve_h5(rid, tier_roots, cfg['date_prefix'])
            except AssertionError as exc:
                why = str(exc)
                break
            loaded[role] = load_case(path, expected, '%s_%s' % (tier, role))
            provenance['%s_%s' % (tier, role)] = dict(run_id=rid, path=path.as_posix(),
                                                      resolution=how, tried=tried)
            two_d_cases['%s_%s' % (tier, role)] = case_entry(loaded[role])
            two_d_cases['%s_%s' % (tier, role)]['side'] = '2D'
            two_d_cases['%s_%s' % (tier, role)]['grid_tier'] = tier
        if why is not None:
            unavailable_2d.append(dict(tier=tier, bg_run_id=cfg['bg'], tgt_run_id=cfg['tgt'],
                                       reason=why))
            continue
        d = paired_difference(loaded['TGT'], loaded['BG'], tier)
        d['label'] = tier
        d['bg_run_id'], d['tgt_run_id'] = cfg['bg'], cfg['tgt']
        two_d[tier] = d
    # A missing 2D tier aborts the analysis by default (the cross-dimension rows are the point of
    # this batch); the selftest is a pipeline smoke test, so it only records the gap.
    if unavailable_2d and not (a.allow_missing_2d or a.selftest):
        raise AssertionError('2D tier input unavailable (use --allow-missing-2d to record it as '
                             'unavailable instead): %s' % unavailable_2d)
    cases.update(two_d_cases)
    for u in unavailable_2d:
        print('WARNING: 2D tier %s unavailable, its cross-dimension row is recorded as '
              'computed=false: %s' % (u['tier'], u['reason']))

    # --- results
    diff_3d = difference_block(d3, bg_id, tgt_id, '3D5CM' if not a.selftest else 'selftest-pseudo-3D5CM')
    cross = {tier: cross_row(tier, two_d[tier], d3, cfg, three_d['BG']['grid'])
             for tier, cfg in TIERS_2D.items() if tier in two_d}
    for tier in TIERS_2D:
        if tier not in cross:
            cross[tier] = dict(tier=tier, computed=False,
                               note='2D tier input unavailable; no other tier, family or dimension '
                                    'was substituted (proposal §2.2 items 2-4)')
    # Cross-check against the read-only reference magnitudes (reported, never asserted).
    ref_check = {}
    for tier, d in two_d.items():
        env = np.abs(hilbert(d['diff_t']))
        recomputed = float(d['t'][int(np.argmax(env))] * 1e9)
        reported = REFERENCE_CONTEXT['two_d_a0_anchor']['envelope_peak_time_ns'].get(tier)
        ref_check[tier] = dict(recomputed_envelope_peak_time_ns=recomputed,
                               reference_context_envelope_peak_time_ns=reported,
                               recomputed_minus_reference_ns=None if reported is None
                               else float(recomputed - reported),
                               status='reported only; not asserted and not a pass criterion')
    cross['reporting'] = (
        'each tier is a separate control against the same 3D difference: no cross-grid subtraction, '
        'no cross-tier merging, no cross-tier subtraction, no cross-family differencing '
        '(proposal §2.2 items 2-4). Both rows may be quoted side by side, never combined.')
    cross['reference_context_check'] = ref_check

    inputs = {}
    for obj in list(three_d.values()):
        inputs[obj['path'].as_posix()] = sha256_file(obj['path'])
    for k in two_d_cases:
        inputs[two_d_cases[k]['path']] = two_d_cases[k]['sha256']
    for f in CONTRACT_FILES:
        fp = Path(f)
        if fp.is_file():
            inputs[fp.as_posix()] = sha256_file(fp)

    result = dict(
        batch='a0_3d_v1',
        analysis='3D cross-check of the C3-family A0 anchor (D10m): shape-only comparison of the 3D '
                 'TGT-BG paired difference against the two archived 2D A0 anchor differences',
        selftest=bool(a.selftest),
        convention_source=dict(
            proposal='docs/research/2026-09-26_a0_3d_validation_contract_proposal.md §1-§3, §8',
            proposal_sha256=inputs.get(CONTRACT_FILES[0]),
            proposal_sha256_is_asserted=False,
            proposal_sha256_note='recorded as a self-check only; the contract text may be revised, '
                                 'the structural assertions on the h5 are the enforced part',
            metric_definitions='proposal §2.1 (seven quantities); implementation reused from '
                               'scripts/analyze_dep3d_gold.py',
            hard_limits='proposal §2.2 (six "explicitly not done" items) copied into hard_limits',
            solver_invoked=False),
        frequency_Hz='linspace(20e6,170e6,501)',
        time_window_ns=1200, tail_windows_ns=dict(main=200, robustness=400),
        diagnostic_window_ns=list(DIAG_WINDOW_NS),
        receiver=dict(name='measurement', component=COMPONENT),
        grid_tiers=dict(three_d='3D5CM', two_d=list(TIERS_2D), tiers_merged=False),
        cases=cases,
        provenance=provenance,
        difference_3d_internal=diff_3d,
        cross_dimension_shape_comparison=cross,
        unavailable_2d_tiers=unavailable_2d,
        reference_context=REFERENCE_CONTEXT,
        hard_limits=dict(
            cross_dimension_absolute_amplitude_compared=False,
            grid_convergence_certified=False,
            direction_consistency_threshold=None,
            physical_acceptance_threshold=None,
            training_labels_generated=False,
            clean_truth_generated=False,
            physical_label_eligible=False,
            training_eligible=False,
            cross_grid_subtraction_performed=False,
            cross_tier_merged=False,
            cross_tier_subtraction_performed=False,
            cross_family_differencing=False,
            phase_not_transferable_across_dimensions=True,
            phase_reported_only_as='linear fit of the phase difference vs frequency; the equivalent '
                                   'time offset is a formal quantity, not a physical delay',
            far_field_transform_is_conclusion_basis=False,
            dep3d_d20m_reference_is='mutually quotable magnitude reference measured on the D20m '
                                    'class; not an expected value and not a threshold for this D10m '
                                    'anchor (proposal §2.3/§8)',
            reference_state='numerically_unresolved',
            solver_invoked=False,
            fdtd_solver_invoked=False,
            conclusions_scope='limited to the a0_3d_v1 scenario family (C3 A0 anchor, 3D 5 cm '
                              'isotropic and the two 2D grid tiers analysed); not field performance, '
                              'not a maximum-depth or detectability conclusion'),
        inputs=inputs)
    (output / 'results.json').write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + '\n',
        encoding='utf-8')

    arrays = dict(frequency_Hz=FREQ, time_3d_s=d3['t'], diff_t_3d=d3['diff_t'],
                  diff_3d_200ns=d3['diff'][MAIN], diff_3d_400ns=d3['diff'][ROBUST])
    for tier, d in two_d.items():
        arrays[f'time_{tier}_s'] = d['t']
        arrays[f'diff_t_{tier}'] = d['diff_t']
        arrays[f'diff_{tier}_200ns'] = d['diff'][MAIN]
        arrays[f'diff_{tier}_400ns'] = d['diff'][ROBUST]
    np.savez_compressed(output / 'arrays.npz', **arrays)

    print(json.dumps(dict(batch='a0_3d_v1', selftest=bool(a.selftest),
                          three_d_pair=[bg_id, tgt_id],
                          two_d_tiers_computed=sorted(two_d),
                          n_cases_loaded=len(cases),
                          unavailable_2d_tiers=[u['tier'] for u in unavailable_2d],
                          results_json_sha256=sha256_file(output / 'results.json'),
                          output=output.as_posix()), sort_keys=True))


if __name__ == '__main__':
    main()
