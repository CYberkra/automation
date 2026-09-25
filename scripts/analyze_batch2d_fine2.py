"""Cross-tier direction re-check for the batch2d_v1 fine2 (FINE2) review subset (CPU/NumPy only; no solver).

SFCW machinery and taper conventions are reused unchanged from scripts/analyze_batch2d.py: the
official gprMax.toolboxes.SFCW.processing load_source/load_receiver/direct_frequency_response on
501 frequency points ``linspace(20e6, 170e6, 501)``, a 1200 ns time window, a 200 ns tail taper as
the main result and a 400 ns tail taper as the robustness variant (same
``(round(tail/dt) - 0.25) / n`` taper-fraction convention). Nothing here launches, imports or writes
gprMax inputs and no solver is invoked; this is read-only post-processing of archived runs.

Case inventory, seeds, mother models and groups come from
configs/research/batch2d_v1_fine2/{cases,groups}.json per
docs/research/2026-09-26_batch_2d_spec_v1.md §1.6 (fine2 review segment), §2.1 (frequency table,
Ex component, 200/400 ns tail windows), §4.1-§4.2 (seed rule and same-mother-model grouping) and
§6-P1/P2/P3/P6/P7.

Three jobs:
(a) per-case load/recheck of the 8 fine2 runs, with the same grid-alignment, receiver-position,
    dtype, finite-value and supervision-excerpt checks as the BASE analysis. A missing or incomplete
    run is recorded with a reason and skipped; it is never fatal.
(b) P2 direction re-check: per family, the coarse BASE anchor paired difference
    (``B2D-C<h>m-D10m-W4m-T0.5m-E20-S0.02``, grid_tier BASE, taken from the archived BASE analysis)
    is compared with the physically identical fine2 grid variant. Reported per family only:
    per-frequency sign agreement of the complex difference, shape correlation of |difference|,
    envelope-peak arrival direction and amplitude-ratio direction. Families are never averaged
    together; the BASE and FINE2 tiers are never merged, never differenced against each other as if
    one were truth, and no direction-consistency threshold is introduced.
(c) reproducibility evidence: ``B2D-C3m-BG-F2`` vs the archived ``2026-09-25_DEP_BG_ZFINE2`` run
    (same grid tier, same receiver, input differing only in the ``#title`` line). Their Ex receiver
    waveforms are compared by max absolute difference, relative L2 and correlation, reported as
    measured.

Hard limits follow the spec and the two existing analyses: reference_state stays
``numerically_unresolved``; no physical or training label and no clean truth is produced;
``physical_acceptance_threshold`` stays None; absolute accuracy and grid convergence are not
claimed; conclusions stay within mechanism / operator-evaluation evidence for the four fine2
scenario families. results.json is deterministic (no randomness, no timestamps; inputs are
identified by SHA-256 only), so two runs on identical inputs give byte-identical JSON. The
auxiliary arrays.npz is a zip archive and therefore carries zip entry timestamps, so only
results.json is byte-stable.
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
TIER = 'FINE2'
COARSE_TIER = 'BASE'
FAMILY_ORDER = ('C1', 'C3', 'C5', 'C8')
ROLE_ORDER = ('bg', 'target', 'nc', 'off')
DATE_PREFIX = '2026-09-26_'
ANCHOR_SUFFIX = 'D10m-W4m-T0.5m-E20-S0.02'
FAMILY_HEIGHT = {'C1': '1', 'C3': '3', 'C5': '5', 'C8': '8'}
# Reproducibility pair of (c): the only fine2 BG whose archived same-grid twin exists.
REPRO_FINE2_RUN_ID = 'B2D-C3m-BG-F2'
REPRO_REF_STEM = 'DEP_BG_ZFINE2'
CASES_DEFAULT = Path('configs/research/batch2d_v1_fine2/cases.json')
GROUPS_DEFAULT = Path('configs/research/batch2d_v1_fine2/groups.json')
RUNS_ROOT_DEFAULT = Path('artifacts/simulations')
COARSE_RESULTS_DEFAULT = Path('artifacts/research_checks/2026-09-26_batch2d_analysis_r1/results.json')
COARSE_ARRAYS_DEFAULT = Path('artifacts/research_checks/2026-09-26_batch2d_analysis_r1/arrays.npz')
COARSE_CASES_DEFAULT = Path('configs/research/batch2d_v1/cases.json')
REPRO_REF_DEFAULT = Path('artifacts/simulations/2026-09-25_DEP_BG_ZFINE2')
OUTPUT_DEFAULT = Path('artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1')
NO_THRESHOLD = 'none declared: the numbers are reported as directions/agreement rates only; ' \
               'no pass threshold exists in this analysis or in the frozen contract (spec §6-P2)'


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_supervision(path):
    """Excerpt of the Windows supervisor record: reason, exit code, wall clock and committed Job memory."""
    rec = json.loads(path.read_text(encoding='utf-8'))
    out = dict(reason=rec['reason'], exit_code=rec['exit_code'], wall_s=float(rec['wall_s']),
               peak_job_commit_bytes=int(rec['peak_job_commit_bytes']),
               memory_semantics=rec.get('memory_semantics', 'undeclared'))
    if rec['reason'] != 'completed':
        out['note'] = 'supervision reason is not "completed"; archived as-is, results reported unchanged'
    return out


def case_sort_key(case):
    """Deterministic order: family, then role, then run_id. Unknown family/role is an error."""
    return (FAMILY_ORDER.index(case['family']), ROLE_ORDER.index(case['role']), case['run_id'])


def load_case(run_dir, case, grid):
    """Load one archived fine2 case: source, Ex receiver, grid identity checks, receiver position,
    dtype, finiteness, time base and the supervision excerpt. Returns (None, reason) when the run
    evidence is incomplete; nothing outside exceptional-code territory is fatal here."""
    h5_path = run_dir / (case['run_id'] + '.h5')
    sup_path = run_dir / 'supervision.json'
    if not h5_path.is_file():
        return None, 'h5 not found in %s' % run_dir.as_posix()
    if not sup_path.is_file():
        return None, 'supervision.json not found in %s' % run_dir.as_posix()
    src = load_source(h5_path)
    rx = load_receiver(h5_path, receiver_path='name:measurement',
                       component=case['geometry']['component'])
    with h5py.File(h5_path, 'r') as h:
        assert np.array_equal(h.attrs['nx_ny_nz'], grid['nx_ny_nz']), (case['run_id'], h.attrs['nx_ny_nz'])
        assert np.allclose(h.attrs['dx_dy_dz'], grid['dx_dy_dz'], rtol=0, atol=1e-12), case['run_id']
        assert np.allclose(h[rx.path].parent.attrs['Position'][1:], case['geometry']['rx_m'],
                           rtol=0, atol=1e-12), case['run_id']
        dtype = h[rx.path].dtype
        iterations = int(h.attrs['Iterations'])
    assert dtype == np.float64, case['run_id']
    assert np.isclose(src.dt, rx.dt, rtol=0, atol=0), case['run_id']
    assert np.isclose(rx.dt, case['dt_s'], rtol=1e-9, atol=0), (case['run_id'], rx.dt)
    assert src.time_offset == src.dt / 2 and rx.time_offset == 0, case['run_id']
    assert len(rx.samples) == iterations, (case['run_id'], len(rx.samples))
    assert np.all(np.isfinite(rx.samples)) and np.all(np.isfinite(src.samples)), case['run_id']
    assert iterations >= case['time_steps'], (case['run_id'], iterations)
    n = min(len(rx.samples), int(np.floor(TIME_WINDOW_S / rx.dt)) + 1)
    return dict(path=h5_path, supervision=read_supervision(sup_path), src=src, rx=rx, dt=rx.dt,
                n=n, t=rx.times[:n], iterations=iterations), None


def case_spectra(obj):
    """Official SFCW direct frequency response for both tail windows of one case."""
    spec, tapers = {}, {}
    for tag, ns in TAIL_NS.items():
        taper = (round(ns * 1e-9 / obj['dt']) - .25) / obj['n']
        r = direct_frequency_response(obj['src'], replace(obj['rx'], samples=obj['rx'].samples[:obj['n']]),
                                      FREQ, tail_taper_fraction=taper)
        assert bool(np.all(r.source_valid)) and bool(np.all(np.isfinite(r.response))), tag
        spec[tag] = r.response
        tapers[tag] = float(taper)
    return spec, tapers


def selected_norm_db(mag):
    norm_db = 20 * np.log10(mag / mag.max())
    return {f'{s:g}': float(norm_db[int(np.argmin(abs(FREQ - s * 1e6)))]) for s in SEL_MHZ}


def band_energy(diff, bg):
    """In-band energy ratio sum|DIFF|^2 / sum|BG|^2 on the 501 requested frequency points.
    Diagnostic only: it is neither detectability, nor SNR, nor a physical acceptance measure."""
    e_diff = float(np.sum(np.abs(diff) ** 2))
    e_bg = float(np.sum(np.abs(bg) ** 2))
    assert e_bg > 0.0, 'family BG in-band energy is zero'
    ratio = e_diff / e_bg
    return dict(in_band_energy_ratio_DIFF_over_BG=float(ratio),
                in_band_energy_ratio_DIFF_over_BG_dB=float(10.0 * np.log10(ratio)))


def difference_metrics(diff, bg_main, diff_t, t):
    """Metrics of one family-internal FINE2 TGT-BG paired difference (spec §6-P7)."""
    mag = np.abs(diff[MAIN])
    i = int(np.argmax(mag))
    env = np.abs(hilbert(diff_t))
    out = dict(
        grid_tier=TIER,
        time_diff_peak_V_m=float(np.max(np.abs(diff_t))),
        envelope_peak_time_ns=float(t[int(np.argmax(env))] * 1e9),
        spectral_peak_frequency_MHz=float(FREQ[i] / 1e6),
        spectral_peak_phase_deg=float(np.degrees(np.angle(diff[MAIN][i]))),
        robustness_400ns_spectral_peak_frequency_MHz=float(FREQ[int(np.argmax(np.abs(diff[ROBUST])))] / 1e6),
        robustness_400ns_relative_L2_vs_200ns=float(np.linalg.norm(diff[MAIN] - diff[ROBUST])
                                                   / np.linalg.norm(diff[MAIN])),
        normalized_mag_dB_at_selected_MHz=selected_norm_db(mag),
        note='time_diff_peak_V_m carries the pair-internal excitation scale; no absolute detectability claim')
    out.update(band_energy(diff[MAIN], bg_main))
    return out


def seed_of(cases_doc, case):
    """Seed rule, spec §4.1: sha256('<batch_id>|<group_id>|<case_id>|<variant_tag>')[:8] big-endian."""
    txt = '|'.join([cases_doc['batch_id'], case['group_id'], case['case_id'], case['variant_tag']])
    return int.from_bytes(hashlib.sha256(txt.encode('utf-8')).digest()[:8], 'big')


def load_coarse_anchor(fam, coarse_case, coarse_entry, arrays):
    """Coarse BASE anchor difference taken from the archived BASE analysis arrays. Returns
    (None, reason) when the anchor evidence is not present in that archive."""
    rid = coarse_case['run_id']
    need = [f'diff_{rid}_{MAIN}', f'diff_{rid}_{ROBUST}', f'diff_t_{rid}', f'time_{rid}_s']
    missing = [k for k in need if k not in arrays]
    if missing:
        return None, 'BASE analysis array(s) missing: %s' % ', '.join(missing)
    if coarse_entry is None:
        return None, 'anchor absent from BASE analysis results.json'
    diff = coarse_entry.get('difference', {})
    if not diff.get('computed'):
        return None, ('BASE anchor difference not computed in the archived analysis (%s)'
                      % diff.get('note', 'no note'))
    if coarse_entry.get('grid_tier') != COARSE_TIER:
        return None, 'BASE anchor carries grid_tier %r, expected %r' % (coarse_entry.get('grid_tier'),
                                                                        COARSE_TIER)
    bg_rid = diff.get('bg_run_id')
    if bg_rid != f'B2D-C{FAMILY_HEIGHT[fam]}m-BG':
        return None, 'BASE anchor is differenced against %r, not its own family BG (spec §6-P7)' % bg_rid
    t = np.asarray(arrays[f'time_{rid}_s'])
    diff_t = np.asarray(arrays[f'diff_t_{rid}'])
    env = np.abs(hilbert(diff_t))
    rec = float(t[int(np.argmax(env))] * 1e9)
    rep = diff.get('envelope_peak_time_ns')
    return dict(run_id=rid, case=coarse_case, grid_tier=COARSE_TIER, bg_run_id=bg_rid,
                spec={tag: np.asarray(arrays[f'diff_{rid}_{tag}']) for tag in TAIL_NS},
                diff_t=diff_t, t=t, reported_envelope_peak_time_ns=rep,
                recomputed_envelope_peak_time_ns=rec,
                envelope_peak_time_matches_reported=bool(rep is not None
                                                         and np.isclose(rec, float(rep), rtol=1e-12, atol=0))), None


def direction_label(value, eps=0.0):
    if value > eps:
        return 'FINE2 larger'
    if value < -eps:
        return 'FINE2 smaller'
    return 'equal at the reported precision'


def time_direction_label(delta_ns):
    if delta_ns > 0:
        return 'FINE2 envelope peak later'
    if delta_ns < 0:
        return 'FINE2 envelope peak earlier'
    return 'identical envelope peak index time'


def direction_row(fam, coarse, fine):
    """One family x tier-stratified direction re-check row (b). Nothing is averaged here and the two
    tiers are not differenced against each other as if one were the truth."""
    c_main, f_main = coarse['spec'][MAIN], fine['diff'][MAIN]
    c_rob, f_rob = coarse['spec'][ROBUST], fine['diff'][ROBUST]
    rates = {}
    for part in ('real', 'imag'):
        a = np.sign(getattr(c_main, part))
        b = np.sign(getattr(f_main, part))
        rates[f'per_frequency_sign_agreement_{part}'] = float(np.mean(a == b))
        rates[f'n_points_with_zero_{part}'] = int(np.sum((a == 0) | (b == 0)))
    c_mag, f_mag = np.abs(c_main), np.abs(f_main)
    c_env, f_env = np.abs(hilbert(coarse['diff_t'])), np.abs(hilbert(fine['diff_t']))
    c_peak_t = float(coarse['t'][int(np.argmax(c_env))] * 1e9)
    f_peak_t = float(fine['t'][int(np.argmax(f_env))] * 1e9)
    spec_ratio = float(np.linalg.norm(f_main) / np.linalg.norm(c_main))
    time_ratio = float(np.max(np.abs(fine['diff_t'])) / np.max(np.abs(coarse['diff_t'])))
    row = dict(
        family=fam, group_id=fine['group_id'], cover_thickness_m=float(fine['cover_thickness_m']),
        coarse_grid_tier=COARSE_TIER, coarse_anchor_run_id=coarse['run_id'],
        fine2_grid_tier=TIER, fine2_anchor_run_id=fine['run_id'],
        coarse_anchor_bg_run_id=coarse['bg_run_id'], fine2_anchor_bg_run_id=fine['bg_run_id'],
        computed=True,
        main_tail_ns=200, robustness_tail_ns=400,
        per_frequency_sign_agreement=rates,
        n_frequency_points=int(FREQ.size),
        zero_masking='none: every requested frequency point is counted, including zero parts',
        difference_spectrum_shape_correlation=dict(
            main_200ns=shape_corr(c_mag, f_mag), robustness_400ns=shape_corr(np.abs(c_rob), np.abs(f_rob))),
        envelope_peak_direction=dict(
            coarse_base_ns=c_peak_t, fine2_ns=f_peak_t,
            fine2_minus_base_ns=float(f_peak_t - c_peak_t),
            direction=time_direction_label(f_peak_t - c_peak_t),
            native_dt_s=dict(coarse_base=float(coarse['t'][1] - coarse['t'][0]),
                             fine2=float(fine['t'][1] - fine['t'][0])),
            note='each tier is read on its own native time axis; the two dt differ, so the comparison '
                 'is of the peak index time, not of a common grid'),
        amplitude_ratio_direction=dict(
            spectral_L2_ratio_FINE2_over_BASE=spec_ratio,
            spectral_L2_direction=direction_label(spec_ratio - 1.0, 0.0),
            time_peak_abs_ratio_FINE2_over_BASE=time_ratio,
            time_peak_abs_direction=direction_label(time_ratio - 1.0, 0.0),
            note='reported as a direction only: a cross-tier ratio mixes grid dispersion, source '
                 'discretisation and the taper sample count, so it is not an accuracy statement and '
                 'not a grid-error estimate'),
        spectral_peak_frequency_MHz=dict(coarse_base=float(FREQ[int(np.argmax(c_mag))] / 1e6),
                                         fine2=float(FREQ[int(np.argmax(f_mag))] / 1e6)),
        threshold=NO_THRESHOLD,
        interpretation='direction/agreement evidence only; not cross-tier averaging, not tier merging, '
                       'not a convergence or absolute-accuracy conclusion')
    return row


def shape_corr(mag_a, mag_b):
    """Correlation of the two |difference| spectra after each is normalised by its own maximum
    (same convention as scripts/analyze_dep3d_gold.py)."""
    return float(np.corrcoef(mag_a / mag_a.max(), mag_b / mag_b.max())[0, 1])


def load_repro_reference(ref_dir, grid):
    """The archived DEP_BG_ZFINE2 run: same FINE2 grid, same receiver, #title-only input difference."""
    h5_path = ref_dir / (REPRO_REF_STEM + '.h5')
    if not h5_path.is_file():
        return None, 'reference h5 not found in %s' % ref_dir.as_posix()
    src = load_source(h5_path)
    rx = load_receiver(h5_path, receiver_path='name:measurement', component='Ex')
    with h5py.File(h5_path, 'r') as h:
        assert np.array_equal(h.attrs['nx_ny_nz'], grid['nx_ny_nz']), h.attrs['nx_ny_nz']
        assert np.allclose(h.attrs['dx_dy_dz'], grid['dx_dy_dz'], rtol=0, atol=1e-12), REPRO_REF_STEM
        assert np.allclose(h[rx.path].parent.attrs['Position'][1:], [16.65, 45.0],
                           rtol=0, atol=1e-12), REPRO_REF_STEM
        assert h[rx.path].dtype == np.float64, REPRO_REF_STEM
        iterations = int(h.attrs['Iterations'])
    assert np.all(np.isfinite(rx.samples)) and np.all(np.isfinite(src.samples)), REPRO_REF_STEM
    out = dict(path=h5_path, sha256=sha256_file(h5_path), grid_tier=TIER,
               nx_ny_nz=list(grid['nx_ny_nz']), dx_dy_dz=list(grid['dx_dy_dz']),
               receiver_component='Ex', receiver_position_m=[16.65, 45.0], dtype='float64',
               dt_s=float(rx.dt), iterations=iterations, samples=int(len(rx.samples)),
               all_finite=True, src=src, rx=rx)
    sup = ref_dir / 'supervision.json'
    if sup.is_file():
        out['supervision'] = read_supervision(sup)
    else:
        out['supervision'] = None
        out['supervision_note'] = 'no supervision.json beside the reference h5'
    return out, None


def waveform_agreement(a, b):
    """Waveform-level agreement of two same-dt traces trimmed to their common length: max absolute
    difference, relative L2 (against b) and Pearson correlation."""
    m = min(len(a), len(b))
    x, y = np.asarray(a[:m]), np.asarray(b[:m])
    denom = float(np.linalg.norm(y))
    return dict(samples_compared=int(m),
                max_abs_difference=float(np.max(np.abs(x - y))),
                relative_L2_difference=float(np.linalg.norm(x - y) / denom) if denom > 0 else None,
                correlation_coefficient=float(np.corrcoef(x, y)[0, 1]),
                reference_of_relative_L2='the archived DEP_BG_ZFINE2 trace')


def input_diff_evidence(a_path, b_path):
    """Compare two gprMax input files. Only the #title line is allowed to differ."""
    la = a_path.read_text(encoding='utf-8').splitlines()
    lb = b_path.read_text(encoding='utf-8').splitlines()
    title_idx = [i for i, line in enumerate(la) if line.startswith('#title:')]
    differing = [i for i in range(min(len(la), len(lb))) if la[i] != lb[i] and i not in title_idx]
    differing += list(range(min(len(la), len(lb)), max(len(la), len(lb))))
    return dict(a=dict(path=a_path.as_posix(), sha256=sha256_file(a_path), lines=len(la)),
                b=dict(path=b_path.as_posix(), sha256=sha256_file(b_path), lines=len(lb)),
                title_line_index=title_idx,
                differing_line_indices_excluding_title=differing,
                identical_except_title=not differing and len(la) == len(lb))


def repro_block(ref, ref_reason, loaded_fine2, grid, contract_cases_dir):
    """(c) Reproducibility evidence for the DEP ZFINE2 pair; computed only when both sides exist."""
    out = dict(
        pair=dict(fine2_run_id=REPRO_FINE2_RUN_ID, reference_run_id=REPRO_REF_STEM,
                  reference_dir='artifacts/simulations/2026-09-25_DEP_BG_ZFINE2'),
        grid_tier=TIER,
        claim='the two runs differ only in the #title line (same grid, same receiver, same materials '
              'and same source); any non-zero difference is therefore run-to-run reproducibility '
              'evidence for this grid tier, reported as measured',
        note='not a convergence result, not an accuracy statement; reference_state stays '
             'numerically_unresolved (spec §6-P2/P6)')
    inp_ref = Path('artifacts/simulations/2026-09-25_DEP_BG_ZFINE2') / (REPRO_REF_STEM + '.in')
    inp_f2 = contract_cases_dir / (REPRO_FINE2_RUN_ID + '.in')
    if loaded_fine2 is not None:
        run_dir_in = loaded_fine2['path'].parent / (REPRO_FINE2_RUN_ID + '.in')
        if run_dir_in.is_file():
            inp_f2 = run_dir_in
    if inp_ref.is_file() and inp_f2.is_file():
        out['input_diff'] = input_diff_evidence(inp_f2, inp_ref)
    else:
        out['input_diff'] = dict(computed=False,
                                 note='one or both input files not found',
                                 looked_for=[inp_f2.as_posix(), inp_ref.as_posix()])
    if ref is None:
        out.update(computed=False, reference_status='unavailable: %s' % ref_reason,
                   grid=dict(nx_ny_nz=list(grid['nx_ny_nz']), dx_dy_dz=list(grid['dx_dy_dz'])))
    else:
        out['reference'] = dict(sha256=ref['sha256'], grid_tier=ref['grid_tier'],
                                nx_ny_nz=list(ref['nx_ny_nz']), dx_dy_dz=list(ref['dx_dy_dz']),
                                receiver_component=ref['receiver_component'],
                                receiver_position_m=list(ref['receiver_position_m']),
                                dtype=ref['dtype'], dt_s=ref['dt_s'], iterations=ref['iterations'],
                                samples=ref['samples'], all_finite=ref['all_finite'],
                                supervision=ref['supervision'],
                                path=ref['path'].as_posix())
        if ref.get('supervision_note'):
            out['reference']['supervision_note'] = ref['supervision_note']
    if loaded_fine2 is None:
        out['fine2_status'] = 'not_analyzed: fine2 run evidence missing for %s' % REPRO_FINE2_RUN_ID
        out['computed'] = False
        out['waveform_comparison'] = dict(computed=False, note='missing fine2 side; no partial metric invented')
        return out
    f2 = dict(path=loaded_fine2['path'].as_posix(), sha256=sha256_file(loaded_fine2['path']),
              grid_tier=TIER, dt_s=float(loaded_fine2['dt']), iterations=loaded_fine2['iterations'],
              samples=int(len(loaded_fine2['rx'].samples)), dtype='float64', all_finite=True)
    out['fine2'] = f2
    out['fine2_status'] = 'analyzed'
    if ref is None:
        out['waveform_comparison'] = dict(computed=False, note='missing reference side')
        return out
    assert np.isclose(loaded_fine2['dt'], ref['dt_s'], rtol=1e-9, atol=0), 'pair dt mismatch'
    out['computed'] = True
    out['waveform_comparison'] = dict(
        computed=True, component='Ex',
        receiver=waveform_agreement(loaded_fine2['rx'].samples, ref['rx'].samples),
        source=waveform_agreement(loaded_fine2['src'].samples, ref['src'].samples),
        dt_identical=bool(np.isclose(loaded_fine2['dt'], ref['dt_s'], rtol=0, atol=0)),
        sample_count=dict(fine2=int(len(loaded_fine2['rx'].samples)), reference=int(len(ref['rx'].samples))))
    return out


def draw_figure(path, rows, plot_data, families, cover_by_family, block_c):
    """Comparison figure on physical axes: |TGT-BG| against frequency [MHz] and the normalised
    envelope of the time-domain difference against time [ns], one row per family."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(families), 2, figsize=(11.0, 2.6 * len(families)), layout='constrained')
    ax = np.atleast_2d(ax)
    for i, fam in enumerate(families):
        row = rows.get(fam)
        data = plot_data.get(fam)
        cov = cover_by_family.get(fam)
        skip = row is None or not row.get('computed') or data is None
        for j in range(2):
            a = ax[i][j]
            a.set_title('%s (cover %.1f m)' % (fam, cov) if cov is not None else fam,
                        fontsize=9, loc='left')
            if skip:
                a.set_axis_off()
                a.text(0.5, 0.5, 'no FINE2 anchor evidence yet', ha='center', va='center',
                       fontsize=8, color='0.35')
                continue
            c = data['coarse']
            f = data['fine']
            if j == 0:
                a.plot(FREQ / 1e6, np.abs(c['spec'][MAIN]), '--', lw=1.0,
                       label='BASE, dy=dz=25 mm')
                a.plot(FREQ / 1e6, np.abs(f['diff'][MAIN]), '-', lw=1.0,
                       label='FINE2, dy 12.5 / dz 3.125 mm')
                a.set(xlabel='Frequency [MHz]', ylabel='|TGT-BG difference| [(V/m)/A]')
                txt = ('sign agree Re %.3f / Im %.3f\nshape corr %.3f (400 ns %.3f)'
                       % (row['per_frequency_sign_agreement']['per_frequency_sign_agreement_real'],
                          row['per_frequency_sign_agreement']['per_frequency_sign_agreement_imag'],
                          row['difference_spectrum_shape_correlation']['main_200ns'],
                          row['difference_spectrum_shape_correlation']['robustness_400ns']))
                a.text(0.02, 0.97, txt, transform=a.transAxes, fontsize=6.5, va='top', ha='left')
                a.legend(fontsize=6.5, loc='upper right')
            else:
                c_env = np.abs(hilbert(c['diff_t'])) / np.max(np.abs(hilbert(c['diff_t'])))
                f_env = np.abs(hilbert(f['diff_t'])) / np.max(np.abs(hilbert(f['diff_t'])))
                a.plot(c['t'] * 1e9, c_env, '--', lw=1.0, label='BASE envelope')
                a.plot(f['t'] * 1e9, f_env, '-', lw=1.0, label='FINE2 envelope')
                a.axvline(row['envelope_peak_direction']['coarse_base_ns'], color='C0',
                          ls=':', lw=0.8)
                a.axvline(row['envelope_peak_direction']['fine2_ns'], color='C1', ls=':', lw=0.8)
                a.set(xlabel='Time [ns]', ylabel='Normalised envelope of TGT-BG', xlim=(0, 700))
                txt = ('env peak BASE %.1f ns -> FINE2 %.1f ns\n%s\nspectral L2 ratio FINE2/BASE %.3f'
                       % (row['envelope_peak_direction']['coarse_base_ns'],
                          row['envelope_peak_direction']['fine2_ns'],
                          row['envelope_peak_direction']['direction'],
                          row['amplitude_ratio_direction']['spectral_L2_ratio_FINE2_over_BASE']))
                a.text(0.02, 0.97, txt, transform=a.transAxes, fontsize=6.5, va='top', ha='left')
                a.legend(fontsize=6.5, loc='center right')
            a.grid(alpha=.2)
    repro = ''
    if block_c.get('computed'):
        wc = block_c['waveform_comparison']['receiver']
        repro = ('\n(c) %s vs %s Ex: max|d| %.3e, rel L2 %.3e, corr %.6f'
                 % (REPRO_FINE2_RUN_ID, REPRO_REF_STEM, wc['max_abs_difference'],
                    wc['relative_L2_difference'], wc['correlation_coefficient']))
    elif block_c.get('fine2_status'):
        repro = '\n(c) pending: %s' % block_c['fine2_status']
    fig.suptitle('batch2d_v1 fine2 review: per-family BASE-anchor vs FINE2-variant TGT-BG differences\n'
                 '(200 ns tail; family rows only, no cross-family average, no cross-tier merging)'
                 + repro, fontsize=8.5)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=OUTPUT_DEFAULT)
    p.add_argument('--runs-root', type=Path, default=RUNS_ROOT_DEFAULT,
                   help='parent of the 2026-09-26_<run_id>/ directories '
                        '(artifacts/simulations before archiving, artifacts/research_checks after)')
    p.add_argument('--cases', type=Path, default=CASES_DEFAULT)
    p.add_argument('--groups', type=Path, default=GROUPS_DEFAULT)
    p.add_argument('--coarse-cases', type=Path, default=COARSE_CASES_DEFAULT)
    p.add_argument('--coarse-results', type=Path, default=COARSE_RESULTS_DEFAULT)
    p.add_argument('--coarse-arrays', type=Path, default=COARSE_ARRAYS_DEFAULT)
    p.add_argument('--repro-reference-dir', type=Path, default=REPRO_REF_DEFAULT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)

    cases_doc = json.loads(a.cases.read_text(encoding='utf-8'))
    groups_doc = json.loads(a.groups.read_text(encoding='utf-8'))
    assert cases_doc['batch_id'] == groups_doc['batch_id'] == 'batch2d_v1', 'wrong batch'
    assert cases_doc['grid_tier'] == groups_doc['grid_tier'] == TIER, 'this script only handles FINE2'
    assert cases_doc['segment'] == groups_doc['segment'] == 'fine2', 'wrong segment'
    assert len(cases_doc['cases']) == cases_doc['n_cases'], 'case count mismatch'
    # x-infinite TMx models carry one x cell; dx inherits dy (verified against every h5 in load_case).
    g = cases_doc['grid']
    grid = dict(nx_ny_nz=[1, g['ny'], g['nz']], dx_dy_dz=[g['dy_m'], g['dy_m'], g['dz_m']])
    # FINE2 runs sit on the same grid and time window as the archived DEP ZFINE2 reference chain.
    assert np.isclose(g['dt_s'], 1.0112647039820335e-11, rtol=0, atol=0), 'FINE2 dt changed'
    groups_by_id = {x['group_id']: x for x in groups_doc['groups']}

    cases = sorted(cases_doc['cases'], key=case_sort_key)
    for c in cases:
        assert c['grid_tier'] == TIER and c['segment'] == 'fine2', c['run_id']
        grp = groups_by_id[c['group_id']]
        assert grp['family'] == c['family'] and grp['cover_thickness_m'] == c['cover_thickness_m'], c['run_id']
        assert c['seed'] == seed_of(cases_doc, c), c['run_id']
        assert c['variant_tag'] == 'fine2', c['run_id']

    bg_by_family = {}
    for c in cases:
        if c['role'] == 'bg':
            assert c['family'] not in bg_by_family, c['family']
            bg_by_family[c['family']] = c['run_id']
    for c in cases:
        if c['role'] != 'bg':
            assert c['family'] in bg_by_family, f"family {c['family']} has no BG in the frozen contract"

    loaded, unavailable = {}, []
    for c in cases:
        obj, reason = load_case(a.runs_root / f"{DATE_PREFIX}{c['run_id']}", c, grid)
        if obj is None:
            run_dir = a.runs_root / f"{DATE_PREFIX}{c['run_id']}"
            unavailable.append(dict(run_id=c['run_id'], family=c['family'], role=c['role'],
                                    grid_tier=TIER, reason=reason,
                                    run_dir_present=run_dir.is_dir(),
                                    stdout_log_present=(run_dir / 'stdout.log').is_file()))
            continue
        loaded[c['run_id']] = obj
    for rid, obj in loaded.items():
        spec, tapers = case_spectra(obj)
        obj['spec'], obj['tapers'] = spec, tapers

    per_case, family_diff = {}, {}
    for c in cases:
        rid = c['run_id']
        entry = dict(run_id=rid, family=c['family'], role=c['role'], group_id=c['group_id'],
                     seed=c['seed'], variant_tag=c['variant_tag'], grid_tier=TIER,
                     cover_thickness_m=float(c['cover_thickness_m']),
                     mother_model_id=c['mother_model_id'],
                     scan_axis={k: v for k, v in c['scan_axis'].items() if v is not None},
                     exceptions=c['exceptions'])
        if rid not in loaded:
            entry['analysis_status'] = 'not_analyzed: %s' % [u['reason'] for u in unavailable
                                                             if u['run_id'] == rid][0]
            entry['difference'] = dict(computed=False)
            per_case[rid] = entry
            continue
        obj = loaded[rid]
        entry['analysis_status'] = 'analyzed'
        entry['h5'] = dict(path=obj['path'].as_posix(), sha256=sha256_file(obj['path']),
                           nx_ny_nz=list(grid['nx_ny_nz']), dx_dy_dz=list(grid['dx_dy_dz']),
                           dt_s=float(obj['dt']), iterations=obj['iterations'],
                           contract_time_steps=c['time_steps'],
                           iterations_minus_contract_time_steps=obj['iterations'] - int(c['time_steps']),
                           time_steps_note='the solver iteration count follows the ceil(window/dt)+1 '
                                          'convention seen on the BASE tier (20352) and on the archived '
                                          'DEP ZFINE2 reference (118665); the frozen contract field is a '
                                          'floor-based planning figure',
                           transformed_samples=int(obj['n']),
                           receiver_component=c['geometry']['component'],
                           receiver_position_m=list(c['geometry']['rx_m']), dtype='float64',
                           all_finite=True)
        entry['supervision'] = obj['supervision']
        entry['tail_taper_fraction'] = obj['tapers']
        bg_rid = bg_by_family[c['family']]
        if c['role'] == 'bg':
            entry['difference'] = dict(computed=False, reference_role='bg',
                                       note='family BG is the differencing reference itself (spec §6-P7)')
        elif bg_rid not in loaded:
            entry['difference'] = dict(computed=False, bg_run_id=bg_rid,
                                       note='family BG unavailable; no cross-family and no cross-tier '
                                            'substitution (spec §6-P7)')
        else:
            b = loaded[bg_rid]
            assert np.array_equal(obj['src'].samples, b['src'].samples), rid
            assert np.isclose(obj['dt'], b['dt'], rtol=0, atol=0) and obj['n'] == b['n'], rid
            dt_pair = {tag: obj['spec'][tag] - b['spec'][tag] for tag in TAIL_NS}
            dt_time = obj['rx'].samples[:obj['n']] - b['rx'].samples[:obj['n']]
            entry['difference'] = dict(computed=True, bg_run_id=bg_rid, grid_tier=TIER,
                                       source_samples_identical_to_family_bg=True)
            entry['difference'].update(difference_metrics(dt_pair, b['spec'][MAIN], dt_time, obj['t']))
            family_diff[c['family']] = dict(run_id=rid, group_id=c['group_id'], bg_run_id=bg_rid,
                                            case_id=c['case_id'], mother_model_id=c['mother_model_id'],
                                            cover_thickness_m=float(c['cover_thickness_m']),
                                            diff=dt_pair, diff_t=dt_time, t=obj['t'])
        per_case[rid] = entry

    # (b) P2 direction re-check, family rows only, against the archived BASE analysis.
    coarse_cases_doc = json.loads(a.coarse_cases.read_text(encoding='utf-8'))
    assert coarse_cases_doc['grid_tier'] == COARSE_TIER, 'coarse contract is not BASE'
    coarse_results = json.loads(a.coarse_results.read_text(encoding='utf-8'))
    assert coarse_results['batch'] == 'batch2d_v1' and coarse_results['grid_tier'] == COARSE_TIER, \
        'coarse results are not the batch2d_v1 BASE analysis'
    coarse_cases = coarse_results['cases']
    with np.load(a.coarse_arrays) as z:
        assert np.array_equal(z['frequency_Hz'], FREQ), 'coarse analysis frequency grid mismatch'
        arrays = {k: z[k] for k in z.files}
    p2, p2_private, unavailable_anchors = {}, {}, []
    for fam in FAMILY_ORDER:
        if fam not in family_diff:
            p2[fam] = dict(family=fam, grid_tier_pair=[COARSE_TIER, TIER], computed=False,
                           note='no fine2 anchor difference available for this family; nothing '
                                'substituted from another family or from the coarse tier alone')
            continue
        fine2 = family_diff[fam]
        mother = fine2['mother_model_id']
        coarse_case = next((x for x in coarse_cases_doc['cases'] if x['case_id'] == mother), None)
        if coarse_case is None:
            p2[fam] = dict(family=fam, grid_tier_pair=[COARSE_TIER, TIER], computed=False,
                           fine2_anchor_run_id=fine2['run_id'],
                           note='BASE contract has no case_id %r' % mother)
            continue
        assert coarse_case['run_id'] == f'B2D-C{FAMILY_HEIGHT[fam]}m-{ANCHOR_SUFFIX}', coarse_case['run_id']
        coarse_obj, reason = load_coarse_anchor(fam, coarse_case, coarse_cases.get(coarse_case['run_id']),
                                                arrays)
        if coarse_obj is None:
            unavailable_anchors.append(dict(family=fam, run_id=coarse_case['run_id'],
                                            grid_tier=COARSE_TIER, fine2_run_id=fine2['run_id'],
                                            reason=reason))
            p2[fam] = dict(family=fam, grid_tier_pair=[COARSE_TIER, TIER], computed=False,
                           coarse_anchor_run_id=coarse_case['run_id'], fine2_anchor_run_id=fine2['run_id'],
                           note='coarse anchor unavailable: %s' % reason)
            continue
        # Arrays stay out of results.json: kept separately and used only by the figure.
        p2_private[fam] = dict(coarse=coarse_obj, fine=fine2)
        p2[fam] = direction_row(fam, coarse_obj, fine2)
    p2['reporting'] = (
        'rows are family x grid-tier-stratified: no cross-family average, no cross-family '
        'differencing, and the BASE and FINE2 results are never merged, subtracted as error, or used '
        'to certify grid convergence (spec §6-P1/P2, §7.8); each tier keeps its own native dt, its own '
        'tail window samples and its own grid_tier label')
    p2['threshold'] = NO_THRESHOLD
    p2['reference_state'] = 'numerically_unresolved'

    # (c) Reproducibility evidence for the DEP ZFINE2 twin pair.
    ref, ref_reason = load_repro_reference(a.repro_reference_dir, grid)
    block_c = repro_block(ref, ref_reason, loaded.get(REPRO_FINE2_RUN_ID), grid,
                          Path('configs/research/batch2d_v1_fine2'))

    stratified = {}
    for fam in FAMILY_ORDER:
        fam_cases = [c for c in cases if c['family'] == fam]
        if not fam_cases:
            continue
        roles = {role: [c['run_id'] for c in fam_cases if c['role'] == role] for role in ROLE_ORDER}
        bg_rid = bg_by_family.get(fam)
        status = 'analyzed' if bg_rid in loaded else 'not_analyzed: family BG unavailable (spec §6-P7)'
        stratified[fam] = dict(
            grid_tier=TIER, cover_thickness_m=float(fam_cases[0]['cover_thickness_m']),
            group_id=fam_cases[0]['group_id'],
            mother_model_hash=groups_by_id[fam_cases[0]['group_id']]['mother_model_hash'],
            bg_run_id=bg_rid, differencing_family_status=status, roles=roles, n_cases=len(fam_cases),
            n_cases_analyzed=sum(1 for c in fam_cases if c['run_id'] in loaded),
            n_differences_computed=sum(1 for c in fam_cases if per_case[c['run_id']]['difference'].get('computed')),
            p2_direction_row_computed=bool(p2.get(fam, {}).get('computed')),
            reporting='family rows only; no cross-family average, no cross-family and no cross-tier '
                      'differencing (spec §6-P2/P7)')

    inputs = {loaded[c['run_id']]['path'].as_posix(): sha256_file(loaded[c['run_id']]['path'])
              for c in cases if c['run_id'] in loaded}
    inputs[a.cases.as_posix()] = sha256_file(a.cases)
    inputs[a.groups.as_posix()] = sha256_file(a.groups)
    inputs[a.coarse_cases.as_posix()] = sha256_file(a.coarse_cases)
    inputs[a.coarse_results.as_posix()] = sha256_file(a.coarse_results)
    inputs[a.coarse_arrays.as_posix()] = sha256_file(a.coarse_arrays)
    if ref is not None:
        inputs[ref['path'].as_posix()] = ref['sha256']
        ref_in = Path('artifacts/simulations/2026-09-25_DEP_BG_ZFINE2') / (REPRO_REF_STEM + '.in')
        if ref_in.is_file():
            inputs[ref_in.as_posix()] = sha256_file(ref_in)
    for extra in block_c.get('input_diff', {}).get('looked_for', []):
        if Path(extra).is_file():
            inputs[Path(extra).as_posix()] = sha256_file(Path(extra))

    result = dict(
        batch='batch2d_v1', segment='fine2', grid_tier=TIER,
        spec='docs/research/2026-09-26_batch_2d_spec_v1.md §1.6, §2.1, §4.1-§4.2, §6-P1/P2/P3/P6/P7',
        frequency_Hz='linspace(20e6,170e6,501)',
        time_window_ns=1200, tail_windows_ns=dict(main=200, robustness=400),
        n_cases_contract=int(cases_doc['n_cases']), n_cases_analyzed=len(loaded),
        pairing_policy=dict(
            rule='every non-BG fine2 case is differenced against the same-family same-tier BG only',
            cross_family_differencing=False, cross_tier_differencing=False,
            family_bg={fam: bg_by_family.get(fam) for fam in FAMILY_ORDER}),
        cases=per_case,
        stratified_by_family=stratified,
        p2_direction_consistency=p2,
        reproduction_dep_zfine2=block_c,
        unavailable_cases=unavailable,
        unavailable_coarse_anchors=unavailable_anchors,
        hard_limits=dict(
            grid_tier=TIER, grid_tiers_merged=False,
            cross_family_differencing=False, cross_tier_differencing=False,
            cross_family_average=False, family_bg_required=True,
            reference_state='numerically_unresolved',
            physical_acceptance_threshold=None, direction_consistency_threshold=None,
            physical_label_eligible=False, training_eligible=False,
            training_labels_generated=False, clean_truth_generated=False,
            grid_convergence_certified=False, absolute_accuracy_claimed=False,
            in_band_energy_ratio_is='diagnostic energy ratio on the 501-point analysis grid; '
                                    'not detectability, not SNR, not a physical threshold',
            amplitude_ratio_direction_is='cross-tier ratio reported as a direction indicator only; '
                                         'it mixes grid dispersion, source discretisation and taper '
                                         'sample count, so it is not an accuracy or grid-error measure',
            reproduction_pair_is='run-to-run reproducibility evidence on one grid tier, not a '
                                 'convergence result and not an accuracy statement',
            conclusions_scope='mechanism and operator-evaluation evidence for the four batch2d_v1 '
                              'fine2 scenario families only (spec §6-P2/P6); not field performance, '
                              'not 3D, not absolute accuracy, no maximum-depth claim',
            fdtd_solver_invoked=False),
        inputs=inputs)
    (a.output / 'results.json').write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + '\n', encoding='utf-8')

    arrays_out = dict(frequency_Hz=FREQ)
    for rid, obj in loaded.items():
        for tag in TAIL_NS:
            arrays_out[f'resp_{rid}_{tag}'] = obj['spec'][tag]
        arrays_out[f'time_{rid}_s'] = obj['t']
    for fam, d in family_diff.items():
        for tag in TAIL_NS:
            arrays_out[f'diff_{fam}_{tag}'] = d['diff'][tag]
        arrays_out[f'diff_t_{fam}'] = d['diff_t']
    np.savez_compressed(a.output / 'arrays.npz', **arrays_out)

    draw_figure(a.output / 'comparison.png', p2, p2_private, FAMILY_ORDER,
                {fam: float(stratified[fam]['cover_thickness_m']) for fam in stratified}, block_c)
    print(json.dumps(dict(batch='batch2d_v1', segment='fine2', grid_tier=TIER,
                          n_cases_analyzed=len(loaded),
                          n_family_differences_computed=len(family_diff),
                          n_unavailable=len(unavailable),
                          n_p2_rows_computed=sum(1 for fam in FAMILY_ORDER
                                                 if p2.get(fam, {}).get('computed')),
                          reproduction_computed=block_c.get('computed', False),
                          output=a.output.as_posix()), sort_keys=True))


if __name__ == '__main__':
    main()
