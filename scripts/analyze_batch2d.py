"""Stratified paired-difference analysis for the batch2d_v1 BASE 2D batch (CPU/NumPy only; no solver).

SFCW machinery reused from scripts/analyze_dep3d_gold.py: the official direct frequency response on
501 points ``linspace(20e6, 170e6, 501)``, a 1200 ns time window, 200 ns tail taper (main result) and
400 ns tail taper (robustness), with the same ``(round(tail/dt) - 0.25) / n`` taper-fraction
convention. Nothing here launches, imports or writes gprMax inputs; this is read-only post-processing.

Case inventory, roles, families, seeds and groups come from
configs/research/batch2d_v1/{cases,groups}.json per docs/research/2026-09-26_batch_2d_spec_v1.md
§1.6 (case roles), §4.1-§4.2 (seed rule and same-mother-model grouping), §6-P2 (grid-tier-stratified
reporting, no absolute-accuracy claim), §6-P6 (reference_state stays ``numerically_unresolved``),
§6-P7 (every target/NC/OFF case is differenced only against its own family BG) and §7 (prohibitions).

Reporting rules enforced here: results are stratified by family x role and never averaged across
families; each result carries ``grid_tier='BASE'``; negative controls stay in their own rows (NC with
the near-zero expectation, OFF with its recorded side-PML margin exception) and are never merged into
the target rows.

Hard limits follow the dep3d_gold analysis and the spec: reference_state stays
``numerically_unresolved``, physical_acceptance_threshold is None, no training labels and no clean
truths are produced, grid tiers are not merged, and conclusions stay within mechanism /
operator-evaluation evidence. Running twice on identical inputs yields byte-identical results.json
(no randomness, no timestamps; inputs are identified by SHA-256 only).
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
TAIL_NS = {'200ns': 200.0, '400ns': 400.0}
MAIN = '200ns'
ROBUST = '400ns'
FAMILY_ORDER = ('C1', 'C3', 'C5', 'C8')
ROLE_ORDER = ('bg', 'target', 'nc', 'off')
NC_EXPECTATION = '预期近零，偏离指示几何/平滑残差'
CASES_DEFAULT = Path('configs/research/batch2d_v1/cases.json')
GROUPS_DEFAULT = Path('configs/research/batch2d_v1/groups.json')
RUNS_ROOT_DEFAULT = Path('artifacts/simulations')


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
    """Load one archived case: source, Ex receiver, grid identity checks (spec §1.5 alignment), finite
    checks and the supervision excerpt. Returns (None, reason) when the run evidence is incomplete."""
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
    assert len(rx.samples) == iterations == case['time_steps'], (case['run_id'], len(rx.samples))
    assert np.all(np.isfinite(rx.samples)) and np.all(np.isfinite(src.samples)), case['run_id']
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
    """Metrics of one TGT-BG / NC-BG / OFF-BG paired difference (spec §6-P7)."""
    mag = np.abs(diff[MAIN])
    if float(mag.max()) == 0.0 and float(np.abs(diff[ROBUST]).max()) == 0.0 \
            and float(np.abs(diff_t).max()) == 0.0:
        return dict(
            grid_tier='BASE',
            time_diff_peak_V_m=0.0,
            envelope_peak_time_ns=None,
            spectral_peak_frequency_MHz=None,
            spectral_peak_phase_deg=None,
            robustness_400ns_spectral_peak_frequency_MHz=None,
            robustness_400ns_relative_L2_vs_200ns=None,
            normalized_mag_dB_at_selected_MHz=None,
            zero_difference=True,
            in_band_energy_ratio_DIFF_over_BG=0.0,
            in_band_energy_ratio_DIFF_over_BG_dB=None,
            note='paired difference is identically zero on both tail windows and in time (expected for '
                 'the null-contrast control); peak/phase/dB metrics not applicable, reported as null')
    i = int(np.argmax(mag))
    env = np.abs(hilbert(diff_t))
    out = dict(
        grid_tier='BASE',
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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--runs-root', type=Path, default=RUNS_ROOT_DEFAULT,
                   help='parent of the 2026-09-26_<run_id>/ directories '
                        '(artifacts/simulations before archiving, artifacts/research_checks after)')
    p.add_argument('--cases', type=Path, default=CASES_DEFAULT)
    p.add_argument('--groups', type=Path, default=GROUPS_DEFAULT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)

    cases_doc = json.loads(a.cases.read_text(encoding='utf-8'))
    groups_doc = json.loads(a.groups.read_text(encoding='utf-8'))
    assert cases_doc['batch_id'] == groups_doc['batch_id'] == 'batch2d_v1', 'wrong batch'
    assert cases_doc['grid_tier'] == groups_doc['grid_tier'] == 'BASE', 'this script only handles BASE'
    assert len(cases_doc['cases']) == cases_doc['n_cases'], 'case count mismatch'
    # x-infinite TMx models carry one x cell; dx inherits dy (verified against every h5 in load_case).
    grid = dict(nx_ny_nz=[1, cases_doc['grid']['ny'], cases_doc['grid']['nz']],
                dx_dy_dz=[cases_doc['grid']['dy_m'], cases_doc['grid']['dy_m'], cases_doc['grid']['dz_m']])
    groups_by_id = {g['group_id']: g for g in groups_doc['groups']}
    exceptions_by_run = {e['run_id']: e for e in cases_doc['exceptions_summary']}

    cases = sorted(cases_doc['cases'], key=case_sort_key)
    for c in cases:
        assert c['grid_tier'] == 'BASE', c['run_id']
        g = groups_by_id[c['group_id']]
        assert g['family'] == c['family'] and g['cover_thickness_m'] == c['cover_thickness_m'], c['run_id']
        # Seed rule, spec §4.1: sha256('<batch_id>|<group_id>|<case_id>|<variant_tag>')[:8] big-endian.
        seed_txt = '|'.join([cases_doc['batch_id'], c['group_id'], c['case_id'], c['variant_tag']])
        assert c['seed'] == int.from_bytes(hashlib.sha256(seed_txt.encode('utf-8')).digest()[:8], 'big'), c['run_id']

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
        obj, reason = load_case(a.runs_root / f"2026-09-26_{c['run_id']}", c, grid)
        if obj is None:
            unavailable.append(dict(run_id=c['run_id'], family=c['family'], role=c['role'],
                                    grid_tier='BASE', reason=reason))
            continue
        loaded[c['run_id']] = obj
    for rid, obj in loaded.items():
        spec, tapers = case_spectra(obj)
        obj['spec'], obj['tapers'] = spec, tapers

    per_case, diffs, diff_t, diff_order = {}, {}, {}, []
    for c in cases:
        rid = c['run_id']
        entry = dict(run_id=rid, family=c['family'], role=c['role'], group_id=c['group_id'],
                     seed=c['seed'], variant_tag=c['variant_tag'], grid_tier=c['grid_tier'],
                     cover_thickness_m=float(c['cover_thickness_m']),
                     scan_axis={k: v for k, v in c['scan_axis'].items()},
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
                           transformed_samples=int(obj['n']), receiver_component=c['geometry']['component'],
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
                                       note='family BG unavailable; no cross-family substitution (spec §6-P7)')
        else:
            b = loaded[bg_rid]
            assert np.array_equal(obj['src'].samples, b['src'].samples), rid
            assert np.isclose(obj['dt'], b['dt'], rtol=0, atol=0) and obj['n'] == b['n'], rid
            dt_pair = {tag: obj['spec'][tag] - b['spec'][tag] for tag in TAIL_NS}
            dt_time = obj['rx'].samples[:obj['n']] - b['rx'].samples[:obj['n']]
            entry['difference'] = dict(computed=True, bg_run_id=bg_rid, grid_tier='BASE',
                                       source_samples_identical_to_family_bg=True)
            entry['difference'].update(difference_metrics(dt_pair, b['spec'][MAIN], dt_time, obj['t']))
            diffs[rid], diff_t[rid], diff_order = dt_pair, dt_time, diff_order + [rid]
        per_case[rid] = entry

    # Negative controls stay single-listed: a paired difference can cancel a common residue that still
    # produces a false response in one control (spec §1.4, §7.9, background_mismatch_results.md:20-24).
    nc_rows, off_rows = [], []
    for c in cases:
        rid = c['run_id']
        row = dict(run_id=rid, family=c['family'], group_id=c['group_id'], grid_tier=c['grid_tier'])
        diff = per_case[rid]['difference']
        if diff.get('computed'):
            row.update(in_band_energy_ratio_DIFF_over_BG=diff['in_band_energy_ratio_DIFF_over_BG'],
                       in_band_energy_ratio_DIFF_over_BG_dB=diff['in_band_energy_ratio_DIFF_over_BG_dB'],
                       envelope_peak_time_ns=diff['envelope_peak_time_ns'],
                       spectral_peak_frequency_MHz=diff['spectral_peak_frequency_MHz'],
                       computed=True)
        else:
            row.update(computed=False, in_band_energy_ratio_DIFF_over_BG=None, note=diff.get('note'))
        if c['role'] == 'nc':
            row['expectation'] = NC_EXPECTATION
            nc_rows.append(row)
        elif c['role'] == 'off':
            exc = exceptions_by_run.get(rid)
            row['side_pml_margin_exception'] = dict(
                exceptions=exc['exceptions'] if exc else c['exceptions'],
                static_check={k: exc['static_check'][k] for k in ('side_pml_margin_m',
                                                                  'target_fully_in_sandstone')} if exc else None,
                freeze_disposition=exc['freeze_disposition'] if exc else None,
                note='OFF keeps a target contrast but sits 5 m off the Tx-Rx baseline; it is a negative '
                     'control for false positives, not evidence of detectability. Its lateral PML margin '
                     'is below the family 13 m reference, so it is single-listed (spec §1.4/§1.5).')
            off_rows.append(row)

    stratified = {}
    for fam in FAMILY_ORDER:
        fam_cases = [c for c in cases if c['family'] == fam]
        if not fam_cases:
            continue
        roles = {role: [c['run_id'] for c in fam_cases if c['role'] == role] for role in ROLE_ORDER}
        bg_rid = bg_by_family.get(fam)
        status = 'analyzed' if bg_rid in loaded else 'not_analyzed: family BG unavailable (spec §6-P7)'
        stratified[fam] = dict(
            grid_tier='BASE', cover_thickness_m=float(fam_cases[0]['cover_thickness_m']),
            group_id=fam_cases[0]['group_id'], mother_model_hash=groups_by_id[fam_cases[0]['group_id']]['mother_model_hash'],
            bg_run_id=bg_rid, differencing_family_status=status,
            roles=roles, n_cases=len(fam_cases),
            n_cases_analyzed=sum(1 for c in fam_cases if c['run_id'] in loaded),
            n_differences_computed=sum(1 for c in fam_cases if per_case[c['run_id']]['difference'].get('computed')),
            reporting='family rows only; no cross-family average, no cross-family differencing (spec §6-P2/P7)')

    inputs = {loaded[c['run_id']]['path'].as_posix(): sha256_file(loaded[c['run_id']]['path'])
              for c in cases if c['run_id'] in loaded}
    inputs[a.cases.as_posix()] = sha256_file(a.cases)
    inputs[a.groups.as_posix()] = sha256_file(a.groups)

    result = dict(
        batch='batch2d_v1', segment='base', grid_tier='BASE',
        spec='docs/research/2026-09-26_batch_2d_spec_v1.md §1.6, §4.1-§4.2, §6-P2/P6/P7',
        frequency_Hz='linspace(20e6,170e6,501)',
        time_window_ns=1200, tail_windows_ns=dict(main=200, robustness=400),
        n_cases_contract=int(cases_doc['n_cases']), n_cases_analyzed=len(loaded),
        pairing_policy=dict(
            rule='every non-BG case is differenced against the same-family same-grid BG only',
            cross_family_differencing=False,
            family_bg={fam: bg_by_family.get(fam) for fam in FAMILY_ORDER}),
        cases=per_case,
        stratified_by_family=stratified,
        negative_controls=dict(
            note='negative controls are single-listed and never merged into the target rows',
            nc=nc_rows, off=off_rows),
        unavailable_cases=unavailable,
        hard_limits=dict(
            grid_tier='BASE', grid_tiers_merged=False,
            cross_family_differencing=False, family_bg_required=True,
            negative_controls_merged_into_targets=False,
            reference_state='numerically_unresolved',
            physical_acceptance_threshold=None,
            physical_label_eligible=False, training_eligible=False,
            training_labels_generated=False, clean_truth_generated=False,
            grid_convergence_certified=False, absolute_accuracy_claimed=False,
            in_band_energy_ratio_is='diagnostic energy ratio on the 501-point analysis grid; '
                                    'not detectability, not SNR, not a physical threshold',
            conclusions_scope='mechanism and operator-evaluation evidence for the batch2d_v1 BASE '
                              'scenario families only (spec §6-P2/P6); not field performance, not 3D, '
                              'not absolute accuracy, no maximum-depth claim',
            fdtd_solver_invoked=False),
        inputs=inputs)
    (a.output / 'results.json').write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + '\n', encoding='utf-8')

    arrays = {'frequency_Hz': FREQ}
    for rid in diff_order:
        for tag in TAIL_NS:
            arrays[f'diff_{rid}_{tag}'] = diffs[rid][tag]
        arrays[f'diff_t_{rid}'] = diff_t[rid]
    for rid, obj in loaded.items():
        for tag in TAIL_NS:
            arrays[f'resp_{rid}_{tag}'] = obj['spec'][tag]
        arrays[f'time_{rid}_s'] = obj['t']
    np.savez_compressed(a.output / 'arrays.npz', **arrays)

    with (a.output / 'difference_spectra.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['frequency_Hz'] + [f'{rid}_{tag}_{part}'
                                       for rid in diff_order for tag in TAIL_NS for part in ('re', 'im')])
        w.writerows(zip(FREQ, *[c for rid in diff_order for tag in TAIL_NS
                                for c in (diffs[rid][tag].real, diffs[rid][tag].imag)]))
    print(json.dumps(dict(batch='batch2d_v1', grid_tier='BASE', n_cases_analyzed=len(loaded),
                          n_differences_computed=len(diff_order), n_unavailable=len(unavailable),
                          output=a.output.as_posix()), sort_keys=True))


if __name__ == '__main__':
    main()
