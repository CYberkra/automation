"""Paired-difference analysis for the batch2d_slope_t2 batch (CPU/NumPy only; no solver).

SFCW machinery identical to scripts/analyze_batch2d.py (official direct frequency response on
501 points linspace(20e6, 170e6, 501), 1200 ns window, 200 ns main / 400 ns robustness tail taper
with the same taper-fraction convention). Read-only post-processing over archived runs.

Analysis plan (first slope batch: C3 family x tier T2 x {noTZ, TZ} x {BG, NC, TGT}):
  1. Within-variant paired differences (spec 2026-09-26_batch_2d_spec_v1.md section 6-P7):
     NC-BG and TGT-BG for S2 and for S2TZ separately. Each slope variant is its own
     mother-model family with its own BG; no differencing against the flat batch2d_v1 BG.
  2. Cross-variant background difference BG_S2TZ - BG_S2: isolates the transition-zone effect
     on the interface response (same mother except the 1 m tzone strips).
  3. Cross-geometry diagnostic vs the archived flat anchor (batch2d_v1, runs of 2026-09-26):
     S2 BG vs flat C3 BG spectra and interface-echo envelope-peak shift. Explicitly NOT a paired
     difference: the T2 tier shifts cover center 3.0 -> 4.2 m (naming caveat in cases.json), so
     timing/energy differences confound slope angle with cover thickness. Reported descriptively.

Hard limits identical to analyze_batch2d.py: reference_state stays numerically_unresolved,
physical_acceptance_threshold is None, no training labels, no clean truths, no cross-tier merging,
negative controls single-listed. Running twice on identical inputs yields byte-identical
results.json (no randomness, no timestamps; inputs identified by SHA-256 only).
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
BATCH = 'batch2d_slope_t2'
RUN_DATE = '2026-09-27'
VARIANT_ORDER = ('S2', 'S2TZ')
ROLE_ORDER = ('bg', 'nc', 'tgt')
GRID = dict(nx_ny_nz=[1, 1280, 2000], dx_dy_dz=[0.025, 0.025, 0.025])
FLAT_ANCHOR = {  # cross-geometry diagnostic only (batch2d_v1, archived runs of 2026-09-26)
    'BG': 'B2D-C3m-BG',
    'TGT': 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02',
}
# NC attempt 1 was killed by a harness-level task timeout after the solver had finished
# (h5 written, supervision.json never finalized); attempt 2 ran the identical hash-locked
# input under the same authorization (gate reactivation_note 2026-09-27). Analysis uses att2.
RUN_DIR_OVERRIDES = {'B2D-C3mS2-NC': '2026-09-27_B2D-C3mS2-NC_att2'}
NC_ATT1_DIR = '2026-09-27_B2D-C3mS2-NC'  # determinism comparison only (Ex samples vs att2)
NC_EXPECTATION = '预期近零，偏离指示几何/平滑残差'


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_supervision(path):
    rec = json.loads(path.read_text(encoding='utf-8'))
    out = dict(reason=rec['reason'], exit_code=rec['exit_code'], wall_s=float(rec['wall_s']),
               peak_job_commit_bytes=int(rec['peak_job_commit_bytes']),
               memory_semantics=rec.get('memory_semantics', 'undeclared'))
    if rec['reason'] != 'completed':
        out['note'] = 'supervision reason is not "completed"; archived as-is, results reported unchanged'
    return out


def load_case(run_dir, run_id, component, rx_m, dt_s, time_steps):
    h5_path = run_dir / (run_id + '.h5')
    sup_path = run_dir / 'supervision.json'
    if not h5_path.is_file():
        return None, 'h5 not found in %s' % run_dir.as_posix()
    if not sup_path.is_file():
        return None, 'supervision.json not found in %s' % run_dir.as_posix()
    src = load_source(h5_path)
    rx = load_receiver(h5_path, receiver_path='name:measurement', component=component)
    with h5py.File(h5_path, 'r') as h:
        assert np.array_equal(h.attrs['nx_ny_nz'], GRID['nx_ny_nz']), (run_id, h.attrs['nx_ny_nz'])
        assert np.allclose(h.attrs['dx_dy_dz'], GRID['dx_dy_dz'], rtol=0, atol=1e-12), run_id
        assert np.allclose(h[rx.path].parent.attrs['Position'][1:], rx_m, rtol=0, atol=1e-12), run_id
        dtype = h[rx.path].dtype
        iterations = int(h.attrs['Iterations'])
    assert dtype == np.float64, run_id
    assert np.isclose(src.dt, rx.dt, rtol=0, atol=0), run_id
    assert np.isclose(rx.dt, dt_s, rtol=1e-9, atol=0), (run_id, rx.dt)
    assert src.time_offset == src.dt / 2 and rx.time_offset == 0, run_id
    assert len(rx.samples) == iterations == time_steps, (run_id, len(rx.samples))
    assert np.all(np.isfinite(rx.samples)) and np.all(np.isfinite(src.samples)), run_id
    n = min(len(rx.samples), int(np.floor(TIME_WINDOW_S / rx.dt)) + 1)
    return dict(path=h5_path, supervision=read_supervision(sup_path), src=src, rx=rx, dt=rx.dt,
                n=n, t=rx.times[:n], iterations=iterations), None


def case_spectra(obj):
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
    e_diff = float(np.sum(np.abs(diff) ** 2))
    e_bg = float(np.sum(np.abs(bg) ** 2))
    assert e_bg > 0.0, 'reference in-band energy is zero'
    ratio = e_diff / e_bg
    return dict(in_band_energy_ratio_DIFF_over_REF=float(ratio),
                in_band_energy_ratio_DIFF_over_REF_dB=float(10.0 * np.log10(ratio)))


def difference_metrics(diff, ref_main, diff_t, t, ref_label):
    mag = np.abs(diff[MAIN])
    if float(mag.max()) == 0.0 and float(np.abs(diff[ROBUST]).max()) == 0.0 \
            and float(np.abs(diff_t).max()) == 0.0:
        return dict(
            grid_tier='BASE', time_diff_peak_V_m=0.0, envelope_peak_time_ns=None,
            spectral_peak_frequency_MHz=None, spectral_peak_phase_deg=None,
            robustness_400ns_spectral_peak_frequency_MHz=None,
            robustness_400ns_relative_L2_vs_200ns=None,
            normalized_mag_dB_at_selected_MHz=None, zero_difference=True,
            in_band_energy_ratio_DIFF_over_REF=0.0, in_band_energy_ratio_DIFF_over_REF_dB=None,
            note='paired difference is identically zero on both tail windows and in time; '
                 'peak/phase/dB metrics not applicable, reported as null')
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
        reference_label=ref_label,
        note='time_diff_peak_V_m carries the pair-internal excitation scale; no absolute detectability claim')
    out.update(band_energy(diff[MAIN], ref_main))
    return out


def envelope_peak(obj):
    env = np.abs(hilbert(obj['rx'].samples[:obj['n']]))
    i = int(np.argmax(env))
    return dict(envelope_peak_time_ns=float(obj['t'][i] * 1e9),
                envelope_peak_V_m=float(env[i]),
                time_peak_abs_V_m=float(np.max(np.abs(obj['rx'].samples[:obj['n']]))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--runs-root', type=Path, default=Path('artifacts/simulations'))
    p.add_argument('--flat-runs-root', type=Path, default=Path('artifacts/simulations'),
                   help='parent of the archived 2026-09-26 flat-anchor run directories')
    p.add_argument('--cases', type=Path, default=Path('configs/research/batch2d_slope_t2/cases.json'))
    p.add_argument('--groups', type=Path, default=Path('configs/research/batch2d_slope_t2/groups.json'))
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)

    cases_doc = json.loads(a.cases.read_text(encoding='utf-8'))
    groups_doc = json.loads(a.groups.read_text(encoding='utf-8'))
    assert cases_doc['batch_id'] == groups_doc['batch_id'] == BATCH, 'wrong batch'
    assert cases_doc['grid_tier'] == groups_doc['grid_tier'] == 'BASE'
    assert len(cases_doc['cases']) == cases_doc['n_cases'] == 6

    cases = sorted(cases_doc['cases'],
                   key=lambda c: (VARIANT_ORDER.index(c['family'].replace('C3', '')),
                                  ROLE_ORDER.index(c['role'])))
    for c in cases:
        assert c['grid_tier'] == 'BASE', c['run_id']
        seed_txt = '|'.join([BATCH, c['group_id'], c['case_id'], c['variant_tag']])
        assert c['seed'] == int.from_bytes(hashlib.sha256(seed_txt.encode('utf-8')).digest()[:8], 'big'), c['run_id']

    loaded, unavailable = {}, []
    for c in cases:
        run_dir = a.runs_root / RUN_DIR_OVERRIDES.get(c['run_id'], f"{RUN_DATE}_{c['run_id']}")
        obj, reason = load_case(run_dir, c['run_id'],
                                c['geometry']['component'], c['geometry']['rx_m'],
                                c['dt_s'], c['time_steps'])
        if obj is None:
            unavailable.append(dict(run_id=c['run_id'], family=c['family'], role=c['role'],
                                    grid_tier='BASE', reason=reason))
            continue
        loaded[c['run_id']] = obj
    for rid, obj in loaded.items():
        spec, tapers = case_spectra(obj)
        obj['spec'], obj['tapers'] = spec, tapers

    # flat anchor (cross-geometry diagnostic only; never a pairing reference for slope cases)
    flat = {}
    for role, rid in FLAT_ANCHOR.items():
        obj, reason = load_case(a.flat_runs_root / f'2026-09-26_{rid}', rid, 'Ex', [16.65, 45.0],
                                5.896635841874211e-11, 20352)
        assert obj is not None, f'flat anchor {rid}: {reason}'
        spec, tapers = case_spectra(obj)
        obj['spec'], obj['tapers'] = spec, tapers
        flat[role] = obj

    per_case, diffs, diff_t, diff_order = {}, {}, {}, []
    bg_by_variant = {}
    for c in cases:
        if c['role'] == 'bg':
            bg_by_variant[c['family']] = c['run_id']
    for c in cases:
        rid = c['run_id']
        entry = dict(run_id=rid, family=c['family'], role=c['role'], group_id=c['group_id'],
                     seed=c['seed'], variant_tag=c['variant_tag'], grid_tier='BASE',
                     cover_thickness_center_m=c['cover_thickness_center_m'],
                     theta_eff_deg=c['scan_axis']['theta_eff_deg'],
                     transition_zone=c['scan_axis']['transition_zone'],
                     scan_axis={k: v for k, v in c['scan_axis'].items()})
        if rid not in loaded:
            entry['analysis_status'] = 'not_analyzed: %s' % [u['reason'] for u in unavailable
                                                             if u['run_id'] == rid][0]
            entry['difference'] = dict(computed=False)
            per_case[rid] = entry
            continue
        obj = loaded[rid]
        entry['analysis_status'] = 'analyzed'
        entry['h5'] = dict(path=obj['path'].as_posix(), sha256=sha256_file(obj['path']),
                           nx_ny_nz=GRID['nx_ny_nz'], dx_dy_dz=GRID['dx_dy_dz'],
                           dt_s=float(obj['dt']), iterations=obj['iterations'],
                           transformed_samples=int(obj['n']),
                           receiver_component=c['geometry']['component'],
                           receiver_position_m=list(c['geometry']['rx_m']), dtype='float64',
                           all_finite=True)
        entry['supervision'] = obj['supervision']
        entry['tail_taper_fraction'] = obj['tapers']
        entry['direct_wave_region'] = envelope_peak(obj)
        bg_rid = bg_by_variant[c['family']]
        if c['role'] == 'bg':
            entry['difference'] = dict(computed=False, reference_role='bg',
                                       note='variant BG is the differencing reference itself (spec 6-P7)')
        elif bg_rid not in loaded:
            entry['difference'] = dict(computed=False, bg_run_id=bg_rid,
                                       note='variant BG unavailable; no cross-variant substitution')
        else:
            b = loaded[bg_rid]
            assert np.array_equal(obj['src'].samples, b['src'].samples), rid
            assert np.isclose(obj['dt'], b['dt'], rtol=0, atol=0) and obj['n'] == b['n'], rid
            dt_pair = {tag: obj['spec'][tag] - b['spec'][tag] for tag in TAIL_NS}
            dt_time = obj['rx'].samples[:obj['n']] - b['rx'].samples[:obj['n']]
            entry['difference'] = dict(computed=True, bg_run_id=bg_rid, grid_tier='BASE',
                                       source_samples_identical_to_variant_bg=True)
            entry['difference'].update(difference_metrics(dt_pair, b['spec'][MAIN], dt_time, obj['t'],
                                                          ref_label='same-variant BG'))
            diffs[rid], diff_t[rid] = dt_pair, dt_time
            diff_order.append(rid)
        per_case[rid] = entry

    # ---- NC attempt-1 vs attempt-2 determinism check (raw Ex samples) ----
    nc_determinism = None
    att1_h5 = a.runs_root / NC_ATT1_DIR / 'B2D-C3mS2-NC.h5'
    if att1_h5.is_file() and 'B2D-C3mS2-NC' in loaded:
        rx1 = load_receiver(att1_h5, receiver_path='name:measurement', component='Ex')
        att2_samples = loaded['B2D-C3mS2-NC']['rx'].samples
        same = bool(np.array_equal(rx1.samples, att2_samples))
        nc_determinism = dict(
            att1_h5=att1_h5.as_posix(), att1_h5_sha256=sha256_file(att1_h5),
            att2_h5=loaded['B2D-C3mS2-NC']['path'].as_posix(),
            ex_samples_identical=same,
            max_abs_diff_V_m=0.0 if same else float(np.max(np.abs(rx1.samples - att2_samples))),
            note='attempt 1 supervision record lost (harness task timeout, not a solver failure); '
                 'attempt 2 re-ran the identical hash-locked input; identical Ex samples would be '
                 'consistent with a deterministic solver build, not a proof of correctness')

    # ---- cross-variant background difference: BG_S2TZ - BG_S2 (transition-zone effect) ----
    cross_variant = None
    rid_a, rid_b = bg_by_variant.get('C3S2'), bg_by_variant.get('C3S2TZ')
    if rid_a in loaded and rid_b in loaded:
        oa, ob = loaded[rid_a], loaded[rid_b]
        assert np.array_equal(oa['src'].samples, ob['src'].samples)
        d = {tag: ob['spec'][tag] - oa['spec'][tag] for tag in TAIL_NS}
        dt_time = ob['rx'].samples[:ob['n']] - oa['rx'].samples[:oa['n']]
        cross_variant = dict(
            computed=True, minus=rid_b, plus_reference=rid_a,
            semantics='BG_S2TZ - BG_S2: transition-zone effect on the background/interface response; '
                      'same mother except the 1 m tzone strips',
            **difference_metrics(d, oa['spec'][MAIN], dt_time, oa['t'], ref_label='BG_S2'))
        diffs['XSVTZ_BG_S2TZ_minus_BG_S2'] = d
        diff_t['XSVTZ_BG_S2TZ_minus_BG_S2'] = dt_time
        diff_order.append('XSVTZ_BG_S2TZ_minus_BG_S2')

    # ---- cross-geometry diagnostic vs flat anchor (NOT a paired difference) ----
    cross_geometry = {}
    for tag, rid in (('S2_BG', bg_by_variant.get('C3S2')), ('S2TZ_BG', bg_by_variant.get('C3S2TZ'))):
        if rid not in loaded:
            continue
        obj = loaded[rid]
        fb = flat['BG']
        rel_l2 = float(np.linalg.norm(obj['spec'][MAIN] - fb['spec'][MAIN])
                       / np.linalg.norm(fb['spec'][MAIN]))
        cross_geometry[tag] = dict(
            slope_run_id=rid, flat_run_id='B2D-C3m-BG',
            relative_L2_vs_flat_BG_200ns=rel_l2,
            envelope_peak_time_ns_slope=envelope_peak(obj)['envelope_peak_time_ns'],
            envelope_peak_time_ns_flat=envelope_peak(fb)['envelope_peak_time_ns'],
            note='cross-geometry diagnostic, NOT a paired difference: cover center differs '
                 '(flat 3.0 m vs T2 4.2 m) in addition to the interface tilt; differences confound '
                 'slope angle with cover thickness (cases.json naming_caveat)')
    flat_tgt = flat['TGT']
    cross_geometry['flat_anchor_reference'] = dict(
        flat_tgt_envelope_peak_time_ns=envelope_peak(flat_tgt)['envelope_peak_time_ns'],
        note='flat C3 anchor TGT direct-region envelope peak, archived for orientation only')

    nc_rows = []
    for c in cases:
        if c['role'] != 'nc':
            continue
        rid = c['run_id']
        row = dict(run_id=rid, family=c['family'], group_id=c['group_id'], grid_tier='BASE')
        diff = per_case[rid]['difference']
        if diff.get('computed'):
            row.update(in_band_energy_ratio_DIFF_over_REF=diff['in_band_energy_ratio_DIFF_over_REF'],
                       in_band_energy_ratio_DIFF_over_REF_dB=diff['in_band_energy_ratio_DIFF_over_REF_dB'],
                       envelope_peak_time_ns=diff['envelope_peak_time_ns'],
                       spectral_peak_frequency_MHz=diff['spectral_peak_frequency_MHz'],
                       computed=True)
        else:
            row.update(computed=False, note=diff.get('note'))
        row['expectation'] = NC_EXPECTATION
        nc_rows.append(row)

    inputs = {loaded[c['run_id']]['path'].as_posix(): sha256_file(loaded[c['run_id']]['path'])
              for c in cases if c['run_id'] in loaded}
    for role, obj in flat.items():
        inputs[obj['path'].as_posix()] = sha256_file(obj['path'])
    inputs[a.cases.as_posix()] = sha256_file(a.cases)
    inputs[a.groups.as_posix()] = sha256_file(a.groups)

    result = dict(
        batch=BATCH, segment='base', grid_tier='BASE',
        spec=('docs/research/2026-09-26_batch_2d_spec_v1.md 6-P2/P6/P7 applied per slope variant; '
              'docs/research/2026-09-27_slope_family_design_draft.md (kimi-accepted draft)'),
        frequency_Hz='linspace(20e6,170e6,501)',
        time_window_ns=1200, tail_windows_ns=dict(main=200, robustness=400),
        n_cases_contract=6, n_cases_analyzed=len(loaded),
        pairing_policy=dict(
            rule='every non-BG case is differenced against the same-variant same-grid BG only',
            cross_variant_differencing=False,
            cross_geometry_vs_flat='diagnostic only, never a paired difference (cover-center confound)',
            variant_bg={v: bg_by_variant.get(f'C3{v}') for v in VARIANT_ORDER}),
        cases=per_case,
        cross_variant_background_difference=cross_variant,
        cross_geometry_diagnostic_vs_flat=cross_geometry,
        nc_attempt_determinism=nc_determinism,
        negative_controls=dict(
            note='negative controls are single-listed and never merged into the target rows',
            nc=nc_rows),
        unavailable_cases=unavailable,
        hard_limits=dict(
            grid_tier='BASE', grid_tiers_merged=False,
            cross_variant_differencing=False, variant_bg_required=True,
            negative_controls_merged_into_targets=False,
            reference_state='numerically_unresolved',
            physical_acceptance_threshold=None,
            physical_label_eligible=False, training_eligible=False,
            training_labels_generated=False, clean_truth_generated=False,
            grid_convergence_certified=False, absolute_accuracy_claimed=False,
            in_band_energy_ratio_is='diagnostic energy ratio on the 501-point analysis grid; '
                                    'not detectability, not SNR, not a physical threshold',
            conclusions_scope='mechanism and operator-evaluation evidence for the batch2d_slope_t2 '
                              'tier-T2 slope variants only; not field performance, not 3D, '
                              'not absolute accuracy, no cross-tier slope-angle conclusion',
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
    for role, obj in flat.items():
        arrays[f'resp_flat_{role}_200ns'] = obj['spec'][MAIN]
    np.savez_compressed(a.output / 'arrays.npz', **arrays)

    with (a.output / 'difference_spectra.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['frequency_Hz'] + [f'{rid}_{tag}_{part}'
                                       for rid in diff_order for tag in TAIL_NS for part in ('re', 'im')])
        w.writerows(zip(FREQ, *[c for rid in diff_order for tag in TAIL_NS
                                for c in (diffs[rid][tag].real, diffs[rid][tag].imag)]))
    print(json.dumps(dict(batch=BATCH, grid_tier='BASE', n_cases_analyzed=len(loaded),
                          n_differences_computed=len(diff_order), n_unavailable=len(unavailable),
                          cross_variant_computed=bool(cross_variant and cross_variant.get('computed')),
                          output=a.output.as_posix()), sort_keys=True))


if __name__ == '__main__':
    main()
