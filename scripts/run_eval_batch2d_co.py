"""Evaluation runner (S2/S3, common-offset variant): run the 27-entry operator
catalogue on the batch2d_v1_co common-offset gathers against the frozen event table.

Implements the common-offset B-scan semantics (the core question this batch exists
to answer, vs the variable-offset MT gather):
  * The 11 traces are 11 independent single-trace simulations at Rx y = 12.65..22.65 m
    (spacing 1.0 m), each with Tx = Rx - 1.30 m: a CONSTANT Tx-Rx offset of 1.30 m.
    Unlike batch2d_v1_mt, the trace axis IS a spatial profile axis.
  * Event windows are the SAME frozen sample-axis intervals as the MT evaluation
    (anchor-geometry two-way times +/- 80 ns, inward rounding) and are applied
    UNshifted to every trace. No per-trace alignment is performed: windows were
    frozen before any candidate ran, and per-trace shifting would redefine them.
  * With constant offset, the direct-coupled wave is aligned across traces; the
    background is a layered half-space, so the gather is y-invariant up to the
    target box (analysis: pairwise bit-equality within sets; anchor t05 bit-for-bit
    equal to the archived single-trace mothers, 2/2).
  * Layering/stratigraphy is laterally invariant in these models: the mean/SVD
    operators therefore act on a spatially coherent B-scan, which is exactly the
    operating regime the operator catalogue was designed for; this complements
    the MT gather (unaligned direct wave) rather than replacing it.

Gate discipline (design v0.1 §A.3), identical to the MT runner:
  1. Event table must be status=frozen (verified by file SHA-256).
  2. Input must satisfy min(shape) >= 2 (11 traces; single-trace gate now
     passes; recorded as single_trace_gate_blocked=false).
  3. BG-order gain configs persist the true pre-gain result; GB-order only
     persists G^-1 Y, labeled as a non-true intermediate. Negative-control
     rows are reported separately and never merged into target rows.
  4. Gather assembly is verified per trace: rx Position y == 12.65+k,
     src y == rx-1.30, 20352 finite samples; anchor trace t05 is asserted
     bit-for-bit equal to the archived single-trace mother h5 (2/2 mothers).

Metrics: with reference_state=numerically_unresolved, waveform metrics D/A/H
and rho are NOT computible (recorded as unavailable with the verbatim gate
reason); R_c is never computed (no legitimate pure-clutter window). The only
v0.2 §4 metric computed is the negative-control residual N_b (energy_metrics);
gain risk diagnostics and resource usage are recorded. Feasibility labels use
configuration_labels with all-null limits -> undetermined/partial_only.

Heavy arrays (pre-gain results and GB G^-1 Y audits) persist as .npy under
artifacts/simulations/2026-09-26_eval_batch2d_co_arrays/ (git-ignored); records
and manifest under artifacts/research_checks/<run_dir>/ carry path + SHA-256.

Read-only inputs: archived CO h5 files, archived single-trace mother h5 files,
frozen event table, operator contract. The runner never writes into input
directories and never invokes a solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import tracemalloc

import h5py
import numpy as np

from research_operator_contract import apply_configuration, catalogue, VERSION as OP_VERSION
from research_evaluation_contract import configuration_labels, energy_metrics, waveform_metrics

ROOT = Path(__file__).resolve().parents[1]
EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'
CO_ANALYSIS = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_co_analysis/results.json'
CO_ARCHIVE = 'artifacts/research_checks/2026-09-26_{rid}/'  # evidence dirs hold the h5 copies
ARRAYS_DIR = ROOT / 'artifacts/simulations/2026-09-26_eval_batch2d_co_arrays'
N_TRACES = 11
N_SAMPLES = 20352
RX_Y = [12.65 + k for k in range(N_TRACES)]  # t01..t11
OFFSET_M = 1.30
GATHER_SEMANTICS_VERSION = 'common_offset_bscan_v1'

# The CO batch simulated exactly two mother models (22 runs = 2 x 11 traces).
# Every event-table entry whose mother is not one of these has no common-offset
# case and is recorded unavailable (no_common_offset_case_in_batch2d_v1_co).
CO_CASES = {
    'B2D-C3m-BG': 'BG',
    'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02': 'TGT',
}

HARD_LIMITS = {
    "absolute_accuracy_claimed": False,
    "clean_truth_generated": False,
    "conclusions_scope": "operator-evaluation evidence for the batch2d_v1_co common-offset subset (C3 family only) only; not field performance, not 3D, not absolute accuracy, no maximum-depth claim",
    "cross_family_differencing": False,
    "cross_tier_differencing": False,
    "family_bg_required": True,
    "fdtd_solver_invoked": False,
    "grid_convergence_certified": False,
    "grid_tier": "BASE",
    "grid_tiers_merged": False,
    "in_band_energy_ratio_is": "diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold",
    "negative_controls_merged_into_targets": False,
    "physical_acceptance_threshold": None,
    "physical_label_eligible": False,
    "reference_state": "numerically_unresolved",
    "training_eligible": False,
    "training_labels_generated": False,
    "n_traces": 11,
    "trace_axis_semantics": "common-offset b-scan: 11 independent single-trace runs, Rx y 12.65->22.65 m spacing 1.0 m, Tx = Rx-1.30 m (constant offset); trace axis IS a spatial profile axis",
    "anchor_trace": "t05 (y16.65 m) bit-for-bit equal to archived single-trace mother (2/2 mothers)",
    "nc_multitrace_control": "no NC/OFF mother simulated in batch2d_v1_co; all nc_zero/off_path rows unavailable; layered BG makes BG-NEG windows the only in-batch negative rows",
    "gather_assembly": "11 separate single-trace h5 stacked to [20352, 11]; per-trace layout asserted (rx y, src y=rx-1.30, finite); anchor t05 asserted bit-for-bit vs archived mother h5",
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_case(case_key):
    """Assemble and verify one [20352, 11] common-offset gather from 11 archived
    single-trace h5 files. Returns (gather, {run_id: (rel_path, sha)}, mother_check)."""
    if case_key == 'BG':
        mother = 'B2D-C3m-BG'
    else:
        mother = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02'
    cols, files = [], {}
    for k in range(1, N_TRACES + 1):
        rid = f'{mother}-CO11-t{k:02d}'
        h5_path = ROOT / CO_ARCHIVE.format(rid=rid) / f'{rid}.h5'
        with h5py.File(h5_path, 'r') as f:
            g = f['rxs/rx1']
            pos = np.asarray(g.attrs.get('Position'), dtype=float)
            src_pos = np.asarray(f['srcs/src1'].attrs['Position'], dtype=float)
            x = g['Ex'][:].astype(np.float64)
        # per-trace layout gate (same checks as analyze_batch2d_co.py)
        assert x.shape == (N_SAMPLES,), f'{rid}: unexpected sample count {x.shape}'
        assert abs(pos[1] - RX_Y[k - 1]) <= 1e-9, f'{rid}: rx y {pos[1]} != {RX_Y[k - 1]}'
        assert abs(src_pos[1] - (RX_Y[k - 1] - OFFSET_M)) <= 1e-9, f'{rid}: src y mispositioned'
        assert np.isfinite(x).all(), f'{rid}: non-finite samples'
        cols.append(x)
        files[rid] = (str(h5_path.relative_to(ROOT)), sha256_file(h5_path))
    gather = np.stack(cols, axis=1)
    # anchor gate: t05 column bit-for-bit equal to archived single-trace mother
    mother_h5 = ROOT / CO_ARCHIVE.format(rid=mother) / f'{mother}.h5'
    with h5py.File(mother_h5, 'r') as f:
        x_mother = f['rxs/rx1/Ex'][:].astype(np.float64)
    anchor_max_diff = float(np.max(np.abs(gather[:, 4] - x_mother)))
    assert anchor_max_diff == 0.0, f'{mother}: anchor t05 not bit-for-bit (max {anchor_max_diff})'
    mother_check = {
        'mother_h5': str(mother_h5.relative_to(ROOT)),
        'mother_h5_sha256': sha256_file(mother_h5),
        'anchor_trace': 't05',
        'anchor_max_abs_difference': anchor_max_diff,
    }
    return gather, files, mother_check


def persist_array(tag, array):
    ARRAYS_DIR.mkdir(parents=True, exist_ok=True)
    path = ARRAYS_DIR / f'{tag}.npy'
    np.save(path, array)
    return str(path.relative_to(ROOT)), sha256_file(path)


def run_eval(run_tag):
    started = time.perf_counter()
    # --- gate 1: frozen event table, verified by whole-file SHA-256 ---------
    table_sha = sha256_file(EVENT_TABLE)
    assert table_sha == EVENT_TABLE_SHA256, 'event table hash mismatch (gate 1)'
    table = json.loads(EVENT_TABLE.read_text(encoding='utf-8'))
    assert table['status'] == 'frozen' and table['n_entries'] == len(table['entries'])
    dt_s = table['time_axis']['dt_s']

    co = json.loads(CO_ANALYSIS.read_text(encoding='utf-8'))
    assert co['batch'] == 'batch2d_v1_co' and co['n_traces'] == N_TRACES
    assert co['all_layout_ok'] is True
    assert all(co['anchor_t05_bit_for_bit'].values()), 'CO analysis anchor gate failed'

    configs = catalogue()
    assert len(configs) == 27, 'catalogue must hold the 27 frozen candidates'
    op_script_sha = sha256_file(ROOT / 'scripts/research_operator_contract.py')
    eval_script_sha = sha256_file(ROOT / 'scripts/run_eval_batch2d_co.py')

    # --- assemble the two gathers (BG and TGT) with layout + anchor gates ----
    cases = {}
    for case_key in ('BG', 'TGT'):
        cases[case_key] = load_case(case_key)

    # --- expand event rows to (event, CO case) evaluation rows --------------
    rows = []  # (event, case_key or None)
    for e in table['entries']:
        if e.get('applies_to') == 'all_cases_of_family':
            fam_prefix = 'B2D-' + e['family'].capitalize() + 'm'  # c3 -> B2D-C3m
            matched = [m for m in sorted(CO_CASES) if m.startswith(fam_prefix + '-')]
            if matched:
                for m in matched:
                    rows.append((e, CO_CASES[m]))
            else:
                rows.append((e, None))
        else:
            rows.append((e, CO_CASES.get(e['mother_model_id'])))

    records = []
    availability = {}
    reason_hist = {}
    tracemalloc.start()
    mem_base = tracemalloc.get_traced_memory()[1]
    for e, case_key in rows:
        evaluable = case_key is not None
        x_window = None
        case_files = case_sha = None
        if evaluable:
            x, case_files, mother_check = cases[case_key]
            lo, hi = int(e['index_lo']), int(e['index_hi'])
            assert 0 <= lo <= hi and hi + 1 <= x.shape[0], f"window {e['event_id']} outside gather"
            x_window = x[lo:hi + 1, :]
            assert np.isfinite(x_window).all() and x_window.shape[1] == N_TRACES
            case_sha = {rid: sha for rid, (_, sha) in case_files.items()}
        for cfg in configs:
            rec = {
                'record_id': f"{e['event_id']}::{case_key or 'NO_CO_CASE'}::{cfg['id']}",
                'event_id': e['event_id'], 'event_role': e['role'],
                'case': case_key, 'family': e['family'],
                'mother_model_id': e['mother_model_id'],
                'candidate_id': cfg['id'], 'background': cfg['background'],
                'parameter': cfg['parameter'], 'end_gain': cfg['end_gain'],
                'requested_order': cfg['order'],
                'row_type': 'negative_control' if e['role'] in
                            ('nc_zero', 'bg_absent', 'off_path') else 'event',
                'window': {'index_lo': int(e['index_lo']) if evaluable else None,
                           'index_hi': int(e['index_hi']) if evaluable else None,
                           'n_window_samples': e['n_window_samples'],
                           'window_hash': e['window_hash'], 'mask_hash': e['mask_hash'],
                           'time_axis_convention_version': e['time_axis_convention_version']},
                'gather_semantics': {
                    'version': GATHER_SEMANTICS_VERSION,
                    'n_traces': N_TRACES if evaluable else None,
                    'trace_axis_is_spatial_profile': True,
                    'window_applied_unshifted_to_all_traces': True,
                    'trace_offset_range_m': [OFFSET_M, OFFSET_M],
                    'trace_spacing_m': 1.0,
                    'rx_y_range_m': [RX_Y[0], RX_Y[-1]],
                    'anchor_trace': 't05',
                },
                'provenance': {
                    'case_h5_files': case_files,
                    'case_h5_sha256': case_sha,
                    'anchor_check': None if not evaluable else mother_check,
                    'event_table_sha256': table_sha,
                    'operator_script_sha256': op_script_sha,
                    'eval_script_sha256': eval_script_sha,
                    'operator_version': OP_VERSION,
                    'dt_s': dt_s,
                },
            }
            if not evaluable:
                rec.update(availability='unavailable', failure_reason='no_common_offset_case_in_batch2d_v1_co',
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] = availability.get('unavailable', 0) + 1
                reason_hist['no_common_offset_case_in_batch2d_v1_co'] = reason_hist.get('no_common_offset_case_in_batch2d_v1_co', 0) + 1
                records.append(rec)
                continue
            t0 = time.perf_counter()
            try:
                res = apply_configuration(x_window, cfg['id'])
                failure = None
            except ValueError as exc:
                rec.update(availability='unavailable', failure_reason=str(exc),
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] = availability.get('unavailable', 0) + 1
                reason_hist[str(exc)] = reason_hist.get(str(exc), 0) + 1
                records.append(rec)
                continue
            except Exception as exc:  # ConfigUnavailable family and numeric failures
                rec.update(availability='unavailable', failure_reason=str(exc),
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] = availability.get('unavailable', 0) + 1
                reason_hist[str(exc)] = reason_hist.get(str(exc), 0) + 1
                records.append(rec)
                continue
            wall = time.perf_counter() - t0
            peak = tracemalloc.get_traced_memory()[1] - mem_base
            availability['ran'] = availability.get('ran', 0) + 1

            y = res['output']
            mask = np.ones(y.shape, dtype=bool)
            nb = energy_metrics(y, x_window, mask)
            wf = waveform_metrics(y, x_window, mask, reference_kind='paired_contrast',
                                  state='numerically_unresolved', scope='contrast')
            svd_diag = next((s['diagnostics'] for s in res['steps']
                             if s['diagnostics'].get('kind') == 'svd'), None)
            gain_diag = None
            if cfg['end_gain'] > 1:
                pre = res['pre_gain'] if cfg['order'] == 'BG' else res['audit_before_gain']
                kind = 'pre_gain_true_intermediate' if cfg['order'] == 'BG' else 'ginv_audit_not_true_intermediate'
                path, fsha = persist_array(f"{e['event_id']}__{case_key}__{cfg['id']}", pre)
                curve = np.asarray(res['gain_curve'])
                cpath, csha = persist_array(f"{e['event_id']}__{case_key}__{cfg['id']}__gaincurve", curve)
                gain_diag = {
                    'persisted_array': path, 'persisted_array_sha256': fsha,
                    'persisted_kind': kind,
                    'max_gain': cfg['end_gain'],
                    'max_abs_input': float(np.max(np.abs(x_window))),
                    'max_abs_pre_gain': float(np.max(np.abs(pre))),
                    'max_abs_output': float(np.max(np.abs(y))),
                    'gain_energy_ratio_output_over_pre': float(np.mean(y * y) / np.mean(pre * pre)),
                    'clipped_samples': int(np.count_nonzero(np.abs(y) > cfg['end_gain'] * np.max(np.abs(pre)) + 1e-300)),
                    'overflow': bool(not np.isfinite(y).all()),
                    'gain_curve_ref': cpath, 'gain_curve_sha256': csha,
                    'gain_curve_first': float(curve[0]), 'gain_curve_end': float(curve[-1]),
                }
            rec.update(
                availability='ran', failure_reason=failure,
                executed_order=''.join(s['operator'] for s in res['steps']),
                fell_back=False, fallback_reason=None,
                metrics={
                    'N_b': nb,
                    'waveform': {'available': wf['available'], 'reason': wf['reason'],
                                 'metrics': None if wf['metrics'] is None else
                                            {k: v for k, v in wf['metrics'].items()}},
                },
                svd_diagnostics=svd_diag,
                gain_diagnostics=gain_diag,
                resource={'wall_s': wall, 'peak_tracemalloc_bytes': int(peak),
                          'compute_dtype': 'float64', 'input_dtype': res['input_dtype'],
                          'n_samples': int(y.shape[0]), 'n_traces': int(y.shape[1])},
            )
            records.append(rec)
    tracemalloc.stop()

    # --- feasibility with all-null limits (G4 unresolved) -------------------
    cand_records = [{'id': r['record_id'], 'available': r['availability'] == 'ran',
                     'metrics': {'nrmse': None, 'amplitude_error': None,
                                 'negative_control_residual': None},
                     'failure_reason': r['failure_reason']} for r in records]
    labels = configuration_labels(cand_records, limits={'nrmse': None, 'amplitude_error': None,
                                                        'negative_control_residual': None},
                                  objectives=['nrmse', 'amplitude_error'])

    manifest = {
        'run_tag': run_tag,
        'runner': 'scripts/run_eval_batch2d_co.py',
        'event_table': str(EVENT_TABLE.relative_to(ROOT)),
        'event_table_sha256': table_sha,
        'time_axis_convention_version': table['time_axis']['convention_version'],
        'grid_tier': 'BASE',
        'gather_semantics_version': GATHER_SEMANTICS_VERSION,
        'trace_axis_semantics': HARD_LIMITS['trace_axis_semantics'],
        'window_mapping': 'frozen sample-axis window applied unshifted to all 11 traces; no per-trace alignment',
        'candidates': {'n': len(configs), 'version': OP_VERSION,
                       'operator_script_sha256': op_script_sha,
                       'eval_script_sha256': eval_script_sha},
        'co_cases': {ck: {'mother_check': cs[2],
                          'h5': {rid: {'path': p, 'sha256': s}
                                 for rid, (p, s) in sorted(cs[1].items())}}
                     for ck, cs in sorted(cases.items())},
        'event_rows': len(rows),
        'n_records': len(records),
        'availability': availability,
        'failure_reason_histogram': reason_hist,
        'gates': {'event_table_frozen': True, 'min_shape_ge_2': True,
                  'single_trace_gate_blocked': False,
                  'blocked_candidate_count': 0,
                  'pre_gain_persisted': True,
                  'negative_controls_separate': True,
                  'gather_layout_verified': True,
                  'anchor_bit_for_bit': True},
        'runner_limits': {'r_c_computed': False, 'isolated_event_granted': False,
                          'cross_family_metrics': False, 'gain_reference_absent': True,
                          'waveform_metrics_computed': False,
                          'waveform_metrics_reason': 'reference_state numerically_unresolved + paired_contrast (design A.4)'},
        'feasibility': {'status': labels['status'],
                        'note': 'all limits null (G4 unresolved): every gate undetermined; no ranking, no unique label'},
        'negative_control_events': sorted({r['event_id'] for r in records if r['row_type'] == 'negative_control'}),
        'hard_limits': HARD_LIMITS,
        'total_wall_s': time.perf_counter() - started,
        'stop_reason': None,
        'arrays_dir': str(ARRAYS_DIR.relative_to(ROOT)),
    }
    return records, manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True,
                   help='must not exist (no overwriting of evidence)')
    p.add_argument('--run-tag', default='r1')
    a = p.parse_args()
    out = (ROOT / a.output_dir).resolve() if not a.output_dir.is_absolute() else a.output_dir
    assert not out.exists(), f'refusing to overwrite {out}'
    records, manifest = run_eval(a.run_tag)
    out.mkdir(parents=True)
    (out / 'records.json').write_text(json.dumps({'records': records}, indent=1) + '\n', encoding='utf-8')
    (out / 'run_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({'records': len(records), 'availability': manifest['availability'],
                      'feasibility': manifest['feasibility']['status'], 'output': str(out)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
