"""Evaluation runner (S2/S3): run the 27-entry operator catalogue on the
batch2d_v1_mt variable-offset multitrace gathers against the frozen event table.

Implements the variable-offset gather semantics:
  * The 33 traces are a single-shot multi-receiver gather (Tx fixed y15.35,
    Tx-Rx offset 2.70 -> 1.30 -> 5.30 m along trace index). This is NOT a
    common-offset B-scan; the trace axis is not a spatial profile axis.
  * Event windows are frozen sample-axis intervals (anchor-geometry two-way
    times +/- 80 ns, inward rounding) and are applied UNshifted to every
    trace. No per-trace alignment is performed: windows were frozen before
    any candidate ran, and per-trace shifting would redefine the windows.
  * Consequence recorded per record: background/event coherence across the
    trace axis is weakened by the offset variation (direct-coupled wave
    differential ~147 samples end-to-end, proposal §2.4 self-computed).
    Mean/SVD operators therefore see a gather, not a flattened B-scan.

Gate discipline (design v0.1 §A.3):
  1. Event table must be status=frozen (verified by file SHA-256).
  2. Input must satisfy min(shape) >= 2 (33 traces; single-trace gate now
     passes; recorded as single_trace_gate_blocked=false).
  3. BG-order gain configs persist the true pre-gain result; GB-order only
     persists G^-1 Y, labeled as a non-true intermediate. Negative-control
     rows are reported separately and never merged into target rows.

Metrics: with reference_state=numerically_unresolved, waveform metrics D/A/H
and rho are NOT computible (recorded as unavailable with the verbatim gate
reason); R_c is never computed (no legitimate pure-clutter window). The only
v0.2 §4 metric computed is the negative-control residual N_b (energy_metrics);
gain risk diagnostics and resource usage are recorded. Feasibility labels use
configuration_labels with all-null limits -> undetermined/partial_only.

Heavy arrays (pre-gain results and GB G^-1 Y audits) persist as .npy under
artifacts/simulations/2026-09-26_eval_batch2d_mt_arrays/ (git-ignored); records
and manifest under artifacts/research_checks/<run_dir>/ carry path + SHA-256.

Read-only inputs: archived MT h5 files, frozen event table, operator contract.
The runner never writes into input directories and never invokes a solver.
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
MT_ANALYSIS = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_mt_analysis/results_r1.json'
MT_ARCHIVE = 'artifacts/research_checks/2026-09-26_{rid}/'  # evidence dirs hold the h5 copies
ARRAYS_DIR = ROOT / 'artifacts/simulations/2026-09-26_eval_batch2d_mt_arrays'
N_TRACES = 33
GATHER_SEMANTICS_VERSION = 'variable_offset_csg_not_common_offset_bscan_v1'
DIRECT_WAVE_DIFF_SAMPLES = 147  # proposal §2.4 self-computed (mt01 vs mt33)

HARD_LIMITS = {
    "absolute_accuracy_claimed": False,
    "clean_truth_generated": False,
    "conclusions_scope": "operator-evaluation evidence for the batch2d_v1_mt BASE multitrace scenario families only; not field performance, not 3D, not absolute accuracy, no maximum-depth claim",
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
    "n_traces": 33,
    "trace_axis_semantics": "variable-offset single-shot multi-receiver gather (Tx fixed y15.35, offset 2.70->1.30->5.30 m along trace index); NOT a common-offset B-scan",
    "anchor_trace": "mt17 (y16.65 m, grid y666) bit-for-bit equal to archived single-trace mother (22/22)",
    "nc_multitrace_control": "NC-BG paired difference identically zero on all 33 traces, all four families (4/4)",
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_gather(rid):
    """[n_samples, 33] float64 gather from an archived MT evidence h5."""
    h5_path = ROOT / MT_ARCHIVE.format(rid=rid) / f'{rid}.h5'
    with h5py.File(h5_path, 'r') as f:
        assert int(f.attrs['nrx']) == N_TRACES
        x = np.stack([f[f'rxs/rx{k}/Ex'][:] for k in range(1, N_TRACES + 1)], axis=1)
    return x, str(h5_path.relative_to(ROOT)), sha256_file(h5_path)


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

    mt = json.loads(MT_ANALYSIS.read_text(encoding='utf-8'))
    assert mt['all_cases_pass'] and mt['n_cases'] == 22
    mt_mothers = {c['mother']: c['run_id'] for c in mt['cases']}

    configs = catalogue()
    assert len(configs) == 27, 'catalogue must hold the 27 frozen candidates'
    op_script_sha = sha256_file(ROOT / 'scripts/research_operator_contract.py')
    eval_script_sha = sha256_file(ROOT / 'scripts/run_eval_batch2d_mt.py')

    # --- expand event rows to (event, MT case) evaluation rows --------------
    rows = []  # (event, case_run_id or None)
    for e in table['entries']:
        if e.get('applies_to') == 'all_cases_of_family':
            fam_prefix = 'B2D-' + e['family'].capitalize() + 'm'  # c1 -> B2D-C1m
            for mother, rid in sorted(mt_mothers.items()):
                if mother.startswith(fam_prefix + '-'):
                    rows.append((e, rid))
        else:
            rows.append((e, mt_mothers.get(e['mother_model_id'])))

    gathers = {}   # rid -> (x, rel_path, sha)
    def get_gather(rid):
        if rid not in gathers:
            gathers[rid] = load_gather(rid)
        return gathers[rid]

    records = []
    availability = {}
    reason_hist = {}
    tracemalloc.start()
    mem_base = tracemalloc.get_traced_memory()[1]
    for e, rid in rows:
        evaluable = rid is not None
        x_window = None
        case_file = case_sha = None
        nc_diff_max = None
        if evaluable:
            x, case_file, case_sha = get_gather(rid)
            lo, hi = int(e['index_lo']), int(e['index_hi'])
            assert 0 <= lo <= hi and hi + 1 <= x.shape[0], f"window {e['event_id']} outside gather"
            x_window = x[lo:hi + 1, :]
            assert np.isfinite(x_window).all() and x_window.shape[1] == N_TRACES
            if e['role'] == 'nc_zero':
                bg_rid = mt_mothers[f'B2D-{rid.split("-")[1]}-BG']
                xb, _, _ = get_gather(bg_rid)
                nc_diff_max = float(np.max(np.abs(x_window - xb[lo:hi + 1, :])))
        for cfg in configs:
            rec = {
                'record_id': f"{e['event_id']}::{rid or 'NO_MT_CASE'}::{cfg['id']}",
                'event_id': e['event_id'], 'event_role': e['role'],
                'case_run_id': rid, 'family': e['family'],
                'candidate_id': cfg['id'], 'background': cfg['background'],
                'parameter': cfg['parameter'], 'end_gain': cfg['end_gain'],
                'requested_order': cfg['order'],
                'row_type': 'negative_control' if e['role'] in
                            ('nc_zero', 'bg_absent', 'off_path') else 'event',
                'window': {'index_lo': lo if evaluable else None,
                           'index_hi': hi if evaluable else None,
                           'n_window_samples': e['n_window_samples'],
                           'window_hash': e['window_hash'], 'mask_hash': e['mask_hash'],
                           'time_axis_convention_version': e['time_axis_convention_version']},
                'gather_semantics': {
                    'version': GATHER_SEMANTICS_VERSION,
                    'n_traces': N_TRACES if evaluable else None,
                    'trace_axis_is_spatial_profile': False,
                    'window_applied_unshifted_to_all_traces': True,
                    'direct_wave_differential_samples_end_to_end': DIRECT_WAVE_DIFF_SAMPLES,
                    'trace_offset_range_m': [2.70, 5.30],
                    'anchor_trace': 'mt17',
                },
                'provenance': {
                    'case_h5': case_file, 'case_h5_sha256': case_sha,
                    'event_table_sha256': table_sha,
                    'operator_script_sha256': op_script_sha,
                    'eval_script_sha256': eval_script_sha,
                    'operator_version': OP_VERSION,
                    'dt_s': dt_s,
                },
            }
            if not evaluable:
                rec.update(availability='unavailable', failure_reason='no_multitrace_case_in_batch2d_v1_mt',
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] = availability.get('unavailable', 0) + 1
                reason_hist['no_multitrace_case_in_batch2d_v1_mt'] = reason_hist.get('no_multitrace_case_in_batch2d_v1_mt', 0) + 1
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
                path, fsha = persist_array(f"{e['event_id']}__{rid}__{cfg['id']}", pre)
                curve = np.asarray(res['gain_curve'])
                cpath, csha = persist_array(f"{e['event_id']}__{rid}__{cfg['id']}__gaincurve", curve)
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
                nc_minus_bg_window_max_abs=nc_diff_max if e['role'] == 'nc_zero' else None,
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
        'runner': 'scripts/run_eval_batch2d_mt.py',
        'event_table': str(EVENT_TABLE.relative_to(ROOT)),
        'event_table_sha256': table_sha,
        'time_axis_convention_version': table['time_axis']['convention_version'],
        'grid_tier': 'BASE',
        'gather_semantics_version': GATHER_SEMANTICS_VERSION,
        'trace_axis_semantics': HARD_LIMITS['trace_axis_semantics'],
        'window_mapping': 'frozen sample-axis window applied unshifted to all 33 traces; no per-trace alignment',
        'candidates': {'n': len(configs), 'version': OP_VERSION,
                       'operator_script_sha256': op_script_sha,
                       'eval_script_sha256': eval_script_sha},
        'event_rows': len(rows),
        'n_records': len(records),
        'availability': availability,
        'failure_reason_histogram': reason_hist,
        'gates': {'event_table_frozen': True, 'min_shape_ge_2': True,
                  'single_trace_gate_blocked': False,
                  'blocked_candidate_count': 0,
                  'pre_gain_persisted': True,
                  'negative_controls_separate': True},
        'runner_limits': {'r_c_computed': False, 'isolated_event_granted': False,
                          'cross_family_metrics': False, 'gain_reference_absent': True,
                          'waveform_metrics_computed': False,
                          'waveform_metrics_reason': 'reference_state numerically_unresolved + paired_contrast (design A.4)'},
        'feasibility': {'status': labels['status'],
                        'note': 'all limits null (G4 unresolved): every gate undetermined; no ranking, no unique label'},
        'negative_control_events': sorted({r['event_id'] for r in records if r['row_type'] == 'negative_control'}),
        'hard_limits': HARD_LIMITS,
        'cases': {rid: {'h5': g[1], 'h5_sha256': g[2]} for rid, g in sorted(gathers.items())},
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
