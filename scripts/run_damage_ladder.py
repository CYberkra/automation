"""G4 step-2 damage ladder runner: apply closed-form damages to archived MT/CO
gather windows (dev groups {B2D-C1m, B2D-C3m} per the frozen S4 split) and score
the 27-entry catalogue + frozen S4 L/R baselines against identity(damaged input)
as the constructed reference.

Epistemics (fixed): the ladder reference is the damaged input itself, exact by
construction. Absolute metrics (D/A/H) are legal ONLY within this ladder via
reference_kind='complete_clean', state='mixed', scope='complete' on the identity
output of the damaged input. This is constructed-truth evidence about operator
behaviour under known perturbations — NOT physical clean truth, NOT field
detectability. Hard limits keep reference_state=numerically_unresolved,
physical thresholds null, training off; no ranking or labels are produced.

Damage types (levels are pre-freeze design starting points):
  * amplitude_scale {0.9, 0.5, 0.1}  on a predefined local rectangle
  * polarity_flip                    on the same rectangle
  * sample_shift {+1,-1,+4,-4,+16,-16} samples, zero-filled, inside the rectangle
  * trace_deletion {1, 4, 8} traces of the predefined central subset
  * nc_amplify {2, 4}                negative-control rows only (gain risk E4)
Rectangle: central half of window samples x central min(8, n_traces) traces.

Metrics: N_b (energy vs damaged input) + waveform D/A/H (constructed reference)
+ svd/gain diagnostics (scalar only; gain arrays are NOT persisted — ladder
outputs are deterministically reproducible from archived inputs, and array
volume would be excessive; recorded in runner_limits).

Deterministic; supports --chunk/--n-chunks for session tool time limits.
Read-only inputs; never invokes a solver.

Split mode (--split {dev,test}, default 'dev'):
  dev  -> dev families {c1, c3}: the frozen G4 S1-S3 behaviour, unchanged by
          this runner revision (same families, same work-item order, same record
          content, same manifest fields and values). Note: every record carries
          provenance.runner_script_sha256 = SHA-256 of THIS file, which is not
          stripped by the merge, so a dev re-run after any edit to this file
          changes that one field (and only that field) in the merged bytes.
  test -> test families {c5, c8}: G4 step 5 / S4 test-family ladder. This is
          data generation ONLY -- damage instances, rectangle, metrics, chunk
          protocol and determinism are identical; only the family set and the
          split-declaring manifest fields change.
No threshold is applied in either mode; the frozen dev-derived threshold is
applied once, later, in S5.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import tracemalloc

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from research_operator_contract import (apply_configuration, catalogue, local_mean_control,
                                        VERSION as OP_VERSION)
from research_evaluation_contract import configuration_labels, energy_metrics, waveform_metrics
from run_eval_batch2d_mt import (EVENT_TABLE, EVENT_TABLE_SHA256, MT_ANALYSIS,
                                 GATHER_SEMANTICS_VERSION as MT_SEMANTICS,
                                 load_gather as load_mt_gather, N_TRACES as MT_N_TRACES)
from run_eval_batch2d_co import (GATHER_SEMANTICS_VERSION as CO_SEMANTICS,
                                 CO_CASES, load_case, N_TRACES as CO_N_TRACES)

ROOT = Path(__file__).resolve().parents[1]
CO_ANALYSIS = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_co_analysis/results.json'
DEV_FAMILIES = ('c1', 'c3')
TEST_FAMILIES = ('c5', 'c8')
DEV_GROUPS = ['B2D-C1m', 'B2D-C3m']

AMP_LEVELS = [0.9, 0.5, 0.1]
SHIFT_LEVELS = [1, -1, 4, -4, 16, -16]
DELETE_LEVELS = [1, 4, 8]
NC_GAIN_LEVELS = [2, 4]
L_LAMBDAS = [0.25, 0.5, 1.0]
L_WIDTHS = [5, 11, 21]
L_GAINS = [1, 2, 4]
R_THRESHOLDS = [0.05, 0.1, 0.2]

MT_SEMANTICS_DICT = {
    'version': MT_SEMANTICS,
    'trace_axis_is_spatial_profile': False,
    'window_applied_unshifted_to_all_traces': True,
    'trace_offset_range_m': [2.70, 5.30],
    'anchor_trace': 'mt17',
}
CO_SEMANTICS_DICT = {
    'version': CO_SEMANTICS,
    'trace_axis_is_spatial_profile': True,
    'window_applied_unshifted_to_all_traces': True,
    'trace_offset_range_m': [1.30, 1.30],
    'trace_spacing_m': 1.0,
    'rx_y_range_m': [12.65, 22.65],
    'anchor_trace': 't05',
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixed_gain_curve(n_samples, end_gain):
    return np.power(float(end_gain), np.arange(n_samples) / (n_samples - 1))[:, None]


def damage_rectangle(shape):
    """Central half of samples x central min(8, n_traces) traces."""
    n_s, n_t = shape
    s0, s1 = n_s // 4, (3 * n_s) // 4
    w_t = min(8, n_t)
    t0 = (n_t - w_t) // 2
    return s0, s1, t0, w_t


def apply_damage(x, damage_type, level):
    """Return (damaged, damage_provenance). Exact by construction."""
    s0, s1, t0, w_t = damage_rectangle(x.shape)
    if damage_type == 'amplitude_scale':
        y = x.copy()
        y[s0:s1, t0:t0 + w_t] *= level
        return y, {'region': [s0, s1, t0, w_t]}
    if damage_type == 'polarity_flip':
        y = x.copy()
        y[s0:s1, t0:t0 + w_t] *= -1.0
        return y, {'region': [s0, s1, t0, w_t]}
    if damage_type == 'sample_shift':
        y = x.copy()
        band = x[s0:s1, t0:t0 + w_t].copy()
        shifted = np.zeros_like(band)
        d = int(level)
        if d > 0:
            shifted[d:] = band[:-d] if d < band.shape[0] else 0.0
        elif d < 0:
            shifted[:d] = band[-d:] if -d < band.shape[0] else 0.0
        y[s0:s1, t0:t0 + w_t] = shifted
        return y, {'region': [s0, s1, t0, w_t], 'zero_filled': True}
    if damage_type == 'trace_deletion':
        cols = list(range(t0, min(t0 + int(level), t0 + w_t)))
        y = np.delete(x, cols, axis=1)
        return y, {'deleted_columns': cols}
    if damage_type == 'nc_amplify':
        y = fixed_gain_curve(x.shape[0], level) * x
        return y, {'gain': level, 'curve': 'fixed_shared_gain_end'}
    raise ValueError(f'unknown damage type {damage_type}')


def damage_instances(row_type):
    insts = ([('amplitude_scale', v) for v in AMP_LEVELS]
             + [('polarity_flip', 1.0)]
             + [('sample_shift', v) for v in SHIFT_LEVELS]
             + [('trace_deletion', v) for v in DELETE_LEVELS])
    if row_type == 'negative_control':
        insts += [('nc_amplify', v) for v in NC_GAIN_LEVELS]
    return insts


def expand_dev_rows(table, mt_mothers, families):
    rows = []
    for e in table['entries']:
        if e['family'] not in families:
            continue
        cases = []
        if e.get('applies_to') == 'all_cases_of_family':
            fam_prefix = 'B2D-' + e['family'].capitalize() + 'm'
            cases += [('mt', mt_mothers[m]) for m in sorted(mt_mothers)
                      if m.startswith(fam_prefix + '-')]
            cases += [('co', CO_CASES[m]) for m in sorted(CO_CASES)
                      if m.startswith(fam_prefix + '-')]
        else:
            if e['mother_model_id'] in mt_mothers:
                cases.append(('mt', mt_mothers[e['mother_model_id']]))
            if e['mother_model_id'] in CO_CASES:
                cases.append(('co', CO_CASES[e['mother_model_id']]))
        for geom, case in cases:
            rows.append((e, geom, case))
    return rows


def get_window(geom, case, mt_gathers, co_cases):
    if geom == 'mt':
        if case not in mt_gathers:
            mt_gathers[case] = load_mt_gather(case)
        return mt_gathers[case][0], MT_N_TRACES, MT_SEMANTICS_DICT
    return co_cases[case][0], CO_N_TRACES, CO_SEMANTICS_DICT


def run_ladder(run_tag, chunk, n_chunks, split='dev'):
    if split not in ('dev', 'test'):
        raise ValueError(f'unknown split {split}')
    families = DEV_FAMILIES if split == 'dev' else TEST_FAMILIES
    started = time.perf_counter()
    table_sha = sha256_file(EVENT_TABLE)
    assert table_sha == EVENT_TABLE_SHA256, 'event table hash mismatch'
    table = json.loads(EVENT_TABLE.read_text(encoding='utf-8'))
    assert table['status'] == 'frozen'
    dt_s = table['time_axis']['dt_s']

    mt = json.loads(MT_ANALYSIS.read_text(encoding='utf-8'))
    assert mt['all_cases_pass'] and mt['n_cases'] == 22
    mt_mothers = {c['mother']: c['run_id'] for c in mt['cases']}
    co = json.loads(CO_ANALYSIS.read_text(encoding='utf-8'))
    assert co['batch'] == 'batch2d_v1_co' and co['all_layout_ok'] is True
    assert all(co['anchor_t05_bit_for_bit'].values())

    configs = catalogue()
    assert len(configs) == 27
    op_script_sha = sha256_file(ROOT / 'scripts/research_operator_contract.py')
    eval_contract_sha = sha256_file(ROOT / 'scripts/research_evaluation_contract.py')
    runner_sha = sha256_file(ROOT / 'scripts/run_damage_ladder.py')

    rows = expand_dev_rows(table, mt_mothers, families)
    co_cases = {ck: load_case(ck) for ck in ('BG', 'TGT')}
    mt_gathers = {}

    # deterministic work-item order: row-major, then damage instances
    work = []
    for e, geom, case in rows:
        row_type = 'negative_control' if e['role'] in ('nc_zero', 'bg_absent', 'off_path') else 'event'
        for dtype, level in damage_instances(row_type):
            work.append((e, geom, case, row_type, dtype, level))
    selected = [w for i, w in enumerate(work) if i % n_chunks == chunk]

    records = []
    availability = {'ran': 0, 'unavailable': 0}
    reason_hist = {}
    tracemalloc.start()
    mem_base = tracemalloc.get_traced_memory()[1]
    for e, geom, case, row_type, dtype, level in selected:
        x, n_traces, sem = get_window(geom, case, mt_gathers, co_cases)
        lo, hi = int(e['index_lo']), int(e['index_hi'])
        assert hi + 1 <= x.shape[0]
        xw = x[lo:hi + 1, :]
        damaged, dprov = apply_damage(xw, dtype, level)
        assert np.isfinite(damaged).all()
        mask = np.ones(damaged.shape, dtype=bool)
        ref = damaged  # constructed reference: identity(damaged input)
        base = {
            'record_prefix': f"{e['event_id']}::{geom}:{case}::{dtype}:{level}",
            'event_id': e['event_id'], 'event_role': e['role'],
            'family': e['family'], 'group_id': e['group_id'],
            'geometry': geom, 'case': case, 'row_type': row_type,
            'damage': {'type': dtype, 'level': level, **dprov,
                       'reference': 'identity_of_damaged_input_constructed_not_physical_truth'},
            'window': {'index_lo': lo, 'index_hi': hi,
                       'n_window_samples': e['n_window_samples'],
                       'window_hash': e['window_hash'], 'mask_hash': e['mask_hash'],
                       'time_axis_convention_version': e['time_axis_convention_version']},
            'gather_semantics': dict(sem, n_traces=int(damaged.shape[1])),
            'provenance': {
                'event_table_sha256': table_sha,
                'operator_script_sha256': op_script_sha,
                'evaluation_contract_sha256': eval_contract_sha,
                'runner_script_sha256': runner_sha,
                'operator_version': OP_VERSION,
                'dt_s': dt_s,
            },
        }

        def score(y, t0):
            wall = time.perf_counter() - t0
            peak = tracemalloc.get_traced_memory()[1] - mem_base
            nb = energy_metrics(y, ref, mask)
            wf = waveform_metrics(y, ref, mask, reference_kind='complete_clean',
                                  state='mixed', scope='complete')
            return {
                'availability': 'ran', 'failure_reason': None,
                'metrics': {
                    'N_b': nb,
                    'waveform': {'available': wf['available'],
                                 'reason': wf['reason'],
                                 'absolute_preservation_eligible': wf['absolute_preservation_eligible'],
                                 'reference_kind': 'complete_clean', 'reference_state': 'mixed',
                                 'reference_provenance': 'identity_of_damaged_input_constructed',
                                 'metrics': None if wf['metrics'] is None else dict(wf['metrics'])},
                },
                'resource': {'wall_s': wall, 'peak_tracemalloc_bytes': int(peak),
                             'compute_dtype': 'float64',
                             'n_samples': int(y.shape[0]), 'n_traces': int(y.shape[1])},
            }

        # --- 27-entry catalogue ------------------------------------------
        for cfg in configs:
            rec = dict(base, record_id=f"{base['record_prefix']}::{cfg['id']}",
                       baseline='catalogue', candidate_id=cfg['id'],
                       background=cfg['background'], parameter=cfg['parameter'],
                       end_gain=cfg['end_gain'], requested_order=cfg['order'])
            if min(damaged.shape) < 2:
                rec.update(availability='unavailable',
                           failure_reason='min_shape_lt_2_after_damage',
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] += 1
                reason_hist['min_shape_lt_2_after_damage'] = reason_hist.get('min_shape_lt_2_after_damage', 0) + 1
                records.append(rec)
                continue
            t0 = time.perf_counter()
            try:
                res = apply_configuration(damaged, cfg['id'])
            except Exception as exc:
                rec.update(availability='unavailable', failure_reason=str(exc),
                           metrics=None, gain_diagnostics=None, resource=None)
                availability['unavailable'] += 1
                reason_hist[str(exc)] = reason_hist.get(str(exc), 0) + 1
                records.append(rec)
                continue
            y = res['output']
            rec.update(score(y, t0),
                       executed_order=''.join(s['operator'] for s in res['steps']),
                       svd_diagnostics=next((s['diagnostics'] for s in res['steps']
                                             if s['diagnostics'].get('kind') == 'svd'), None))
            if cfg['end_gain'] > 1:
                pre = res['pre_gain'] if cfg['order'] == 'BG' else res['audit_before_gain']
                rec['gain_diagnostics'] = {
                    'persisted_array': None,
                    'persisted_kind': 'pre_gain_true_intermediate' if cfg['order'] == 'BG'
                                      else 'ginv_audit_not_true_intermediate',
                    'arrays_persisted': False,
                    'max_gain': cfg['end_gain'],
                    'max_abs_pre_gain': float(np.max(np.abs(pre))),
                    'max_abs_output': float(np.max(np.abs(y))),
                    'gain_energy_ratio_output_over_pre': float(np.mean(y * y) / np.mean(pre * pre)),
                    'clipped_samples': int(np.count_nonzero(
                        np.abs(y) > cfg['end_gain'] * np.max(np.abs(pre)) + 1e-300)),
                    'overflow': bool(not np.isfinite(y).all()),
                }
            else:
                rec['gain_diagnostics'] = None
            availability['ran'] += 1
            records.append(rec)

        # --- L baseline grid ----------------------------------------------
        for lam in L_LAMBDAS:
            for w in L_WIDTHS:
                for q in L_GAINS:
                    cid = f"L_lam{lam}_w{w}_q{q}"
                    rec = dict(base, record_id=f"{base['record_prefix']}::{cid}",
                               baseline='L', candidate_id=cid,
                               strength_lambda=lam, window_width=w, end_gain=q,
                               requested_order='LM_then_G')
                    if w > damaged.shape[1] or min(damaged.shape) < 2:
                        rec.update(availability='unavailable',
                                   failure_reason='width_exceeds_n_traces_after_damage' if w > damaged.shape[1]
                                                  else 'min_shape_lt_2_after_damage',
                                   metrics=None, gain_diagnostics=None, resource=None)
                        availability['unavailable'] += 1
                        reason_hist[rec['failure_reason']] = reason_hist.get(rec['failure_reason'], 0) + 1
                        records.append(rec)
                        continue
                    t0 = time.perf_counter()
                    pre = local_mean_control(damaged, w, strength=lam)
                    y = pre if q == 1 else fixed_gain_curve(pre.shape[0], q) * pre
                    rec.update(score(y, t0), executed_order='LM' + ('G' if q > 1 else ''))
                    if q > 1:
                        rec['gain_diagnostics'] = {
                            'persisted_array': None, 'persisted_kind': 'pre_gain_true_intermediate',
                            'arrays_persisted': False, 'max_gain': q,
                            'max_abs_pre_gain': float(np.max(np.abs(pre))),
                            'max_abs_output': float(np.max(np.abs(y))),
                            'gain_energy_ratio_output_over_pre': float(np.mean(y * y) / np.mean(pre * pre)),
                            'clipped_samples': int(np.count_nonzero(
                                np.abs(y) > q * np.max(np.abs(pre)) + 1e-300)),
                            'overflow': bool(not np.isfinite(y).all()),
                        }
                    else:
                        rec['gain_diagnostics'] = None
                    availability['ran'] += 1
                    records.append(rec)

        # --- R baseline rule ------------------------------------------------
        if min(damaged.shape) >= 2:
            scale = float(np.max(np.abs(damaged)))
            svd_vals = np.linalg.svd(damaged / scale if scale > 0 else damaged,
                                     compute_uv=False)
            relative = svd_vals / svd_vals[0]
            rank_floor = max(damaged.shape) * np.finfo(np.float64).eps
            numerical_rank = int(np.count_nonzero(relative > rank_floor))
            ks = [k for k in (1, 2, 3) if k < relative.shape[0]]
            gaps = [float(relative[k - 1] - relative[k]) for k in ks]
            k_star = int(np.argmax(gaps)) + 1
            rejected = bool((k_star < numerical_rank and gaps[k_star - 1] <= 1e-8)
                            or relative[k_star - 1] <= rank_floor)
        else:
            gaps, k_star, rejected, numerical_rank, rank_floor = [], None, True, 0, None
        for tau in R_THRESHOLDS:
            cid = f"R_tau{tau}_q1"
            rec = dict(base, record_id=f"{base['record_prefix']}::{cid}",
                       baseline='R', candidate_id=cid, threshold_tau=tau, q=1)
            if min(damaged.shape) < 2:
                rec.update(availability='unavailable', failure_reason='min_shape_lt_2_after_damage',
                           metrics=None, resource=None)
                availability['unavailable'] += 1
                records.append(rec)
                continue
            select_k = None if (not gaps or max(gaps) < tau or rejected) else k_star
            t0 = time.perf_counter()
            if select_k is None:
                y = damaged.copy()
                executed = 'identity'
            else:
                res = apply_configuration(damaged, f'B{3 + select_k}_G1_BG')
                y = res['output']
                executed = f'svd_k{select_k}_BG'
            rec.update(score(y, t0), executed_order=executed,
                       gap_spectrum_d1_d2_d3=gaps, k_star=k_star, max_gap=max(gaps),
                       numerical_rank=numerical_rank, numerical_rank_floor=rank_floor,
                       selection_rejected_numerically=rejected,
                       selection='identity' if select_k is None else f'svd_k{select_k}')
            availability['ran'] += 1
            records.append(rec)
    tracemalloc.stop()

    cands = [{'id': r['record_id'], 'available': r['availability'] == 'ran',
              'metrics': {'nrmse': None, 'amplitude_error': None,
                          'negative_control_residual': None},
              'failure_reason': r['failure_reason']} for r in records]
    labels = configuration_labels(cands, limits={'nrmse': None, 'amplitude_error': None,
                                                   'negative_control_residual': None},
                                  objectives=['nrmse', 'amplitude_error'])

    manifest = {
        'run_tag': run_tag,
        'runner': 'scripts/run_damage_ladder.py',
        'chunk': chunk,
        'n_chunks': n_chunks,
        'event_table': str(EVENT_TABLE.relative_to(ROOT)),
        'event_table_sha256': table_sha,
        'time_axis_convention_version': table['time_axis']['convention_version'],
        'grid_tier': 'BASE',
        'dev_families': list(DEV_FAMILIES),
        'n_rows': len(rows),
        'n_work_items': len(work),
        'n_selected_work_items': len(selected),
        'n_records': len(records),
        'availability': availability,
        'failure_reason_histogram': reason_hist,
        'damage_types': {
            'amplitude_scale': AMP_LEVELS,
            'polarity_flip': [1.0],
            'sample_shift': SHIFT_LEVELS,
            'trace_deletion': DELETE_LEVELS,
            'nc_amplify': {'levels': NC_GAIN_LEVELS, 'rows': 'negative_control only'},
        },
        'damage_rectangle': 'central half of window samples x central min(8, n_traces) traces',
        'operator_set': '27-entry catalogue + S4 L grid (27) + R rule (3 tau); identity = B0_G1_BG anchor',
        'gates': {'event_table_frozen': True, 'mt_analysis_pass': True,
                  'co_analysis_pass': True, 'dev_groups_only': True},
        'runner_limits': {
            'r_c_computed': False, 'isolated_event_granted': False,
            'cross_family_metrics': False, 'gain_reference_absent': True,
            'gain_arrays_persisted': False,
            'gain_arrays_reason': 'ladder volumes; outputs deterministically reproducible from archived inputs + frozen script; scalar diagnostics only',
            'waveform_metrics_computed_within_ladder': True,
            'waveform_reference': 'complete_clean/mixed/complete on identity(damaged input); constructed only',
        },
        'feasibility': {'status': labels['status'],
                        'note': 'all limits null (G4 unresolved): every gate undetermined; no ranking, no unique label'},
        'hard_limits': {
            'absolute_accuracy_claimed': False,
            'clean_truth_generated': False,
            'conclusions_scope': 'operator behaviour under closed-form damages on archived dev-group windows only; not field performance, not detectability, not physical thresholds',
            'constructed_reference_is_physical_truth': False,
            'cross_family_differencing': False,
            'cross_tier_differencing': False,
            'fdtd_solver_invoked': False,
            'grid_convergence_certified': False,
            'grid_tier': 'BASE',
            'physical_acceptance_threshold': None,
            'physical_label_eligible': False,
            'reference_state': 'numerically_unresolved',
            'training_eligible': False,
            'training_labels_generated': False,
            'negative_controls_merged_into_targets': False,
            'test_groups_excluded': ['B2D-C5m', 'B2D-C8m'],
        },
        'total_wall_s': time.perf_counter() - started,
        'stop_reason': None,
    }
    if split == 'test':
        manifest['split'] = 'test'
        manifest['families_used'] = list(families)
        # key kept (empty) so existing readers that assume 'dev_families' exists
        # still find it; the split's real family set is 'families_used'.
        manifest['dev_families'] = []
        manifest['dev_families_note'] = ('test split: no dev-family rows expanded; '
                                         'dev_families retained empty for schema compatibility')
        manifest['gates']['dev_groups_only'] = False
        manifest['gates']['test_groups_only'] = True
        hl = {}
        for k, v in manifest['hard_limits'].items():
            if k == 'test_groups_excluded':
                hl['dev_groups_excluded'] = list(DEV_GROUPS)
                continue
            if k == 'conclusions_scope':
                hl[k] = v.replace('dev-group', 'test-group')
                continue
            hl[k] = v
        manifest['hard_limits'] = hl
    return records, manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True,
                   help='must not exist (no overwriting of evidence)')
    p.add_argument('--run-tag', default='r1')
    p.add_argument('--chunk', type=int, default=0)
    p.add_argument('--n-chunks', type=int, default=1)
    p.add_argument('--split', choices=['dev', 'test'], default='dev',
                   help="dev={c1,c3} (frozen G4 S1-S3 behaviour, unchanged); "
                        "test={c5,c8} (G4 step 5 / S4 test-family data generation)")
    a = p.parse_args()
    out = (ROOT / a.output_dir).resolve() if not a.output_dir.is_absolute() else a.output_dir
    assert not out.exists(), f'refusing to overwrite {out}'
    assert 0 <= a.chunk < a.n_chunks
    records, manifest = run_ladder(a.run_tag, a.chunk, a.n_chunks, a.split)
    out.mkdir(parents=True)
    (out / 'records.json').write_text(json.dumps({'records': records}, indent=1) + '\n',
                                      encoding='utf-8')
    (out / 'run_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n',
                                           encoding='utf-8')
    print(json.dumps({'records': len(records), 'availability': manifest['availability'],
                      'chunk': a.chunk, 'n_chunks': a.n_chunks,
                      'feasibility': manifest['feasibility']['status'], 'output': str(out)},
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
