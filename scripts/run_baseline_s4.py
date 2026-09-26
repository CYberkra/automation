"""S4 baseline group runner: evaluate the frozen L/R baselines on the same frozen
event windows against the same metrics contract as the main eval runners (MT + CO
geometries, layered separately, never merged or cross-paired), with F as a
no-computation dev-group view over the existing 27-candidate eval records and O
declared unavailable (no clean-truth reference, no approximation).

Implements configs/research/s4_baseline_group_v0.1.json:
  * L: local_mean_control external grid lambda x width x gain, BG-like order
    (local mean then fixed shared gain, no GB variant); width > n_traces is
    unavailable for that window (no auto-shortening).
  * R: input-spectrum rule d_k=(sigma_k-sigma_{k+1})/sigma_1, k in {1,2,3};
    argmax d_k if >= tau else identity; three threshold candidates run side by
    side as mechanism diagnostics; q=1 ablation only. Threshold choice itself is
    undetermined (G4) and is NOT made here.
  * F: selection protocol exists but is not executed (G4); dev-group projection
    of existing eval records only, no new computation.
  * O: unavailable_no_computable_reference per row, no approximation.

Discipline: reference_state=numerically_unresolved; D/A/H/rho not computible;
R_c never computed; feasibility partial_only; no ranking, no labels, no training.
Read-only inputs; never invokes a solver. Heavy arrays persist to a git-ignored
dir with SHA-256 in records.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import tracemalloc

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from research_operator_contract import apply_configuration, local_mean_control, VERSION as OP_VERSION
from research_evaluation_contract import configuration_labels, energy_metrics
from run_eval_batch2d_mt import (EVENT_TABLE, EVENT_TABLE_SHA256, MT_ANALYSIS,
                                 GATHER_SEMANTICS_VERSION as MT_SEMANTICS,
                                 DIRECT_WAVE_DIFF_SAMPLES, load_gather as load_mt_gather,
                                 N_TRACES as MT_N_TRACES)
from run_eval_batch2d_co import (GATHER_SEMANTICS_VERSION as CO_SEMANTICS,
                                 CO_CASES, load_case, N_TRACES as CO_N_TRACES)

ROOT = Path(__file__).resolve().parents[1]
S4_CONTRACT = ROOT / 'configs/research/s4_baseline_group_v0.1.json'
S4_CONTRACT_SHA256 = '8676f37d6aee603b7f2481779ffe361d2878735316863da0bed0db82e16b5593'
CO_ANALYSIS = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_co_analysis/results.json'
ARRAYS_DIR = ROOT / 'artifacts/simulations/2026-09-26_s4_baseline_arrays'
MT_EVAL_R1 = 'artifacts/research_checks/2026-09-26_eval_batch2d_mt_r1'
CO_EVAL_R1 = 'artifacts/research_checks/2026-09-26_eval_batch2d_co_r1'

L_LAMBDAS = [0.25, 0.5, 1.0]
L_WIDTHS = [5, 11, 21]
L_GAINS = [1, 2, 4]
R_THRESHOLDS = [0.05, 0.1, 0.2]
R_KS = [1, 2, 3]
DEV_GROUPS = ['B2D-C1m', 'B2D-C3m']
TEST_GROUPS = ['B2D-C5m', 'B2D-C8m']

MT_SEMANTICS_DICT = {
    'version': MT_SEMANTICS,
    'trace_axis_is_spatial_profile': False,
    'window_applied_unshifted_to_all_traces': True,
    'direct_wave_differential_samples_end_to_end': DIRECT_WAVE_DIFF_SAMPLES,
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


def persist_array(tag, array):
    ARRAYS_DIR.mkdir(parents=True, exist_ok=True)
    path = ARRAYS_DIR / f'{tag}.npy'
    np.save(path, array)
    return str(path.relative_to(ROOT)), sha256_file(path)


def expand_rows(table, mt_mothers):
    """(event, geometry, case_key or rid or None) rows, mirroring both eval runners."""
    rows = []
    for e in table['entries']:
        if e.get('applies_to') == 'all_cases_of_family':
            fam_prefix = 'B2D-' + e['family'].capitalize() + 'm'
            mt_matches = sorted(m for m in mt_mothers if m.startswith(fam_prefix + '-'))
            co_matches = sorted(m for m in CO_CASES if m.startswith(fam_prefix + '-'))
            if mt_matches:
                for m in mt_matches:
                    rows.append((e, 'mt', mt_mothers[m]))
            else:
                rows.append((e, 'mt', None))
            if co_matches:
                for m in co_matches:
                    rows.append((e, 'co', CO_CASES[m]))
            else:
                rows.append((e, 'co', None))
        else:
            rows.append((e, 'mt', mt_mothers.get(e['mother_model_id'])))
            rows.append((e, 'co', CO_CASES.get(e['mother_model_id'])))
    return rows


def run_baselines(run_tag):
    started = time.perf_counter()
    # --- gates: frozen S4 contract, frozen event table, analysis gates ---------
    contract_sha = sha256_file(S4_CONTRACT)
    assert contract_sha == S4_CONTRACT_SHA256, 'S4 contract hash mismatch'
    contract = json.loads(S4_CONTRACT.read_text(encoding='utf-8'))
    assert contract['status'] == 'frozen'
    for rel, fsha in contract['inputs'].items():
        assert sha256_file(ROOT / rel) == fsha, f'input hash mismatch: {rel}'

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

    op_script_sha = sha256_file(ROOT / 'scripts/research_operator_contract.py')
    eval_contract_sha = sha256_file(ROOT / 'scripts/research_evaluation_contract.py')
    runner_sha = sha256_file(ROOT / 'scripts/run_baseline_s4.py')
    freeze_script_sha = sha256_file(ROOT / 'scripts/freeze_s4_baseline_group.py')

    rows = expand_rows(table, mt_mothers)
    co_cases = {ck: load_case(ck) for ck in ('BG', 'TGT')}
    mt_gathers = {}

    def get_window(geom, case):
        if geom == 'mt':
            if case not in mt_gathers:
                mt_gathers[case] = load_mt_gather(case)
            x = mt_gathers[case][0]
            n_traces = MT_N_TRACES
            sem = MT_SEMANTICS_DICT
        else:
            x = co_cases[case][0]
            n_traces = CO_N_TRACES
            sem = CO_SEMANTICS_DICT
        return x, n_traces, sem

    L_records, R_records = [], []
    availability = {'L': {}, 'R': {}}
    reason_hist = {}
    tracemalloc.start()
    mem_base = tracemalloc.get_traced_memory()[1]
    for e, geom, case in rows:
        evaluable = case is not None
        base = {
            'event_id': e['event_id'], 'event_role': e['role'],
            'family': e['family'], 'group_id': e['group_id'],
            'geometry': geom, 'case': case,
            'mother_model_id': e['mother_model_id'],
            'row_type': 'negative_control' if e['role'] in
                        ('nc_zero', 'bg_absent', 'off_path') else 'event',
            'dev_group': e['group_id'] in DEV_GROUPS,
            'window': {'index_lo': int(e['index_lo']) if evaluable else None,
                       'index_hi': int(e['index_hi']) if evaluable else None,
                       'n_window_samples': e['n_window_samples'],
                       'window_hash': e['window_hash'], 'mask_hash': e['mask_hash'],
                       'time_axis_convention_version': e['time_axis_convention_version']},
            'provenance': {
                'event_table_sha256': table_sha,
                's4_contract_sha256': contract_sha,
                'operator_script_sha256': op_script_sha,
                'evaluation_contract_sha256': eval_contract_sha,
                'runner_script_sha256': runner_sha,
                'freeze_script_sha256': freeze_script_sha,
                'operator_version': OP_VERSION,
                'dt_s': dt_s,
            },
        }
        if not evaluable:
            rec = dict(base, record_id=f"{e['event_id']}::{geom}:{case or 'NO_CASE'}::L",
                       baseline='L', availability='unavailable',
                       failure_reason='no_case_in_geometry')
            L_records.append(rec)
            availability['L']['unavailable'] = availability['L'].get('unavailable', 0) + 1
            reason_hist['no_case_in_geometry'] = reason_hist.get('no_case_in_geometry', 0) + 1
            rec = dict(base, record_id=f"{e['event_id']}::{geom}:{case or 'NO_CASE'}::R",
                       baseline='R', availability='unavailable',
                       failure_reason='no_case_in_geometry')
            R_records.append(rec)
            availability['R']['unavailable'] = availability['R'].get('unavailable', 0) + 1
            continue

        x, n_traces, sem = get_window(geom, case)
        lo, hi = int(e['index_lo']), int(e['index_hi'])
        assert 0 <= lo <= hi and hi + 1 <= x.shape[0]
        xw = x[lo:hi + 1, :]
        assert np.isfinite(xw).all()

        # --- L: local mean external grid ------------------------------------
        for lam in L_LAMBDAS:
            for w in L_WIDTHS:
                for q in L_GAINS:
                    cid = f"L_lam{lam}_w{w}_q{q}"
                    rec = dict(base, record_id=f"{e['event_id']}::{geom}:{case}::{cid}",
                               baseline='L', candidate_id=cid,
                               gather_semantics=dict(sem, n_traces=n_traces),
                               strength_lambda=lam, window_width=w, end_gain=q,
                               requested_order='LM_then_G')
                    if w > n_traces:
                        rec.update(availability='unavailable',
                                   failure_reason='width_exceeds_n_traces_no_auto_shortening',
                                   metrics=None, gain_diagnostics=None, resource=None)
                        availability['L']['unavailable'] = availability['L'].get('unavailable', 0) + 1
                        reason_hist['width_exceeds_n_traces_no_auto_shortening'] = \
                            reason_hist.get('width_exceeds_n_traces_no_auto_shortening', 0) + 1
                        L_records.append(rec)
                        continue
                    t0 = time.perf_counter()
                    pre = local_mean_control(xw, w, strength=lam)
                    y = pre
                    if q > 1:
                        y = fixed_gain_curve(y.shape[0], q) * pre
                    assert np.isfinite(y).all()
                    wall = time.perf_counter() - t0
                    peak = tracemalloc.get_traced_memory()[1] - mem_base
                    availability['L']['ran'] = availability['L'].get('ran', 0) + 1
                    mask = np.ones(y.shape, dtype=bool)
                    nb = energy_metrics(y, xw, mask)
                    gdiag = None
                    if q > 1:
                        path, fsha = persist_array(
                            f"{e['event_id']}__{geom}_{case}__{cid}__pregain", pre)
                        curve = fixed_gain_curve(y.shape[0], q)[:, 0]
                        cpath, csha = persist_array(
                            f"{e['event_id']}__{geom}_{case}__{cid}__gaincurve", curve)
                        gdiag = {
                            'persisted_array': path, 'persisted_array_sha256': fsha,
                            'persisted_kind': 'pre_gain_true_intermediate',
                            'max_gain': q,
                            'max_abs_input': float(np.max(np.abs(xw))),
                            'max_abs_pre_gain': float(np.max(np.abs(pre))),
                            'max_abs_output': float(np.max(np.abs(y))),
                            'gain_energy_ratio_output_over_pre': float(np.mean(y * y) / np.mean(pre * pre)),
                            'clipped_samples': int(np.count_nonzero(
                                np.abs(y) > q * np.max(np.abs(pre)) + 1e-300)),
                            'overflow': bool(not np.isfinite(y).all()),
                            'gain_curve_ref': cpath, 'gain_curve_sha256': csha,
                            'gain_curve_first': float(curve[0]), 'gain_curve_end': float(curve[-1]),
                        }
                    rec.update(
                        availability='ran', failure_reason=None,
                        executed_order='LM' + ('G' if q > 1 else ''),
                        metrics={'N_b': nb},
                        gain_diagnostics=gdiag,
                        resource={'wall_s': wall, 'peak_tracemalloc_bytes': int(peak),
                                  'compute_dtype': 'float64',
                                  'n_samples': int(y.shape[0]), 'n_traces': int(y.shape[1])},
                    )
                    L_records.append(rec)

        # --- R: input-spectrum rule diagnostics ------------------------------
        scale = float(np.max(np.abs(xw)))
        svd_vals = np.linalg.svd(xw / scale if scale > 0 else xw,
                                 compute_uv=False)
        relative = svd_vals / svd_vals[0]
        rank_floor = max(xw.shape) * np.finfo(np.float64).eps
        numerical_rank = int(np.count_nonzero(relative > rank_floor))
        gaps = [float(relative[k - 1] - relative[k]) for k in R_KS]
        k_star = int(np.argmax(gaps)) + 1
        # same rejection口径 as the operator contract (gap_rtol=1e-8): a cutoff
        # whose gap is an unresolved tie below the numerical rank floor, or
        # whose singular value itself is at/below the floor, selects identity.
        rejected = bool(
            (k_star < numerical_rank and gaps[k_star - 1] <= 1e-8)
            or relative[k_star - 1] <= rank_floor)
        for tau in R_THRESHOLDS:
            cid = f"R_tau{tau}_q1"
            select_k = None if (max(gaps) < tau or rejected) else k_star
            chosen = 'identity' if select_k is None else f'svd_k{select_k}'
            t0 = time.perf_counter()
            if select_k is None:
                y = xw.copy()
                executed = 'identity'
            else:
                res = apply_configuration(xw, f'B{3 + select_k}_G1_BG')
                y = res['output']
                executed = f'svd_k{select_k}_BG'
            wall = time.perf_counter() - t0
            peak = tracemalloc.get_traced_memory()[1] - mem_base
            availability['R']['ran'] = availability['R'].get('ran', 0) + 1
            mask = np.ones(y.shape, dtype=bool)
            nb = energy_metrics(y, xw, mask)
            R_records.append(dict(
                base, record_id=f"{e['event_id']}::{geom}:{case}::{cid}",
                baseline='R', candidate_id=cid,
                gather_semantics=dict(sem, n_traces=n_traces),
                threshold_tau=tau, q=1,
                gap_spectrum_d1_d2_d3=gaps,
                k_star=k_star, max_gap=max(gaps),
                numerical_rank_floor=rank_floor,
                numerical_rank=numerical_rank,
                selection_rejected_numerically=rejected,
                selection=chosen,
                availability='ran', failure_reason=None,
                executed_order=executed,
                metrics={'N_b': nb},
                resource={'wall_s': wall, 'peak_tracemalloc_bytes': int(peak),
                          'compute_dtype': 'float64',
                          'n_samples': int(y.shape[0]), 'n_traces': int(y.shape[1])},
            ))
    tracemalloc.stop()

    # --- feasibility with all-null limits (G4 unresolved) -------------------
    def feas(records):
        cands = [{'id': r['record_id'], 'available': r['availability'] == 'ran',
                  'metrics': {'nrmse': None, 'amplitude_error': None,
                              'negative_control_residual': None},
                  'failure_reason': r['failure_reason']} for r in records]
        return configuration_labels(cands, limits={'nrmse': None, 'amplitude_error': None,
                                                   'negative_control_residual': None},
                                    objectives=['nrmse', 'amplitude_error'])['status']
    feasibility = {'L': feas(L_records), 'R': feas(R_records)}

    manifest = {
        'run_tag': run_tag,
        'runner': 'scripts/run_baseline_s4.py',
        's4_contract': str(S4_CONTRACT.relative_to(ROOT)),
        's4_contract_sha256': contract_sha,
        'event_table': str(EVENT_TABLE.relative_to(ROOT)),
        'event_table_sha256': table_sha,
        'time_axis_convention_version': table['time_axis']['convention_version'],
        'grid_tier': 'BASE',
        'selection_status': contract['selection_status'],
        'L_grid': {'lambdas': L_LAMBDAS, 'widths': L_WIDTHS, 'gains': L_GAINS,
                   'legality_rule': 'width <= n_traces else unavailable (no auto-shortening)',
                   'order': 'local_mean_then_fixed_shared_gain (no GB variant)'},
        'R_rule': {'gaps': 'd_k=(sigma_k-sigma_{k+1})/sigma_1, k in {1,2,3}',
                   'threshold_candidates': R_THRESHOLDS,
                   'threshold_status': 'uncalibrated; mechanism diagnostics only; q=1 ablation',
                   'numerical_rejection': 'ties below numerical rank floor -> identity'},
        'F_view': {
            'computation': 'none',
            'reason': 'F single-config selection undetermined (G4); dev-group projection only',
            'eval_records': {'mt': MT_EVAL_R1, 'co': CO_EVAL_R1},
            'dev_groups': DEV_GROUPS, 'test_groups': TEST_GROUPS,
            'protocol': contract['F_protocol'],
        },
        'O_status': 'unavailable_no_computable_reference_no_approximation',
        'event_rows': len(rows),
        'n_L_records': len(L_records),
        'n_R_records': len(R_records),
        'availability': availability,
        'failure_reason_histogram': reason_hist,
        'geometry_layering': 'MT and CO records layered separately; never merged or cross-paired; per-record gather_semantics.version carried',
        'gates': {'s4_contract_frozen': True, 'event_table_frozen': True,
                  'contract_input_hashes_verified': True,
                  'mt_analysis_pass': True, 'co_analysis_pass': True,
                  'pre_gain_persisted': True, 'negative_controls_separate': True},
        'runner_limits': {'r_c_computed': False, 'isolated_event_granted': False,
                          'cross_family_metrics': False, 'gain_reference_absent': True,
                          'waveform_metrics_computed': False,
                          'waveform_metrics_reason': 'reference_state numerically_unresolved + paired_contrast (design A.4)',
                          'd_a_h_rho': 'not_computable',
                          'thresholds_all_null': True},
        'feasibility': {'status': feasibility,
                        'note': 'all limits null (G4 unresolved): every gate undetermined; no ranking, no unique label'},
        'dev_groups': DEV_GROUPS,
        'test_groups': TEST_GROUPS,
        'hard_limits': {
            'absolute_accuracy_claimed': False,
            'clean_truth_generated': False,
            'conclusions_scope': 'S4 baseline mechanism diagnostics on frozen event windows (MT variable-offset gather + CO common-offset profile, C-families) only; not field performance, not 3D, no absolute accuracy',
            'cross_family_differencing': False,
            'cross_tier_differencing': False,
            'fdtd_solver_invoked': False,
            'grid_convergence_certified': False,
            'grid_tier': 'BASE',
            'grid_tiers_merged': False,
            'physical_acceptance_threshold': None,
            'physical_label_eligible': False,
            'reference_state': 'numerically_unresolved',
            'training_eligible': False,
            'training_labels_generated': False,
            'negative_controls_merged_into_targets': False,
            'N_b_is': 'absolute residual diagnostic; not detectability, not SNR, not a physical threshold',
            'selection_made': False,
            'selection_status': 'F config / R tau / L params all undetermined (G4); O unavailable',
        },
        'total_wall_s': time.perf_counter() - started,
        'stop_reason': None,
        'arrays_dir': str(ARRAYS_DIR.relative_to(ROOT)),
    }
    return L_records, R_records, manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True,
                   help='must not exist (no overwriting of evidence)')
    p.add_argument('--run-tag', default='r1')
    a = p.parse_args()
    out = (ROOT / a.output_dir).resolve() if not a.output_dir.is_absolute() else a.output_dir
    assert not out.exists(), f'refusing to overwrite {out}'
    L_records, R_records, manifest = run_baselines(a.run_tag)
    out.mkdir(parents=True)
    (out / 'records.json').write_text(
        json.dumps({'L_records': L_records, 'R_records': R_records}, indent=1) + '\n',
        encoding='utf-8')
    (out / 'run_manifest.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({'L_records': len(L_records), 'R_records': len(R_records),
                      'availability': manifest['availability'],
                      'feasibility': manifest['feasibility']['status'],
                      'output': str(out)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
