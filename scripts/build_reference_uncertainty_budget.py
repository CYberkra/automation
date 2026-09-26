"""G4 step 1: build the reference uncertainty budget from archived analyses.

CPU read-only; no solver, no timestamps, no network. Deterministic output:
run twice (r1/r2) into separate directories and compare results.json bytes
(see scripts/check_reference_uncertainty_budget_acceptance.py).

Design authority: docs/research/2026-09-26_g4_calibration_plan_v0.1.md §2.
The budget aggregates, per frozen event window and per metric class
(arrival-time / shape / amplitude), the reference-uncertainty magnitudes
measured in four archived analyses:

  E6  grid chain      artifacts/research_checks/2026-09-25_deep_zfine2_analysis_r1
  E7  dimension diff  artifacts/research_checks/2026-09-26_dep3d_gold_analysis_r1
  fine2 tier direction artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1
  error budget        artifacts/research_checks/2026-09-25_error_budget

The report does NOT upgrade reference_state; every number here is an
adjacent-difference / direction magnitude, not an error versus an exact
solution (the archived interpretations are carried verbatim).
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

EVENT_TABLE = ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json'
EVENT_TABLE_SHA256 = 'b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c'

E6_RESULTS = ROOT / 'artifacts/research_checks/2026-09-25_deep_zfine2_analysis_r1/results.json'
E7_RESULTS = ROOT / 'artifacts/research_checks/2026-09-26_dep3d_gold_analysis_r1/results.json'
FINE2_RESULTS = ROOT / 'artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1/results.json'
EB_RESTORATION = ROOT / 'artifacts/research_checks/2026-09-25_error_budget/restoration_results.json'
EB_ATTESTATION = ROOT / 'artifacts/research_checks/2026-09-25_error_budget/three_grid_attestation.json'

PLAN = 'docs/research/2026-09-26_g4_calibration_plan_v0.1.md'

JUDGMENT_RULE_ZH = '算子间指标差异 ≤ 对应窗参考不确定性 ⇒ 该差异记 undetermined，不得进入排序'
JUDGMENT_RULE_EN = ('an inter-operator metric difference at or below the reference uncertainty '
                    'of its event window is recorded undetermined and must not enter ranking')

AMPLITUDE_G1_G6_NOTE = ('amplitude-class assertions are G1/G6-limited: quantitative amplitude '
                        'claims require a real antenna model (Warren 2016, gprMax antenna '
                        'modelling guidance); all amplitude magnitudes here are simulation-domain '
                        'diagnostic quantities, not field amplitudes')


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def extract_fine2(p2):
    """Per-family BASE<->FINE2 tier-difference magnitudes (direct evidence for
    the batch2d_v1 c1/c3/c5/c8 event windows)."""
    out = {}
    for fam in ('C1', 'C3', 'C5', 'C8'):
        r = p2[fam]
        env = r['envelope_peak_direction']
        amp = r['amplitude_ratio_direction']
        sign = r['per_frequency_sign_agreement']
        out[fam.lower()] = {
            'grid_tiers_compared': [r['coarse_grid_tier'], r['fine2_grid_tier']],
            'anchor_pair': {
                'coarse': r['coarse_anchor_run_id'],
                'fine2': r['fine2_anchor_run_id'],
            },
            'envelope_peak_shift_ns': {
                'coarse_base_ns': env['coarse_base_ns'],
                'fine2_ns': env['fine2_ns'],
                'fine2_minus_base_ns': env['fine2_minus_base_ns'],
                'direction': env['direction'],
            },
            'difference_spectrum_shape_correlation_main_200ns':
                r['difference_spectrum_shape_correlation']['main_200ns'],
            'difference_spectrum_shape_correlation_deviation_from_unity':
                1.0 - r['difference_spectrum_shape_correlation']['main_200ns'],
            'per_frequency_sign_agreement_real': sign['per_frequency_sign_agreement_real'],
            'per_frequency_sign_agreement_imag': sign['per_frequency_sign_agreement_imag'],
            'spectral_L2_ratio_FINE2_over_BASE': amp['spectral_L2_ratio_FINE2_over_BASE'],
            'time_peak_abs_ratio_FINE2_over_BASE': amp['time_peak_abs_ratio_FINE2_over_BASE'],
            'spectral_peak_frequency_MHz': r['spectral_peak_frequency_MHz'],
            'interpretation': r['interpretation'],
            'threshold': r['threshold'],
        }
    return out


def build_budget(et, e6, e7, fine2, eb_rest, eb_attest):
    fam = extract_fine2(fine2['p2_direction_consistency'])

    c2o = e6['comparisons_to_original_2p5cm']
    adj1 = e6['zfine_vs_uniform_fine']
    adj2 = e6['zfine2_vs_zfine']
    grid_chain = {
        'source': str(E6_RESULTS.relative_to(ROOT)),
        'family_scope': 'deep DEP scenario family (DEP_BG/DEP_20 chain, 2D)',
        'chain_full_band_relative_L2': {
            'original_2p5cm_to_FINE': c2o['FINE']['all_relative_L2_difference'],
            'FINE_to_ZFINE': adj1['relative_L2_difference'],
            'ZFINE_to_ZFINE2': adj2['relative_L2_difference'],
        },
        'chain_max_abs_phase_deg': {
            'original_2p5cm_to_FINE': c2o['FINE']['max_abs_phase_error_deg'],
            'FINE_to_ZFINE': adj1['max_abs_phase_error_deg'],
            'ZFINE_to_ZFINE2': adj2['max_abs_phase_error_deg'],
        },
        'chain_max_abs_amplitude_error_dB': {
            'original_2p5cm_to_FINE': c2o['FINE']['max_abs_amplitude_error_dB'],
            'FINE_to_ZFINE': adj1['max_abs_amplitude_error_dB'],
            'ZFINE_to_ZFINE2': adj2['max_abs_amplitude_error_dB'],
        },
        'vertical_second_order_ratios': {
            'relative_L2_FINE_to_ZFINE_over_ZFINE_to_ZFINE2':
                adj1['relative_L2_difference'] / adj2['relative_L2_difference'],
            'phase_FINE_to_ZFINE_over_ZFINE_to_ZFINE2':
                adj1['max_abs_phase_error_deg'] / adj2['max_abs_phase_error_deg'],
        },
        'tail_window_relative_L2_by_tier': e6['tail_window_relative_L2'],
        'interpretation_carried_verbatim': e6['interpretation'],
        'grid_convergence_certified': False,
        'richardson_extrapolation_performed': False,
    }

    cd = e7['cross_dimension_same_spacing']
    cross_dim = {
        'source': str(E7_RESULTS.relative_to(ROOT)),
        'family_scope': 'dep3d_gold_v1 scenario family (deep scenario, 2D vs 3D same grid spacing)',
        'confounds_carried_verbatim': cd['confounds'],
        'pairs': {
            name: {
                'envelope_peak_time_difference_ns': v['envelope_peak_time_difference_ns'],
                'spectral_peak_frequency_difference_MHz': v['spectral_peak_frequency_difference_MHz'],
                'normalized_spectrum_shape_correlation': v['normalized_spectrum_shape_correlation'],
                'waveform_max_normalized_correlation': v['waveform_max_normalized_correlation'],
                'waveform_optimal_lag_ns': v['waveform_optimal_lag_ns'],
                'normalized_mag_trend_difference_dB': v['normalized_mag_trend_difference_dB'],
                'phase_difference_residual_std_deg':
                    v['phase_difference_linear_fit']['residual_std_deg'],
                'far_field_exploratory_shape_correlation':
                    v['far_field_3d_to_2d_exploratory']['transformed_vs_2D_shape_correlation'],
            }
            for name, v in cd.items() if isinstance(v, dict) and 'envelope_peak_time_difference_ns' in v
        },
        'convergence_3d_dz_refinement': {
            'comparison': e7['convergence_3d_dz_refinement']['comparison'],
            'envelope_peak_time_difference_ns':
                e7['convergence_3d_dz_refinement']['envelope_peak_time_difference_ns'],
            'normalized_spectrum_shape_correlation':
                e7['convergence_3d_dz_refinement']['normalized_spectrum_shape_correlation'],
            'waveform_max_normalized_correlation':
                e7['convergence_3d_dz_refinement']['waveform_max_normalized_correlation'],
            'threshold_commitment':
                e7['convergence_3d_dz_refinement']['threshold_commitment'],
        },
        'vs_dep2d_zfine2_chain': e7['vs_dep2d_zfine2_chain'],
        'hard_limits_carried': e7['hard_limits'],
    }

    constructed_array = {
        'restoration_source': str(EB_RESTORATION.relative_to(ROOT)),
        'attestation_source': str(EB_ATTESTATION.relative_to(ROOT)),
        'evidence_level': eb_rest['evidence_level'],
        'statistics': eb_rest['statistics'],
        'budget_max': eb_rest['statistics']['budget_max'],
        'limitations_carried_verbatim': eb_rest['limitations'],
        'attestation_reason': eb_attest['refined_reason'],
        'attestation_note_carried_verbatim': eb_attest['note'],
        'applicability': ('synthetic constructed-array processing certificates only; the budget '
                          'certifies numerical identifiability of array algebra, not physical '
                          'accuracy; it does not bound FDTD trace uncertainty'),
    }

    scenario_level = {
        'grid_chain_E6': grid_chain,
        'cross_dimension_E7': cross_dim,
        'constructed_array_identifiability': constructed_array,
        'applicability_note': ('E6/E7 are deep-scenario (DEP/dep3d) family evidence; for batch2d_v1 '
                               'c1-c8 windows they are cross-scenario reference magnitudes, not '
                               'measurements on the window itself. The per-window direct evidence '
                               'is the fine2 BASE<->FINE2 tier difference of the same family.'),
    }

    windows = []
    for e in et['entries']:
        f = fam[e['family']]
        windows.append({
            'event_id': e['event_id'],
            'family': e['family'],
            'role': e['role'],
            'group_id': e['group_id'],
            'mother_model_id': e['mother_model_id'],
            't_lo_ns': e['t_lo_ns'],
            't_hi_ns': e['t_hi_ns'],
            'half_width_ns': e['half_width_ns'],
            'arrival_time_class': {
                'direct_evidence': {
                    'grid_tier_shift_ns': f['envelope_peak_shift_ns']['fine2_minus_base_ns'],
                    'direction': f['envelope_peak_shift_ns']['direction'],
                    'provenance': f['anchor_pair'],
                    'comparison_axis_note': ('each tier read on its own native time axis; the two dt '
                                             'differ, so the comparison is of peak index time, not '
                                             'of a common grid'),
                },
                'scenario_reference': 'scenario_level.grid_chain_E6.chain_max_abs_phase_deg '
                                      '(phase-derived; formal time offset only) and '
                                      'scenario_level.cross_dimension_E7 (tens of ns; deep scenario)',
                'magnitude_characterization': ('systematic inter-tier envelope-peak shift of the '
                                               'family anchor pair; differences between operators '
                                               'below this magnitude for the same family are '
                                               'recorded undetermined'),
            },
            'shape_class': {
                'direct_evidence': {
                    'difference_spectrum_shape_correlation_main_200ns':
                        f['difference_spectrum_shape_correlation_main_200ns'],
                    'deviation_from_unity':
                        f['difference_spectrum_shape_correlation_deviation_from_unity'],
                    'per_frequency_sign_agreement_real': f['per_frequency_sign_agreement_real'],
                    'per_frequency_sign_agreement_imag': f['per_frequency_sign_agreement_imag'],
                    'provenance': f['anchor_pair'],
                },
                'scenario_reference': 'scenario_level.cross_dimension_E7.pairs '
                                      '(normalized_spectrum_shape_correlation 0.935-0.958 '
                                      'same-spacing 2D/3D; deep scenario)',
                'magnitude_characterization': ('BASE<->FINE2 difference-spectrum shape correlation '
                                               '0.99813-0.99955 and per-frequency sign agreement '
                                               '88.0%-93.8% across the four families'),
            },
            'amplitude_class': {
                'g1_g6_limited': True,
                'g1_g6_note': AMPLITUDE_G1_G6_NOTE,
                'direct_evidence': {
                    'spectral_L2_ratio_FINE2_over_BASE': f['spectral_L2_ratio_FINE2_over_BASE'],
                    'time_peak_abs_ratio_FINE2_over_BASE': f['time_peak_abs_ratio_FINE2_over_BASE'],
                    'provenance': f['anchor_pair'],
                    'ratio_note': ('cross-tier ratio is a direction indicator only; it mixes grid '
                                   'dispersion, source discretisation and taper sample count, so '
                                   'it is not an accuracy or grid-error measure'),
                },
                'scenario_reference': 'scenario_level.grid_chain_E6.chain_max_abs_amplitude_error_dB '
                                      '(adjacent-grid 0.035-0.142 dB; vs original grid up to 0.77 dB; '
                                      'deep scenario) and scenario_level.cross_dimension_E7 '
                                      '(normalized magnitude trend differences up to ~17.4 dB at low '
                                      'frequency; cross-dimension, confounded)',
            },
        })

    budget = {
        'schema': 'reference-uncertainty-budget/1',
        'batch_id': et['batch_id'],
        'grid_tier_of_event_table': et['grid_tier'],
        'plan': PLAN,
        'plan_section': 'G4 step 1 (§2 reference calibration report v0.1)',
        'event_table': str(EVENT_TABLE.relative_to(ROOT)),
        'event_table_sha256': EVENT_TABLE_SHA256,
        'judgment_rule_preregistered': {
            'zh': JUDGMENT_RULE_ZH,
            'en': JUDGMENT_RULE_EN,
            'same_ruler_as': 'S4/B.4 stable-difference judgment (evaluation protocol v0.2 §7)',
            'operationalization': ('same-metric comparison only: an operator difference expressed '
                                   'in a given metric is compared against the reference-uncertainty '
                                   'magnitude of that same metric class for the same family/window'),
        },
        'reference_state': et['declarations']['reference_state'],
        'reference_state_unchanged_by_this_budget': True,
        'physical_acceptance_threshold': None,
        'n_event_windows': len(windows),
        'per_family_direct_evidence': fam,
        'scenario_level_uncertainty': scenario_level,
        'per_event_window': windows,
        'family_level_notes': {
            'c8_resolution_attention': ('FINE2 C8 anchor difference in-band energy ratio 9.37e-11: '
                                        'the family sits in the numerical-resolution attention zone; '
                                        'no detectability statement is made for C8 windows'),
            'fine2_family_direction': ('all four families are direction-consistent across tiers '
                                       '(fine2 analysis 8/8); nothing is downgraded to pending-recheck'),
        },
        'hard_limits': {
            'absolute_accuracy_claimed': False,
            'clean_truth_generated': False,
            'grid_convergence_certified': False,
            'richardson_extrapolation_performed': False,
            'conclusions_scope': ('reference-uncertainty magnitudes for the frozen batch2d_v1 event '
                                  'windows plus deep-scenario archived evidence; mechanism and '
                                  'calibration input only; not field performance, not absolute '
                                  'accuracy, no threshold values'),
            'negative_controls_merged_into_targets': False,
            'physical_acceptance_threshold': None,
            'physical_label_eligible': False,
            'reference_state': 'numerically_unresolved',
            'solver_invoked': False,
            'training_eligible': False,
            'training_labels_generated': False,
            'cross_family_differencing': False,
            'cross_tier_differencing': False,
        },
        'provenance': {
            str(p.relative_to(ROOT)): sha256_file(p)
            for p in (E6_RESULTS, E7_RESULTS, FINE2_RESULTS, EB_RESTORATION, EB_ATTESTATION)
        },
    }
    return budget


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-dir', required=True)
    args = ap.parse_args()
    out = Path(args.output_dir)
    assert not out.exists(), f'output dir must not exist: {out}'

    actual_sha = sha256_file(EVENT_TABLE)
    assert actual_sha == EVENT_TABLE_SHA256, \
        f'frozen event table SHA mismatch: {actual_sha} != {EVENT_TABLE_SHA256}'

    et = load(EVENT_TABLE)
    budget = build_budget(et, load(E6_RESULTS), load(E7_RESULTS), load(FINE2_RESULTS),
                          load(EB_RESTORATION), load(EB_ATTESTATION))

    out.mkdir(parents=True)
    text = json.dumps(budget, indent=1, sort_keys=True, ensure_ascii=False)
    target = out / 'results.json'
    target.write_bytes(text.encode('utf-8'))
    print(f'wrote {target}')
    print('sha256:', sha256_file(target))
    print('n_event_windows:', budget['n_event_windows'])


if __name__ == '__main__':
    main()
