"""Known-error direction controls; no FDTD, field data, or physical thresholds."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from research_evaluation_contract import waveform_metrics


def run():
    checks = []

    def check(name, ok):
        if not ok:
            raise AssertionError(name)
        checks.append({'id': name, 'passed': True})

    s = np.array([[1., -1.], [0., 0.]])
    mask = np.ones(s.shape, dtype=bool)

    def evaluate(y, ref=s, bound=None):
        return waveform_metrics(y, ref, mask, reference_kind='complete_clean',
                                state='isolated', direction_error_bound=bound)

    missing = evaluate(s)
    check('undeclared_budget_abstains', missing['metrics']['signed_cosine'] is None
          and missing['metric_reasons']['signed_cosine'] == 'direction_error_bound_not_declared')
    check('missing_direction_does_not_drop_fidelity', missing['available']
          and missing['metrics']['nrmse'] == 0 and missing['metrics']['amplitude_error'] == 0)
    check('exact_fixture_identity', evaluate(s, bound=0.)['metrics']['signed_cosine'] == 1.)
    check('exact_fixture_polarity', evaluate(-s, bound=0.)['metrics']['signed_cosine'] == -1.)
    check('zero_output_missing', evaluate(0*s, bound=0.)['metric_reasons']['signed_cosine'] == 'zero_output')
    for sign in (-1, 1):
        r = evaluate(sign*1e-16*s, bound=1e-14)
        check(f'roundoff_sign_{sign}_missing', r['metrics']['signed_cosine'] is None
              and r['metric_reasons']['signed_cosine'] == 'output_direction_numerically_unresolved')
        check(f'roundoff_sign_{sign}_loss_retained', np.isclose(r['metrics']['nrmse'], 1.)
              and np.isclose(r['metrics']['amplitude_error'], 1.))
    small = evaluate(1e-8*s, bound=1e-12)
    check('weak_resolved_output_not_zeroed', np.isclose(small['metrics']['signed_cosine'], 1.)
          and small['metrics']['amplitude_factor'] > 0)
    baseline = evaluate(.5*s, bound=.01)
    for factor in (1e-100, 1e-12, 1e12, 1e100):
        scaled = evaluate(.5*s*factor, ref=s*factor, bound=.01*factor)
        check(f'unit_scale_invariant_{factor}', all(np.isclose(scaled['metrics'][k], value)
              for k, value in baseline['metrics'].items()))
    tiny = evaluate(s*1e-100, ref=s*1e-100, bound=1e-110)
    check('tiny_absolute_signal_resolved', tiny['metrics']['nrmse'] == 0
          and np.isclose(tiny['metrics']['signed_cosine'], 1.))
    norm = float(np.linalg.norm(s))
    check('boundary_is_unresolved', evaluate(s, bound=norm)['metrics']['signed_cosine'] is None)
    check('above_boundary_is_defined', evaluate(s, bound=norm*.99)['metrics']['signed_cosine'] is not None)
    for value in (-1., float('inf'), float('nan')):
        try:
            evaluate(s, bound=value)
        except ValueError:
            rejected = True
        else:
            rejected = False
        check(f'invalid_budget_{value}', rejected)
    bad_state = waveform_metrics(s, s, mask, reference_kind='isolated_event',
                                 state='numerically_unresolved', direction_error_bound=0.)
    check('budget_does_not_promote_physical_reference', not bad_state['available'])
    root = Path(__file__).parent
    return {'schema': 'direction-diagnostic-checks/1', 'checks': checks,
            'evidence_level': 'constructed_array_error_budgets_only', 'numpy': np.__version__,
            'physical_thresholds_frozen': False, 'fdtd_executed': False, 'training_executed': False,
            'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                             [Path(__file__), root/'research_evaluation_contract.py']}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'passed': len(result['checks']), 'output': str(args.output)}))
