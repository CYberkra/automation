"""Small regression checks; metadata/synthetic arrays only, no held-out evaluation."""
import copy
import json
import numpy as np
from run_g4_mission_capability import expected_record_layout, validate_records, TEST_FAMILIES
from check_damage_ladder_acceptance import comparison_bytes


def main():
    rows = []
    for key, (family, geometry, case, cid, role, dtype, level) in expected_record_layout(TEST_FAMILIES).items():
        rows.append(dict(record_id=key, family=family, geometry=geometry, case=case,
                         candidate_id=cid, row_type=role, damage=dict(type=dtype, level=level),
                         availability='ran' if cid == 'B0_G1_BG' else 'unavailable',
                         metrics={'waveform': {'metrics': {'nrmse': 0.0}}}))
    validate_records(rows, TEST_FAMILIES)
    bad = {
        'empty': [], 'missing_row': rows[:-1], 'duplicate': rows + [rows[0]],
        'missing_family': [r for r in rows if r['family'] != 'c8'],
        'missing_identity': [r for r in rows if r['candidate_id'] != 'B0_G1_BG'],
    }
    wrong = copy.deepcopy(rows)
    wrong[0]['family'] = 'c1'
    bad['wrong_family'] = wrong
    unavailable = copy.deepcopy(rows)
    next(r for r in unavailable if r['candidate_id'] == 'B0_G1_BG')['availability'] = 'unavailable'
    bad['unavailable_identity'] = unavailable
    for name, data in bad.items():
        try:
            validate_records(data, TEST_FAMILIES)
        except ValueError:
            continue
        raise AssertionError(f'{name} was accepted')
    old = [dict(metrics={'nrmse': 0.25}, resource={'wall_s': 1},
                provenance={'runner_script_sha256': 'a'*64, 'input_sha256': 'c'*64})]
    new = copy.deepcopy(old)
    new[0]['provenance']['runner_script_sha256'] = 'b'*64
    new[0]['resource']['wall_s'] = 2
    assert comparison_bytes(old) != comparison_bytes(new)
    assert comparison_bytes(old, cross_version=True) == comparison_bytes(new, cross_version=True)
    assert old[0]['provenance']['runner_script_sha256'] == 'a'*64
    shared = [dict(old[0]), dict(old[0])]
    assert comparison_bytes(shared, cross_version=True)
    assert shared[0]['provenance']['runner_script_sha256'] == 'a'*64
    for section, key, value in [('metrics', 'nrmse', .26), ('provenance', 'input_sha256', 'd'*64)]:
        changed = copy.deepcopy(new)
        changed[0][section][key] = value
        assert comparison_bytes(old, cross_version=True) != comparison_bytes(changed, cross_version=True)
    x = np.sin(np.arange(128)*.3) + .2*np.cos(np.arange(128)*.7)
    assert np.isclose(np.corrcoef(abs(np.fft.rfft(x)), abs(np.fft.rfft(-x)))[0,1], 1)
    assert np.linalg.norm(-x-x)/np.linalg.norm(x) == 2
    print(json.dumps(dict(expected_test_records=len(rows), rejected_cases=list(bad),
                         regression_exclusions_checked=True, spectrum_counterexample_checked=True,
                         solver_invoked=False, held_out_waveforms_evaluated=False)))


if __name__ == '__main__':
    main()
