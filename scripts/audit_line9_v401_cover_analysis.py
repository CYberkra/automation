"""Recompute four native spectra and fixed-window metrics independently of analysis."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from audit_line9_v401_delivery import native_response, read, sha


def main(a):
    assert not a.out.exists()
    c = read(a.source / 'execution_contract.json')
    v = read(a.source / 'completed_verification.json')
    p = read(a.public / 'analysis.json')
    assert v['completed'] and p['contract_sha256'] == v['contract_sha256'] == sha(a.source / 'execution_contract.json')
    assert p['verification_sha256'] == sha(a.source / 'completed_verification.json')
    assert p['numerical_sha256'] == sha(a.numerical)
    ids = ['base_snapshot_H0', 'flat_snapshot_H0', 'far_cover_removed_H0', 'near_cover_removed_H0']
    assert [g['id'] for g in c['groups']] == [g['id'] for g in v['groups']] == ids
    assert v['groups'][0]['observer_receiver_and_source_bitwise_equal']
    f = 20e6 + np.arange(501) * 300000.
    time_ns = np.arange(4008) / (4008 * 300000.) * 1e9
    direct = []
    errors = []
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:], f)
        assert h['ids'].asstr()[:].tolist() == ids
        stored = h['response'][:]
    for j, row in enumerate(v['groups']):
        path = a.source / (row['id'] + '.h5')
        q = native_response(path, .025, row['native_sha256'])
        direct.append(q)
        error = float(np.linalg.norm(q - stored[:, j]) / np.linalg.norm(q))
        assert error < 1e-9
        errors.append(error)
    z = np.column_stack(direct)
    gates = dict(early=[0., 120.], wide=[300., 450.], inherited_basal=[345.6180971390552, 374.6081170991351], late=[450., 1100.])
    checks = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w = w / w.mean()
        count, worst = 0, 0.
        assert set(p['metrics'][window]) == set(gates)
        for name, bounds in gates.items():
            record = p['metrics'][window][name]
            assert record['gate_ns'] == bounds and [r['id'] for r in record['rows']] == ids
            keep = (time_ns >= bounds[0]) & (time_ns <= bounds[1])
            x = np.exp(2j * np.pi * time_ns[keep, None] * 1e-9 * f) @ (z * w[:, None]) / 501
            ref = x[:, 0]
            norm = np.linalg.norm
            for j, row in enumerate(record['rows']):
                q = x[:, j]
                values = dict(norm=float(norm(q)), norm_over_base=float(norm(q) / norm(ref)),
                    difference_over_base=float(norm(ref - q) / norm(ref)),
                    correlation_with_base=float(abs(np.vdot(ref, q)) / (norm(ref) * norm(q))),
                    peak_ns=float(time_ns[keep][np.argmax(abs(q))]))
                for key, value in values.items():
                    error = abs(value - row[key]) / max(1., abs(value))
                    assert error < 1e-8, (window, name, row['id'], key, value, row[key])
                    worst = max(worst, error)
                    count += 1
        checks[window] = dict(scalar_checks=count, max_scaled_error=worst)
    result = dict(status='PASS_INDEPENDENT_RAW_DFT_DIRECT_INVERSE_FOUR_CONTROLS',
        auditor_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        contract_sha256=sha(a.source / 'execution_contract.json'), independent_DFT_relative_L2=errors,
        checks=checks, limits='Audits native trace identity and arithmetic only. Full snapshot fields audited on remote host, not re-read here. No pure-path identification, real material calibration or field equivalence.')
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
