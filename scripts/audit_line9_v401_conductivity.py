"""Independent native, single-parameter, direct-DFT and inverse audit.

Uses the independent delivery auditor, never the production SFCW helpers.
"""
import argparse
import copy
import json
from pathlib import Path

import h5py
import numpy as np
from audit_line9_v401_delivery import native_response, read, sha


def main(a):
    assert not a.out.exists()
    pub = read(a.public / 'analysis.json')
    contracts = [read(p / 'execution_contract.json') for p in [a.baseline, a.source]]
    prepared = [read(p / 'manifest.json') for p in [a.baseline_package, a.package]]
    f = 20e6 + np.arange(501) * 300000.
    raw = []
    identities = []
    dft_errors = []
    for root, package, contract, manifest in zip(
            [a.baseline, a.source], [a.baseline_package, a.package], contracts, prepared):
        ver = read(root / 'completed_verification.json')
        assert ver['completed'] and ver['contract_sha256'] == sha(root / 'execution_contract.json')
        responses = []
        for g, record, original in zip(contract['groups'], ver['groups'], manifest['groups']):
            assert g['id'] == record['id'] == original['id']
            for key in ['input', 'geometry', 'material']:
                assert sha(package / original[key]) == g[key + '_sha256'] == original[key + '_sha256']
            p = root / (g['id'] + '.h5')
            q = native_response(p, .025, record['native_sha256'])
            with h5py.File(p) as h:
                assert h.attrs['gprMax'] == '4.0.1'
                assert h.attrs['Iterations'] == g['expected_samples']
                np.testing.assert_allclose(h.attrs['dt'], g['dt_s'], rtol=1e-14, atol=0)
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'], g['native_shape'])
                # Physical decimal positions are stored after integer-grid multiplication.
                np.testing.assert_allclose(h['srcs/src1'].attrs['Position'], g['tx_m'], rtol=0, atol=1e-12)
                np.testing.assert_allclose(h['rxs/rx1'].attrs['Position'], g['rx_m'], rtol=0, atol=1e-12)
                np.testing.assert_array_equal(h['srcs/src1'].attrs['GridPosition'], np.rint(np.array(g['tx_m'])/.025))
                np.testing.assert_array_equal(h['rxs/rx1'].attrs['GridPosition'], np.rint(np.array(g['rx_m'])/.025))
            responses.append(q)
            identities.append(dict(id=g['id'], native_sha256=sha(p)))
        z = np.column_stack([*responses, responses[1] - responses[0]])
        with h5py.File(root / 'sfcw.h5') as h:
            np.testing.assert_array_equal(h['frequency_Hz'][:], f)
            error = float(np.linalg.norm(z - h['response'][:]) / np.linalg.norm(z))
        assert error < 1e-9
        dft_errors.append(error)
        raw.append(z)
    assert pub['numerical_sha256'] == sha(a.source / 'sfcw.h5')
    assert pub['baseline_numerical_sha256'] == sha(a.baseline / 'sfcw.h5')

    for old, new in zip(prepared[0]['groups'], prepared[1]['groups']):
        for key in ['input_sha256', 'geometry_sha256', 'tx_m', 'rx_m', 'native_shape', 'dt_s']:
            assert old[key] == new[key]
        before = read(a.baseline_package / old['material'])
        after = read(a.package / new['material'])
        mud_old = before['materials']['material_002_mudstone']
        mud_new = after['materials']['material_002_mudstone']
        assert mud_old['base']['electric_conductivity_s_per_m'] == .003
        assert mud_new['base']['electric_conductivity_s_per_m'] == .0003
        restored = copy.deepcopy(after)
        restored['database'] = copy.deepcopy(before['database'])
        restored['materials']['material_002_mudstone']['metadata'] = copy.deepcopy(mud_old['metadata'])
        restored['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m'] = .003
        assert restored == before
        def spectrum(m):
            base = m['base']
            value = np.full(501, base['relative_permittivity'], dtype=complex)
            for pole in m['poles']:
                value += pole['relative_permittivity_difference'] / (1 + 2j*np.pi*f*pole['relaxation_time_s'])
            return value - 1j*base['electric_conductivity_s_per_m']/(2*np.pi*f*8.8541878128e-12)
        np.testing.assert_array_equal(spectrum(mud_old).real, spectrum(mud_new).real)
        for name in ['planar_H0', 'planar_H1']:
            with h5py.File(a.baseline/(name+'.h5')) as h0, h5py.File(a.source/(name+'.h5')) as h1:
                np.testing.assert_array_equal(h0['srcs/src1/excitation/samples'][:], h1['srcs/src1/excitation/samples'][:])

    t = np.arange(4008)/(4008*300000.)
    lo, hi = contracts[1]['study_manifest']['basal_gate_ns']
    keep = (t*1e9 >= lo) & (t*1e9 <= hi)
    e = np.exp(2j*np.pi*t[keep, None]*f)
    metrics = {}
    for name, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w = w/w.mean()
        old, new = [e@(z*w[:, None])/501 for z in raw]
        norm = np.linalg.norm
        corr = lambda x, y: float(abs(np.vdot(x, y))/(norm(x)*norm(y)))
        values = dict(baseline_H0_over_delta=float(norm(old[:, 0])/norm(old[:, 2])),
                      changed_H0_over_delta=float(norm(new[:, 0])/norm(new[:, 2])),
                      changed_delta_over_baseline_delta=float(norm(new[:, 2])/norm(old[:, 2])),
                      changed_H0_over_baseline_H0=float(norm(new[:, 0])/norm(old[:, 0])),
                      changed_H1_vs_delta_correlation=corr(new[:, 2], new[:, 1]),
                      baseline_H1_vs_delta_correlation=corr(old[:, 2], old[:, 1]))
        for key, value in values.items():
            assert abs(value - pub['metrics'][name][key]) < 1e-9
        metrics[name] = values
    result = dict(status='PASS_INDEPENDENT_SINGLE_FACTOR_NATIVE_DFT_AND_INVERSE',
                  script_sha256=sha(__file__), analysis_sha256=sha(a.public/'analysis.json'),
                  physical_change_only='mudstone sigmaDC .003 -> .0003 S/m',
                  full_501_real_permittivity_unchanged=True, native_identities=identities,
                  independent_DFT_relative_L2=dft_errors, direct_inverse_metrics=metrics,
                  limits='Planar material sensitivity only; not site calibration, full-line or 3D certification.')
    a.out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'baseline', 'package', 'baseline-package', 'public', 'out']:
        p.add_argument('--'+key, type=Path, required=True)
    main(p.parse_args())
