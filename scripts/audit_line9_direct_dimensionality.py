"""Independent dyadic/free-space integral/native DFT/direct inverse CPU checks."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import quad
from scipy.special import hankel2
from audit_line9_v401_delivery import native_response, read, sha


def main(a):
    assert not a.out.exists()
    p = read(a.public / 'analysis.json')
    assert p['numerical_sha256'] == sha(a.numerical)
    f = 20e6 + np.arange(501) * 300000.
    c, eps = 299792458., 8.8541878128e-12
    mu = 1 / (eps * c * c)
    r = float(np.linalg.norm(np.array(p['receiver_m']) - p['source_tx_m']))
    assert abs(r - p['air_distance_m']) < 1e-14
    k = 2 * np.pi * f / c
    with h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:], f)
        stored = h['response'][:]
        assert stored.shape == (501, 4)
        assert h['ids'].asstr()[:].tolist() == p['response_columns']
        integrated = h['integrated_3D_line_response_E_over_I'][:]
    # Transverse component of (I+gradgrad/k^2) exp(-ikr)/(4pi r).
    dipole = -1j * (2 * np.pi * f) * mu * np.exp(-1j * k * r) / (4 * np.pi * r) * (1 - 1j / (k * r) - 1 / (k * r)**2)
    line = -(2 * np.pi * f) * mu / 4 * hankel2(0, k * r)
    error_3d = float(np.linalg.norm(dipole - stored[:, 3]) / np.linalg.norm(dipole))
    assert error_3d < 1e-12
    np.testing.assert_allclose(stored[:, 2], line / .025, rtol=1e-13, atol=0)
    integral_errors = np.linalg.norm(integrated - line[:, None], axis=0) / np.linalg.norm(line)
    np.testing.assert_allclose(integral_errors, p['integrated_3D_vs_Hankel_relative_L2'], rtol=0, atol=1e-12)
    assert max(integral_errors) < 1e-8
    # A different scalar3D Green integral; derivative boundary terms vanish at
    # infinity. Adaptive real/imag integration independently checks Hankel sign.
    rotation = np.exp(-1j * np.pi / 4)
    scalar_errors = []
    for frequency in [20e6, 95e6, 170e6]:
        kk = 2 * np.pi * frequency / c
        def integrand(u):
            radius = np.sqrt(r * r + (rotation * u / kk)**2)
            return 2 * rotation / kk * np.exp(-1j * kk * radius) / (4 * np.pi * radius)
        value = quad(lambda u: integrand(u).real, 0, 80, epsabs=1e-12, epsrel=1e-12)[0]
        value += 1j * quad(lambda u: integrand(u).imag, 0, 80, epsabs=1e-12, epsrel=1e-12)[0]
        expected = -1j / 4 * hankel2(0, kk * r)
        error = float(abs(value - expected) / abs(expected))
        assert error < 1e-9
        scalar_errors.append(error)
    errors = []
    for j, sid in enumerate(['high_x19000_H0', 'low_x19000_H0']):
        z = native_response(a.source / (sid + '.h5'), .025, p['native'][sid])
        error = float(np.linalg.norm(z - stored[:, j]) / np.linalg.norm(z))
        assert error < 1e-9
        errors.append(error)
    tn = np.arange(4008) / (4008 * 300000.) * 1e9
    gates = dict(early=[0., 20.], wide=[300., 450.],
                 inherited_basal=[345.6180971390552, 374.6081170991351], late=[450., 1100.])
    count = 0
    worst = 0.
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean()
        early = (tn >= 0) & (tn <= 20)
        x0 = np.exp(2j * np.pi * tn[early, None] * 1e-9 * f) @ (stored * w[:, None]) / 501
        norm0 = np.linalg.norm(x0, axis=0)
        for j in [0, 1]:
            value = float(np.linalg.norm(x0[:, j] - x0[:, 2]) / norm0[j])
            recorded = p['native_H0_vs_2D_direct_unfitted_early_relative_L2'][window][j]
            assert abs(value - recorded) < 1e-9
            count += 1
        for name, bounds in gates.items():
            record = p['metrics'][window][name]
            assert record['gate_ns'] == bounds
            keep = (tn >= bounds[0]) & (tn <= bounds[1])
            x = np.exp(2j * np.pi * tn[keep, None] * 1e-9 * f) @ (stored * w[:, None]) / 501
            values = np.linalg.norm(x, axis=0) / norm0
            for j, value in enumerate(values):
                error = float(abs(value - record['norm_over_own_early'][j]) / max(1, abs(value)))
                assert error < 1e-9
                assert abs(tn[keep][np.argmax(abs(x[:, j]))] - record['peak_ns'][j]) < 1e-10
                worst = max(worst, error)
                count += 2
    result = dict(status='PASS_INDEPENDENT_DYADIC_SCALAR_GREEN_NATIVE_DFT_AND_DIRECT_INVERSE',
        auditor_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        transverse_dyadic_vs_radiation_expression_relative_L2=error_3d,
        adaptive_scalar3D_integral_vs_line_Green_relative_errors=scalar_errors,
        native_DFT_relative_L2=errors, direct_inverse_scalar_checks=count, max_scaled_metric_error=worst,
        calls_solver=False, calls_CUDA=False,
        limits='Mathematical/source/data consistency checks only. No3Dnative geology or antenna/port response. Does not attribute actual late H0 solely to direct wave or validate field calibration.')
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
