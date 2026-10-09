"""Analytic Maxwell plane-wave and interpolation controls, no solver."""
import json
import numpy as np
from analyze_line9_snapshot_directions import align_h, transport


def main():
    checks = []
    impedance = 376.730313668
    e = np.array([.2, -.4, 1.3])
    for angle in (0, 45, 90, 180, 270):
        radians = np.deg2rad(angle)
        hx, hy = e * np.sin(radians)/impedance, -e*np.cos(radians)/impedance
        sx, sy = transport(e, hx, hy)
        expected = e**2 / impedance
        np.testing.assert_allclose(sx, expected*np.cos(radians), rtol=1e-14, atol=1e-17)
        np.testing.assert_allclose(sy, expected*np.sin(radians), rtol=1e-14, atol=1e-17)
        checks.append(f'Maxwell_progressive_wave_{angle}_degrees')
    # Equal counterpropagating waves make a standing field: E is not zero,
    # yet the whole-cycle mean transport vanishes. Field brightness is no path proof.
    phase = 2*np.pi*np.arange(4096)/4096
    x = .37
    upward = np.cos(phase-x); downward = np.cos(phase+x)
    ez = upward+downward; hx = (upward-downward)/impedance
    _, sy = transport(ez, hx, np.zeros_like(hx))
    assert np.max(abs(ez)) > 1 and abs(np.mean(sy)) < 1e-17
    assert np.max(abs(sy)) > 1e-4
    checks.append('standing_wave_bright_field_zero_cycle_mean_nonzero_instantaneous_flux')
    alpha = .5/34
    f = lambda t: 2+3*t
    for method in ('forward_linear', 'backward_linear', 'quadratic'):
        np.testing.assert_allclose(align_h(f(-1),f(0),f(1),alpha,method), f(alpha), atol=1e-15, rtol=0)
        checks.append(f'linear_time_polynomial_{method}')
    f = lambda t: 2+3*t+5*t*t
    np.testing.assert_allclose(align_h(f(-1),f(0),f(1),alpha,'quadratic'), f(alpha), atol=1e-15, rtol=0)
    checks.append('quadratic_time_polynomial_exact')
    for alpha, method in [(-.01,'forward_linear'),(1,'quadratic'),(.01,'unknown')]:
        try:
            align_h(1.,2.,3.,alpha,method)
        except ValueError:
            checks.append(f'reject_invalid_alignment_{alpha}_{method}')
        else:
            raise AssertionError('Invalid interpolation accepted')
    print(json.dumps({'status':'PASS_ANALYTIC_DIRECTION_CONTROLS','count':len(checks),'checks':checks,
                      'limits':'Analytic controls verify arithmetic and sign conventions only, not snapshot physical accuracy.'}))


if __name__ == '__main__':
    main()
