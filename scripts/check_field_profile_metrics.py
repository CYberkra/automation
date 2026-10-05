"""Counterexamples for endpoint contamination and signed-profile peak metrics."""
import json
import numpy as np
from scipy.signal import hilbert
from field_profile_metrics import late_peak_metrics


def main():
    time = np.linspace(0, 700, 501)
    y = np.zeros_like(time)
    y[10] = 2.; y[300] = .002; y[-1] = .5
    result = late_peak_metrics(abs(y), time)
    assert abs(result['median_dB']+60) < 1e-12
    assert result['median_peak_time_ns'] == time[300]
    assert late_peak_metrics(abs(-10*y), time)['median_dB'] == result['median_dB']
    for bad, axis, stop in [(abs(y), time, 700), (abs(y), time[::-1], 650),
                            (np.zeros_like(y), time, 650), (abs(y)[:-1], time, 650)]:
        try:
            late_peak_metrics(bad, axis, late_stop_ns=stop)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid metric input accepted')
    # An early-only record has exactly zero signed samples in the late window,
    # yet a periodic FFT Hilbert envelope creates a large endpoint response.
    early_only = np.zeros_like(time); early_only[0] = 1.; early_only[10] = 10.
    env = abs(hilbert(early_only))
    assert np.max(abs(early_only)[time >= 300]) == 0
    assert np.max(env[time > 650]) > 10*np.max(env[(time >= 300)&(time <= 600)])
    print(json.dumps({'status': 'PASS', 'known_late_pulse_dB': result['median_dB'],
                      'endpoint_outlier_excluded': True, 'sign_scale_invariant': True,
                      'invalid_inputs_rejected': 4, 'early_only_fft_endpoint_counterexample': True,
                      'solver_called': False}))


if __name__ == '__main__':
    main()
