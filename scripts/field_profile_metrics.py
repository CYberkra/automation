"""Peak ratios in explicitly bounded interiors of exported real profiles.

FFT-envelope endpoints are excluded from late-event metrics. This is a
diagnostic measurement-window rule, not a geological event classifier.
"""
import numpy as np


def late_peak_metrics(envelope, time_ns, *, early_stop_ns=150., late_start_ns=300., late_stop_ns=650.):
    values = np.asarray(envelope, dtype=float)
    time = np.asarray(time_ns, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if (values.ndim != 2 or time.ndim != 1 or values.shape[0] != len(time)
            or not np.isfinite(values).all() or not np.isfinite(time).all()
            or np.any(values < 0) or np.any(np.diff(time) <= 0)):
        raise ValueError('finite nonnegative [sample,trace] envelope and increasing time required')
    if not time[0] < early_stop_ns <= late_start_ns < late_stop_ns < time[-1]:
        raise ValueError('separate early/late windows and a strictly interior late endpoint required')
    early = time < early_stop_ns
    late = (time >= late_start_ns)&(time <= late_stop_ns)
    if not early.any() or not late.any():
        raise ValueError('empty measurement window')
    denominator = values[early].max(axis=0)
    numerator = values[late].max(axis=0)
    if np.any(denominator <= 0) or np.any(numerator <= 0):
        raise ValueError('zero peak denominator/numerator has no finite dB ratio')
    ratios = 20*np.log10(numerator/denominator)
    return {'median_dB': float(np.median(ratios)),
            'p05_p95_dB': np.percentile(ratios, [5, 95]).tolist(),
            'median_peak_time_ns': float(np.median(time[late][np.argmax(values[late], axis=0)])),
            'early_window_ns': [float(time[0]), early_stop_ns],
            'late_window_ns': [late_start_ns, late_stop_ns],
            'endpoint_guard_ns': float(time[-1]-late_stop_ns)}
