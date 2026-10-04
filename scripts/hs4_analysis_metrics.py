"""CPU diagnostics; no physical resolution or numerical-error certification."""
import numpy as np


def air_propagation_mask(frequency_hz, kx_rad_per_m, speed_m_per_s=299792458.0):
    f = np.asarray(frequency_hz, dtype=float)
    k = np.asarray(kx_rad_per_m, dtype=float)
    if (f.ndim != 1 or k.ndim != 1 or not f.size or not k.size
            or not np.isfinite(f).all() or not np.isfinite(k).all()
            or np.any(f <= 0) or not np.isfinite(speed_m_per_s)
            or speed_m_per_s <= 0):
        raise ValueError('finite positive frequencies/speed and finite 1D wave numbers required')
    return np.abs(k)[None, :] <= 2 * np.pi * f[:, None] / speed_m_per_s


def split_profile(profile, spacing_m, split_scale_m):
    p = np.asarray(profile, dtype=float)
    if (p.ndim != 1 or len(p) < 3 or not np.isfinite(p).all()
            or not np.isfinite(spacing_m) or spacing_m <= 0
            or not np.isfinite(split_scale_m) or split_scale_m <= 0):
        raise ValueError('finite profile and positive spatial scales required')
    frequency = np.fft.fftfreq(len(p), d=spacing_m)
    centered = p - p.mean()
    low = np.fft.ifft(np.fft.fft(centered)
                      * (np.abs(frequency) <= 1 / split_scale_m)).real
    return low, centered - low


def correlation(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a - a.mean(), b - b.mean()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return None if na == 0 or nb == 0 else float(np.dot(a / na, b / nb))


def matched_profile_metrics(ridge, truth, spacing_m, split_scale_m):
    r, t = np.asarray(ridge, dtype=float), np.asarray(truth, dtype=float)
    if r.shape != t.shape:
        raise ValueError('ridge and truth must share samples')
    rl, rh = split_profile(r, spacing_m, split_scale_m)
    tl, th = split_profile(t, spacing_m, split_scale_m)

    def component(a, b):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        return {'corr': correlation(a, b),
                'rmse_m': float(np.sqrt(np.mean((a - b) ** 2))),
                'relative_L2_error': None if nb == 0 else float(np.linalg.norm(a - b) / nb),
                'amplitude_norm_ratio': None if nb == 0 else float(na / nb),
                'signed_projection_gain': None if nb == 0 else float(np.dot(a, b) / nb ** 2),
                'numerical_error_budget': None,
                'interpretation': 'Diagnostic only; small components require a separate error budget.'}

    return {'corr_full': correlation(r, t),
            'rmse_m': float(np.sqrt(np.mean((r - t) ** 2))),
            'mean_bias_m': float(np.mean(r - t)),
            'lowpass_matched': component(rl, tl),
            'highpass_matched': component(rh, th),
            'diagnostic_split_scale_m': float(split_scale_m),
            'split_rule': 'Same periodic FFT mask on truth and ridge; not a measured resolution bound.',
            'window_boundary': 'Finite interval treated as periodic; split sensitive to interval edges.',
            'recoverability_certified': False}
