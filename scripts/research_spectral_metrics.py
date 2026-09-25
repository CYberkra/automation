"""Existing spectral comparison definitions, independent of solver imports."""
import numpy as np


def metrics(h, ref):
    error = abs(h-ref)/abs(ref)
    return dict(max_relative_complex_error=float(max(error)),
                median_relative_complex_error=float(np.median(error)),
                max_abs_amplitude_error_dB=float(max(abs(20*np.log10(abs(h)/abs(ref))))),
                max_abs_phase_error_deg=float(max(abs(np.angle(h/ref, deg=True)))))
