"""Official Waveform(user_func=...) callback for frozen Yee-lattice replay.

This is defined only at the source evaluation lattice. It never invents an
off-lattice interpolation or changes the saved incident field samples.
"""
import numpy as np


def lattice_waveform(values, dt, offset):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('Expected finite 1D replay history with at least two samples')
    if not np.isfinite(dt) or dt <= 0 or offset not in (0., .5*dt):
        raise ValueError('Expected positive dt and integer/half-step source offset')

    def sample(t):
        # Official Waveform checks the function at t=0 at build time. Electric
        # half-step source history starts later, so return zero before its start.
        if t < offset or t > offset+(len(values)-1)*dt+dt*1e-8:
            return 0.
        location = (t-offset)/dt
        index = round(location)
        if abs(location-index) > 1e-8:
            raise ValueError('Replay evaluated off its frozen Yee time lattice')
        if not 0 <= index < len(values):
            return 0.
        return float(values[index])

    return sample
