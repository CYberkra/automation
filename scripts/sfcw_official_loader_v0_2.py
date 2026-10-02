"""Official-carrier signed SFCW loader (v0.2, dev-side fix for the 75 MHz shift).

Bug being fixed (see docs/research/2026-10-02_model_design_review.md §2):
legacy loaders (study_t3_damage_ladder.load_bscan,
run_reward_protocol_b2_pilot_v0_1.load_bscan,
freeze_s1s3_reference_window_v0_1.load_trace) rebuilt the signed trace as

    |env| * cos(2*pi*FC*t + angle(env))   with FC = 95 MHz (band centre)

which equals Re(complex_envelope * exp(2j*pi*95MHz*t)). The official gprMax
chain (gprMax.toolboxes.SFCW.processing.reconstruct_time_response, SHA-256
adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b) modulates
the complex envelope at frequencies[0] = 20 MHz and returns
real_bandpass = 2*Re(complex_envelope * exp(2j*pi*20MHz*t)). The legacy form
therefore (a) shifts the physical band up by +75 MHz to 95-245 MHz and
(b) carries half the official amplitude (amplitude convention may be
re-declared; the frequency shift may not).

This module loads each trace ONCE through the official chain (identical
preprocessing to the legacy loaders: 1200 ns truncation, same tail taper,
same 501-point 20-170 MHz grid, zero_pad_factor=8, Hann) and exposes both
representations derived from the SAME complex envelope, so old-vs-new diffs
never conflate loader drift with the carrier fix.

Scope: dev-side comparison only. Frozen scripts/outputs stay untouched;
no solver runs; no test-family (C5/C8) data.
"""

from dataclasses import replace
import hashlib
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import (
    load_source, load_receiver, direct_frequency_response,
    reconstruct_time_response)
import gprMax.toolboxes.SFCW.processing as _processing

OFFICIAL_PROCESSING_SHA256 = 'adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b'


def verify_official_runtime():
    digest = hashlib.sha256(Path(_processing.__file__).read_bytes()).hexdigest()
    if digest != OFFICIAL_PROCESSING_SHA256:
        raise ValueError('official SFCW processing SHA-256 differs from reviewed version')

ROOT = Path(__file__).resolve().parents[1]
SIMS = ROOT / 'artifacts/simulations'

FREQ = np.linspace(20e6, 170e6, 501)
F0_OFFICIAL = float(FREQ[0])                 # 20 MHz — official carrier
FC_LEGACY = 0.5 * (FREQ[0] + FREQ[-1])       # 95 MHz — legacy (wrong) carrier


def trace_time_response(h5_path):
    """One archived CO trace through the official SFCW chain (Hann).

    Preprocessing is ladder-exact: 1200 ns window, same tail-taper fraction
    formula, 501 x 0.3 MHz grid, zero_pad_factor=8, Hann window.
    """
    verify_official_runtime()
    h5_path = Path(h5_path)
    with h5py.File(h5_path, 'r') as h:
        dt = float(h.attrs['dt'])
        items = list(h['rxs'].items())
        if len(items) != 1:
            raise ValueError('carrier loader requires exactly one receiver')
        rx_name = items[0][1].attrs['Name']
        raw = items[0][1]['Ex'][:]
    src = load_source(h5_path)
    rxt = load_receiver(h5_path, receiver_path='name:' + rx_name,
                        component='Ex')
    n = min(len(raw), int(np.floor(1200e-9 / dt)) + 1)
    taper = (round(200.0e-9 / dt) - .25) / n
    rx = replace(rxt, samples=raw[:n])
    r = direct_frequency_response(src, rx, FREQ, tail_taper_fraction=taper)
    return reconstruct_time_response(r, zero_pad_factor=8, window='hann')


def signed_official(tr):
    """Official signed bandpass: 2*Re(env * exp(2j*pi*20MHz*t))."""
    return np.asarray(tr.real_bandpass, dtype=np.float64)


def signed_legacy(tr):
    """Legacy (buggy) signed trace rebuilt from the SAME envelope.

    Kept for dev-side diffs only; never use for new evidence.
    """
    env_c = np.asarray(tr.complex_envelope, dtype=np.complex128)
    t = np.asarray(tr.time, dtype=np.float64)
    return np.real(env_c * np.exp(2j * np.pi * FC_LEGACY * t))


def time_ns(tr):
    return np.asarray(tr.time, dtype=np.float64) * 1e9


def dominant_freq_mhz(sig, dt_s):
    """Dominant Fourier frequency (MHz) of a real signed trace."""
    sp = np.abs(np.fft.rfft(np.asarray(sig, dtype=np.float64)))
    fr = np.fft.rfftfreq(sig.size, d=dt_s)
    return float(fr[int(np.argmax(sp))] * 1e-6)


def load_bscan_both(mother, date, n_traces=33, *, input_records=None,
                    audited_inputs=None):
    """Load one mother model once; return (t_ns, sig_official, sig_legacy).

    Both signed arrays derive from the same per-trace complex envelopes.
    """
    sig_off, sig_leg, t = [], [], None
    for k in range(n_traces):
        rid = f'{mother}-CO33-t{k + 1:02d}'
        h5 = SIMS / f'{date}_{rid}' / f'{rid}.h5'
        if input_records is not None:
            from sfcw_carrierfix_acceptance import verify_input
            record = input_records[rid]
            h5 = SIMS / record['source_dir'] / f'{rid}.h5'
            checked = verify_input(h5, record)
            if audited_inputs is not None:
                audited_inputs.append(checked)
        tr = trace_time_response(h5)
        if input_records is not None:
            verify_input(h5, record)  # Reject a concurrent change during loading.
        sig_off.append(signed_official(tr))
        sig_leg.append(signed_legacy(tr))
        current_t = time_ns(tr)
        if t is not None and not np.array_equal(t, current_t):
            raise ValueError('reconstructed axes differ between traces')
        t = current_t
    return t, np.stack(sig_off), np.stack(sig_leg)
