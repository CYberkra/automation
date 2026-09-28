"""Damage ladder on the t3 BG CO B-scans: metric-response calibration for reward draft §4.

Applies known, metered damages to the archived batch2d_slope_t3_co signed radargrams
(official SFCW chain, Hann window) and records the v0.2 metric battery (D/a/A/H/rho
on the interface event window, N_b ratio on the negative-control window, arrival
drift diagnostic). The ladder is the calibration input for task tolerances (step 3);
it uses development-family data only, no test family, no solver runs.

Damages (event = interface band, Fermat line er=18 +/- 10 ns per trace):
  amp_db      scale event band by -1/-3/-10 dB
  polarity    negate event band
  shift_smp   roll event band by 1/2/4 samples of the reconstructed time axis
              (axis step ~0.832 ns => 0.83/1.66/3.33 ns; the reconstruction grid
              IS the device data product's resolution — sub-sample shifts are not
              representable; supersedes the draft's raw-dt wording)
  delete      zero the event band (interface deletion)
  nc_amp      add seeded gaussian noise at +6/+20 dB over NC-window RMS (NC window
              250-400 ns, per reward draft §3.2)
  identity    zero damage (sanity: D=0, A=0, H=0, rho=1)

Deterministic (fixed seeds); run writes r1 then reruns in-process and asserts the
r2 records are identical, saving both.
"""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import (
    load_source, load_receiver, direct_frequency_response, reconstruct_time_response)

ROOT = Path(__file__).resolve().parents[1]
SIMS = ROOT / 'artifacts/simulations'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-28_t3_damage_ladder_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-28_t3_damage_ladder_r2.json'

FREQ = np.linspace(20e6, 170e6, 501)
FC = 0.5 * (FREQ[0] + FREQ[-1])
C_M_NS = 0.299792458
OFFSET = 1.30
N_TRACES = 33
Y_FIRST = 12.65
DY_RX = 0.25
ER_FERMAT = 18.0
EV_HALF_NS = 10.0          # reward draft §3.2: interface window half-width
NC_WIN = (250.0, 400.0)    # reward draft §3.2: negative-control window
AMP_DB = [-1.0, -3.0, -10.0]
SHIFT_SMP = [1, 2, 4]
NC_AMP_DB = [6.0, 20.0]
SEED = 20260928

MOTHERS = [
    ('S2X', 'B2D-C3mS2X-BG', 'slope'),
    ('S2TZX', 'B2D-C3mS2TZX-BG', 'slope'),
    ('C3mX', 'B2D-C3mX-BG', 'flat'),
]


def fermat_times(ifz):
    v2 = C_M_NS / np.sqrt(ER_FERMAT)
    out = []
    for k in range(N_TRACES):
        rx_y = Y_FIRST + k * DY_RX
        tx = rx_y - OFFSET
        y1g = np.linspace(rx_y - 8, rx_y + 14, 50)
        yqg = np.linspace(12.0, 36.9, 250)
        y2g = np.linspace(rx_y - 6, rx_y + 16, 50)
        Y1, YQ, Y2 = np.meshgrid(y1g, yqg, y2g, indexing='ij')
        ZQ = ifz(YQ)
        t = ((np.hypot(Y1 - tx, 15.0) + np.hypot(Y2 - rx_y, 15.0)) / C_M_NS
             + (np.hypot(YQ - Y1, 30.0 - ZQ) + np.hypot(Y2 - YQ, 30.0 - ZQ)) / v2)
        out.append(float(t.min()))
    return np.array(out)


def load_bscan(mother):
    sigs, t = [], None
    for k in range(N_TRACES):
        rid = f'{mother}-CO33-t{k + 1:02d}'
        h5 = SIMS / f'2026-09-28_{rid}' / f'{rid}.h5'
        with h5py.File(h5, 'r') as h:
            dt = float(h.attrs['dt'])
            items = list(h['rxs'].items())
            assert len(items) == 1
            rx_name = items[0][1].attrs['Name']
            raw = items[0][1]['Ex'][:]
        src = load_source(h5)
        rxt = load_receiver(h5, receiver_path='name:' + rx_name, component='Ex')
        n = min(len(raw), int(np.floor(1200e-9 / dt)) + 1)
        taper = (round(200.0e-9 / dt) - .25) / n
        rx = replace(rxt, samples=raw[:n])
        r = direct_frequency_response(src, rx, FREQ, tail_taper_fraction=taper)
        tr = reconstruct_time_response(r, zero_pad_factor=8, window='hann')
        tt = np.asarray(tr.time, dtype=float) * 1e9
        env = np.abs(np.asarray(tr.complex_envelope, dtype=np.complex128))
        sig = env * np.cos(2 * np.pi * FC * tr.time + np.angle(tr.complex_envelope))
        sigs.append(sig)
        t = tt
    return np.stack(sigs), t


def event_mask(t, t_ev):
    m = np.zeros((N_TRACES, len(t)), dtype=bool)
    for i in range(N_TRACES):
        m[i] = (t >= t_ev[i] - EV_HALF_NS) & (t <= t_ev[i] + EV_HALF_NS)
    return m


def metrics(z, s, m_ev, m_nc, t, t_ev):
    """v0.2 battery on the event window + NC ratio + arrival drift diagnostic."""
    zs, ss = z[m_ev], s[m_ev]
    denom = float(np.linalg.norm(ss))
    assert denom > 0
    d = float(np.linalg.norm(zs - ss) / denom)
    a = float(np.dot(zs, ss) / np.dot(ss, ss))
    h = float(np.linalg.norm(zs - a * ss) / denom)
    zn = float(np.linalg.norm(zs))
    rho = float(np.dot(zs, ss) / (zn * denom)) if zn > 0 else None
    nb_in = float(np.mean(s[m_nc] ** 2))
    nb_out = float(np.mean(z[m_nc] ** 2))
    drift = []
    for i in range(N_TRACES):
        seg_z = z[i][m_ev[i]]
        drift.append(round(float(t[m_ev[i]][int(np.argmax(np.abs(seg_z)))] - t_ev[i]), 4))
    return {
        'D': round(d, 6), 'a': round(a, 6), 'A': round(abs(a - 1.0), 6),
        'H': round(h, 6), 'rho': (round(rho, 6) if rho is not None else None),
        'Nb_ratio': round(nb_out / nb_in, 6) if nb_in > 0 else None,
        'arrival_drift_ns_median': round(float(np.median(drift)), 4),
        'arrival_drift_ns_max_abs': round(float(np.max(np.abs(drift))), 4),
    }


def apply_damage(s, m_ev, m_nc, t, kind, level, rng):
    z = s.copy()
    if kind == 'identity':
        return z
    if kind == 'amp_db':
        z[m_ev] *= 10.0 ** (level / 20.0)
        return z
    if kind == 'polarity':
        z[m_ev] *= -1.0
        return z
    if kind == 'delete':
        z[m_ev] = 0.0
        return z
    if kind == 'shift_smp':
        k = int(level)
        assert k >= 1
        for i in range(N_TRACES):
            band = z[i][m_ev[i]]
            z[i][m_ev[i]] = np.roll(band, k) if k < len(band) else 0.0
            if k < len(band):
                z[i][m_ev[i]][:k] = 0.0
        return z
    if kind == 'nc_amp':
        rms = float(np.sqrt(np.mean(s[m_nc] ** 2)))
        noise = rng.standard_normal(s.shape) * (rms * 10.0 ** (level / 20.0))
        z = z + noise * m_nc
        return z
    raise ValueError(kind)


def build_records():
    records = []
    for fam, mother, geo in MOTHERS:
        ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
            else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
        s, t = load_bscan(mother)
        t_ev = fermat_times(ifz)
        m_ev = event_mask(t, t_ev)
        m_nc = (t >= NC_WIN[0]) & (t <= NC_WIN[1])
        m_nc = np.broadcast_to(m_nc, s.shape)
        rng = np.random.default_rng(SEED)
        damages = [('identity', 0.0)]
        damages += [('amp_db', v) for v in AMP_DB]
        damages += [('polarity', 0.0)]
        damages += [('shift_smp', v) for v in SHIFT_SMP]
        damages += [('delete', 0.0)]
        damages += [('nc_amp', v) for v in NC_AMP_DB]
        for kind, level in damages:
            z = apply_damage(s, m_ev, m_nc, t, kind, level, rng)
            rec = {'family': fam, 'mother': mother, 'damage': kind, 'level': level,
                   **metrics(z, s, m_ev, m_nc, t, t_ev)}
            records.append(rec)
        print(fam, 'done:', len(damages), 'damage instances')
    return records


def main():
    doc = {
        'schema': 't3_damage_ladder/1',
        'date': '2026-09-28',
        'basis': 'reward metrics draft v0.1 §4 step 2 (damage ladder), user "确认，启动吧"',
        'data': 'batch2d_slope_t3_co archived h5 (BG only); official SFCW chain, Hann window, '
                'signed radargram at 95 MHz band centre',
        'event_window': {'construction': 'Fermat line (er=18, smooth t3 line clamped at '
                         'outcrop; flat family z=27) +/- %.1f ns per trace' % EV_HALF_NS},
        'reconstruction_axis_note': 'reconstructed time axis step ~0.832 ns (501 pts x 0.3 MHz '
                                    'grid); all peak-based readouts are quantised to it; '
                                    'shift_smp levels are integer samples of this axis',
        'nc_window_ns': list(NC_WIN),
        'seed': SEED,
        'note': 'metric-response calibration only; thresholds are NOT set here (draft §4 '
                'step 3); development family; no solver runs; no test-family data',
        'records': build_records(),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    rec2 = build_records()  # in-process rerun for byte-identity (same seed path)
    doc2 = dict(doc, records=rec2)
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1) + '\n'
    assert text1 == text2, 'r1/r2 mismatch'
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text1, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('records:', len(doc['records']),
          '| sha256', hashlib.sha256(text1.encode('utf-8')).hexdigest()[:16])
    print('wrote', OUT_R1.relative_to(ROOT), 'and', OUT_R2.relative_to(ROOT))


if __name__ == '__main__':
    main()
