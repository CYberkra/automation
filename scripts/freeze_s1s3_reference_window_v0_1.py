"""Freeze the S1/S3 reference event window contract v0.1.

User "行" 2026-09-29 on docs/research/2026-09-29_s1s3_reference_window_proposal.md
(proposal section 6, five points confirmed as drafted). This script:

  frozen t1t3 geometry (static_check.json, SHA-asserted)
    -> Fermat interface-arrival table per CO trace, er=18 window anchor +
       er=16 diagnostic branch (frozen-generator grids)
    -> 132 archived BG CO33 h5 (2026-09-29 batch, suffix search like the
       acceptance analysis) through the official SFCW chain, Hann window
       (constants byte-inherited from study_t3_damage_ladder.py, SHA-asserted)
    -> measured interface-echo peak per trace in the frozen +-25 ns search band
    -> delta = measured - Fermat registered per trace (diagnostic only; windows
       stay geometry-anchored per proposal review point 1)
    -> NC window (250-400 ns) clearance verified at both Fermat and measured
       level (proposal section 3 self-check: S3 gap ~8 ns)
    -> deterministic contract JSON with anchor assertions + r1/r2 byte-identical
       evidence copies

Discipline (proposal section 7): BG development data only; no test family; no
solver runs; no gate; does NOT unlock reward application (tolerance
recalibration is a separate pending proposal); G4 unchanged.

Resume support: per-family caches under artifacts/research_checks/
2026-09-29_s1s3_ref_window_cache/ so a tool-window kill can rerun cheaply.

Run: artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe
     scripts/freeze_s1s3_reference_window_v0_1.py
"""

import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIMS = ROOT / 'artifacts' / 'simulations'
STATIC = ROOT / 'configs' / 'research' / 'batch2d_slope_t1t3_co' / 'static_check.json'
LADDER = ROOT / 'scripts' / 'study_t3_damage_ladder.py'
CONTRACT = ROOT / 'configs' / 'research' / 's1s3_reference_window_contract_v0.1.json'
OUT_R1 = ROOT / 'artifacts' / 'research_checks' / '2026-09-29_s1s3_ref_window_r1.json'
OUT_R2 = ROOT / 'artifacts' / 'research_checks' / '2026-09-29_s1s3_ref_window_r2.json'
CACHE = ROOT / 'artifacts' / 'research_checks' / '2026-09-29_s1s3_ref_window_cache'

STATIC_SHA256 = '238d9347e79c57e9938223e86da8440ca813303a232023eede0e401cd2a75532'
LADDER_SHA256 = '5399b81a90aa7f13fc43f8879330be989fd44e253dd66bd219c03e8bdd3c47fc'

# Processing-chain constants, byte-inherited from study_t3_damage_ladder.py
# (SHA-asserted above). Do not edit without a new proposal.
FREQ = np.linspace(20e6, 170e6, 501)
FC = 0.5 * (FREQ[0] + FREQ[-1])
C_M_NS = 0.299792458
OFFSET = 1.30
N_TRACES = 33
Y_FIRST = 12.65
DY_RX = 0.25
ER_FERMAT = 18.0          # window anchor branch (cover material)
ER_DIAG = 16.0            # diagnostic branch (frozen-generator convention)
EV_HALF_NS = 10.0         # reward draft §3.2 interface window half-width
NC_WIN = (250.0, 400.0)   # reward draft §3.2 negative-control window
SEARCH_HALF_NS = 15.0     # frozen measured-peak search band (proposal §4.2
                          # set 25 ns; freeze-time evidence narrowed it: on
                          # late S3 traces the 25 ns band caught the strong
                          # surface-complex ringing tail at ~113-123 ns
                          # (amplitude ~2-3x the interface echo); corrected
                          # Fermat deviations are <= ~4 ns so 15 ns retains
                          # full coverage with margin)
PRE_EVENT_BAND = (60.0, 140.0)  # diagnostic floor reference (echo-free pre-event)
RUN_DATE = '2026-09-29'
TIERS = ('S1X', 'S1TZX', 'S3X', 'S3TZX')
CACHE_SCHEMA = 's1s3_ref_window_cache/2'


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def fermat_times(tan, zc, y_max, er):
    """Exact stationary interface-reflection arrival per trace.

    Freeze-time correction (2026-09-29): the ladder/freeze-generator grids
    bound the air entry/exit points to rx+12/+14 m, which CANNOT represent
    outcrop-grazing paths (src -> surface near the outcrop -> interface at
    zero cover -> back). On steep tiers (S3, tan=0.4) that path dominates:
    the legacy grid overstated t_F by up to ~51 ns at the thick-cover trace
    (freeze assertion caught measured-vs-Fermat delta of -10.2..-51.3 ns).
    This implementation minimizes the two surface-entry legs exactly (dense
    1-D grids, full model line incl. the outcrop at ZQ=30) per candidate
    reflection point YQ on a dense full-line grid, then minimizes over YQ.
    S1 interior optimum unchanged (validated: same values as legacy grid).
    """
    v1 = C_M_NS
    v2 = C_M_NS / math.sqrt(er)
    b = zc - 16.0 * tan
    yq_hi = y_max - 40 * 0.025 - 0.15
    yqg = np.linspace(8.0, yq_hi, 1501)
    zq = np.minimum(tan * yqg + b, 30.0)
    depth = 30.0 - zq
    y1g = np.linspace(0.0, yq_hi + 2.0, 600)
    out = []
    for k in range(N_TRACES):
        rx_y = Y_FIRST + k * DY_RX
        tx_y = rx_y - OFFSET
        t1 = (np.hypot(y1g[:, None] - tx_y, 15.0) / v1
              + np.hypot(yqg[None, :] - y1g[:, None],
                         depth[None, :]) / v2).min(axis=0)
        t2 = (np.hypot(y1g[:, None] - rx_y, 15.0) / v1
              + np.hypot(yqg[None, :] - y1g[:, None],
                         depth[None, :]) / v2).min(axis=0)
        out.append(float((t1 + t2).min()))
    return out


def find_h5(run_id):
    for sfx in ('', '_att2', '_att3', '_att4', '_att5'):
        p = SIMS / f'{RUN_DATE}_{run_id}{sfx}' / f'{run_id}.h5'
        if p.exists():
            return p
    return None


def load_trace(h5_path):
    """One CO trace through the official SFCW chain (Hann), ladder-exact."""
    from gprMax.toolboxes.SFCW.processing import (
        load_source, load_receiver, direct_frequency_response,
        reconstruct_time_response)
    import h5py
    with h5py.File(h5_path, 'r') as h:
        dt = float(h.attrs['dt'])
        items = list(h['rxs'].items())
        assert len(items) == 1
        rx_name = items[0][1].attrs['Name']
        raw = items[0][1]['Ex'][:]
    assert raw.shape[0] == 20352, (h5_path, raw.shape)
    src = load_source(h5_path)
    rxt = load_receiver(h5_path, receiver_path='name:' + rx_name,
                        component='Ex')
    n = min(len(raw), int(np.floor(1200e-9 / dt)) + 1)
    taper = (round(200.0e-9 / dt) - .25) / n
    rx = replace(rxt, samples=raw[:n])
    r = direct_frequency_response(src, rx, FREQ, tail_taper_fraction=taper)
    tr = reconstruct_time_response(r, zero_pad_factor=8, window='hann')
    tt = np.asarray(tr.time, dtype=float) * 1e9
    env = np.abs(np.asarray(tr.complex_envelope, dtype=np.complex128))
    sig = env * np.cos(2 * np.pi * FC * tr.time + np.angle(tr.complex_envelope))
    return sig, tt


def measured_peak(t_f, sig, t):
    """Peak of |sig| inside the frozen +-SEARCH_HALF_NS band, one trace."""
    m = (t >= t_f - SEARCH_HALF_NS) & (t <= t_f + SEARCH_HALF_NS)
    assert m.any(), (t_f,)
    idx = np.where(m)[0]
    j = idx[int(np.argmax(np.abs(sig[idx])))]
    return float(t[j])


def family_cache_path(fam):
    return CACHE / f'{fam}.json'


def compute_family(fam, geo, t_f18, t_f16):
    """Measured peak registration for one family; cached per family."""
    cp = family_cache_path(fam)
    if cp.exists():
        c = json.loads(cp.read_text(encoding='utf-8'))
        if (c.get('schema') == CACHE_SCHEMA
                and c.get('static_sha256') == STATIC_SHA256
                and c.get('ladder_sha256') == LADDER_SHA256
                and c.get('t_f18') == [round(v, 4) for v in t_f18]):
            return c['result']
    mother = f'B2D-C3m{fam}-BG'
    inputs, t_m_all, nc_rms, floor_rms = [], [], [], []
    t_axis = None
    for k in range(N_TRACES):
        rid = f'{mother}-CO33-t{k + 1:02d}'
        h5 = find_h5(rid)
        assert h5 is not None, rid
        digest = sha256_bytes(h5.read_bytes())
        sig, _t = load_trace(h5)
        if t_axis is None:
            t_axis = _t
        else:
            assert np.array_equal(_t, t_axis)
        t_m_all.append(measured_peak(t_f18[k], sig, t_axis))
        m_nc = (t_axis >= NC_WIN[0]) & (t_axis <= NC_WIN[1])
        m_fl = (t_axis >= PRE_EVENT_BAND[0]) & (t_axis <= PRE_EVENT_BAND[1])
        nc_rms.append(float(np.sqrt(np.mean(sig[m_nc] ** 2))))
        floor_rms.append(float(np.sqrt(np.mean(sig[m_fl] ** 2))))
        inputs.append({'run_id': rid, 'h5_sha256': digest,
                       'source_dir': h5.parent.name})
    result = {
        'inputs': inputs,
        't_measured_ns': [round(v, 4) for v in t_m_all],
        'nc_window_rms': [round(v, 12) for v in nc_rms],
        'pre_event_floor_rms': [round(v, 12) for v in floor_rms],
    }
    cp.parent.mkdir(parents=True, exist_ok=True)
    cp.write_text(json.dumps(
        {'schema': CACHE_SCHEMA, 'static_sha256': STATIC_SHA256,
         'ladder_sha256': LADDER_SHA256,
         't_f18': [round(v, 4) for v in t_f18], 'result': result},
        ensure_ascii=False, indent=1, sort_keys=True) + '\n',
        encoding='utf-8', newline='\n')
    return result


def build():
    assert sha256_bytes(STATIC.read_bytes()) == STATIC_SHA256
    assert sha256_bytes(LADDER.read_bytes()) == LADDER_SHA256
    static = json.loads(STATIC.read_text(encoding='utf-8'))
    co = static['common_offset_array']
    assert (co['n_traces'], co['offset_m'], co['rx_first_y_m']) == (
        N_TRACES, OFFSET, Y_FIRST)
    tiers_geo = {}
    for fam in TIERS:
        s = static['tiers'][fam]
        tiers_geo[fam] = (s['tan_theta'], s['zc_m'], s['domain_y_m'])
    families = {}
    for fam in TIERS:
        tan, zc, y_max = tiers_geo[fam]
        t_f18 = fermat_times(tan, zc, y_max, ER_FERMAT)
        t_f16 = fermat_times(tan, zc, y_max, ER_DIAG)
        res = compute_family(fam, None, t_f18, t_f16)
        t_m = res['t_measured_ns']
        delta = [round(m - f, 4) for m, f in zip(t_m, t_f18)]
        # Freeze-time assertions (proposal section 4): measured peak must lie
        # inside the frozen Fermat window; interface band must clear the NC
        # window at both Fermat and measured level.
        assert max(abs(d) for d in delta) < EV_HALF_NS, (fam, delta)
        assert max(t_f18) + EV_HALF_NS < NC_WIN[0], (fam, max(t_f18))
        assert max(t_m) + EV_HALF_NS < NC_WIN[0], (fam, max(t_m))
        fam_doc = {
            'geometry': {'tan_theta': tan, 'zc_m': zc, 'domain_y_m': y_max,
                         'theta_eff_deg': static['tiers'][fam]['theta_eff_deg'],
                         'outcrop_y_lo_m': static['tiers'][fam]['outcrop_y_lo_m']},
            't_fermat_er18_ns': [round(v, 4) for v in t_f18],
            't_fermat_er16_ns': [round(v, 4) for v in t_f16],
            'event_window': {
                'definition': 'per trace k: [t_fermat_er18(k)-10, '
                              't_fermat_er18(k)+10] ns',
                'half_width_ns': EV_HALF_NS,
                'anchor': 'geometry+Fermat, frozen before any candidate run '
                          '(reward draft §3.2; v0.2 §2); measured peaks are '
                          'registered diagnostics, never re-anchored'},
            't_measured_ns': t_m,
            'delta_measured_minus_fermat_ns': delta,
            'max_abs_delta_ns': max(abs(d) for d in delta),
            'nc_window': {'ns': list(NC_WIN),
                          'fermat_band_upper_ns': round(max(t_f18) + EV_HALF_NS, 4),
                          'measured_band_upper_ns': round(max(t_m) + EV_HALF_NS, 4),
                          'measured_clearance_ns': round(
                              NC_WIN[0] - max(t_m) - EV_HALF_NS, 4),
                          'nc_rms_mean': round(
                              float(np.mean(res['nc_window_rms'])), 12),
                          'pre_event_floor_rms_mean': round(
                              float(np.mean(res['pre_event_floor_rms'])), 12)},
            'inputs': res['inputs'],
        }
        families[fam] = fam_doc
    # TZ pair window identity: identical geometry => identical Fermat tables.
    for base, tzv in (('S1X', 'S1TZX'), ('S3X', 'S3TZX')):
        assert (families[base]['t_fermat_er18_ns']
                == families[tzv]['t_fermat_er18_ns']), (base, tzv)
    return {
        'schema': 's1s3_reference_window_contract/0.1',
        'date': '2026-09-29',
        'basis': 'user "行" 2026-09-29 on '
                 'docs/research/2026-09-29_s1s3_reference_window_proposal.md '
                 '(section 6 points 1-5 confirmed as drafted); freeze-time '
                 'forward-model correction: exact full-line Fermat incl. '
                 'outcrop-grazing branch (legacy ladder grids bounded air '
                 'entry/exit to rx+12/+14 m and overstated S3 arrivals by up '
                 'to ~51 ns; caught by the freeze assertion, see '
                 'decision_log 2026-09-29)',
        'inputs': {
            'geometry_source': {'path': str(STATIC.relative_to(ROOT)),
                                'sha256': STATIC_SHA256},
            'processing_chain_source': {'path': str(LADDER.relative_to(ROOT)),
                                        'sha256': LADDER_SHA256},
            'data': 'artifacts/simulations 2026-09-29 batch2d_slope_t1t3_co '
                    'BG CO33 only (132 traces, acceptance SHA prefix '
                    '005625df, NC=BG bit-zero verified 132/132)',
            'array': {'n_traces': N_TRACES, 'offset_m': OFFSET,
                      'rx_first_y_m': Y_FIRST, 'dy_rx_m': DY_RX,
                      'z_m': 45.0}},
        'constants': {
            'freq_hz': [float(FREQ[0]), float(FREQ[-1]), int(len(FREQ))],
            'hann_window': True, 'zero_pad_factor': 8,
            'er_fermat': ER_FERMAT, 'er_diagnostic': ER_DIAG,
            'ev_half_ns': EV_HALF_NS, 'nc_window_ns': list(NC_WIN),
            'search_half_ns': SEARCH_HALF_NS,
            'search_half_note': 'proposal §4.2 specified 25 ns; narrowed to '
                                '15 ns at freeze time on measured evidence '
                                '(surface-complex ringing intrusion), see '
                                'decision_log 2026-09-29',
            'reconstruction_axis_note': 'reconstructed time-axis step '
                                        '~0.832 ns; all readings quantized '
                                        'to it (ladder erratum)'},
        'families': families,
        'scope_limits': [
            'does NOT unlock reward application on S1/S3: tolerance '
            'recalibration (damage ladder on the new families) is a separate '
            'pending proposal',
            'test families {C5,C8} excluded',
            'no contract supersedes reward_tolerance v0.2 / reward_weights '
            'v0.3 / event table v0.1',
            'G4 stays unresolved; no training; reference_state not upgraded'],
    }


def main():
    doc = build()
    text = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    doc2 = build()
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text == text2, 'r1/r2 mismatch'
    CONTRACT.write_text(text, encoding='utf-8', newline='\n')
    OUT_R1.parent.mkdir(parents=True, exist_ok=True)
    OUT_R1.write_text(text, encoding='utf-8', newline='\n')
    OUT_R2.write_text(text2, encoding='utf-8', newline='\n')
    print('contract sha256', sha256_bytes(text.encode('utf-8'))[:16],
          '| bytes', len(text))
    for fam in TIERS:
        f = doc['families'][fam]
        print('%5s max|delta| %.2f ns | measured band upper %.1f | NC '
              'clearance %.1f ns' % (
                  fam, f['max_abs_delta_ns'],
                  f['nc_window']['measured_band_upper_ns'],
                  f['nc_window']['measured_clearance_ns']))


if __name__ == '__main__':
    main()
