"""Gain-class operator effects table v0.1 (design step 3, user "可以" 2026-09-29).

Applies the frozen gain_control_cases_v0.1 contract attenuations to the t3 development
families and runs the pre-registered candidate set (catalogue identity/G2/G4 plus
mechanism references M1-M9 grounded in docs/research/2026-09-29_gain_operator_research.md)
as RECOVERY operators on the attenuated radargrams. Reference = original unattenuated
event (constructed_reference), same loading chain as the gain ladder step 1.

Full-scale convention (design doc §2.3 "simulated_full_scale", reading B): the full
scale FS = max|s| of the ORIGINAL signature is fixed per family; the operator input is
x_in = (attenuated array)/FS in [-1,1]-scale units; clipping is measured on the operator
OUTPUT against |y|<=1. Per-case input renormalization is NOT used: it would cancel
time_global attenuation in the normalized input and make the pre-registered nc/structure
contrasts (F1/F3) unobservable. Reference for the battery is x_ref = s/FS so that a
perfect recovery reads a=1, D_e=0.

Hard gates reuse the FROZEN tolerance contract (SHA asserted): tau_A=0.20, tau_D=0.95,
eps_Nb=1.5. The CLIP threshold is NOT gated here (no existing number may be reused);
clip_ratio is recorded as calibration data for step 4 (threshold proposal -> user).

Pre-registered expectations (research note §3-4, checked in the results doc, not here):
  F1 catalogue G2/G4 under-recover by tens of dB at deep levels;
  F2 M2/M3 (+40/+60 dB fixed exponential) fire the clip scale on strong early events;
  F3 nc_attenuated "uniform recovery" feasible beats "event-directed amplification";
  F4 AGC/envelope rows fail the a-gate despite excellent visual recovery;
  F5 AGC fires N_b on nc_untouched cells;
  F6 smart-interp gain (no-clip design) does NOT fire the clip scale.

Deterministic; r1/r2 byte-identity asserted. No solver; no test families {C5,C8};
no threshold is derived from candidate performance in this script.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from study_t3_damage_ladder import (MOTHERS, NC_WIN, load_bscan, fermat_times,
                                    event_mask)
from study_gain_ladder_v0_1 import factor_map, battery

CONTRACT = ROOT / 'configs/research/gain_control_cases_v0.1.json'
TOL = ROOT / 'configs/research/reward_tolerance_contract_v0.1.json'
OUT_R1 = ROOT / 'artifacts/research_checks/2026-09-29_gain_effects_r1.json'
OUT_R2 = ROOT / 'artifacts/research_checks/2026-09-29_gain_effects_r2.json'
CONTRACT_SHA256 = 'c027798904a4506bf3bd74d12fef00a4dc3df94b3cf6a938af170960b6689a18'
TOL_SHA256 = '9e7a0551e90e0bc8cd025f0d70b1bb434441ef980e2901bb61f4bfa82c023bf3'

# Cover-layer material of the t3 mothers (B2D-C3mS2X-BG-CO33-t01.in: "#material:
# 18.017 0.003 1 0 cover"); used ONLY to derive the SEC absorption rate alpha as a
# model-derived value (not a field claim; research note §2).
ER_COVER = 18.017
SIGMA_COVER = 0.003          # S/m
FC_HZ = 95.0e6               # SFCW band centre (20-170 MHz)
ETA0 = 376.73031346177066    # ohm
C_M_NS = 0.299792458

T_REF_NS = 10.0              # SEC/power-gain reference time (pre-registered)
T_GAIN_END_NS = 400.0        # gain-curve coordinate end = repo metric window
                             # (event ~150-200 ns, NC 250-400 ns; all metric
                             # windows of this repo live inside 0-400 ns; the
                             # reconstruction axis extends to ~3332 ns but every
                             # battery window is inside 400 ns). Curves reach
                             # their end value at 400 ns and HOLD it beyond.
G_MAX_LIN = 1000.0           # +60 dB industry-style Maximum Gain cap (pre-registered)
AGC_WIN_NS = 50.0            # RMS-AGC / envelope window (pre-registered)
AGC_FLOOR = 1e-4             # RMS/envelope floor in FS units (-80 dB), pre-registered
SMART_N_PTS = 16             # smart-interp control points (GPR-Slice style)
SMART_TARGET = 0.95          # maps local peak to 0.95 FS (never-clip design)

LEVELS = [-3.0, -6.0, -10.0, -20.0, -40.0]
STRUCTURES = ['event_only', 'time_global']
STRUCTURES_BAND = ['event_band']           # sampled levels only (research note §4)
BAND_LEVELS = [-10.0, -40.0]
NC_STATES = ['nc_untouched', 'nc_attenuated']

# Material-derived SEC absorption slope (two-way), dB/ns:
# low-loss check sigma/(omega*eps0*epsr) ~ 0.031 << 1 (documented in results doc).
_W = 2.0 * np.pi * FC_HZ
_EPS = ER_COVER * 8.8541878128e-12
LOW_LOSS_RATIO = SIGMA_COVER / (_W * _EPS)
ALPHA_DB_PER_M = 8.686 * 0.5 * SIGMA_COVER * ETA0 / np.sqrt(ER_COVER)
V_M_NS = C_M_NS / np.sqrt(ER_COVER)
A1_DB_PER_NS = ALPHA_DB_PER_M * V_M_NS          # two-way absorption slope
assert LOW_LOSS_RATIO < 0.1, 'low-loss formula precondition failed'


def load_contracts():
    assert hashlib.sha256(CONTRACT.read_bytes()).hexdigest() == CONTRACT_SHA256
    assert hashlib.sha256(TOL.read_bytes()).hexdigest() == TOL_SHA256
    tol = json.loads(TOL.read_text(encoding='utf-8'))['tolerances']
    return (float(tol['tau_A']['value']), float(tol['tau_D']['value']),
            float(tol['eps_Nb']['value']))


def sec_gain_db(t, a1):
    """SEC curve in dB on the metric-window coordinate: spherical spreading
    20log10(t/t_ref) + two-way absorption, capped at +60 dB Maximum Gain."""
    tc = np.minimum(t, T_GAIN_END_NS)
    g_db = np.zeros_like(t)
    m = tc > T_REF_NS
    g_db[m] = 20.0 * np.log10(tc[m] / T_REF_NS) + a1 * (tc[m] - T_REF_NS)
    return np.minimum(g_db, 60.0)   # Maximum Gain cap (+60 dB, pre-registered)


def fixed_exp_curve(t, end_gain):
    """Fixed shared exponential on the metric-window coordinate (end value
    reached at T_GAIN_END_NS and held beyond). Used by mechanism references."""
    u = np.minimum(t, T_GAIN_END_NS) / T_GAIN_END_NS
    return np.power(float(end_gain), u)


def fixed_exp_curve_full_axis(t, end_gain):
    """Frozen catalogue coordinate u=(t-t0)/(t_end-t0) over the FULL record
    (catalogue gain_coordinate 'u=i/(n_samples-1)'). Catalogue baselines only."""
    u = (t - t[0]) / (t[-1] - t[0])
    return np.power(float(end_gain), u)


def power_curve(t):
    tc = np.minimum(t, T_GAIN_END_NS)
    return np.maximum(tc / T_REF_NS, 1.0)


def _hilbert_env(x):
    """Per-trace Hilbert envelope via numpy FFT (axis 1 = time; deterministic)."""
    n = x.shape[1]
    X = np.fft.fft(x, axis=1)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = 1.0
        h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.abs(np.fft.ifft(X * h[None, :], axis=1))


def _smooth_same(x, half_w):
    """Centred moving mean along axis 1 (time), edges shrink to available samples."""
    csum = np.cumsum(x, axis=1)
    out = np.empty_like(x)
    n = x.shape[1]
    for i in range(n):
        lo = max(0, i - half_w)
        hi = min(n, i + half_w + 1)
        out[:, i] = (csum[:, hi - 1] - (csum[:, lo - 1] if lo > 0 else 0.0)) / (hi - lo)
    return out


def rms_agc(x, t, dt_ns):
    half = max(1, int(round(0.5 * AGC_WIN_NS / dt_ns)))
    rms = np.sqrt(_smooth_same(x ** 2, half))
    g = 1.0 / np.maximum(rms, AGC_FLOOR)
    return np.minimum(g, G_MAX_LIN)


def env_gain(x, t, dt_ns):
    half = max(1, int(round(0.5 * AGC_WIN_NS / dt_ns)))
    env = _smooth_same(_hilbert_env(x), half)
    g = 1.0 / np.maximum(env, AGC_FLOOR)
    return np.minimum(g, G_MAX_LIN)


def smart_interp_gain(x, t):
    """Shared per-line curve: local peak over all traces at 16 control times,
    linearly interpolated; maps each local peak to 0.95 FS (never-clip design)."""
    tc = np.linspace(t[0], t[-1], SMART_N_PTS)
    half_seg = 0.5 * (t[-1] - t[0]) / SMART_N_PTS
    g_pts = []
    for tc_i in tc:
        m = (t >= tc_i - half_seg) & (t <= tc_i + half_seg)
        local_max = float(np.max(np.abs(x[:, m]))) if m.any() else AGC_FLOOR
        g_pts.append(SMART_TARGET / max(local_max, AGC_FLOOR))
    g = np.interp(t, tc, g_pts)
    return np.minimum(g, G_MAX_LIN)


def apply_candidate(name, x, t, k, db, f):
    """Returns (gained_array, gain_diag_dict). x is [trace, sample] (time = axis 1).
    Shared curves broadcast along axis 0; per-trace adaptive gains are full arrays.
    f = the imposed attenuation factor map of this cell (needed by M1 oracle)."""
    dt_ns = float(t[1] - t[0])
    if name == 'identity':
        return x.copy(), {'kind': 'identity'}
    if name == 'cat_G2':
        return x * fixed_exp_curve_full_axis(t, 2.0)[None, :], {'kind': 'fixed_shared_exp_catalogue_coordinate', 'end_db': round(float(20 * np.log10(2.0)), 4)}
    if name == 'cat_G4':
        return x * fixed_exp_curve_full_axis(t, 4.0)[None, :], {'kind': 'fixed_shared_exp_catalogue_coordinate', 'end_db': round(float(20 * np.log10(4.0)), 4)}
    if name == 'M1_matched_inverse':
        # Oracle anchor: exact inverse of the imposed contract factor map
        # (contract-derived, NOT fitted to candidate behaviour). Recovers the
        # original signature exactly: a=1, D_e=0; NC ratio restored to 1 on
        # nc_attenuated cells and untouched NC stays at 1 on nc_untouched.
        g = 1.0 / f
        return x * g, {'kind': 'exact_inverse_of_imposed_factor_map'}
    if name == 'M2_exp_end_40dB':
        return x * fixed_exp_curve(t, 100.0)[None, :], {'kind': 'fixed_shared_exp', 'end_db': 40.0}
    if name == 'M3_exp_end_60dB':
        return x * fixed_exp_curve(t, 1000.0)[None, :], {'kind': 'fixed_shared_exp', 'end_db': 60.0}
    if name == 'M4_sec_model':
        g = np.power(10.0, sec_gain_db(t, A1_DB_PER_NS) / 20.0)
        return x * g[None, :], {'kind': 'sec_model', 'a1_db_per_ns': round(float(A1_DB_PER_NS), 6), 'end_db': round(float(20 * np.log10(g[-1])), 3)}
    if name == 'M5_sec_2a':
        g = np.power(10.0, sec_gain_db(t, 2.0 * A1_DB_PER_NS) / 20.0)
        return x * g[None, :], {'kind': 'sec_2a', 'a1_db_per_ns': round(float(2 * A1_DB_PER_NS), 6), 'end_db': round(float(20 * np.log10(g[-1])), 3)}
    if name == 'M6_power_t1':
        g = power_curve(t)
        return x * g[None, :], {'kind': 'power', 'alpha': 1.0, 'end_db': round(float(20 * np.log10(g[-1])), 3)}
    if name == 'M7_rms_agc_50ns':
        g = rms_agc(x, t, dt_ns)
        return g * x, {'kind': 'rms_agc', 'window_ns': AGC_WIN_NS, 'g_max_db': 60.0}
    if name == 'M8_env_gain_50ns':
        g = env_gain(x, t, dt_ns)
        return g * x, {'kind': 'env_gain', 'window_ns': AGC_WIN_NS, 'g_max_db': 60.0}
    if name == 'M9_smart_interp_16pt':
        g = smart_interp_gain(x, t)
        return x * g[None, :], {'kind': 'smart_interp', 'n_pts': SMART_N_PTS, 'target': SMART_TARGET}
    raise ValueError(name)


CANDIDATES = ['identity', 'cat_G2', 'cat_G4', 'M1_matched_inverse', 'M2_exp_end_40dB',
              'M3_exp_end_60dB', 'M4_sec_model', 'M5_sec_2a', 'M6_power_t1',
              'M7_rms_agc_50ns', 'M8_env_gain_50ns', 'M9_smart_interp_16pt']


def run_family(fam, mother, geo, tau_A, tau_D, eps_Nb):
    ifz = (lambda y: np.minimum(0.2 * y + 22.625, 30.0)) if geo == 'slope' \
        else (lambda y: np.zeros_like(np.asarray(y, float)) + 27.0)
    s, t = load_bscan(mother)
    t_ev = fermat_times(ifz)
    m_ev = event_mask(t, t_ev)
    m_nc = np.broadcast_to((t >= NC_WIN[0]) & (t <= NC_WIN[1]), s.shape)
    fs = float(np.max(np.abs(s)))
    x_ref = s / fs
    rows = []
    for db in LEVELS:
        for st in STRUCTURES + [b for b in STRUCTURES_BAND
                                if b == 'event_band' and db in BAND_LEVELS]:
            for nc in NC_STATES:
                k = 10.0 ** (db / 20.0)
                f = factor_map(t, t_ev, m_ev, k, st, nc)
                z0 = f * s
                x_in = z0 / fs
                # Construction sanity anchor (unrounded): identity must read a=k
                # exactly on every cell, as in the step-1 ladder.
                zs, ss = x_in[m_ev], x_ref[m_ev]
                a_in = float(np.dot(zs, ss) / np.dot(ss, ss))
                assert abs(a_in - k) <= 1e-12 * max(k, 1e-12), \
                    f'input anchor failed: {fam} {db} {st} {nc}'
                for cand in CANDIDATES:
                    y, diag = apply_candidate(cand, x_in, t, k, db, f)
                    if cand == 'M1_matched_inverse':
                        # Oracle anchor (unrounded): exact inverse must recover
                        # a=1, D_e=0 on every cell, as identity did in the ladder.
                        zs_o, ss_o = y[m_ev], x_ref[m_ev]
                        a_o = float(np.dot(zs_o, ss_o) / np.dot(ss_o, ss_o))
                        d_o = float(np.linalg.norm(zs_o - ss_o)
                                     / np.linalg.norm(ss_o))
                        assert abs(a_o - 1.0) <= 1e-12 and d_o <= 1e-12, \
                            f'oracle anchor failed: {fam} {db} {st} {nc}'
                    b = battery(y, x_ref, m_ev)
                    nc_ratio = float(np.mean(y[m_nc] ** 2) / np.mean(x_ref[m_nc] ** 2))
                    over = np.abs(y) > 1.0
                    clip_ratio = float(np.sum(y[over] ** 2) / np.sum(y ** 2)) \
                        if over.any() else 0.0
                    a_pass = abs(b['a'] - 1.0) <= tau_A
                    d_pass = b['D_e'] <= tau_D
                    nb_pass = nc_ratio <= eps_Nb
                    rows.append({
                        'family': fam, 'level_db': db, 'structure': st,
                        'nc_state': nc, 'candidate': cand, **b,
                        'rec_dB': round(20.0 * np.log10(max(b['a'], 1e-12) / k), 4),
                        'nc_ratio': round(nc_ratio, 9),
                        'clip_ratio': round(clip_ratio, 9),
                        'clip_frac': round(float(np.mean(over)), 9),
                        'a_pass': a_pass, 'D_pass': d_pass, 'nb_pass': nb_pass,
                        'feasible': bool(a_pass and d_pass and nb_pass),
                        'gain_diag': diag,
                    })
    print(fam, 'done:', len(rows), 'effect rows')
    return rows


def build_records():
    tau_A, tau_D, eps_Nb = load_contracts()
    records = []
    for fam, mother, geo in MOTHERS:
        records.extend(run_family(fam, mother, geo, tau_A, tau_D, eps_Nb))
    return records


def main():
    doc = {
        'schema': 'gain_effects/0.1',
        'date': '2026-09-29',
        'basis': 'gain class design step 3 (user "可以" to candidate set M1-M9, '
                 '2026-09-29); contract gain_control_cases_v0.1.json (SHA asserted) '
                 'drives every pre-registered value; tolerance gates reuse frozen '
                 'reward_tolerance_contract_v0.1.json (SHA asserted); candidate '
                 'mechanisms per 2026-09-29_gain_operator_research.md',
        'full_scale': 'FS = max|original signature| per family (fixed); operator '
                      'input = attenuated/FS; clip measured on output vs |y|<=1 '
                      '(simulated_full_scale, reading B; per-case input '
                      'renormalization would cancel time_global contrasts F1/F3)',
        'gain_curve_coordinate': 'mechanism-reference curves (M2-M6) are defined '
                                 'on the repo metric window [0, 400 ns] (u = '
                                 'min(t,400)/400; end value held beyond 400 ns); '
                                 'the catalogue baselines cat_G2/cat_G4 keep the '
                                 'FROZEN catalogue coordinate over the full '
                                 'record (gain_coordinate u=i/(n_samples-1)). '
                                 'Justification: every metric window of this repo '
                                 '(event +-10 ns ~150-200 ns, NC 250-400 ns) lives '
                                 'inside 0-400 ns; the reconstruction axis extends '
                                 'to ~3332 ns; anchoring curves to the metric '
                                 'window is an a-priori convention choice (repo '
                                 'norm), NOT fitted to candidate results. '
                                 'Recorded-axis note: with g(0)=1, monotone rising '
                                 'curves cannot clip the strong EARLY events '
                                 '(they sit at gain ~1); the pre-registered F2 '
                                 'mechanism wording (strong early events) is '
                                 'amended BEFORE analysis to: clipping under fixed '
                                 'curves arises where curve x local amplitude '
                                 'exceeds FS, primarily the amplified deep-time '
                                 'tail under end-value-hold extrapolation.',
        'reference': 'original unattenuated event, x_ref = s/FS (constructed_reference)',
        'gates': {'tau_A': 0.20, 'tau_D': 0.95, 'eps_Nb': 1.5,
                  'note': 'from frozen tolerance contract; clip NOT gated here'},
        'sec_material': {'er_cover': ER_COVER, 'sigma_cover': SIGMA_COVER,
                         'fc_hz': FC_HZ, 'low_loss_ratio': round(float(LOW_LOSS_RATIO), 6),
                         'alpha_db_per_m': round(float(ALPHA_DB_PER_M), 6),
                         'v_m_ns': round(float(V_M_NS), 6),
                         'a1_db_per_ns_two_way': round(float(A1_DB_PER_NS), 6),
                         'source': 'B2D-C3mS2X-BG-CO33-t01.in cover material, '
                                   'model-derived (not a field claim)'},
        'pre_registered_expectations': {
            'F1': 'cat_G2/cat_G4 under-recover deep levels by tens of dB '
                  '(catalogue full-axis coordinate delivers only a fraction of '
                  'end gain at event times)',
            'F2_amended': 'M2/M3 fire the clip scale (mechanism per '
                          'gain_curve_coordinate note: end-value-hold tail '
                          'amplification; rising curves do not clip early events)',
            'F3': 'nc_ratio separates uniform vs event-directed treatment: global '
                  'curves on nc_untouched cells push NC ratio over eps; M1 oracle '
                  'is feasible on all cells (contract-derived exact inverse)',
            'F4': 'M7/M8 fail a-gate despite large rec_dB (normalisation destroys '
                  'absolute amplitude; cap pins a)',
            'F5': 'M7 fires nb gate on nc_untouched cells (quiet-window noise '
                  'boosting)',
            'F6': 'M9 (0.95-FS smart-interp design target) fires little or no '
                  'clip relative to M2/M3'},
        'records': build_records(),
    }
    text1 = json.dumps(doc, ensure_ascii=False, indent=1) + '\n'
    rec2 = build_records()
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
