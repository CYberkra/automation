"""Material-scan analysis: echo/direct, decay slope, relief fidelity vs measured Line9.

Reads the completed 2026-10-05_hs4_material_scan_r1 capsule (20 groups x 13 traces,
8 m altitude) and the measured Line9(36) CSV (read-only), then:
- per material: echo/direct and echo/surface dB (gated contrast peak), relief-pick
  correlation across the 13 stations, mean total-envelope decay curve;
- real-data anchors recomputed from the CSV: direct peak 15.4 ns, ground bounce at
  ~60 ns, interface band 400-500 ns, late-peak/early-peak ratio;
- depth extrapolation (DECLARED MODEL): conductive two-way loss alpha=sigma*eta/2 at
  95 MHz applied to assumed cover thickness candidates {13, 18.5} m. This is not a
  full-depth forward run or field material calibration. The historical ~-27 dB
  late peak was an FFT endpoint artifact; the guarded diagnostic is ~-52.6 dB.
  Field stripe means and simulated contrast peaks are different observables.
Outputs summary.json + figures into a fresh result directory.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from hs_capsule_identity import sha256
from analyze_hs4_height_wavefield import response
from check_hs4_v4_factor_evidence import direct_response
from plot_hs4_permittivity_bscans import C_AIR, SURFACE_Z, interface_relief
from sfcw_official_loader_v0_2 import verify_official_runtime
from field_profile_metrics import late_peak_metrics
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

CHECKS = ROOT / 'artifacts/research_checks'
CAPSULE = CHECKS / '2026-10-05_hs4_material_scan_r1'
REAL_CSV = ROOT / 'real_data_yingshan/营山测线数据/Line9origin(36).csv'
TAKE = np.array([0, 83, 167, 250, 333, 417, 500])
STEP, RX_OFFSET, ANT_Z = 0.5, 1.3, 20.0
F_BAND_HZ = 95e6
REAL_DEPTHS_M = [13.0, 18.5]  # measured interface band 370 ns two-way at v=c/3 (eps 9) and v=c/sqrt(18)


def relief_twoway(ant_z, xmids, fx, fz, v_cover):
    out = []
    for xmid in xmids:
        z = float(np.interp(xmid, fx, fz))
        out.append(2*np.hypot(0.65, ant_z - SURFACE_Z)/C_AIR + 2*(SURFACE_Z - z)/v_cover)
    return np.array(out)


def ratio_db(num, den):
    if num <= 0.0 or den <= 0.0:
        return None  # degenerate: zero interface contrast (eps9_s0.001 cover == rock)
    return float(20*np.log10(num/den))


def real_anchors(csv_path=None):
    csv_path = REAL_CSV if csv_path is None else Path(csv_path)
    with open(csv_path, encoding='utf-8') as f:
        hdr = [f.readline() for _ in range(4)]
    ns = int(hdr[0].split('=')[1].split(',')[0]); T = float(hdr[1].split('=')[1].split(',')[0])
    nt = int(hdr[2].split('=')[1].split(',')[0])
    d = np.loadtxt(csv_path, delimiter=',', skiprows=4)
    if d.shape != (nt*ns, 5) or not np.isfinite(d).all():
        raise ValueError('finite five-column stacked field payload required')
    amp = d[:, 3].reshape(nt, ns).T
    t = np.linspace(0, T, ns)
    from scipy.signal import hilbert
    env = np.abs(hilbert(amp, axis=0))
    mean_env = env.mean(axis=1)
    mdb = 20*np.log10(mean_env/mean_env.max())
    def at(tt):
        return float(mdb[np.argmin(np.abs(t-tt))])
    # Whole-record FFT Hilbert envelopes wrap strong early samples into the
    # end of the record. Historical >=300ns max picked that endpoint artifact.
    late = late_peak_metrics(env, t, late_stop_ns=T-50.)
    return {'t_ns': t, 'mean_db': mdb,
            'direct_peak_ns': float(t[np.argmax(mean_env)]),
            'ground_bounce_db': at(60.0),
            'interface_band_db_mean': float(np.mean(mdb[(t >= 400) & (t <= 500)])),
            'late_peak_over_early_peak_db_median': late['median_dB'],
            'late_peak_window_ns': late['late_window_ns'],
            'late_peak_endpoint_guard_ns': late['endpoint_guard_ns'],
            'anchor_scope': 'Unpadded Hilbert interior-window diagnostic; no calibrated scene/instrument equivalence. Historical unrestricted endpoint ratio invalid.'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    contract = json.loads((CAPSULE/'execution_contract.json').read_text('utf-8'))
    ver = json.loads((CAPSULE/'completed_verification.json').read_text('utf-8'))
    if ver['status'] != 'PASS' or ver['contract_sha256'] != sha256(CAPSULE/'execution_contract.json'):
        raise ValueError('completed verification PASS required')

    fx, fz = interface_relief()
    station_x = np.array([14.6 + STEP*k for k in range(13)])
    xmids = station_x + RX_OFFSET/2

    groups = {g['id']: g for g in contract['groups']}
    tags = sorted({g['material_tag'] for g in contract['groups']},
                  key=lambda t: (next(m[1] for m in [(g['material_tag'], g['cover_eps_r']) for g in contract['groups']] if m[0] == t)))
    spectra, dft_err = {}, {}
    for g in contract['groups']:
        audit = json.loads((Path(g['input']).parent/'audit.json').read_text('utf-8'))
        for row in audit:
            raw = Path(g['input']).parent/row['file']
            if sha256(raw) != row['sha256']:
                raise ValueError('raw identity differs: '+raw.name)
            r = response(raw)
            spectra[g['material_tag'], g['role'], row['station_index']] = r
            exact = direct_response(raw, r.frequency[TAKE])
            err = float(np.linalg.norm(exact - r.response[TAKE])/np.linalg.norm(exact))
            if err > 1e-9:
                raise ValueError('independent DFT failed: '+raw.name)
            dft_err[f"{g['id']}/{row['file']}"] = err

    template = next(iter(spectra.values()))
    products = {}
    for (tag, role, si), r in spectra.items():
        if role == 'rough':
            contrast_spec = r.response - spectra[tag, 'halfspace', si].response
            products[tag, si, 'total'] = reconstruct_time_response(r, window='hann', zero_pad_factor=8)
            products[tag, si, 'contrast'] = reconstruct_time_response(
                replace(template, response=contrast_spec), window='hann', zero_pad_factor=8)
    time = next(iter(products.values())).time * 1e9

    surface_t = 2*np.hypot(0.65, ANT_Z-SURFACE_Z)/C_AIR  # ~53.5 ns
    per_material = {}
    for g in contract['groups']:
        tag = g['material_tag']
        if tag in per_material:
            continue
        eps = g['cover_eps_r']
        v_cover = C_AIR/np.sqrt(eps)
        relief = relief_twoway(ANT_Z, xmids, fx, fz, v_cover)
        rows = []
        for i, si in enumerate(range(1, 122, 10)):
            total = products[tag, si, 'total']
            contrast = products[tag, si, 'contrast']
            mgate = np.abs(time - relief[i]) <= 25.0
            msurf = np.abs(time - surface_t) <= 12.0
            a_echo = float(np.max(np.abs(contrast.complex_envelope[mgate])))
            a_direct = float(np.max(np.abs(total.complex_envelope)))
            a_surface = float(np.max(np.abs(total.complex_envelope[msurf])))
            # contrast floor in the surface->echo gap (cross-scene subtraction floor probe)
            mfloor = (time >= surface_t + 15) & (time <= relief[i] - 30)
            floor = float(np.max(np.abs(contrast.complex_envelope[mfloor])))
            rows.append({'station_index': si,
                         'echo_over_direct_dB': ratio_db(a_echo, a_direct),
                         'echo_over_surface_dB': ratio_db(a_echo, a_surface),
                         'contrast_floor_over_echo_dB': ratio_db(floor, a_echo),
                         'pick_ns': float(time[mgate][np.argmax(np.abs(contrast.complex_envelope[mgate]))])})
        degenerate = all(r['echo_over_direct_dB'] is None for r in rows)
        picks = np.array([r['pick_ns'] for r in rows])
        corr = None if degenerate else float(np.corrcoef(picks, relief)[0, 1])
        # mean total envelope decay
        envs = np.array([np.abs(products[tag, si, 'total'].complex_envelope) for si in range(1, 122, 10)])
        mean_env = envs.mean(axis=0)
        mdb = 20*np.log10(mean_env/mean_env.max())
        eta = 376.730313668/np.sqrt(eps)
        alpha = g['cover_sigma']*eta/2  # Np/m, band centre
        model_cover_m = float(np.mean(SURFACE_Z - np.interp(xmids, fx, fz)))
        tw_model_dB = 20*np.log10(np.exp(1))*2*alpha*model_cover_m
        med_eod = None if degenerate else float(np.median([r['echo_over_direct_dB'] for r in rows]))
        extrap = None if degenerate else {
            f'{d:g}m': med_eod + 20*np.log10(np.exp(1))*2*alpha*(model_cover_m - d)
            for d in REAL_DEPTHS_M}
        per_material[tag] = {
            'cover_eps_r': eps, 'cover_sigma': g['cover_sigma'], 'debye': g['debye'],
            'degenerate_zero_contrast': degenerate,
            'echo_over_direct_dB_median': med_eod,
            'echo_over_surface_dB_median': None if degenerate else float(np.median([r['echo_over_surface_dB'] for r in rows])),
            'contrast_floor_over_echo_dB_worst': None if degenerate else float(max(r['contrast_floor_over_echo_dB'] for r in rows)),
            'relief_pick_corr': corr,
            'alpha_Np_per_m': alpha, 'two_way_conductive_dB_at_model_cover': tw_model_dB,
            'extrapolated_echo_over_direct_dB_at_real_depth': extrap,
            'decay_db': {'surface': float(mdb[np.argmin(np.abs(time-surface_t))]),
                         'echo': float(np.mean([mdb[np.argmin(np.abs(time-rt))] for rt in relief]))},
            'per_station': rows}

    real = real_anchors()
    real_summary = {k: v for k, v in real.items() if k not in ('t_ns', 'mean_db')}

    # ---- figure 1: B-scan panels, total (top) / contrast (bottom), grey envelope dB self-normalized
    nmat = len(per_material)
    fig, axes = plt.subplots(2, nmat, figsize=(2.6*nmat, 7.5), sharey='row', layout='constrained')
    order = list(per_material)
    for col, tag in enumerate(order):
        for row, kind in enumerate(('total', 'contrast')):
            ax = axes[row, col]
            mat = np.stack([np.abs(products[tag, si, kind].complex_envelope) for si in range(1, 122, 10)], axis=1)
            if mat.max() <= 0.0:
                mat = np.full_like(mat, 1e-300)  # degenerate zero-contrast control: show floor
            mdb = 20*np.log10(np.maximum(mat, mat.max()*1e-9)/mat.max())
            ax.imshow(mdb, aspect='auto', cmap='gray', vmin=-60, vmax=0,
                      extent=[station_x[0]-STEP/2, station_x[-1]+RX_OFFSET+STEP/2, time[-1], time[0]])
            eps = per_material[tag]['cover_eps_r']
            relief = relief_twoway(ANT_Z, xmids, fx, fz, C_AIR/np.sqrt(eps))
            ax.plot(xmids, relief, 'r-', lw=0.8)
            ax.set_xlim(station_x[0]-STEP/2, station_x[-1]+RX_OFFSET+STEP/2)
            ax.set_ylim(260, 0)
            ax.set_title(f'{tag}\n{kind}', fontsize=7)
            if row == 1:
                ax.set_xlabel('src x (m)', fontsize=7)
            if col == 0:
                ax.set_ylabel('time (ns)', fontsize=8)
    fig.suptitle('Material scan B-scans at 8 m (self-normalized envelope dB; red = true relief two-way per material)')
    a.out.mkdir(parents=True)
    fig.savefig(a.out/'material_scan_bscans.png', dpi=140)
    plt.close(fig)

    # ---- figure 2: decay curves + real
    fig2, ax2 = plt.subplots(figsize=(9, 5.5), layout='constrained')
    for tag in order:
        envs = np.array([np.abs(products[tag, si, 'total'].complex_envelope) for si in range(1, 122, 10)])
        mdb = 20*np.log10(envs.mean(axis=0)/envs.mean(axis=0).max())
        ax2.plot(mdb, time, lw=1.0, label=tag)
    mreal = real['t_ns'] <= 650  # crop Hilbert edge artefact at the 700 ns endpoint
    ax2.plot(real['mean_db'][mreal], real['t_ns'][mreal], 'k-', lw=2.0, label='REAL Line9(36) mean')
    ax2.invert_yaxis(); ax2.set_xlabel('mean envelope dB'); ax2.set_ylabel('time (ns)')
    ax2.grid(alpha=0.3); ax2.legend(fontsize=6, ncol=2); ax2.set_xlim(-90, 0)
    ax2.set_ylim(700, 0)
    ax2.set_title('Mean trace decay: sim per material vs measured (own time axes; real interface ~430 ns)')
    fig2.savefig(a.out/'decay_vs_real.png', dpi=140)
    plt.close(fig2)

    # ---- figure 3: echo levels bar chart vs real anchor
    fig3, ax3 = plt.subplots(figsize=(10, 4.5), layout='constrained')
    vals = [per_material[t]['echo_over_surface_dB_median'] if per_material[t]['echo_over_surface_dB_median'] is not None else np.nan for t in order]
    ax3.bar(range(nmat), vals)
    ax3.axhline(real['interface_band_db_mean']-real['ground_bounce_db'], color='k', ls='--',
                label=f"real interface-band/ground-bounce ≈ {real['interface_band_db_mean']-real['ground_bounce_db']:.1f} dB")
    ax3.set_xticks(range(nmat)); ax3.set_xticklabels(order, rotation=45, ha='right', fontsize=7)
    ax3.set_ylabel('echo/surface dB (median, 13 stations)'); ax3.legend(fontsize=8)
    ax3.set_title('Interface echo relative to surface reflection per material vs measured anchor')
    fig3.savefig(a.out/'echo_over_surface_vs_real.png', dpi=140)
    plt.close(fig3)

    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'code_sha256': sha256(__file__), 'capsule': CAPSULE.name,
               'independent_dft_max_rel_l2': max(dft_err.values()),
               'real_anchors': real_summary,
               'per_material': per_material,
               'limits': '2D line source; contrast floor probed in surface->echo gap; extrapolation uses band-centre conductive alpha only (no Debye tail, 2D spreading); real anchors from a single CSV; site-adaptation diagnostic, not blind-test performance.'}
    with open(a.out/'summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps({'real': real_summary,
                      'per_material': {t: {'echo_over_surface_dB': per_material[t]['echo_over_surface_dB_median'],
                                           'echo_over_direct_dB': per_material[t]['echo_over_direct_dB_median'],
                                           'relief_corr': per_material[t]['relief_pick_corr']}
                                       for t in order}}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
