"""Background suppression on the P2 3D antenna B-scan (user 2026-10-05: 做个背景抑制呢).

Applies the frozen operator-catalogue background operators (mean trace, SVD
rank 1/2/3 removal, RPCA lambda=1.0; scripts/research_operator_contract.py,
configs/research/operator_catalogue_v0.1.json + v0.2 RPCA axis) to the total
rough B-scan (the realistic data product; no halfspace available in the
field), on the SFCW-reconstructed real bandpass in a DECLARED window
0-300 ns, axes [sample, trace].

Evaluation (mechanism diagnostic, no physical thresholds):
- background residual ratio: RMS(0-140 ns) after / before suppression;
- relief-window preservation: per-station signed correlation and relative L2
  against the gated rough-minus-halfspace contrast (reference carries the
  cross-scene numerical floor caveat, see
  docs/research/2026-10-05_cross_scene_subtraction_floor.md);
- negative control: SVD rank-2 on the halfspace B-scan (no interface present)
  to expose what suppression leaves/invents without an interface.
Per AGENTS.md the removed component is reported as 'removed component', not
'noise'; horizontal/continuity is not treated as truth.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from hs_capsule_identity import sha256
from analyze_antenna_bscan_power import (load_traces, gate, station_x,
                                         N_STATIONS, STEP, PANEL_D,  # noqa: F401
                                         ANT_Z_ABS)
from analyze_hs4_power_budget import relief_twoway
from plot_hs4_permittivity_bscans import C_AIR, interface_relief
from research_operator_contract import apply_configuration, ConfigUnavailable
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

WINDOW_NS = (0.0, 300.0)
RELIEF_HALF_WINDOW_NS = 25.0
OPERATORS = [('mean x1.0', 'B3_G1_BG', '0.1'),
             ('svd rank-1', 'B4_G1_BG', '0.1'),
             ('svd rank-2', 'B5_G1_BG', '0.1'),
             ('svd rank-3', 'B6_G1_BG', '0.1'),
             ('rpca lam=1.0', 'B8_G1_BG', '0.2')]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    spectra, midx = load_traces()
    template = spectra['rough_k06']

    def product(spec):
        return reconstruct_time_response(replace(template, response=spec),
                                         window='hann', zero_pad_factor=8)

    # B-scan matrices on the declared window, axes [sample, trace]
    probe = product(spectra['rough_k06'].response)
    t_all = probe.time * 1e9
    mt = (t_all >= WINDOW_NS[0]) & (t_all <= WINDOW_NS[1])
    time = t_all[mt]
    rough = np.stack([product(spectra[f'rough_k{k:02d}'].response).real_bandpass[mt]
                      for k in range(N_STATIONS)], axis=1)
    halfspace = np.stack([product(spectra[f'halfspace_k{k:02d}'].response).real_bandpass[mt]
                          for k in range(N_STATIONS)], axis=1)
    ref = np.stack([gate(product(spectra[f'rough_k{k:02d}'].response
                                 - spectra[f'halfspace_k{k:02d}'].response)).real_bandpass[mt]
                    for k in range(N_STATIONS)], axis=1)

    xmids = np.array([midx[f'rough_k{k:02d}'] for k in range(N_STATIONS)])
    fx, fz = interface_relief()
    relief = relief_twoway(ANT_Z_ABS, xmids, fx, fz, C_AIR / np.sqrt(18.017))

    # relief window mask [sample, trace]
    wmask = np.zeros_like(rough, dtype=bool)
    for i in range(N_STATIONS):
        wmask[:, i] = np.abs(time - relief[i]) <= RELIEF_HALF_WINDOW_NS
    early = time <= 140.0

    results = {}
    suppressed = {}
    for label, cid, ver in OPERATORS:
        try:
            r = apply_configuration(rough, cid, catalogue_version=ver)
        except ConfigUnavailable as exc:
            results[label] = {'config_id': cid, 'status': f'unavailable: {exc}'}
            continue
        y = r['output']
        suppressed[label] = y
        resid = float(np.sqrt(np.mean(y[early] ** 2)) / np.sqrt(np.mean(rough[early] ** 2)))
        per_station = []
        for i in range(N_STATIONS):
            m = wmask[:, i]
            a_, b_ = y[m, i], ref[m, i]
            denom = float(np.linalg.norm(a_) * np.linalg.norm(b_))
            corr = float(np.dot(a_, b_) / denom) if denom > 0 else None
            rl2 = float(np.linalg.norm(a_ - b_) / np.linalg.norm(b_)) if np.linalg.norm(b_) > 0 else None
            per_station.append({'corr': corr, 'rel_l2': rl2})
        diag = r['steps'][0]['diagnostics']
        results[label] = {
            'config_id': cid, 'status': 'ok',
            'removed_component_kind': diag.get('kind'),
            'operator_diagnostics': {k: v for k, v in diag.items()
                                     if k in ('k', 'effective_k', 'lambda', 'lam_factor',
                                              'cutoff_gap_relative_to_s1', 'status',
                                              'rank_L', 'rel_residual', 'iterations')},
            'relative_singular_values': diag.get('relative_singular_values'),
            'background_residual_ratio_0_140ns': resid,
            'relief_window_vs_gated_contrast': per_station,
            'relief_window_corr_median': float(np.nanmedian([p['corr'] for p in per_station
                                                             if p['corr'] is not None]))}

    # negative control: svd rank-2 on halfspace (no interface)
    ctrl = apply_configuration(halfspace, 'B5_G1_BG')
    suppressed['svd rank-2 ON HALFSPACE (control)'] = ctrl['output']
    ctrl_resid = float(np.sqrt(np.mean(ctrl['output'][early] ** 2))
                       / np.sqrt(np.mean(halfspace[early] ** 2)))
    ctrl_late = float(np.sqrt(np.mean(ctrl['output'][wmask] ** 2)))

    # figure
    panels = [('total rough (input)', rough), ('gated contrast (reference, floor caveat)', ref)]
    panels += [(f'suppressed: {l}', suppressed[l]) for l in suppressed]
    vmax = max(float(np.max(np.abs(p))) for _, p in panels[:2])
    vmax_s = max(float(np.max(np.abs(p))) for _, p in panels[2:])
    for cmap, suffix in (('RdBu_r', ''), ('gray', '_gray')):
        fig, axes = plt.subplots(len(panels), 1, figsize=(13, 2.4 * len(panels)),
                                 sharex=True, layout='constrained')
        for ax, (label, p) in zip(axes, panels):
            v = vmax if label.startswith(('total', 'gated')) else vmax_s
            t_edge = np.concatenate(([time[0] - (time[1] - time[0]) / 2],
                                     (time[:-1] + time[1:]) / 2,
                                     [time[-1] + (time[1] - time[0]) / 2]))
            for k in range(N_STATIONS):
                ax.pcolormesh([station_x(k) - STEP / 2, station_x(k) + STEP / 2], t_edge,
                              p[:, [k]], cmap=cmap, shading='flat',
                              norm=matplotlib.colors.SymLogNorm(linthresh=1e-4 * v,
                                                                vmin=-v, vmax=v))
            ax.plot(xmids, relief, 'k-', lw=1.0)
            ax.axvspan(120, 140, color='orange', alpha=0.2, lw=0)
            ax.set_xlim(station_x(0) - STEP / 2, station_x(N_STATIONS - 1) + STEP / 2)
            ax.set_ylim(260, 60)
            ax.set_ylabel('ns'); ax.set_title(label, fontsize=9)
        axes[-1].set_xlabel('source x (m)')
        fig.suptitle('Background suppression, 3D x4 antenna B-scan (declared window 0-300 ns; '
                     'removed component is not certified noise)', fontsize=10)
        a.out.mkdir(parents=True, exist_ok=True)
        fig.savefig(a.out / f'background_suppression{suffix}.png', dpi=150)
        plt.close(fig)

    summary = {'status': 'COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE',
               'code_sha256': sha256(__file__),
               'input_window_ns': WINDOW_NS,
               'input_product': 'SFCW real bandpass (hann, zero_pad 8), axes [sample, trace]',
               'operators': results,
               'negative_control': {'config': 'svd rank-2 on halfspace',
                                    'background_residual_ratio_0_140ns': ctrl_resid,
                                    'relief_window_rms_remaining': ctrl_late},
               'reference_caveat': ('gated contrast carries the cross-scene numerical floor; '
                                    'correlation targets are mechanism-level'),
               'removed_component_statement': ('mean/SVD/RPCA components are coherent-horizontal '
                                             'content (direct+surface waves), not certified noise; '
                                             'a flat real interface would also be removed by these '
                                             'operators'),
               'physical_attribution_certified': False,
               'scope': 'array-level suppression on simulation products; no field data.'}
    (a.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                                   allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({l: (r.get('background_residual_ratio_0_140ns'),
                          r.get('relief_window_corr_median'))
                      for l, r in results.items()}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
