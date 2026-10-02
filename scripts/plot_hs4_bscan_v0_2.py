"""Reproducible CPU diagnostic plots of four archived HS4 groups.

Rebuild all signed arrays using the pinned official chain, verify input
manifests before/after, and write a NEW capsule. No solver or physical gate.
Raw amplitude, rank-1 residual, and display-only AGC have separate scales.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

from hs_capsule_identity import read_manifest, verify_file, sha256
from review_bscan_hs4_20261003 import GROUPS, analyze
from sfcw_official_loader_v0_2 import verify_official_runtime, OFFICIAL_PROCESSING_SHA256


def diagnostic_arrays(t, signed):
    t, signed = np.asarray(t), np.asarray(signed)
    if (t.ndim != 1 or len(t) < 2 or signed.ndim != 2 or signed.shape[0] != len(t)
            or signed.shape[1] < 2 or np.iscomplexobj(signed)
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(signed))
            or not np.all(np.diff(t) > 0)
            or not np.allclose(np.diff(t), t[1] - t[0], rtol=1e-10, atol=1e-12)):
        raise ValueError('requires finite signed [sample, trace] data on a uniform ns axis')
    mask = (t >= 0) & (t <= 250)
    if np.count_nonzero(mask) < 3:
        raise ValueError('0-250 ns diagnostic window empty/too short')
    tm, raw = t[mask], signed[mask].copy()
    centered = raw - raw.mean(axis=0, keepdims=True)
    u, singular, vh = np.linalg.svd(centered, full_matrices=False)
    residual = centered - singular[0] * np.outer(u[:, 0], vh[0])
    # AGC is display only: 20 ns centered RMS, reflected edges, shared 1% floor.
    samples = max(1, int(round(20 / (tm[1] - tm[0]))))
    rms = np.sqrt(np.maximum(0, uniform_filter1d(residual ** 2, samples, axis=0, mode='reflect')))
    floor = .01 * float(rms.max())
    gain = np.divide(1., np.maximum(rms, floor), out=np.ones_like(rms), where=rms > 0)
    return tm, raw, centered, residual, gain, residual * gain, samples, floor


def plot(path, name, tm, midpoints, raw, residual, agc):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    limits = [float(np.quantile(np.abs(v), .999)) or 1. for v in (raw, residual)] + [3.]
    columns = [('Raw signed response', raw, 'signed field/source response'),
               ('Temporal mean + rank-1 removed', residual, 'signed field/source response'),
               ('Residual + AGC (display only)', agc, 'normalized; amplitude not preserved')]
    for row, (lo, hi) in enumerate(((0, 250), (160, 220))):
        mask = (tm >= lo) & (tm <= hi)
        for col, (title, matrix, unit) in enumerate(columns):
            ax = axes[row, col]
            im = ax.pcolormesh(midpoints, tm[mask], matrix[mask], shading='nearest', cmap='gray',
                               vmin=-limits[col], vmax=limits[col], rasterized=True)
            ax.set_ylim(hi, lo)
            ax.set_ylabel('Reconstructed time (ns)')
            ax.set_xlabel('Tx/Rx midpoint along profile (m)')
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            ax.set_title(title)
            fig.colorbar(im, ax=ax, label=unit)
    fig.suptitle(f'{name}: diagnostic only; fixed scale across both time windows\n'
                 'SVD can remove common geological responses; AGC brightness is not SNR')
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return {'full_window_ns': [0, 250], 'detail_window_ns': [160, 220],
            'raw_residual_clip_quantile': .999, 'symmetric_color_limits': limits,
            'time_units': 'ns', 'coordinate': 'actual H5 Tx/Rx midpoint (m)'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, help='NEW directory; historical capsules are read-only')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise ValueError('output exists; refuse to overwrite history')
    if any(out.resolve().is_relative_to(group[1].resolve()) for group in GROUPS):
        raise ValueError('output must be outside read-only input capsules')
    verify_official_runtime()
    manifests = {}
    identities = {}
    for _, folder, _, _, _, _ in GROUPS:
        if folder not in manifests:
            manifests[folder] = (sha256(folder / 'manifest.json'), read_manifest(folder))
            identities[str(folder)] = [verify_file(folder, manifests[folder][1], name)
                                       for name in manifests[folder][1]]
    results = {}
    data = {}
    for group in GROUPS:
        name = group[0]
        metrics, (t, signed, _, _) = analyze(group)
        arrays = diagnostic_arrays(t, signed)
        axis = 1 if name in ('co13', 'col6') else 0
        tx = np.asarray(metrics['source_positions'])[:, axis]
        rx = np.asarray(metrics['receiver_positions'])[:, axis]
        midpoint = (tx + rx) / 2
        if not np.all(np.diff(midpoint) > 0):
            raise ValueError('profile coordinates are not strictly increasing')
        results[name] = metrics
        data[name] = (arrays, tx, rx, midpoint)
    # Reject concurrent capsule changes before creating any result.
    for folder, (digest, records) in manifests.items():
        if sha256(folder / 'manifest.json') != digest:
            raise ValueError('input manifest changed during analysis')
        for name in records:
            verify_file(folder, records, name)
    out.mkdir(parents=True)
    for name, (arrays, tx, rx, midpoint) in data.items():
        tm, raw, centered, residual, gain, agc, samples, floor = arrays
        np.savez_compressed(out / f'{name}_diagnostic.npz', t_ns=tm, raw=raw, centered=centered,
                            rank1_residual=residual, display_agc_gain=gain, display_agc=agc,
                            source_profile_m=tx, receiver_profile_m=rx, midpoint_m=midpoint)
        results[name]['plot'] = plot(out / f'{name}_bscan_v0_2.png', name, tm, midpoint, raw, residual, agc)
        results[name]['agc'] = {'role': 'display_only', 'rms_window_ns': 20, 'window_samples': samples,
                                'edge_mode': 'reflect', 'shared_floor_fraction': .01, 'floor': floor}
    report = {'version': 'hs4-bscan-diagnostic/0.2', 'solver_executed': False,
              'official_processing_sha256': OFFICIAL_PROCESSING_SHA256,
              'physical_validation': 'pending_same_condition_reference_and_numerical_controls',
              'flat_2d_reference_status': 'missing; no step shifts or 7x-smoothing acceptance generated',
              'input_manifests': {str(p): d for p, (d, _) in manifests.items()},
              'input_identities': identities, 'groups': results,
              'code_identities': {p.name: sha256(p) for p in
                  [Path(__file__), Path(__file__).with_name('review_bscan_hs4_20261003.py'),
                   Path(__file__).with_name('hs_capsule_identity.py'),
                   Path(__file__).with_name('sfcw_official_loader_v0_2.py')]}}
    (out / 'diagnostic_metrics.json').write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding='utf-8')
    (out / 'manifest.json').write_text(json.dumps([
        {'file': p.name, 'sha256': sha256(p), 'bytes': p.stat().st_size}
        for p in sorted(out.iterdir()) if p.is_file()], indent=1), encoding='utf-8')
    print(json.dumps({k: {'rebuild_error': v['official_rebuild_max_abs_error'],
                         'ground_interface_energy_ratio': v['ground_to_interface_residual_energy_ratio']}
                      for k, v in results.items()}, indent=1))


if __name__ == '__main__':
    main()
