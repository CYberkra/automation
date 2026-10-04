"""Quantify PML-only effects in raw and fixed-band matched contrasts."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
from analyze_hs4t2d_cause_controls import load_response, weighted_comparison
from analyze_hs4t2d_source_controls import array_comparison
from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import verify_official_runtime, OFFICIAL_PROCESSING_SHA256

ROOT = Path(__file__).resolve().parents[1]


def measured_comparison(reference, candidate):
    result = array_comparison(reference, candidate)
    result.update({'reference_L2': float(np.linalg.norm(reference)),
                   'difference_L2': float(np.linalg.norm(candidate-reference)),
                   'reference_peak_abs': float(np.max(np.abs(reference))),
                   'difference_peak_abs': float(np.max(np.abs(candidate-reference)))})
    return result


def raw_matrix(folder):
    rows = json.loads((folder / 'audit.json').read_text('utf-8'))
    values = []
    for row in rows:
        path = folder / row['file']
        if sha256(path) != row['sha256']:
            raise ValueError('raw identity changed')
        with h5py.File(path) as h:
            values.append(h['rxs/rx1/Ey'][:])
            dt = float(h.attrs['dt'])
    return np.stack(values, axis=1), np.arange(len(values[0])) * dt * 1e9


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new analysis directory required')
    verify_official_runtime()
    c = json.loads(args.contract.read_text('utf-8'))
    events = [json.loads(s) for s in (args.contract.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('status') != 'COMPLETED' or events[-1].get('traces') != 32 or events[0]['contract_sha256'] != sha256(args.contract):
        raise ValueError('completed frozen32-trace study required')
    references = {role: ROOT / ('artifacts/research_checks/2026-10-03_hs4t2d_cause_controls/wide36_' + role)
                  for role in ('rough', 'halfspace')}
    response, baseline_raw = {}, {}
    source_audit = {}
    for role, folder in references.items():
        response[role], _ = load_response(folder)
        baseline_raw[role], raw_time = raw_matrix(folder)
    base = response['rough']
    baseline = base.response - response['halfspace'].response
    if not base.source_valid.all() or not response['halfspace'].source_valid.all():
        raise ValueError('source normalization incomplete')
    spectra = {'baseline_x20': baseline}
    raw_contrast = {'baseline_x20': baseline_raw['rough'] - baseline_raw['halfspace']}
    for g in c['groups']:
        folder = args.contract.parent / g['id']
        r, _ = load_response(folder)
        if not r.source_valid.all() or r.source.dt != base.source.dt or not np.array_equal(r.source.samples, base.source.samples) or r.source.spatial_scale != base.source.spatial_scale or r.source.quantity != base.source.quantity or r.source.units != base.source.units:
            raise ValueError('source conditions differ')
        response[g['id']] = r
        matrix, t = raw_matrix(folder)
        if not np.array_equal(t, raw_time):
            raise ValueError('raw time axis differs')
        baseline_raw[g['id']] = matrix
        source_audit[g['id']] = {'quantity': r.source.quantity, 'units': r.source.units,
                                'spatial_scale': r.source.spatial_scale, 'all501bins_valid': True}
    factors = ['x40wide', 'x60wide', 'bottom40wide', 'top40wide']
    for name in factors:
        spectra[name] = response[name + '_rough'].response - response[name + '_halfspace'].response
        raw_contrast[name] = baseline_raw[name + '_rough'] - baseline_raw[name + '_halfspace']
    products = {name: reconstruct_time_response(replace(base, response=s), window='hann', zero_pad_factor=8)
                for name, s in spectra.items()}
    time = products['baseline_x20'].time * 1e9
    windows = {'ground': (85, 120), 'underground': (160, 220), 'bottom_proxy': (300, 450)}
    metrics = {}
    for name in factors:
        reference = slice(None) if name == 'x40wide' else slice(6, 7)
        item = {'complex_spectral': weighted_comparison(baseline[:, reference], spectra[name]),
                'raw_full_record': measured_comparison(raw_contrast['baseline_x20'][:, reference], raw_contrast[name])}
        for window_name, (lo, hi) in windows.items():
            band_mask = (time >= lo) & (time <= hi)
            raw_mask = (raw_time >= lo) & (raw_time <= hi)
            item[window_name] = {
                'signed_band': measured_comparison(products['baseline_x20'].real_bandpass[band_mask, reference], products[name].real_bandpass[band_mask]),
                'complex_envelope': measured_comparison(products['baseline_x20'].complex_envelope[band_mask, reference], products[name].complex_envelope[band_mask]),
                'raw_FDTD': measured_comparison(raw_contrast['baseline_x20'][raw_mask, reference], raw_contrast[name][raw_mask])}
        metrics[name] = item
    band_mask = (time >= 160) & (time <= 220)
    metrics['x40_to_x60_center'] = {
        'complex_spectral': weighted_comparison(spectra['x40wide'][:,6:7], spectra['x60wide']),
        'underground': {'signed_band': measured_comparison(products['x40wide'].real_bandpass[band_mask,6:7], products['x60wide'].real_bandpass[band_mask]),
                        'complex_envelope': measured_comparison(products['x40wide'].complex_envelope[band_mask,6:7], products['x60wide'].complex_envelope[band_mask])}}
    args.out.mkdir(parents=True)
    arrays = {'time_ns': time, 'raw_time_ns': raw_time, 'midpoint_m': np.arange(13) * .5 + 3.25,
              'frequency_Hz': base.frequency}
    for name in spectra:
        arrays[name + '_frequency_complex'] = spectra[name]
        arrays[name + '_signed'] = products[name].real_bandpass
        arrays[name + '_complex_envelope'] = products[name].complex_envelope
        arrays[name + '_raw_FDTD_contrast'] = raw_contrast[name]
    np.savez_compressed(args.out / 'boundary_arrays.npz', **arrays)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Microsoft YaHei'
    selected = (time >= 160) & (time <= 220)
    x = arrays['midpoint_m']
    before = products['baseline_x20'].real_bandpass[selected]
    after = products['x40wide'].real_bandpass[selected]
    limit = float(np.max(np.abs(np.stack([before, after])))) or 1
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    for ax, a, label in zip(axes, [before, after, after - before], ['左右PML20单元', '左右PML40单元', '两者差值（同一色标）']):
        im = ax.pcolormesh(x, time[selected], a, cmap='gray', vmin=-limit, vmax=limit, shading='nearest')
        ax.set_ylim(220, 160); ax.set_title(label); ax.set_xlabel('原模型收发中点X（m）'); ax.set_ylabel('时间（ns）')
        fig.colorbar(im, ax=ax, label='带符号场/电流响应，无归一化')
    fig.suptitle('36m域 / 15m离地高度 / 相同实际材料网格 / 仅左右PML厚度改变\n全覆盖层匹配复谱差分 / 20–170MHz / Hann / 无SVD或增益')
    fig.savefig(args.out / 'side_pml_bscan.png', dpi=145); plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    for name in spectra:
        k = 6 if products[name].real_bandpass.shape[1] == 13 else 0
        axes[0].plot(time[selected], products[name].real_bandpass[selected, k], label=name)
        axes[1].plot(raw_time, raw_contrast[name][:, k], label=name)
    axes[0].set_xlim(160, 220); axes[0].set_ylabel('带限场/电流响应')
    axes[1].set_xlim(150, 450); axes[1].set_ylabel('原始FDTD Ey差分（V/m）')
    for ax in axes:
        ax.legend(fontsize=8); ax.set_xlabel('时间（ns）')
    fig.suptitle('中心道：分别只加厚侧向/底部/顶部PML；所有输入几何一致\n两行是不同数据域，不比较纵轴幅度；不由峰形强判事件身份')
    fig.savefig(args.out / 'boundary_center_domains.png', dpi=145); plt.close(fig)
    report = {'status': 'COMPLETED', 'new_traces': 32, 'contract_sha256': sha256(args.contract),
              'official_processing_sha256': OFFICIAL_PROCESSING_SHA256, 'code_sha256': sha256(__file__),
              'source_audit': source_audit, 'comparisons': metrics, 'shared_side_bscan_limit': limit,
              'fixed_windows_ns': windows, 'scope': 'PML-only sensitivity in specified2D case/windows. No absolute boundary-error bound, spatial convergence, field claim or event-identity acceptance.'}
    (args.out / 'results.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    (args.out / 'manifest.json').write_text(json.dumps([{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha256(p)}
        for p in sorted(args.out.iterdir())], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
