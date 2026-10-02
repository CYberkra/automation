"""Independent CPU review of archived HS4 B-scans; never invokes a solver.

Default output is a fresh ignored directory. --out must be a NEW directory.
Input capsules are checked before and after analysis; no historical writes.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime
from gprMax.toolboxes.SFCW.processing import (
    load_source, load_receiver, direct_frequency_response, reconstruct_time_response)

ROOT = Path(__file__).resolve().parents[1]
HS = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'
T2 = ROOT / 'artifacts/research_checks/2026-10-02_hs4t2d_transect'
GROUPS = [('co13', HS, 'hs4_co13', 13, 'hs4_co13_bscan_official.npz', 'Ex'),
          ('col6', HS, 'hs4_col6', 25, 'hs4_col6_bscan_official.npz', 'Ex'),
          ('arch2d', T2, 'hs4t2d', 121, 'hs4t2d_bscan_official.npz', 'Ey'),
          ('step2d', T2, 'hs4t2ds', 121, 'hs4t2dstep_bscan_official.npz', 'Ey')]
WINDOWS = {'direct': (0, 20), 'pre_ground': (20, 80), 'ground': (85, 120),
           'intermediate': (120, 160), 'interface': (160, 220)}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(folder):
    return {p.name: sha(p) for p in folder.iterdir() if p.is_file()}


def manifest_check(folder):
    manifest = json.loads((folder / 'manifest.json').read_text(encoding='utf-8'))
    records = manifest if isinstance(manifest, list) else [
        {'file': k, 'sha256': v} for k, v in manifest.items()]
    require(all(sha(folder / r['file']) == r['sha256'] for r in records), 'Manifest mismatch')
    return len(records)


def run_acceptance(out, capsule, tag):
    dest = out / (tag + '_output')
    result = subprocess.run(
        [sys.executable, str(ROOT / 'scripts/run_hs_acceptance_v0_2.py'),
         '--dir', str(capsule), '--out', str(dest), '--fig', str(out / (tag + '.png'))],
        capture_output=True, text=True, encoding='utf-8', cwd=ROOT)
    (out / (tag + '_stderr.txt')).write_text(result.stderr, encoding='utf-8')
    metrics = dest / 'hs_acceptance_metrics_v0_2.json'
    return {'exit_code': result.returncode, 'stderr_tail': result.stderr[-800:],
            'metrics': json.loads(metrics.read_text(encoding='utf-8')) if metrics.exists() else None}


def analyze(group):
    name, folder, prefix, count, npz_name, component = group
    paths = [folder / f'{prefix}_t{k:02d}.h5' for k in range(1, count + 1)]
    first_src = load_source(paths[0])
    first_rx = load_receiver(paths[0], receiver_path='/rxs/rx1', component=component)
    raw, positions, sources = [], [], []
    for path in paths:
        src = load_source(path)
        rx = load_receiver(path, receiver_path='/rxs/rx1', component=component)
        require(src.dt == first_src.dt and src.time_offset == first_src.time_offset
                and np.array_equal(src.samples, first_src.samples), 'Source differs: ' + str(path))
        require(rx.dt == first_rx.dt and rx.time_offset == first_rx.time_offset
                and rx.samples.shape == first_rx.samples.shape, 'Receiver axes differ')
        with h5py.File(path, 'r') as h:
            positions.append(np.asarray(h['rxs/rx1'].attrs['Position']).tolist())
            sources.append(np.asarray(h['srcs/src1'].attrs['Position']).tolist())
        raw.append(rx.samples)
    raw = np.stack(raw, axis=1)
    n = min(len(raw), int(np.floor(1200e-9 / first_rx.dt)) + 1)
    rx = replace(first_rx, samples=raw[:n])
    response = direct_frequency_response(
        first_src, rx, FREQ, tail_taper_fraction=(round(200e-9 / first_rx.dt) - .25) / n)
    reconstructed = reconstruct_time_response(response, zero_pad_factor=8, window='hann')
    t = reconstructed.time * 1e9
    with np.load(folder / npz_name) as z:
        require(np.array_equal(t, z['t']), 'Saved time axis differs')
        S = z['S'].copy()
    error = float(np.max(np.abs(S - reconstructed.real_bandpass)))
    relative = error / float(np.max(np.abs(S)))
    require(relative < 1e-12, 'Saved signed reconstruction differs')
    m = t <= 250
    tm = t[m]
    a = S[m] - S[m].mean(axis=0, keepdims=True)
    u, singular, vh = np.linalg.svd(a, full_matrices=False)
    residual = a - singular[0] * np.outer(u[:, 0], vh[0])
    energy = {}
    for window, (lo, hi) in WINDOWS.items():
        mask = (tm >= lo) & (tm <= hi)
        energy[window] = {'max_abs': float(np.max(np.abs(residual[mask]))),
                          'squared_norm': float(np.sum(residual[mask] ** 2))}
    center = count // 2
    envelope = 2 * np.abs(reconstructed.complex_envelope)
    ground = (t >= 85) & (t <= 115)
    raw_time = first_rx.time_offset * 1e9 + np.arange(len(raw)) * first_rx.dt * 1e9
    diff = raw[:, [0, count - 1]] - raw[:, center, None]
    early = raw_time < 85
    j, k = np.unravel_index(np.argmax(np.abs(diff[early])), diff[early].shape)
    result = {
        'traces': count, 'component': component, 'saved_shape': list(S.shape),
        'official_rebuild_max_abs_error': error, 'official_rebuild_relative_error': relative,
        'source_signals_identical': True, 'first_two_singular_values': singular[:2].tolist(),
        'rank1_residual_peak_ns': float(tm[np.argmax(np.max(np.abs(residual), axis=1))]),
        'residual_windows': energy,
        'ground_to_interface_residual_energy_ratio': energy['ground']['squared_norm'] /
            energy['interface']['squared_norm'],
        'ground_to_interface_residual_peak_ratio': energy['ground']['max_abs'] /
            energy['interface']['max_abs'],
        'official_ground_envelope_peak_ns_center': float(t[ground][np.argmax(envelope[ground, center])]),
        'signed_ground_abs_peak_ns_center': float(t[ground][np.argmax(np.abs(S[ground, center]))]),
        'pre85ns_edge_minus_center_raw_peak_over_raw_peak': float(np.max(np.abs(diff[early])) / np.max(np.abs(raw))),
        'pre85ns_max_difference_time_ns': float(raw_time[early][j]),
        'source_positions': sources, 'receiver_positions': positions,
    }
    return result, (t, S, tm, residual)


def plot(out, groups):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Microsoft YaHei'
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for column, name in enumerate(('col6', 'arch2d')):
        t, S, tm, residual = groups[name]
        x = np.linspace(2.6, 8.6, S.shape[1]) + 0.65
        bound = float(np.quantile(np.abs(residual), .99))
        for row, (lo, hi) in enumerate(((0, 250), (160, 220))):
            ax = axes[row, column]
            mask = (tm >= lo) & (tm <= hi)
            im = ax.pcolormesh(x, tm[mask], residual[mask] / bound, cmap='gray',
                               vmin=-1, vmax=1, shading='nearest', rasterized=True)
            ax.set_ylim(hi, lo)
            ax.set_xlabel('收发中点坐标 (m)；输入 Tx 坐标 +0.65 m')
            ax.set_ylabel('时间 (ns)')
            label = '三维 COL6' if name == 'col6' else '二维 HS4T2D'
            ax.set_title(f'{label}：去道内均值 + 去1阶，{lo}–{hi} ns\n同列两图使用同一色标，无 AGC')
            if row == 0:
                ax.axhspan(85, 120, color='orange', alpha=.12)
                ax.axhspan(160, 220, color='cyan', alpha=.08)
                ax.text(x[0] + .1, 78, '地表窗 85–120 ns', color='darkorange')
                ax.text(x[0] + .1, 154, '地下候选窗 160–220 ns', color='teal')
            fig.colorbar(im, ax=ax, label=f'幅值 / 本列固定尺度 {bound:.3g}；±1截幅')
    fig.suptitle('独立 B-scan 审查：二维最强残差在地表窗；三维与二维分别定色标，幅值不跨维比较')
    fig.savefig(out / 'bscan_residual_diagnostic.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    out = args.out
    if out is None:
        out = Path(tempfile.mkdtemp(prefix='hs4_bscan_review_', dir=ROOT / 'artifacts/local_checks'))
    else:
        out.mkdir(parents=True, exist_ok=False)
    verify_official_runtime()
    before = {str(p): snapshot(p) for p in (HS, T2)}
    manifests = {str(p.relative_to(ROOT)): manifest_check(p) for p in (HS, T2)}
    results, plots = {}, {}
    for group in GROUPS:
        results[group[0]], plots[group[0]] = analyze(group)
        print(f"Rebuilt {group[0]}: {results[group[0]]['traces']} traces", flush=True)
    normal = run_acceptance(out, HS, 'current_manifest')
    require(normal['exit_code'] != 0 and 'AttributeError' in normal['stderr_tail'],
            'Manifest regression no longer reproduces; revisit review')
    capsule = out / 'compatible_manifest_copy'
    capsule.mkdir()
    identity = {}
    for name in ('hs1_flat_halfspace', 'hs2_coveronly_halfspace', 'hs3_domain16_halfspace'):
        shutil.copy2(HS / (name + '.h5'), capsule / (name + '.h5'))
        identity[name + '.h5'] = sha(capsule / (name + '.h5'))
    (capsule / 'manifest.json').write_text(json.dumps(identity), encoding='utf-8')
    compatible = run_acceptance(out, capsule, 'compatible_manifest')
    require(compatible['exit_code'] == 0, 'Compatible manifest control failed')
    p = capsule / 'hs3_domain16_halfspace.h5'
    with h5py.File(p, 'r+') as h:
        h['rxs/rx1/Ex'][:] *= -1
    identity[p.name] = sha(p)
    (capsule / 'manifest.json').write_text(json.dumps(identity), encoding='utf-8')
    polarity = run_acceptance(out, capsule, 'hs3_polarity_flip')
    require(polarity['exit_code'] == 0 and polarity['metrics']['verdict'] == 'PASS',
            'Polarity counterexample no longer passes; revisit review')
    plot(out, plots)
    require(all(before[str(p)] == snapshot(p) for p in (HS, T2)), 'Historical capsule changed')
    result = {'schema': 'hs4_bscan_review/1', 'reviewed_commit': 'd6394b1',
              'solver_invoked': False, 'test_family_read': False,
              'input_manifests_verified': manifests, 'capsules_unchanged': True,
              'groups': results, 'current_acceptance': normal,
              'compatible_manifest_control': compatible, 'polarity_counterexample': polarity,
              'reference_gap': 'No tracked 2D flat BG H5/array or step xcorr/plot script accompanies the claimed TGT-BG validation.',
              'source_sha256': sha(Path(__file__)),
              'versions': {'python': sys.version, 'numpy': np.__version__}}
    with (out / 'review.json').open('x', encoding='utf-8', newline='\n') as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print('Review output: ' + str(out.resolve()), flush=True)


if __name__ == '__main__':
    main()
