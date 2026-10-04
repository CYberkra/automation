"""Interface echo visibility budget for HS4T2D (CPU-only, archived H5 inputs).

Quantifies, in the actual 15 m / 2.5 cm condition, why the bedrock-cover
interface shape is not visible in the unimaged B-scan even though propagation
is fine: per-station level of the interface contrast response (rough -
full-cover) against (a) the no-interface background in the same window and
(b) the total response, plus the common-mode/shape decomposition across the
13 stations. Ratios are windowed waveform norms, not power percentages or
detectability claims. No solver runs; inputs verified against capsule
manifests / independent verification hashes before use.
"""
import argparse
import importlib.machinery
import json
import sys
import types
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REMAINING = ROOT / 'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_remaining'
CENTRE = ROOT / 'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'
CONTINUATION = ROOT / 'artifacts/research_checks/2026-10-04_hs4_height_wavefield_continuation'
PATCH_NPZ = ROOT / 'artifacts/research_checks/2026-10-04_hs4_local_patch_results/patch_arrays.npz'
PROFILE_CSV = ROOT / 'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv'
WINDOWS_NS = {'early': (160.0, 180.0), 'later': (180.0, 220.0), 'full': (160.0, 220.0)}
HIGH_INDICES = list(range(1, 122, 10))


def shim_gprmax(src):
    """Import path shim: use the reviewed processing.py without executing the
    solver-heavy gprMax package __init__ (no cython/terminaltables needed)."""
    try:
        import gprMax.toolboxes.SFCW.processing  # noqa: F401
        return 'installed'
    except ImportError:
        for name, sub in [('gprMax', 'gprMax'), ('gprMax.toolboxes', 'gprMax/toolboxes'),
                          ('gprMax.toolboxes.SFCW', 'gprMax/toolboxes/SFCW'),
                          ('gprMax.toolboxes.Utilities', 'gprMax/toolboxes/Utilities')]:
            mod = types.ModuleType(name)
            mod.__path__ = [str(Path(src) / sub)]
            mod.__spec__ = importlib.machinery.ModuleSpec(name, None, is_package=True)
            sys.modules[name] = mod
        return 'namespace_shim:' + str(src)


def sha256(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest_sha(capsule, rel):
    entries = json.loads((capsule / 'manifest.json').read_text('utf-8'))
    for e in entries:
        if e['file'] == rel:
            return e['sha256']
    raise ValueError(f'{rel} not in {capsule.name} manifest')


def high_path(role, index):
    if index == 61:
        return CENTRE, f'centre_{role}/profile.h5'
    segment, first = ('left', 1) if index < 61 else ('right', 71)
    rel = f'{segment}_{role}/profile{(index - first) // 10 + 1}.h5'
    return REMAINING, rel


def low_identity(role):
    v = json.loads((CONTINUATION / 'independent_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS':
        raise ValueError('continuation independent verification not PASS')
    for g in v['groups']:
        if g['id'] == 'low_' + role:
            return g['raw_sha256']
    raise ValueError('low_' + role + ' not in continuation verification')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--gprmax-src', type=Path, default=Path(r'E:\gprMax-v.4.0.0\gprMax-v.4.0.0'))
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new result directory required')
    import_mode = shim_gprmax(args.gprmax_src)
    sys.path.insert(0, str(ROOT / 'scripts'))
    from sfcw_official_loader_v0_2 import verify_official_runtime
    from analyze_hs4_height_wavefield import response
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    verify_official_runtime()

    # Input identity: every H5 checked against its capsule's pinned hash first.
    inputs = {}
    for role in ('rough', 'halfspace'):
        for idx in HIGH_INDICES:
            capsule, rel = high_path(role, idx)
            actual = sha256(capsule / rel)
            expected = manifest_sha(capsule, rel)
            if actual != expected:
                raise ValueError(f'identity mismatch: {rel}')
            inputs[f'high_{role}_{idx}'] = capsule / rel
    for role in ('rough', 'halfspace'):
        p = CONTINUATION / f'low_{role}' / 'profile.h5'
        if sha256(p) != low_identity(role):
            raise ValueError(f'identity mismatch: low_{role}')
        inputs[f'low_{role}'] = p

    spectra = {k: response(p) for k, p in inputs.items()}
    for k, r in spectra.items():
        if not r.source_valid.all():
            raise ValueError(f'invalid source bins: {k}')

    def product(spec, template):
        return reconstruct_time_response(replace(template, response=spec),
                                         window='hann', zero_pad_factor=8)

    midpoints = 3.25 + 0.5 * np.arange(13)
    rough = np.stack([spectra[f'high_rough_{i}'].response for i in HIGH_INDICES], axis=1)
    halfspace = np.stack([spectra[f'high_halfspace_{i}'].response for i in HIGH_INDICES], axis=1)
    diff = rough - halfspace  # complex pair difference before reconstruction
    template = spectra['high_rough_61']
    prod = {k: product(v, template) for k, v in
            [('rough', rough), ('halfspace', halfspace), ('diff', diff)]}
    t_ns = prod['diff'].time * 1e9

    with np.load(PATCH_NPZ) as old:
        anchor_abs = float(np.max(np.abs(old['high_base_spectrum'] - diff)))
        anchor_rel_l2 = float(np.linalg.norm(old['high_base_spectrum'] - diff)
                              / np.linalg.norm(old['high_base_spectrum']))
    # Chain-consistency gate (not a physical threshold): the archived baseline was
    # produced by the same reviewed processing.py on the same H5 inputs; only
    # environment floating-point noise (numpy/scipy builds) may differ.
    if anchor_rel_l2 > 1e-9:
        raise ValueError('reconstructed pair difference differs from archived local-patch baseline')

    low_diff = spectra['low_rough'].response - spectra['low_halfspace'].response
    prod['low_rough'] = product(spectra['low_rough'].response, spectra['low_rough'])
    prod['low_halfspace'] = product(spectra['low_halfspace'].response, spectra['low_rough'])
    prod['low_diff'] = product(low_diff, spectra['low_rough'])
    air_shift_ns = 2 * 13 / 299792458 * 1e9

    def wnorm(p, mask):
        return np.linalg.norm(p.complex_envelope[mask], axis=0)

    metrics = {'windows_ns': WINDOWS_NS, 'air_shift_ns_low': air_shift_ns}
    for label, (lo, hi) in WINDOWS_NS.items():
        mask = (t_ns >= lo) & (t_ns <= hi)
        n_d, n_hs, n_ro = wnorm(prod['diff'], mask), wnorm(prod['halfspace'], mask), wnorm(prod['rough'], mask)
        entry = {
            'interface_contrast_over_nointerface_background_dB': (20 * np.log10(n_d / n_hs)).tolist(),
            'interface_contrast_over_total_response_dB': (20 * np.log10(n_d / n_ro)).tolist(),
            'nointerface_background_over_total_response_dB': (20 * np.log10(n_hs / n_ro)).tolist(),
        }
        if label == 'full':
            # Common-mode vs shape decomposition on signed traces (linear, windowed).
            for name in ('rough', 'diff'):
                sig = prod[name].real_bandpass[mask]  # [samples, stations]
                common = sig.mean(axis=1, keepdims=True)
                total = float(np.sum(sig ** 2))
                entry[name + '_common_mode_energy_share'] = float(13 * np.sum(common ** 2) / total)
                entry[name + '_shape_energy_share'] = float(np.sum((sig - common) ** 2) / total)
            env = np.abs(prod['diff'].complex_envelope[mask])
            peak_idx = np.argmax(env, axis=0)
            entry['interface_contrast_envelope_peak_ns'] = t_ns[mask][peak_idx].tolist()
            entry['contrast_peak_time_span_ns'] = float(t_ns[mask][peak_idx].max() - t_ns[mask][peak_idx].min())
            full_mask = t_ns <= 250.0
            for name in ('rough',):
                inwin = np.max(np.abs(prod[name].real_bandpass[mask]), axis=0)
                whole = np.max(np.abs(prod[name].real_bandpass[full_mask]), axis=0)
                entry[name + '_signed_window_peak_over_full_peak'] = (inwin / whole).tolist()
        # Low-altitude centre station: fixed air-delay window shift only.
        lmask = (t_ns + air_shift_ns >= lo) & (t_ns + air_shift_ns <= hi)
        ld = float(wnorm(prod['low_diff'], lmask))
        lhs = float(wnorm(prod['low_halfspace'], lmask))
        lro = float(wnorm(prod['low_rough'], lmask))
        entry['low_centre_interface_contrast_over_nointerface_background_dB'] = float(20 * np.log10(ld / lhs))
        entry['low_centre_interface_contrast_over_total_response_dB'] = float(20 * np.log10(ld / lro))
        metrics[label] = entry

    geom = np.genfromtxt(PROFILE_CSV, delimiter=',', names=True)
    bins = np.clip(np.floor(midpoints / 0.25).astype(int), 0, 47)
    cover_depth = geom['cover_depth_m'][bins]
    peaks = np.asarray(metrics['full']['interface_contrast_envelope_peak_ns'])
    corr = float(np.corrcoef(cover_depth, peaks)[0, 1])
    metrics['full']['midpoint_cover_depth_m'] = cover_depth.tolist()
    metrics['full']['contrast_peak_vs_local_cover_depth_correlation'] = corr

    arrays = {'time_ns': t_ns, 'midpoint_m': midpoints, 'frequency_Hz': spectra['high_rough_61'].frequency,
              'rough_spectrum': rough, 'halfspace_spectrum': halfspace, 'diff_spectrum': diff,
              'rough_signed': prod['rough'].real_bandpass, 'halfspace_signed': prod['halfspace'].real_bandpass,
              'diff_signed': prod['diff'].real_bandpass,
              'diff_complex_envelope': prod['diff'].complex_envelope,
              'low_diff_signed': prod['low_diff'].real_bandpass,
              'low_halfspace_signed': prod['low_halfspace'].real_bandpass}
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out / 'budget_arrays.npz', **arrays)
    summary = {
        'status': 'COMPLETED_VISIBILITY_BUDGET_NOT_DETECTABILITY_CLAIM',
        'code_sha256': sha256(__file__),
        'gprmax_import_mode': import_mode,
        'inputs_verified_against': ['joint_grid_remaining/manifest.json',
                                    'joint_grid_centre/manifest.json',
                                    'height_wavefield_continuation/independent_verification.json'],
        'archived_baseline_anchor_max_abs_spectrum_difference': anchor_abs,
        'archived_baseline_anchor_relative_L2': anchor_rel_l2,
        'processing': 'official 501 x 20-170 MHz, 200 ns native tail taper, Hann, zero_pad 8; '
                      'complex pair difference (rough - full-cover) before reconstruction; '
                      'windowed complex-envelope L2 norms; low windows use fixed air-delay shift only',
        'midpoint_m': midpoints.tolist(),
        'metrics': metrics,
        'scope': 'Windowed waveform-norm budget for the actual 15 m / 2.5 cm HS4T2D condition plus the '
                 '2 m centre-station control. Ratios are not power percentages, SNR, or detectability '
                 'claims. Common/shape split is a linear cross-trace decomposition of signed traces. '
                 'No imaging, no per-trace normalization, no new solver runs.',
    }
    (args.out / 'summary.json').write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    plot(args.out, arrays, summary)
    print(json.dumps({'anchor_relative_L2': anchor_rel_l2,
                      'full_contrast_over_background_dB_median':
                          float(np.median(metrics['full']['interface_contrast_over_nointerface_background_dB'])),
                      'diff_common_share': metrics['full']['diff_common_mode_energy_share'],
                      'rough_common_share': metrics['full']['rough_common_mode_energy_share'],
                      'low_centre_contrast_over_background_dB_full':
                          metrics['full']['low_centre_interface_contrast_over_nointerface_background_dB']},
                     indent=2))


def plot(out, a, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Microsoft YaHei'
    t, x = a['time_ns'], a['midpoint_m']
    mask = (t >= 150) & (t <= 230)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    limit = max(float(np.max(np.abs(a[k][mask]))) for k in ('rough_signed', 'halfspace_signed', 'diff_signed'))
    for ax, k, label in zip(axes, ('rough_signed', 'halfspace_signed', 'diff_signed'),
                            ['起伏模型总响应', '全覆盖层（无界面）背景', '界面对比响应（起伏−全覆盖层）']):
        im = ax.pcolormesh(x, t[mask], a[k][mask], cmap='RdBu_r', vmin=-limit, vmax=limit, shading='nearest')
        ax.invert_yaxis()
        ax.set_title(label)
        ax.set_xlabel('收发中点原X（m）')
        ax.set_ylabel('时间（ns）')
    fig.colorbar(im, ax=axes.tolist(), label='带符号响应，三图共享线性色标')
    fig.suptitle('15m / 2.5cm / 20–170MHz：界面回波与同窗无界面背景同尺度对照（非SNR认证）')
    fig.savefig(out / 'echo_vs_background_bscan.png', dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    m = s['metrics']['full']
    ax.plot(x, m['interface_contrast_over_nointerface_background_dB'], 'o-',
            label='界面对比 / 无界面背景')
    ax.plot(x, m['interface_contrast_over_total_response_dB'], 's-',
            label='界面对比 / 总响应')
    ax.plot(x, m['nointerface_background_over_total_response_dB'], '^-',
            label='无界面背景 / 总响应')
    ax.axhline(0, color='k', lw=0.7, ls='--')
    ax.set_xlabel('收发中点原X（m）')
    ax.set_ylabel('160–220 ns 窗包络范数比（dB）')
    ax.legend()
    ax.set_title('界面回波可见性预算（15m实际条件；比值非能量占比）')
    fig.savefig(out / 'visibility_budget_db.png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
