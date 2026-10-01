"""Read-only field comparison for the v0.6 ringing review; no fitting for deployment.

CSV-derived output MUST stay in ignored local_checks, outside the source directory.
Reported MHz/ns use the provisional endpoint-inclusive CSV time convention.
Fixed-template R2 measures agreement with a damped sine, not noise identification.
The synthetic probe executes selected original v0.6 functions without solver imports.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1]
WINS = [(25, 45), (45, 70), (70, 100), (100, 140), (140, 200)]


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def quant(x):
    return {f'p{p}': float(np.percentile(x, p)) for p in (10, 25, 50, 75, 90)}


def fixed_fit(seg, tt):
    B = np.column_stack([np.exp(-tt / 73) * np.sin(2 * np.pi * .107 * tt),
                         np.exp(-tt / 73) * np.cos(2 * np.pi * .107 * tt)])
    coef = seg @ np.linalg.pinv(B).T
    pred = coef @ B.T
    den = np.sum((seg - seg.mean(axis=1, keepdims=True)) ** 2, axis=1)
    r2 = 1 - np.sum((seg - pred) ** 2, axis=1) / den
    # Independent scalar least-squares check on deterministic first/middle/last rows.
    for i in (0, len(seg) // 2, len(seg) - 1):
        c = np.linalg.lstsq(B, seg[i], rcond=None)[0]
        ref = 1 - np.sum((seg[i] - B @ c) ** 2) / den[i]
        assert abs(ref - r2[i]) < 1e-10
    return r2


def inspect_file(p, divisor='N-1'):
    with p.open('r', encoding='utf-8-sig') as f:
        hdr = [f.readline() for _ in range(4)]
    vals = [float(s.split('=')[1].split(',')[0]) for s in hdr]
    ns, tw, nr, dx = int(vals[0]), vals[1], int(vals[2]), vals[3]
    raw = np.loadtxt(p, delimiter=',', skiprows=4, usecols=(3,), dtype=float)
    assert len(raw) == ns * nr, (p.name, len(raw), ns * nr)
    assert np.isfinite(raw).all(), p.name
    xall = raw.reshape(nr, ns)
    dt = tw / (ns - 1 if divisor == 'N-1' else ns)
    tt = np.arange(ns) * dt
    i0all = np.argmax(np.abs(xall[:, tt <= 60]), axis=1)
    pkall = np.abs(xall[np.arange(nr), i0all])
    valid = pkall > 0
    x, i0, pk = xall[valid], i0all[valid], pkall[valid]
    assert len(x) > 1
    t0 = tt[i0]
    normalized = x / pk[:, None]
    sos = butter(4, [.090, .125], btype='band', fs=1 / dt, output='sos')
    env = np.abs(hilbert(sosfiltfilt(sos, normalized, axis=1), axis=1))
    rel = tt[None, :] - t0[:, None]
    prof = np.column_stack([
        np.nanmedian(np.where((rel >= lo) & (rel <= hi), env, np.nan), axis=1)
        for lo, hi in WINS])
    q = prof[:, 1]
    fit = {}
    spectrum = None
    for lo, hi in [(25, 100), (25, 220)]:
        offs = np.arange(int(np.ceil(lo / dt)), int(np.floor(hi / dt)) + 1)
        seg = normalized[np.arange(len(x))[:, None], i0[:, None] + offs]
        fit[f'{lo}_{hi}'] = quant(fixed_fit(seg, offs * dt))
        if hi == 220:
            ft = np.fft.rfft((seg - seg.mean(axis=1, keepdims=True)) * np.hanning(len(offs)), n=2048, axis=1)
            freqs = np.fft.rfftfreq(2048, dt) * 1000
            power = np.abs(ft) ** 2
            band = (freqs >= 20) & (freqs <= 170)
            narrow = (freqs >= 90) & (freqs <= 125)
            dom = freqs[band][np.argmax(power[:, band], axis=1)]
            fraction = power[:, narrow].sum(axis=1) / power[:, band].sum(axis=1)
            spectrum = {'dominant_mhz': quant(dom), 'energy_90_125_over_20_170': quant(fraction)}
    # Full records, same absolute early window as historical full-line texture script.
    seg = xall[:, (tt >= 30) & (tt <= 100)]
    seg -= seg.mean(axis=1, keepdims=True)
    den = np.linalg.norm(seg[:-1], axis=1) * np.linalg.norm(seg[1:], axis=1)
    cc = np.sum(seg[:-1] * seg[1:], axis=1)[den > 0] / den[den > 0]
    # Exact historical lateral proxy: Hilbert of each finite relative 40..100 segment.
    e = []
    for row, z, t in zip(normalized, rel, t0):
        m = (z >= 40) & (z <= 100)
        e.append(np.abs(hilbert(row[m])).mean())
    e = np.asarray(e)
    lat = float(np.sqrt(np.mean(np.diff(e / e.mean()) ** 2)))
    blocks = [float(np.median(q[j:j+33])) for j in range(0, len(q)-32, 33)]
    return {'file': p.name, 'n_traces': nr, 'n_samples': ns, 'time_window_ns': tw,
            'nominal_trace_interval_m': dx, 'dt_ns_assumed': dt, 'nonzero_traces': int(valid.sum()),
            'direct_peak': quant(pk), 'early_peak_t0_ns': quant(t0), 'q': quant(q),
            'q_33trace_block_medians': quant(blocks), 'profile_bandpass_medians': np.median(prof, axis=0).tolist(),
            'fixed_107mhz_73ns_R2': fit, 'spectrum_mixed_window': spectrum,
            'adjacent_corr_30_100ns': quant(cc), 'lateral_diff_rms': lat}, normalized, tt, prof


def original_namespace():
    names = ['augment_sim2real_v0_3.py', 'augment_sim2real_v0_5.py', 'augment_sim2real_v0_6.py']
    wanted = {'smooth_series', 'op2_surface_edit', 't0_index', 'win_mask', '_sos', 'q_of_trace',
              'onset_gate', 'op1_ringing_v06'}
    scope = {'np': np, 'TAU': 73., 'F0': .107, 'BP_LO': .09, 'BP_HI': .125,
             'QW_LO': 45., 'QW_HI': 70., 'ONSET0': 12., 'ONSET1': 22.,
             'JITTER_STD': .16, 'DELTA_NS': 3.5,
             'K_LO': .20, 'K_HI': .64, 'JITTER_NS': 15.}
    hashes = {}
    for name in names:
        p = ROOT / 'scripts' / name
        hashes[name] = digest(p)
        tree = ast.parse(p.read_text(encoding='utf-8-sig'))
        nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), 'exec'), scope)
    assert wanted.issubset(scope)
    return scope, hashes


def probe(scope, qlo, qhi):
    t = np.arange(0, 700, .832)
    # Mechanism fixture, explicitly not replay of archived physical BG simulations.
    a = np.exp(-.5 * ((t - 15) / 5) ** 2) * np.cos(2 * np.pi * .095 * (t - 15))
    late = .02 * np.exp(-.5 * ((t - 240) / 5) ** 2) * np.cos(2 * np.pi * .095 * (t - 240))
    x = np.tile(a, (33, 1))
    y, clip = scope['op1_ringing_v06'](x, t, np.random.default_rng(20261001), qlo, qhi)
    yp, _ = scope['op1_ringing_v06'](x + late, t, np.random.default_rng(20261001), qlo, qhi)
    before = t <= 27
    m = (t >= 220) & (t <= 290)
    diff = yp - y
    assert np.max(np.abs(y[:, before] - x[:, before])) == 0
    surface = .14 * np.exp(-.5 * ((t - 100) / 5) ** 2) * np.cos(2 * np.pi * .095 * (t - 100))
    shallow = .02 * np.exp(-.5 * ((t - 130) / 5) ** 2) * np.cos(2 * np.pi * .095 * (t - 130))
    base = x + surface
    ed0 = scope['op2_surface_edit'](base, t, np.random.default_rng(20261001))
    ed1 = scope['op2_surface_edit'](base + shallow, t, np.random.default_rng(20261001))
    shallow_ref = np.tile(shallow, (33, 1))
    ed_diff = ed1 - ed0
    removed = np.linalg.norm(ed_diff - shallow_ref) / np.linalg.norm(shallow_ref)
    return {'label': 'synthetic mechanism probe, NOT archived simulation replay',
            'clip_count': clip, 'late_event_relative_change_l2': float(np.linalg.norm(diff[:, m] - late[m]) / np.linalg.norm(np.tile(late[m], (33, 1)))),
            'q_output': quant([scope['q_of_trace'](row, t, scope['_sos'](.832)) for row in y]),
            'early_change_max': float(np.max(np.abs(y[:, before] - x[:, before]))),
            'surface_editor_shallow_event_relative_change_l2': float(removed),
            'surface_editor_shallow_event_energy_norm_retained': float(np.linalg.norm(ed_diff) / np.linalg.norm(shallow_ref))}


def make_plot(records, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    lines = ['Line3', 'Line6', 'Line7', 'Line9', 'LineL1', 'LineX1']
    colors = ['#6f8093', '#c37426', '#2666ac']
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.3))
    for j, power in enumerate([20, 30, 36]):
        rr = [next(r for r in records if r['file'] == f'{line}origin({power}).csv') for line in lines]
        xs = np.arange(6) + (j-1) * .23
        yy = np.array([r['q']['p50'] for r in rr])
        err = np.array([[r['q']['p50']-r['q']['p25'] for r in rr], [r['q']['p75']-r['q']['p50'] for r in rr]])
        ax[0].errorbar(xs, yy, yerr=err, fmt='o', capsize=3, color=colors[j], label=f'{power} dBm')
        ax[1].plot(xs, [r['fixed_107mhz_73ns_R2']['25_220']['p50'] for r in rr], 'o', color=colors[j], label=f'{power} dBm')
    l9 = next(r for r in records if r['file'] == 'Line9origin(36).csv')
    ax[0].axhspan(l9['q']['p25'], l9['q']['p75'], color='#2666ac', alpha=.10, label='v0.6 目标采样区间')
    for a in ax[:2]:
        a.set_xticks(np.arange(6), lines)
        a.grid(axis='y', alpha=.25)
        a.legend(fontsize=8)
    ax[0].set_ylabel('q：90–125 MHz 包络 / 早时峰（中位与四分位）')
    ax[0].set_title('单用 Line9(36) 的 q 区间覆盖有限')
    ax[1].set_ylabel('单衰减正弦拟合 R² 中位数')
    ax[1].set_ylim(0, 1)
    ax[1].set_title('固定 107 MHz、73 ns：t0+25–220 ns')
    centers = [(lo+hi)/2 for lo, hi in WINS]
    for line in lines:
        r = next(r for r in records if r['file'] == f'{line}origin(36).csv')
        ax[2].plot(centers, r['profile_bandpass_medians'], 'o-', label=line, lw=1.2)
    ax[2].set_xlabel('相对各道早时峰 t0 的窗中心（ns）')
    ax[2].set_ylabel('90–125 MHz 包络 / 早时峰（全线中位）')
    ax[2].set_title('36 dBm 各线包络形状也不一致')
    ax[2].grid(alpha=.25)
    ax[2].legend(fontsize=8)
    fig.suptitle('营山振铃审查：18 份 CSV / 31,355 道；混合早时窗，未分离地质回波；时间轴暂按 T/(N−1)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(out / 'ringing_field_gap.png', dpi=160)
    fig.savefig(out / 'ringing_field_gap.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    args = ap.parse_args()
    source = args.data_dir.resolve()
    out = args.out_dir.resolve()
    assert (ROOT / 'artifacts/local_checks').resolve() in out.parents
    assert source != out and source not in out.parents
    out.mkdir(parents=True, exist_ok=True)
    records, plots = [], {}
    for p in sorted(source.glob('Line*origin(*).csv')):
        rec, x, t, prof = inspect_file(p)
        rec['sha256'] = digest(p)
        records.append(rec)
        if p.name == 'Line9origin(36).csv':
            plots['x'], plots['t'], plots['profile'] = x, t, prof
        print(f"{p.name}: n={rec['n_traces']} q={rec['q']['p50']:.4f} R2={rec['fixed_107mhz_73ns_R2']['25_220']['p50']:.3f} cc={rec['adjacent_corr_30_100ns']['p50']:.3f}", flush=True)
    scope, hashes = original_namespace()
    l9 = next(r for r in records if r['file'] == 'Line9origin(36).csv')
    for i in (0, len(plots['x']) // 2, len(plots['x']) - 1):
        actual = scope['q_of_trace'](plots['x'][i], plots['t'], scope['_sos'](plots['t'][1]))
        assert abs(actual - plots['profile'][i, 1]) < 1e-10
    sensitivity, *_ = inspect_file(source / l9['file'], divisor='N')
    result = {'source_dir': str(source), 'source_read_only': True,
              'solver_executed': False, 'training_executed': False,
              'status': 'audit only; all observed field datasets are disclosed',
              'time_axis': 'provisional T/(N-1), T/N sensitivity on Line9(36); not physical calibration',
              'windows_ns': WINS, 'source_code_sha256': hashes, 'files': records,
              'audit_script_sha256': digest(Path(__file__)),
              'line9_T_over_N_sensitivity': sensitivity,
              'v06_synthetic_probe': probe(scope, l9['q']['p25'], l9['q']['p75'])}
    (out / 'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    np.savez_compressed(out / 'private_plot_arrays.npz', **plots)
    make_plot(records, out)
    print('saved', out, flush=True)


if __name__ == '__main__':
    main()
