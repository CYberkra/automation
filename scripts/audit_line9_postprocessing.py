"""CPU-only processing audit of three immutable Line9 packages; never invokes FDTD."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from review_line9_result_packages import FREQ, save, sha, rel


def weights(name, count):
    if name == 'cosine_edges_5pct':
        x = np.linspace(0., 1., count)
        w = np.ones(count)
        left, right = x < .05, x > .95
        w[left] = .5 * (1 - np.cos(np.pi * x[left] / .05))
        w[right] = .5 * (1 - np.cos(np.pi * (1 - x[right]) / .05))
        return w / w.mean()
    return sf.spectral_window(name, count)


def inverse(response, frequency, weight, pad=8):
    count = len(frequency)
    padded = np.zeros((count * pad,) + response.shape[1:], complex)
    expand = (-1,) + (1,) * (response.ndim - 1)
    padded[:count] = response * weight.reshape(expand)
    time = np.arange(count * pad) / (count * pad * (frequency[1] - frequency[0]))
    baseband = np.fft.ifft(padded, axis=0) * pad
    carrier = np.exp(2j * np.pi * frequency[0] * time).reshape(expand)
    return baseband * carrier, time


def self_checks():
    results = {}
    delay = 280e-9
    response = np.exp(-2j * np.pi * FREQ * delay)
    for name in ['gaussian', 'hann', 'cosine_edges_5pct']:
        w = weights(name, len(FREQ))
        v, t = inverse(response, FREQ, w)
        direct = np.exp(2j * np.pi * t[::97, None] * FREQ) @ (response * w) / len(FREQ)
        error = rel(v[::97], direct)
        assert error < 1e-10
        peak = float(t[np.argmax(abs(v))])
        assert abs(peak - delay) <= t[1] / 2
        # Zero padding changes sample spacing, not inverse amplitude convention.
        v1, _ = inverse(response, FREQ, w, pad=1)
        assert rel(v[::8], v1) < 1e-12
        results[name] = dict(direct_sum_relative_L2=error, peak_ns=peak * 1e9,
                             padding_shared_sample_relative_L2=rel(v[::8], v1))
    selected = FREQ <= 40e6
    f = FREQ[selected]
    h = response[selected]
    w = weights('gaussian', len(f))
    correct, time = inverse(h, f, w)
    buf = np.zeros(len(f) * 8, complex)
    buf[:len(f)] = h * w
    legacy = np.fft.ifft(buf) * len(buf)
    assert np.allclose(abs(legacy), abs(correct) * len(f), rtol=1e-10, atol=1e-12)
    results['low_band_legacy_amplitude_factor'] = len(f)
    results['low_band_legacy_amplitude_offset_dB'] = float(20 * np.log10(len(f)))
    # A real, perfectly horizontal layer is present at every spatial position.
    flat_layer = np.broadcast_to(abs(correct)[:, None], (len(time), 20)).copy()
    ratio = flat_layer / np.maximum(flat_layer, 1e-30)
    assert np.max(abs(20 * np.log10(ratio))) < 1e-10
    results['horizontal_layer_magnitude_division_display_dB'] = [0., 0.]
    assert np.max(abs(flat_layer - flat_layer.mean(axis=1, keepdims=True))) < 1e-14
    results['horizontal_layer_mean_subtraction_erases_layer'] = True
    # Magnitude subtraction cannot recover a complex primary for both echo phases.
    results['magnitude_demultiple_phase_counterexample'] = {
        'primary': 1., 'multiple_amplitude': .5,
        'opposite_phase_total_magnitude': abs(1 - .5),
        'subtract_multiple_magnitude_then_clip': max(abs(1 - .5) - .5, 0),
        'true_primary_magnitude': 1.}
    return results


def plot_panel(axes, values, time, pos, records, reference, title, signed=False):
    order = np.argsort(pos)
    edges = np.r_[pos[order][0] - .1, (pos[order][:-1] + pos[order][1:]) / 2, pos[order][-1] + .1]
    keep = time <= 600e-9
    t = time[keep] * 1e9
    te = np.r_[t[0] - (t[1] - t[0]) / 2, (t[:-1] + t[1:]) / 2, t[-1] + (t[-1] - t[-2]) / 2]
    if signed:
        # Fixed display clipping only, always preserved separately from complex data.
        data = -2 * values[keep][:, order].real / (2 * reference * 10 ** (-55 / 20))
        pic = axes.pcolormesh(edges, te, data, cmap='gray', vmin=-1, vmax=1, shading='flat', rasterized=True)
    else:
        data = 20 * np.log10(np.maximum(abs(values[keep][:, order]) / reference, 1e-12))
        pic = axes.pcolormesh(edges, te, data, cmap='gray_r', vmin=-110, vmax=-35, shading='flat', rasterized=True)
    for key, colour, label in [('cover_base', '#00a6d6', '覆盖层底'),
                              ('first_sand', '#c25be5', '首个砂岩顶'),
                              ('basal_sand', '#1ec766', '底部砂岩顶')]:
        curve = np.array([r['geometry'][key]['time95_ns'] for r in records])
        axes.plot(pos[order], curve[order], '--', color=colour, lw=.9, label=label)
    axes.set(ylim=(600, 0), xlabel='测线里程 / m（仅完成区）', ylabel='时间 / ns', title=title)
    axes.invert_xaxis()
    axes.legend(fontsize=7, loc='lower right')
    return pic


def main(root, review, out, cache):
    if out.exists() or cache.exists():
        raise ValueError('Use fresh output and private cache directories')
    out.mkdir(parents=True)
    cache.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    report = dict(calls_solver=False, calls_training=False, script_sha256=sha(__file__),
                  official_processing_sha256=sha(sf.__file__), self_checks=self_checks(), packages=[])
    for audit_path in sorted(review.glob('*_audit.json')):
        audit = json.loads(audit_path.read_text(encoding='utf-8'))
        path = root / audit['package']
        records = audit['records']
        raw, stored, native_peak_gates = [], [], []
        for record in records:
            native = path / 'cases' / record['id'] / 'profile.h5'
            if sha(native) != record['native_sha256']:
                raise ValueError('Native input changed')
            with h5py.File(native) as h:
                raw.append(h['rxs/rx1/Ez'][:])
                dt = float(h.attrs['dt'])
                source_offset = float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
                if float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset']) != 0:
                    raise ValueError('Receiver origin needs explicit accounting')
            file = path / 'sfcw' / (record['id'] + '_sfcw.h5')
            if sha(file) != record['sfcw_sha256']:
                raise ValueError('Stored SFCW input changed')
            with h5py.File(file) as h:
                np.testing.assert_array_equal(h['frequency'][:], FREQ)
                stored.append(h['response'][:])
        raw = np.column_stack(raw)
        source = dt * 40 * np.exp(-2j * np.pi * FREQ * source_offset)
        response = sf.engineering_dft(raw, dt, FREQ) / source[:, None] / .025
        full, time = inverse(response, FREQ, weights('gaussian', 501))
        old, _ = inverse(np.column_stack(stored) / .025, FREQ, weights('gaussian', 501))
        hann, _ = inverse(response, FREQ, weights('hann', 501))
        edges, _ = inverse(response, FREQ, weights('cosine_edges_5pct', 501))
        native_time = np.arange(len(raw)) * dt
        early_gate = np.ones(len(raw))
        early_gate[native_time >= 120e-9] = 0
        ramp = (native_time > 100e-9) & (native_time < 120e-9)
        early_gate[ramp] = .5 * (1 + np.cos(np.pi * (native_time[ramp] - 100e-9) / 20e-9))
        late_response = sf.engineering_dft(raw * (1 - early_gate[:, None]), dt, FREQ) / source[:, None] / .025
        late, _ = inverse(late_response, FREQ, weights('gaussian', 501))
        mean_removed, _ = inverse(response - response.mean(axis=1, keepdims=True), FREQ, weights('gaussian', 501))
        variants = [('original_gaussian', old, '原包Gaussian：接收尾200ns渐消'),
                    ('no_tail_gaussian', full, 'Gaussian：取消尾渐消，原始全记录'),
                    ('no_tail_hann', hann, 'Hann：取消尾渐消，压低频带端点'),
                    ('no_tail_cosine_edges_5pct', edges, '探索对照：两端各5%余弦渐消，取消尾渐消'),
                    ('late_only_gaussian', late, '诊断：先平滑移除0–120ns早时场，再Gaussian'),
                    ('complex_mean_removed_gaussian', mean_removed, '风险对照：复频谱均值道扣除，会删除共同地层')]
        reference = float(np.max(abs(full)))
        positions = np.array([r['chainage_m'] for r in records])
        summaries = {}
        for name, v, _ in variants:
            gates = {}
            for key in ['cover_base', 'basal_sand']:
                centers = np.array([r['geometry'][key]['time95_ns'] for r in records])
                gate = abs(time[:, None] * 1e9 - centers[None, :]) <= 10
                gates[key] = dict(complex_relative_change_from_no_tail_gaussian=rel(v[gate], full[gate]),
                                  amplitude_RMS_over_global_reference_dB=float(20 * np.log10(np.linalg.norm(v[gate]) / np.sqrt(gate.sum()) / reference)))
            summaries[name] = gates
        count = len(records) // 5 * 5
        stack = full[:, :count].reshape(len(time), -1, 5).mean(axis=2)
        average_mag = abs(full[:, :count]).reshape(len(time), -1, 5).mean(axis=2)
        stack_metrics = {}
        for key in ['cover_base', 'basal_sand']:
            centers = np.array([r['geometry'][key]['time95_ns'] for r in records[:count]]).reshape(-1, 5).mean(axis=1)
            gate = abs(time[:, None] * 1e9 - centers[None, :]) <= 10
            attenuation = 20 * np.log10(np.maximum(abs(stack[gate]) / np.maximum(average_mag[gate], 1e-30), 1e-15))
            stack_metrics[key] = dict(coherent_over_mean_magnitude_dB_percentile10_50_90=np.percentile(attenuation, [10, 50, 90]).tolist())
        raw_reference = float(np.max(abs(raw)))
        centers = np.array([r['geometry']['basal_sand']['time95_ns'] for r in records])
        raw_gate = abs(native_time[:, None] * 1e9 - centers[None, :]) <= 10
        raw_peak_gates = np.max(np.where(raw_gate, abs(raw), 0), axis=0)
        item = dict(package=audit['package'], audit_sha256=sha(audit_path), actual_traces=len(records),
                    full_frequency_count=501, sample_time_ns=float(time[1] * 1e9), period_ns=float(1e9 / 300000),
                    tail_taper_starts_at_nominal_ns=1000 if 'pkg1' in audit['package'] else 600,
                    processing_variants=summaries, group5_dropped_station_count=len(records) - count,
                    group5_stack_metrics=stack_metrics,
                    native_bipolar_4percent_clip_basal_peak_as_fraction_of_colour_full_scale_percentile10_50_90=np.percentile(raw_peak_gates / raw_reference / .04, [10, 50, 90]).tolist(),
                    private_complex_cache_file=cache.name + '/' + audit['package'] + '.npz')
        np.savez_compressed(cache / (audit['package'] + '.npz'), frequency_Hz=FREQ, response=response,
                            response_early_removed_diagnostic=late_response, time_s=time, chainage_m=positions,
                            ids=np.array([r['id'] for r in records]), global_reference=reference)
        fig, axes = plt.subplots(3, 2, figsize=(16, 13), layout='constrained')
        for ax, (_, value, title) in zip(axes.flat, variants):
            pic = plot_panel(ax, value, time, positions, records, reference, title)
            fig.colorbar(pic, ax=ax, label='复幅度dB / 同一Gaussian全记录峰值')
        fig.suptitle(audit['package'] + '\n20–170MHz/0.3MHz/501点：同一绝对参考、无AGC/叠道/中值滤波；彩线为实际H5的95MHz近似，非检出')
        fig.savefig(out / (audit['package'] + '_processing_ablation.png'), dpi=125)
        plt.close(fig)
        fig, axes = plt.subplots(2, 1, figsize=(13, 8), layout='constrained')
        for ax, value, title in zip(axes, [full, hann], ['Gaussian', 'Hann']):
            pic = plot_panel(ax, value, time, positions, records, reference,
                             title + '：正确实带通2Re[复基带×20MHz载波]，未AGC/叠道', signed=True)
            fig.colorbar(pic, ax=ax, label='显示刻度；零=灰、正=黑、负=白')
        fig.suptitle(audit['package'] + '\n双极显示对照：共同峰值参考，±−55dB幅度饱和，仅显示裁剪；保留原复数数据')
        fig.savefig(out / (audit['package'] + '_signed_sfcw.png'), dpi=130)
        plt.close(fig)
        report['packages'].append(item)
        print(json.dumps(item, ensure_ascii=False), flush=True)
    report['limits'] = 'Processing diagnostics only. No matched target/background, no true SNR or interface recovery certification. Early removal and mean subtraction may remove genuine responses; frequency weights are exploratory, not verified device settings.'
    save(out / 'postprocessing_audit.json', report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'out', 'cache']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.review, args.out, args.cache)
