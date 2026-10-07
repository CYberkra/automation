"""CPU-only attribution to native time intervals and known-echo injection on real spectra."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from audit_line9_postprocessing import FREQ, inverse, weights, plot_panel, sha, save, rel


def partition(time_ns, width):
    transitions = []
    for cut in [120., 240., 400., 600.]:
        x = (time_ns - cut) / width
        s = np.zeros(len(x))
        s[x >= .5] = 1.
        inside = abs(x) < .5
        s[inside] = .5 * (1 + np.sin(np.pi * x[inside]))
        transitions.append(s)
    gates = np.array([1 - transitions[0]] +
                     [transitions[i] - transitions[i + 1] for i in range(3)] +
                     [transitions[-1]])
    assert np.max(abs(gates.sum(axis=0) - 1)) < 1e-15
    assert gates.min() >= 0
    return gates


def correlation(a, b):
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(abs(np.vdot(a, b)) / den) if den > 0 else None


def main(root, review, cache, out):
    if out.exists():
        raise ValueError('Use a fresh output directory')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    labels = ['0–120ns', '120–240ns', '240–400ns', '400–600ns', '600ns–记录末']
    findings = []
    for audit_path in sorted(review.glob('*_audit.json')):
        audit = json.loads(audit_path.read_text(encoding='utf-8'))
        rows = audit['records']
        p = root / audit['package']
        cache_file = cache / (audit['package'] + '.npz')
        with np.load(cache_file) as z:
            response = z['response'].copy()
            np.testing.assert_array_equal(z['frequency_Hz'], FREQ)
            np.testing.assert_array_equal(z['ids'], [r['id'] for r in rows])
        raw = []
        for r in rows:
            f = p / 'cases' / r['id'] / 'profile.h5'
            if sha(f) != r['native_sha256']:
                raise ValueError('Native identity changed')
            with h5py.File(f) as h:
                raw.append(h['rxs/rx1/Ez'][:])
                dt = float(h.attrs['dt'])
                origin = float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
                if float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset']) != 0:
                    raise ValueError('Receiver origin mismatch')
        raw = np.column_stack(raw)
        native_time = np.arange(len(raw)) * dt
        source = dt * 40 * np.exp(-2j * np.pi * FREQ * origin)
        current = sf.engineering_dft(raw, dt, FREQ) / source[:, None] / .025
        assert rel(current, response) < 1e-12
        gaussian, time = inverse(response, FREQ, weights('gaussian', 501))
        hann, _ = inverse(response, FREQ, weights('hann', 501))
        reference = float(np.max(abs(gaussian)))
        positions = np.array([r['chainage_m'] for r in rows])
        base_centers = np.array([r['geometry']['basal_sand']['time95_ns'] for r in rows])
        cover_centers = np.array([r['geometry']['cover_base']['time95_ns'] for r in rows])
        masks = {'cover_pm10ns': abs(time[:, None] * 1e9 - cover_centers[None, :]) <= 10,
                 'basal_pm10ns': abs(time[:, None] * 1e9 - base_centers[None, :]) <= 10,
                 'old_320_497ns': np.broadcast_to(((time >= 320e-9) & (time <= 497e-9))[:, None], gaussian.shape)}
        report = dict(package=audit['package'], native_count=len(rows), audit_sha256=sha(audit_path),
                      private_cache_sha256=sha(cache_file), full_spectrum_cache_relative_L2=rel(current, response), partitions={})
        shown = None
        for width in [10., 20.]:
            gates = partition(native_time * 1e9, width)
            spectra = [sf.engineering_dft(raw * g[:, None], dt, FREQ) / source[:, None] / .025 for g in gates]
            closure = rel(sum(spectra), response)
            assert closure < 1e-12
            result = dict(spectrum_closure_relative_L2=closure, windows={})
            for name, total in [('gaussian', gaussian), ('hann', hann)]:
                pieces = [inverse(s, FREQ, weights(name, 501))[0] for s in spectra]
                closure_time = rel(sum(pieces), total)
                assert closure_time < 1e-12
                numbers = []
                for label, v in zip(labels, pieces):
                    scopes = {}
                    for key, m in masks.items():
                        values = v[m]
                        scopes[key] = dict(component_over_full_complex_L2=float(np.linalg.norm(values) / np.linalg.norm(total[m])),
                                           component_full_complex_correlation=correlation(values, total[m]),
                                           component_RMS_over_gaussian_global_reference_dB=float(20 * np.log10(np.linalg.norm(values) / np.sqrt(m.sum()) / reference)))
                    numbers.append(dict(native_interval=label, scopes=scopes))
                result['windows'][name] = dict(time_closure_relative_L2=closure_time, components=numbers)
                if width == 20:
                    if shown is None:
                        shown = {}
                    shown[name] = pieces
            report['partitions'][str(int(width))] = result
        fig, axes = plt.subplots(3, 2, figsize=(16, 12), layout='constrained')
        for col, name in enumerate(['gaussian', 'hann']):
            total = gaussian if name == 'gaussian' else hann
            for row, value, label in [(0, total, '全部原始记录'),
                                      (1, shown[name][1], '仅原始120–240ns分量'),
                                      (2, shown[name][2], '仅原始240–400ns分量')]:
                pic = plot_panel(axes[row, col], value, time, positions, rows, reference, name + '：' + label)
                fig.colorbar(pic, ax=axes[row, col], label='复幅度dB / 同一Gaussian峰值')
        fig.suptitle(audit['package'] + '\n原始时间来源分割：交界20ns平滑重叠、全部5分量复数相加闭合；分量不是纯地层，非检出证据')
        fig.savefig(out / (audit['package'] + '_native_time_origin.png'), dpi=125)
        plt.close(fig)
        # A known inclined echo is ADDED to actual spectra, with both opposite phases.
        # It is a processing sensitivity fixture, never a claimed recovered target.
        injections = []
        for level in [-70., -80., -90.]:
            amplitude = reference * 10 ** (level / 20)
            injected_response = amplitude * np.exp(-2j * np.pi * FREQ[:, None] * base_centers[None, :] * 1e-9)
            entry = dict(known_echo_amplitude_dB_over_gaussian_reference=level, windows={})
            for name, total in [('gaussian', gaussian), ('hann', hann)]:
                isolated, _ = inverse(injected_response, FREQ, weights(name, 501))
                plus, _ = inverse(response + injected_response, FREQ, weights(name, 501))
                minus, _ = inverse(response - injected_response, FREQ, weights(name, 501))
                recovery = max(rel(plus - total, isolated), rel(total - minus, isolated))
                assert recovery < 1e-9
                peak = time[np.argmax(abs(isolated), axis=0)] * 1e9
                max_delay_error = float(np.max(abs(peak - base_centers)))
                assert max_delay_error <= time[1] * 1e9 / 2
                m = masks['basal_pm10ns']
                entry['windows'][name] = dict(known_difference_relative_L2=recovery, max_known_delay_error_ns=max_delay_error,
                                              actual_total_over_known_echo_complex_L2=float(np.linalg.norm(total[m]) / np.linalg.norm(isolated[m])))
            injections.append(entry)
        report['known_echo_injections'] = injections
        if 'pkg3' in audit['package']:
            amplitude = reference * 1e-4
            inject = amplitude * np.exp(-2j * np.pi * FREQ[:, None] * base_centers[None, :] * 1e-9)
            fig, axes = plt.subplots(3, 2, figsize=(15, 11), layout='constrained')
            for col, name in enumerate(['gaussian', 'hann']):
                for row, r, title in [(0, response, '原始真实数据'),
                                      (1, response + inject, '真实数据＋人为−80dB回波'),
                                      (2, inject, '仅已知人为回波（非仿真检出）')]:
                    value, _ = inverse(r, FREQ, weights(name, 501))
                    pic = plot_panel(axes[row, col], value, time, positions, rows, reference, name + '：' + title)
                    axes[row, col].set_ylim(420, 140)
                    fig.colorbar(pic, ax=axes[row, col], label='复幅度dB / 固定原数据Gaussian峰值')
            fig.suptitle('pkg3实际70道上的已知回波注入检验：共享绝对色标，无AGC/叠道\n沿H5近似层位人为加信号，只检验处理链能否保留形态；中/下行不得当真实砂岩结果')
            fig.savefig(out / 'pkg3_known_echo_injection_on_actual_data.png', dpi=130)
            plt.close(fig)
        findings.append(report)
        print(json.dumps(report, ensure_ascii=False), flush=True)
    save(out / 'native_time_origin_report.json', dict(calls_solver=False, calls_training=False, script_sha256=sha(__file__),
        processing_helper_sha256=sha(Path(__file__).with_name('audit_line9_postprocessing.py')), findings=findings,
        limits='Time partitions are linear diagnostic components with spectral changes at their transitions, not physical target ablations. L2 ratios can exceed1 by coherent cancellation. Injected echoes are known artificial signals, not real geology or FDTD results. True-interface timing uses approximate95MHz phase-ray curves.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'cache', 'out']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.review, args.cache, args.out)
