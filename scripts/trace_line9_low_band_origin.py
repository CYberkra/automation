"""CPU-only native-time origin audit of the exact 20-39.8MHz subset; no FDTD."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from trace_line9_time_origin import partition, correlation
from audit_line9_postprocessing import FREQ, inverse, weights, plot_panel, sha, save, rel, sf


def main(root, review, cache, out):
    if out.exists():
        raise ValueError('Use a new output directory')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    selected = FREQ <= 40e6
    freq = FREQ[selected]
    findings = []
    for ap in sorted(review.glob('*_audit.json')):
        audit = json.loads(ap.read_text(encoding='utf-8'))
        rows = audit['records']
        raw = []
        for r in rows:
            p = root / audit['package'] / 'cases' / r['id'] / 'profile.h5'
            if sha(p) != r['native_sha256']:
                raise ValueError('Native identity changed')
            with h5py.File(p) as h:
                raw.append(h['rxs/rx1/Ez'][:])
                dt = float(h.attrs['dt'])
                offset = float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
                if float(h['rxs/rx1/Ez'].attrs['TimeSampleOffset']) != 0:
                    raise ValueError('Receiver origin mismatch')
        raw = np.column_stack(raw)
        cache_file = cache / (audit['package'] + '.npz')
        with np.load(cache_file) as z:
            np.testing.assert_array_equal(z['ids'], [r['id'] for r in rows])
            full_response = z['response'].copy()
        reference = float(np.max(abs(inverse(full_response, FREQ, weights('gaussian', 501))[0])))
        source = dt * 40 * np.exp(-2j * np.pi * freq * offset)
        row = dict(package=audit['package'], audit_sha256=sha(ap), private_cache_sha256=sha(cache_file), partitions={})
        shown = None
        for width in [10., 20.]:
            g = partition(np.arange(len(raw)) * dt * 1e9, width)
            spectra = [sf.engineering_dft(raw * x[:, None], dt, freq) / source[:, None] / .025
                       for x in [g[0], g[1], g[2:].sum(axis=0)]]
            closure = rel(sum(spectra), full_response[selected])
            assert closure < 1e-12
            result = dict(spectrum_closure_relative_L2=closure, windows={})
            for name in ['gaussian', 'hann']:
                values = [inverse(s, freq, weights(name, len(freq)))[0] for s in spectra]
                time = inverse(spectra[0], freq, weights(name, len(freq)))[1]
                centers = np.array([r['geometry']['basal_sand']['time95_ns'] for r in rows])
                mask = abs(time[:, None] * 1e9 - centers[None, :]) <= 10
                after120 = values[1] + values[2]
                result['windows'][name] = dict(
                    shallow120_240_over_after120_complex_L2=float(np.linalg.norm(values[1][mask]) / np.linalg.norm(after120[mask])),
                    shallow120_240_with_after120_complex_correlation=correlation(values[1][mask], after120[mask]),
                    after240_over_after120_complex_L2=float(np.linalg.norm(values[2][mask]) / np.linalg.norm(after120[mask])))
                if width == 20 and name == 'gaussian':
                    shown = [sum(values), after120, values[1], values[2]]
            row['partitions'][str(int(width))] = result
        positions = np.array([r['chainage_m'] for r in rows])
        fig, axes = plt.subplots(2, 2, figsize=(14, 9), layout='constrained')
        titles = ['完整原始记录重建', '诊断：移除0–120ns后重建', '仅原始120–240ns分量', '仅原始240ns之后分量']
        for ax, v, title in zip(axes.flat, shown, titles):
            pic = plot_panel(ax, v, time, positions, rows, reference, title)
            ax.set_ylim(500, 120)
            fig.colorbar(pic, ax=ax, label='复幅度dB / 固定全频Gaussian峰值')
        fig.suptitle(audit['package'] + '\n20–39.8MHz/67点/Gaussian：时间分量20ns平滑交界；无AGC/叠道\n分量不是纯目标，低频波包展宽可把较早原始响应延伸到更深显示时间')
        fig.savefig(out / (audit['package'] + '_low_band_origin.png'), dpi=130)
        plt.close(fig)
        findings.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
    save(out / 'low_band_origin_report.json', dict(calls_solver=False, calls_training=False,
        script_sha256=sha(__file__), partition_helper_sha256=sha(Path(__file__).with_name('trace_line9_time_origin.py')),
        frequency_count=67, realized_band_Hz=[float(freq[0]), float(freq[-1])], findings=findings,
        limits='Model-informed approximate basal timing gate. Time partitions alter spectra and are not geological ablations. Correlations are not energy fractions; ratios can exceed1. The120-240ns interval may contain shallow reflections and other paths, not a pure material response.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'cache', 'out']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.review, args.cache, args.out)
