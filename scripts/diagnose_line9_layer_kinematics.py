"""Offline dispersive primary-path comparison; model-informed diagnostics, no solver."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_line9_postprocessing import FREQ, inverse, weights, sha, save, rel
from review_line9_result_packages import C0, indices
from trace_line9_time_origin import correlation


def primary_paths(geometry, materials, frequency=FREQ):
    """Normal-incidence Fresnel primaries, with local phase-ray offset correction.

    Engineering convention: n has negative imaginary part; exp(-i*k*distance)
    attenuates. This is a local layered plane-wave proxy, not a 2D Green function.
    It excludes spreading, slopes, lateral paths, internal multiples and PML.
    """
    n = np.array([indices(materials, f) for f in frequency])
    n95 = indices(materials, 95e6)
    thickness = np.zeros(len(materials))
    thickness[0] = geometry['midpoint_agl_m']
    transmission = np.ones(len(frequency), complex)
    previous_depth = 0.
    spectra = []
    for boundary in geometry['boundaries']:
        a, b = boundary['above'], boundary['below']
        thickness[a] += boundary['depth_m'] - previous_depth
        previous_depth = boundary['depth_m']
        correction = boundary['time95_ns'] * 1e-9 - 2 * np.dot(thickness, n95.real) / C0
        reflection = (n[:, a] - n[:, b]) / (n[:, a] + n[:, b])
        propagation = np.exp(-4j * np.pi * frequency * (n @ thickness) / C0)
        spectra.append(transmission * reflection * propagation *
                       np.exp(-2j * np.pi * frequency * correction))
        transmission *= 4 * n[:, a] * n[:, b] / (n[:, a] + n[:, b]) ** 2
    return spectra


def checks():
    materials = {str(i): {'base': {'relative_permittivity': e,
                                  'electric_conductivity_s_per_m': 0.}}
                 for i, e in enumerate([1., 9., 4.])}
    geom = {'midpoint_agl_m': 2., 'boundaries': [
        {'above': 0, 'below': 1, 'depth_m': 0., 'time95_ns': 4 / C0 * 1e9},
        {'above': 1, 'below': 2, 'depth_m': 7., 'time95_ns': 46 / C0 * 1e9}]}
    h = primary_paths(geom, materials)
    expected = .75 * .2 * np.exp(-2j * np.pi * FREQ * 46 / C0)
    err = rel(h[1], expected)
    assert err < 1e-12
    v, t = inverse(h[1], FREQ, weights('hann', len(FREQ)))
    peak_error = abs(t[np.argmax(abs(v))] - 46 / C0) * 1e9
    assert peak_error <= t[1] * 1e9 / 2
    return {'nondispersive_two_layer_complex_relative_L2': err,
            'nondispersive_peak_error_ns': peak_error}


def rms_db(v, reference):
    return float(20 * np.log10(np.linalg.norm(v) / np.sqrt(v.size) / reference))


def main(root, review, cache, out):
    if out.exists():
        raise ValueError('Fresh diagnostic directory required')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    result = dict(calls_solver=False, calls_training=False, self_checks=checks(),
                  script_sha256=sha(__file__), packages=[])
    for ap in sorted(review.glob('*_audit.json')):
        audit = json.loads(ap.read_text('utf-8'))
        rows = audit['records']
        material_file = next((root / audit['package'] / 'geometries').glob('*.json'))
        if sha(material_file) != audit['materials_sha256']:
            raise ValueError('Materials changed')
        materials = json.loads(material_file.read_text('utf-8'))['materials']
        cp = cache / (audit['package'] + '.npz')
        with np.load(cp) as z:
            np.testing.assert_array_equal(z['frequency_Hz'], FREQ)
            np.testing.assert_array_equal(z['ids'], [r['id'] for r in rows])
            response = z['response'].copy()
            positions = z['chainage_m'].copy()
            reference = float(z['global_reference'])
        paths = [primary_paths(r['geometry'], materials) for r in rows]
        item = dict(package=audit['package'], audit_sha256=sha(ap),
                    cache_sha256=sha(cp), materials_sha256=sha(material_file),
                    native_dtype=audit['native_dtype'], interfaces={})
        curves = {}
        # Named interfaces remain meaningful where an interbed pinches out.
        for key in ['cover_base', 'first_sand', 'basal_sand']:
            indices_in_column = [r['geometry']['boundaries'].index(r['geometry'][key]) for r in rows]
            hp = np.column_stack([p[k] for p, k in zip(paths, indices_in_column)])
            phase_centers = np.array([r['geometry'][key]['time95_ns'] for r in rows])
            entry = {}
            for window in ['gaussian', 'hann']:
                w = weights(window, len(FREQ))
                model, time = inverse(hp, FREQ, w)
                actual, _ = inverse(response, FREQ, w)
                peaks = time[np.argmax(abs(model), axis=0)] * 1e9
                gate = abs(time[:, None] * 1e9 - peaks[None, :]) <= 12
                # Controls change the entire time trajectory without refitting the data.
                flat_shift = np.median(peaks) - peaks
                reverse_shift = peaks[::-1] - peaks
                flat, _ = inverse(hp * np.exp(-2j * np.pi * FREQ[:, None] * flat_shift[None, :] * 1e-9), FREQ, w)
                reverse, _ = inverse(hp * np.exp(-2j * np.pi * FREQ[:, None] * reverse_shift[None, :] * 1e-9), FREQ, w)
                a, b = model[gate], actual[gate]
                # One global complex coefficient, no per-trace phase/delay fitting.
                coefficient = np.vdot(a, b) / np.vdot(a, a)
                peak_in_gate = time[np.argmax(np.where(gate, abs(actual), 0.), axis=0)] * 1e9
                metrics = dict(
                    phase95_curve_ns=phase_centers.tolist(), primary_peak_ns=peaks.tolist(),
                    primary_minus_phase95_ns_percentile0_50_100=np.percentile(peaks - phase_centers, [0, 50, 100]).tolist(),
                    complex_correlation_primary=correlation(a, b),
                    complex_correlation_flat_control=correlation(flat[gate], b),
                    complex_correlation_reversed_control=correlation(reverse[gate], b),
                    one_global_fit_relative_L2=rel(coefficient * a, b),
                    one_global_coefficient=[float(coefficient.real), float(coefficient.imag)],
                    actual_peak_in_informed_gate_minus_primary_ns_percentile0_50_100=np.percentile(peak_in_gate - peaks, [0, 50, 100]).tolist(),
                    actual_RMS_dB_over_gaussian_reference=rms_db(b, reference),
                    search='Geometry-informed +-12ns only; not blind detection; primary trajectory not fitted to data')
                # Guard against confusing the basal contact with the interface above it.
                if key == 'basal_sand':
                    above_hp = np.column_stack([p[k - 1] for p, k in zip(paths, indices_in_column)])
                    above, _ = inverse(above_hp, FREQ, w)
                    above_peak = time[np.argmax(abs(above), axis=0)] * 1e9
                    metrics['previous_interface_peak_separation_ns_range'] = [float(min(peaks - above_peak)), float(max(peaks - above_peak))]
                entry[window] = metrics
                if window == 'hann':
                    curves[key] = peaks
            item['interfaces'][key] = entry
        # Every geological contact, including sandstone-to-mudstone, is drawn.
        all_hann_peaks = []
        for p in paths:
            all_hann_peaks.append([float(time[np.argmax(abs(inverse(h, FREQ, weights('hann', len(FREQ)))[0]))] * 1e9) for h in p])
        if 'pkg3' in audit['package']:
            # Compare competing trajectories at the same informed basal gate.
            centers = curves['basal_sand']
            for name in ['right_side_pml_air_mirror_ns', 'top_pml_air_mirror_ns']:
                x = np.array([r[name] for r in rows])
                item[name + '_range'] = [float(min(x)), float(max(x))]
                item[name + '_endpoint_change_ns'] = float(x[-1] - x[0])
            item['basal_primary_endpoint_change_ns'] = float(centers[-1] - centers[0])
        item['all_interface_hann_peak_ns_by_station'] = all_hann_peaks
        gaussian, time = inverse(response, FREQ, weights('gaussian', len(FREQ)))
        hann, _ = inverse(response, FREQ, weights('hann', len(FREQ)))
        lo = max(0., float(min(curves['cover_base'])) - 25)
        hi = float(max(curves['basal_sand'])) + 40
        keep = (time * 1e9 >= lo) & (time * 1e9 <= hi)
        order = np.argsort(positions)
        xe = np.r_[positions[order][0] - .1, (positions[order][:-1] + positions[order][1:]) / 2, positions[order][-1] + .1]
        tt = time[keep] * 1e9
        te = np.r_[tt[0] - (tt[1]-tt[0])/2, (tt[:-1]+tt[1:])/2, tt[-1]+(tt[1]-tt[0])/2]
        fig, axes = plt.subplots(1, 3, figsize=(18, 7), layout='constrained')
        for ax, value, title in zip(axes, [gaussian, hann, hann], [
            'Gaussian总场：无参考线', 'Hann总场：无参考线', '相同Hann总场＋全部界面近似峰位']):
            db = 20*np.log10(np.maximum(abs(value[keep][:, order])/reference, 1e-12))
            pic = ax.pcolormesh(xe, te, db, cmap='gray_r', vmin=-110, vmax=-45, rasterized=True)
            ax.set(xlabel='原测线里程 / m（仅完成区）', ylabel='时间 / ns', ylim=(hi, lo), title=title)
            ax.invert_xaxis()
            fig.colorbar(pic, ax=ax, label='复幅度dB / 同一Gaussian总场峰值')
        for j in range(max(map(len, all_hann_peaks))):
            ys = np.array([p[j] if j < len(p) else np.nan for p in all_hann_peaks])
            # Boundary indices can change when a bed disappears; avoid linking unlike contacts.
            types = [r['geometry']['boundaries'][j]['below'] if j < len(r['geometry']['boundaries']) else None for r in rows]
            if len(set(types)) == 1:
                axes[-1].plot(positions[order], ys[order], '--', color='#e49020', lw=.8)
        for key, color, label in [('cover_base', '#00a6d6', '覆盖层底'),
                                  ('first_sand', '#c25be5', '首个砂岩顶'),
                                  ('basal_sand', '#1ec766', '底部砂岩顶')]:
            axes[-1].plot(positions[order], curves[key][order], '--', color=color, lw=1.1, label=label)
        axes[-1].legend(fontsize=8)
        fig.suptitle(audit['package'] + '\n20–170MHz/501点；同一原始复谱，无去背景/AGC/时间分割；同一色标\n右图曲线来自完整Debye谱局部一次反射近似；橙线为其他接触面，非已认证目标真值')
        fig.savefig(out / (audit['package'] + '_unmarked_and_contacts.png'), dpi=135)
        plt.close(fig)
        result['packages'].append(item)
        print(json.dumps(dict(package=audit['package'], interfaces={k: {w: {n: v for n,v in d.items() if 'percentile' in n or n.startswith('complex_correlation')} for w,d in v.items()} for k,v in item['interfaces'].items()}), ensure_ascii=False), flush=True)
    result['limits'] = ('Local horizontal plane-wave primary approximation using existing H5-derived geometry. '
        'No 2D spreading, slopes, lateral scattering, multiples, PML, antenna or precision certificate. '
        'Geometry-informed correlations, not independent target ablation, blind detection or energy percentages. '
        'Reference curves do not alter actual data. FP32 originals remain FP32 evidence.')
    save(out / 'layer_kinematics.json', result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'cache', 'out']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.review, args.cache, args.out)
