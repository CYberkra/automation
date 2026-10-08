"""Exploratory nonlocal cover-double-trip ray and local planar comparison at190m.

The candidate is selected by minimum optical length among fixed starts, never by
matching an observed peak. No ray amplitude is compared with native field units.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.optimize import minimize
from scipy.special import hankel2
from diagnose_line9_refracted_paths import extract_columns
from diagnose_line9_v5_planar_green import integral, reflect_parts, MU0
from review_line9_result_packages import indices, C0, ray_time
from analyze_line9_v401_version_controls import inverse, FREQ
from hs_capsule_identity import sha256 as sha

SHIFTS = [-10., -5., 0., 5., 10.]
MEDIA = np.array([0, 1, 1, 1, 1, 0])
CONTACTS = [0, 1, 0, 1, 0]


def solve(tx, rx, xx, curves, refractive):
    def nodes(v):
        yy = [np.interp(x, xx, curves[c]) for x, c in zip(v, CONTACTS)]
        return np.vstack([tx[:2], np.column_stack([v, yy]), rx[:2]])
    def objective(v):
        d = np.linalg.norm(np.diff(nodes(v), axis=0), axis=1)
        return float(np.dot(d, refractive[MEDIA].real) / C0 * 1e9)
    midpoint = (tx[0] + rx[0]) / 2; trials = []
    for shift in SHIFTS:
        fit = minimize(objective, np.full(5, np.clip(midpoint + shift, xx[0], xx[-1])),
            method='L-BFGS-B', bounds=[(xx[0], xx[-1])] * 5,
            options=dict(ftol=1e-12, gtol=1e-6, maxiter=1000, maxls=40))
        points = nodes(fit.x); lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
        trials.append(dict(start_shift_m=shift, optimizer_success=bool(fit.success),
            optimizer_message=str(fit.message), phase_time95_ns=float(fit.fun),
            nodes_xy_m=points.tolist(), segment_materials=MEDIA.tolist(),
            path_length_by_material_m=np.bincount(MEDIA, weights=lengths, minlength=4).tolist()))
    assert any(t['optimizer_success'] for t in trials)
    return trials


def support_check(points, xx, curves, data):
    wrong = near = 0.
    for p, q, medium in zip(points[:-1], points[1:], MEDIA):
        length = np.linalg.norm(q - p); count = max(2, int(np.ceil(length / .0125)))
        u = (np.arange(count) + .5) / count; samples = p + u[:, None] * (q - p)
        cells = np.floor(samples / .025).astype(int)
        if (cells < 0).any() or (cells >= data.shape).any():
            return dict(wrong_material_length_outside_contact_band_m=float(length), near_contact_unresolved_length_m=float(near), geometry_check='REJECT_OUTSIDE_DOMAIN')
        distance = np.min(abs(samples[:, 1, None] - np.column_stack([np.interp(samples[:, 0], xx, c) for c in curves])), axis=1)
        contact = distance <= .05; near += length * np.mean(contact)
        wrong += length * np.mean((data[cells[:, 0], cells[:, 1]] != medium) & ~contact)
    return dict(wrong_material_length_outside_contact_band_m=float(wrong), near_contact_unresolved_length_m=float(near),
                geometry_check='NO_OUTSIDE_BAND_MISMATCH' if wrong == 0 else 'REJECT_WRONG_MATERIAL_SEGMENT')


def spectrum(materials, lengths):
    n = np.array([indices(materials, f) for f in FREQ])
    r01 = (n[:, 0] - n[:, 1]) / (n[:, 0] + n[:, 1])
    r12 = (n[:, 1] - n[:, 2]) / (n[:, 1] + n[:, 2])
    return -(1 - r01 * r01) * r01 * r12 * r12 * np.exp(-2j * np.pi * FREQ * (n @ lengths) / C0)


def flat_check():
    xx = np.linspace(-50, 50, 501); curves = np.array([np.zeros(len(xx)), np.full(len(xx), -7.)])
    tx = np.array([0., 10.]); rx = np.array([1.3, 10.]); refractive = np.array([1., 3., 4., 2.])
    trials = solve(tx, rx, xx, curves, refractive)
    expected = ray_time(np.array([10., 14.]), refractive[:2], 1.3)
    error = min(abs(q['phase_time95_ns'] - expected) for q in trials if q['optimizer_success'])
    assert error < 1e-4
    return dict(flat_double_cover_trip_phase_time_error_ns=error)


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    m = read(a.package / 'manifest.json'); pub = read(a.public / 'analysis.json'); audit = read(a.public / 'independent_audit.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(a.public / 'analysis.json')
    assert pub['numerical_sha256'] == sha(a.response)
    native_hashes = {r['id']: r['native_sha256'] for r in pub['native']}
    for variant in ['high', 'low']:
        sid = f'{variant}_x19000_H0'
        assert sha(a.source / (sid + '.h5')) == native_hashes[sid]
    station = next(s for s in m['stations'] if s['chainage_m'] == 190.)
    group = next(g for g in m['groups'] if g['id'] == 'high_x19000_H1')
    gp = a.package / group['geometry']; assert sha(gp) == group['geometry_sha256']
    with h5py.File(gp) as h:
        data = h['data'][:, :, 0]; np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025] * 3)
    columns = extract_columns(data, .025); signature = [(0, 1), (1, 2)]
    valid = np.array([[(a, b) for a, b, y in c[:2]] == signature for c in columns])
    middle = round((station['tx_m'][0] + station['rx_m'][0]) / .05)
    assert valid[middle]; left = right = middle
    while left > 80 and valid[left - 1]: left -= 1
    while right < len(valid) - 81 and valid[right + 1]: right += 1
    dbs = {}; mat_hashes = {}
    for variant in ['high', 'low']:
        g = next(g for g in m['groups'] if g['variant'] == variant)
        db = a.package / g['material']; assert sha(db) == g['material_sha256']
        dbs[variant] = read(db)['materials']; mat_hashes[variant] = sha(db)
    n95 = indices(dbs['high'], 95e6)
    assert indices(dbs['low'], 95e6)[1] == n95[1]
    rays = []; ray_spectra = []; planar = []; convergence = []
    tx = np.array(station['tx_m']); rx = np.array(station['rx_m'])
    for stride in [8, 4]:
        ix = np.unique(np.r_[left, np.arange(left, right + 1, stride), right]); xx = ix * .025
        curves = np.array([[columns[i][j][2] for i in ix] for j in range(2)])
        trials = solve(tx, rx, xx, curves, n95)
        for row in trials:
            row.update(support_check(np.array(row['nodes_xy_m']), xx, curves, data))
            row['spectral_column'] = len(ray_spectra)
            ray_spectra.append(np.column_stack([spectrum(dbs[v], row['path_length_by_material_m']) for v in ['high', 'low']]))
        accepted = [j for j, q in enumerate(trials) if q['optimizer_success'] and q['geometry_check'] == 'NO_OUTSIDE_BAND_MISMATCH']
        assert accepted
        best = min(accepted, key=lambda j: trials[j]['phase_time95_ns'])
        rays.append(dict(stride=stride, interface_sample_spacing_m=stride * .025,
            crossing_interval_m=[float(xx[0]), float(xx[-1])], trials=trials, best_geometry_only_index=best,
            successful_valid_time_range_ns=[min(trials[j]['phase_time95_ns'] for j in accepted), max(trials[j]['phase_time95_ns'] for j in accepted)]))
    for v in ['high', 'low']:
        g = station['geometry'][v]; hs = tx[1] - g['surface_y_m']; hr = rx[1] - g['surface_y_m']
        cover = g['cover_base']['depth_m']; mud = g['basal_sand']['depth_m'] - cover; dx = abs(rx[0] - tx[0])
        n = np.array([indices(dbs[v], f) for f in FREQ]); k = 2 * np.pi * FREQ[:, None] * n / C0
        kernel = lambda q, ky: reflect_parts(q, ky, k[:, 1], k[:, 2], k[:, 3], cover, mud)
        coarse = integral(FREQ, hs + hr, dx, kernel, order=256) * (-FREQ * MU0 / .025)[:, None]
        fine = integral(FREQ, hs + hr, dx, kernel, order=512) * (-FREQ * MU0 / .025)[:, None]
        errors = np.linalg.norm(fine - coarse, axis=0) / np.linalg.norm(fine, axis=0)
        assert max(errors) < 1e-7; convergence.append(errors.tolist())
        direct = -2 * np.pi * FREQ * MU0 / (4 * .025) * hankel2(0, k[:, 0] * np.hypot(dx, hs - hr))
        planar.append(np.column_stack([direct, fine]))
    planar = np.stack(planar, axis=1); ray_spectra = np.stack(ray_spectra, axis=1)
    with h5py.File(a.response) as h:
        ids = list(h['ids'].asstr()[:]); np.testing.assert_array_equal(h['frequency_Hz'][:], FREQ)
        actual = h['response'][:, [ids.index(f'{v}_x19000_H0') for v in ['high', 'low']]]
    metrics = {}; wave = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean(); native, t = inverse(actual, w)
        analytic, _ = inverse(planar.reshape(501, -1), w); analytic = analytic.reshape(4008, 2, 6)
        paths, _ = inverse(ray_spectra.reshape(501, -1), w); paths = paths.reshape(4008, -1, 2)
        gate = (t * 1e9 >= 300) & (t * 1e9 <= 450); rows = {}
        peak = lambda z: float(t[gate][np.argmax(abs(z[gate]))] * 1e9)
        for j, v in enumerate(['high', 'low']):
            total = analytic[:, j, :5].sum(axis=1); q = native[gate, j]; p = total[gate]
            rows[v] = dict(native_H0_peak_ns=peak(native[:, j]), planar_H0_peak_ns=peak(total),
                planar_cover_twice_peak_ns=peak(analytic[:, j, 3]),
                planar_H0_unfitted_relative_L2=float(np.linalg.norm(q - p) / np.linalg.norm(q)),
                planar_H0_complex_correlation=float(abs(np.vdot(q, p)) / (np.linalg.norm(q) * np.linalg.norm(p))),
                ray_peak_ns_by_stride_and_start=[[peak(paths[:, r['spectral_column'], j]) for r in s['trials']] for s in rays],
                geometry_selected_ray_peak_ns=[peak(paths[:, s['trials'][s['best_geometry_only_index']]['spectral_column'], j]) for s in rays])
        metrics[window] = rows; wave[window] = (native, analytic, paths, t)
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']; plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), layout='constrained')
    lo, hi = 150., 175.; ia, ib = round(lo / .025), round(hi / .025)
    axes[0].imshow(data[ia:ib:4, ::4].T, origin='lower', extent=[lo + 20, hi + 20, 0, 42.5], aspect='auto',
        cmap=ListedColormap(['#e8f4fa', '#e0b992', '#9c8c87', '#e6cf70']), norm=BoundaryNorm(np.arange(5) - .5, 4), interpolation='nearest')
    for row, color in zip(rays, ['#2563a5', '#813ca2']):
        best = row['trials'][row['best_geometry_only_index']]; points = np.array(best['nodes_xy_m'])
        axes[0].plot(points[:, 0] + 20, points[:, 1], '.-', color=color, lw=1.5, label=f'界面采样{row["interface_sample_spacing_m"]:g}m / 几何最短候选')
    axes[0].set(xlim=(lo + 20, hi + 20), ylim=(20, 41), title='覆盖层内两次底反射候选，非波场快照', xlabel='剖面里程 / m', ylabel='模型y / m'); axes[0].legend(fontsize=8)
    native, analytic, paths, t = wave['hann']; show = (t * 1e9 >= 300) & (t * 1e9 <= 410)
    for j, v in enumerate(['high', 'low']):
        ax = axes[j + 1]; ax.plot(t[show] * 1e9, abs(native[show, j]), color='black', label='实际无底砂H0')
        ax.plot(t[show] * 1e9, abs(analytic[show, j, :5].sum(axis=1)), color='#d58000', ls='--', label='局部平层完整H0，无拟合')
        ax.plot(t[show] * 1e9, abs(analytic[show, j, 3]), color='#29945e', ls=':', label='平层覆盖层二次往返')
        for row, color in zip(rays, ['#2563a5', '#813ca2']):
            value = metrics['hann'][v]['geometry_selected_ray_peak_ns'][rays.index(row)]
            ax.axvline(value, color=color, ls=':', label=f'离轴候选到时 {value:.2f}ns')
        ax.set(title=('高损耗' if v == 'high' else '低损耗诊断') + ' / Hann / 190m', xlabel='SFCW时间 / ns', ylabel='包络 / (V/m)/(A·m)'); ax.legend(fontsize=8)
    fig.suptitle('190m晚到强波：局部平层与非局部覆盖层候选分开核查\n射线仅预测到时，未拟合幅度/相位/时延；五初值多局部极小，不是唯一物理归因，各幅度面板独立纵轴')
    fig.savefig(a.out / 'nonlocal_cover_multiple_hypothesis.png', dpi=130); plt.close(fig)
    with h5py.File(a.numerical, 'x') as h:
        h.create_dataset('frequency_Hz', data=FREQ); h.create_dataset('planar_components', data=planar)
        h.create_dataset('ray_spectra_dimensionless', data=ray_spectra); h.create_dataset('actual_H0', data=actual)
    result = dict(status='EXPLORATORY_NONLOCAL_COVER_DOUBLE_TRIP_NOT_CAUSAL_CERTIFICATION',
        script_sha256=sha(__file__), geometry_sha256=sha(gp), material_sha256=mat_hashes,
        manifest_sha256=sha(a.package / 'manifest.json'), analysis_sha256=sha(a.public / 'analysis.json'),
        input_audit_sha256=sha(a.public / 'independent_audit.json'), response_sha256=sha(a.response), numerical_sha256=sha(a.numerical),
        helper_sha256={s: sha(Path(__file__).with_name(s)) for s in ['diagnose_line9_refracted_paths.py', 'diagnose_line9_v5_planar_green.py', 'review_line9_result_packages.py']},
        station=station, flat_check=flat_check(), rays=rays, quadrature_256_512_relative_L2=convergence,
        metrics=metrics, native_identities={r['id']: r['native_sha256'] for r in pub['native'] if r['id'].endswith('x19000_H0')},
        new_solves=0, selection='Best phase95 optical length among fixed five starts and valid sampled paths; never closest observed arrival.',
        limits='Posthoc exploratory local planar and real-index ray proxies. Five-start local minima are not global minimum/stationary-branch proof. Two-cell contact ambiguity,linear interface interpolation,normal Fresnel coefficients,no ray spreading/diffraction/frequency-dependent bending. Matching time does not identify an observed path; no ray amplitude calibration,physical error budget or field/whole-line/finite3D validation. No production changes.')
    (a.out / 'analysis.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], metrics=metrics,
        geometry_selected_paths=[r['trials'][r['best_geometry_only_index']] for r in rays])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'package', 'public', 'response', 'out', 'numerical']:
        parser.add_argument('--' + key, type=Path, required=True)
    main(parser.parse_args())
