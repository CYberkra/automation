"""Independent material/length/inverse and adaptive planar checks of a ray hypothesis."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.integrate import quad_vec
from hs_capsule_identity import sha256 as sha

C0 = 299792458.
EPS0 = 8.8541878128e-12


def refractive(db, f):
    values = []
    for medium in db.values():
        b = medium['base']; omega = 2 * np.pi * f
        er = b['relative_permittivity'] + b['electric_conductivity_s_per_m'] / (1j * omega * EPS0)
        for pole in medium.get('poles', []):
            er = er + pole['relative_permittivity_difference'] / (1 + 1j * omega * pole['relaxation_time_s'])
        values.append(np.sqrt(er))
    return np.stack(values, axis=-1)


def main(a):
    assert not a.out.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    r = read(a.public / 'analysis.json'); m = read(a.package / 'manifest.json')
    assert r['manifest_sha256'] == sha(a.package / 'manifest.json') and r['numerical_sha256'] == sha(a.numerical)
    assert r['script_sha256'] == sha(Path(__file__).with_name('diagnose_line9_dense_cover_multiple.py'))
    with h5py.File(a.numerical) as h:
        f = h['frequency_Hz'][:]; spectra = h['ray_spectra_dimensionless'][:]
        planar = h['planar_components'][:]; actual = h['actual_H0'][:]
    np.testing.assert_array_equal(f, 20e6 + np.arange(501) * 300000.)
    dbs = []
    for v in ['high', 'low']:
        group = next(g for g in m['groups'] if g['variant'] == v)
        path = a.package / group['material']; assert sha(path) == r['material_sha256'][v]
        dbs.append(read(path)['materials'])
        native = a.source / f'{v}_x19000_H0.h5'; assert sha(native) == r['native_identities'][native.stem]
    nr = [refractive(db, f) for db in dbs]; time_errors = []; spectrum_errors = []; metric_checks = 0
    for row in r['rays']:
        eligible = []
        for j, trial in enumerate(row['trials']):
            points = np.array(trial['nodes_xy_m']); media = np.array(trial['segment_materials'])
            np.testing.assert_array_equal(media, [0, 1, 1, 1, 1, 0])
            lengths = np.linalg.norm(points[1:] - points[:-1], axis=1)
            by_medium = np.array([sum(lengths[media == k]) for k in range(4)])
            np.testing.assert_allclose(by_medium, trial['path_length_by_material_m'], rtol=0, atol=1e-12)
            phase = float(np.dot(by_medium, nr[0][250].real) / C0 * 1e9)
            err = abs(phase - trial['phase_time95_ns']); assert err < 1e-9; time_errors.append(err)
            for v in range(2):
                n = nr[v]; reflection_air = (n[:, 0] - n[:, 1]) / (n[:, 0] + n[:, 1])
                reflection_cover = (n[:, 1] - n[:, 2]) / (n[:, 1] + n[:, 2])
                coefficient = -(1 - reflection_air ** 2) * reflection_air * reflection_cover ** 2
                expected = coefficient * np.exp(-1j * 2 * np.pi * f * np.sum(n * by_medium, axis=1) / C0)
                stored = spectra[:, trial['spectral_column'], v]
                error = float(np.linalg.norm(stored - expected) / np.linalg.norm(expected)); assert error < 1e-10; spectrum_errors.append(error)
            if trial['optimizer_success'] and trial['geometry_check'] == 'NO_OUTSIDE_BAND_MISMATCH': eligible.append(j)
        assert row['best_geometry_only_index'] == min(eligible, key=lambda j: row['trials'][j]['phase_time95_ns'])
    t = np.arange(4008) / (4008 * 300000.); keep = (t * 1e9 >= 300) & (t * 1e9 <= 450)
    kernel = np.exp(2j * np.pi * t[keep, None] * f) / 501
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w /= w.mean(); natives = kernel @ (actual * w[:, None])
        ps = (kernel @ (planar.reshape(501, -1) * w[:, None])).reshape(keep.sum(), 2, 6)
        rays = (kernel @ (spectra.reshape(501, -1) * w[:, None])).reshape(keep.sum(), -1, 2)
        peak = lambda q: float(t[keep][np.argmax(abs(q))] * 1e9)
        for j, v in enumerate(['high', 'low']):
            recorded = r['metrics'][window][v]; q = natives[:, j]; p = ps[:, j, :5].sum(axis=1)
            values = dict(native_H0_peak_ns=peak(q), planar_H0_peak_ns=peak(p), planar_cover_twice_peak_ns=peak(ps[:, j, 3]),
                planar_H0_unfitted_relative_L2=float(np.linalg.norm(q - p) / np.linalg.norm(q)),
                planar_H0_complex_correlation=float(abs(np.vdot(q, p)) / (np.linalg.norm(q) * np.linalg.norm(p))))
            for name, value in values.items():
                assert abs(value - recorded[name]) / max(1, abs(recorded[name])) < 1e-8; metric_checks += 1
            for i, row in enumerate(r['rays']):
                for k, trial in enumerate(row['trials']):
                    value = peak(rays[:, trial['spectral_column'], j])
                    assert abs(value - recorded['ray_peak_ns_by_stride_and_start'][i][k]) < 1e-8; metric_checks += 1
    adaptive_errors = []; admittance_errors = []
    g = r['station']['geometry']['high']; hs = r['station']['tx_m'][1] - g['surface_y_m']; hr = r['station']['rx_m'][1] - g['surface_y_m']
    dx = abs(r['station']['tx_m'][0] - r['station']['rx_m'][0]); d1 = g['cover_base']['depth_m']; d2 = g['basal_sand']['depth_m'] - d1
    for v, db in enumerate(dbs):
        for j in [0, 250, 500]:
            fj = f[j]; k = 2 * np.pi * fj / C0 * refractive(db, np.array(fj)); k0 = k[0].real
            def parts(q, ky0):
                ky = np.sqrt(k * k - q * q + 0j); ky = np.where(ky.imag > 0, -ky, ky)
                r0 = (ky0 - ky[1]) / (ky0 + ky[1]); r1 = (ky[1] - ky[2]) / (ky[1] + ky[2])
                p1 = np.exp(-2j * ky[1] * d1)
                once = (1 - r0 * r0) * r1 * p1; twice = -(1 - r0 * r0) * r0 * r1 * r1 * p1 * p1
                full = []
                for last in [2, 3]:
                    y = ky[last]
                    for layer in range(last - 1, 0, -1):
                        tangent = np.tan(ky[layer] * [0, d1, d2][layer])
                        y = ky[layer] * (y + 1j * ky[layer] * tangent) / (ky[layer] + 1j * y * tangent)
                    full.append((ky0 - y) / (ky0 + y))
                recurrence = (r0 + r1 * p1) / (1 + r0 * r1 * p1)
                admittance_errors.append(float(abs(full[0] - recurrence)))
                return np.array([r0, once, twice, full[0] - r0 - once - twice, full[1] - full[0]])
            def propagating(theta):
                q = k0 * np.sin(theta); ky = k0 * np.cos(theta)
                return np.cos(q * dx) * np.exp(-1j * ky * (hs + hr)) * parts(q, ky)
            def evanescent(u):
                q = k0 * np.cosh(u); ky = -1j * k0 * np.sinh(u)
                return 1j * np.cos(q * dx) * np.exp(-1j * ky * (hs + hr)) * parts(q, ky)
            pp, _ = quad_vec(propagating, 0, np.pi / 2, epsabs=1e-12, epsrel=1e-11)
            ee, _ = quad_vec(evanescent, 0, np.arcsinh(40 / (k0 * (hs + hr))), epsabs=1e-12, epsrel=1e-11)
            expected = -fj / (EPS0 * C0 * C0 * .025) * (pp + ee)
            errors = abs(expected - planar[j, v, 1:]) / abs(planar[j, v, 1:]); assert max(errors) < 1e-5
            adaptive_errors.append(errors.tolist())
    result = dict(status='PASS_INDEPENDENT_LENGTH_MATERIAL_RAY_INVERSE_ADAPTIVE_PLANAR',
        auditor_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        optical_phase_time_error_ns_max=max(time_errors), independent_ray_spectrum_relative_L2_max=max(spectrum_errors),
        metric_scalar_checks=metric_checks, adaptive_planar_relative_error_max=max(max(x) for x in adaptive_errors),
        input_admittance_reflection_absolute_error_max=max(admittance_errors),
        optimizer_independently_certified=False, physical_path_certified=False,
        limits='Verifies reported arithmetic/material/planar integration, not optimizer uniqueness,global stationary ray,contact interpolation or observed path identity. Independent planar audit uses continuum interfaces; no FDTD agreement implied.')
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8'); print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'package', 'public', 'numerical', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
