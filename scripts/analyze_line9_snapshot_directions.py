"""Explore collocated native TM snapshot transport direction; no solver or SFCW edits."""
import argparse
import hashlib
import json
from pathlib import Path, PureWindowsPath

import h5py
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def align_h(previous, current, following, alpha, method):
    """Evaluate H at electric time, inside the sampled H interval."""
    if not 0 <= alpha < 1:
        raise ValueError('Expected a positive sub-snapshot Yee time offset')
    if method == 'forward_linear':
        return current + alpha * (following - current)
    if method == 'backward_linear':
        return current + alpha * (current - previous)
    if method == 'quadratic':
        return current + .5 * alpha * (following - previous) + .5 * alpha**2 * (following - 2 * current + previous)
    raise ValueError('Unknown alignment method')


def transport(ez, hx, hy):
    # E=(0,0,Ez), H=(Hx,Hy,0). This is total-field E x H, not separated rays.
    return -ez * hy, ez * hx


def intersections(xx, surface, origin, sx, sy):
    if sx <= 0 or sy <= 0:
        return []
    yy = origin[1] + (xx - origin[0]) * sy / sx
    residual = yy - surface
    indices = np.flatnonzero((residual[:-1] * residual[1:] <= 0) & (xx[1:] < origin[0]))
    return [float(xx[i] - residual[i] * (xx[i+1] - xx[i]) / (residual[i+1] - residual[i]))
            for i in indices if residual[i+1] != residual[i]]


def main(a):
    assert not a.out.exists() and not a.cache.exists()
    read = lambda p: json.loads(Path(p).read_text('utf-8'))
    c = read(a.source / 'execution_contract.json')
    v = read(a.source / 'completed_verification.json')
    assert v['completed'] and v['contract_sha256'] == sha(a.source / 'execution_contract.json')
    assert v['groups'][0]['observer_receiver_and_source_bitwise_equal']
    source_ids = read(a.runtime_sources / 'runtime_snapshot_source_identity.json')
    for identity in source_ids:
        key = 'snapshots.py' if identity['name'] == 'snapshots.py' else 'cuda_opencl/knl_snapshots.py'
        assert identity['sha256'] == c['source_identities'][key]
    assert sha(a.geometry) == c['study_manifest']['groups'][0]['geometry_sha256']
    with h5py.File(a.geometry, 'r') as h:
        g = h['data'][:, :, 0]
    native_dl = .025
    yindices = np.arange(g.shape[1]) + 1
    surface_all = np.max(np.where(g != 0, yindices[None, :], 0), axis=1) * native_dl
    bottom_all = np.max(np.where(g == 2, yindices[None, :], 0), axis=1) * native_dl
    xx = 145 + .1 * (np.arange(350) + .5)
    yy = 10 + .1 * (np.arange(315) + .5)
    ix = np.floor(xx / native_dl).astype(int)
    surface, bottom = surface_all[ix], bottom_all[ix]
    desired = np.array([np.full(350, y) for y in (32.05, 34.05, 36.05, 38.05)] +
                       [surface + .35] + [bottom + q * (surface - bottom) for q in (.25, .5, .75)])
    iy = np.argmin(abs(yy[None, None, :] - desired[:, :, None]), axis=2)
    actual = yy[iy]
    assert abs(actual - desired).max() <= .05 + 1e-12
    # Reject an interpolated sample whose underlying native support crosses a contact/PML.
    assert np.all(actual[:5] - surface > .3 - 1e-12)
    assert np.all(actual[5:] > bottom + .3) and np.all(actual[5:] < surface - .3)
    rx = np.array(c['groups'][0]['rx_m'][:2])
    point = (abs(xx - rx[0]).argmin(), abs(yy - rx[1]).argmin())
    aperture_masks = {str(radius): ((abs(xx[:, None] - rx[0]) <= radius + 1e-12) &
                                    (abs(yy[None, :] - rx[1]) <= radius + 1e-12))
                      for radius in (.15, .5, 1.)}
    assert all(np.any(mask) for mask in aperture_masks.values())
    assert all(np.all(yy[np.where(mask)[1]] < 40.5) for mask in aperture_masks.values())
    records = sorted((r for r in v['groups'][0]['snapshots'] if r['roi'] == 'local_view'), key=lambda r: r['iteration'])
    assert len(records) == 250 and [r['iteration'] for r in records] == list(range(0, 8500, 34))
    checked = []
    dt = c['study_manifest']['dt_s']

    def load(row):
        p = a.source / 'profile_snaps' / PureWindowsPath(row['file']).name
        assert sha(p) == row['sha256']
        with h5py.File(p, 'r') as h:
            assert h.attrs['gprMax'] == '4.0.1' and h.attrs['iteration'] == row['iteration']
            assert abs(h.attrs['time'] - row['iteration'] * dt) < 1e-20
            assert abs(h.attrs['magnetic_time'] - (row['iteration'] - .5) * dt) < 1e-20
            np.testing.assert_array_equal(h.attrs['origin'], [145, 10, 0])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.1, .1, .025])
            fields = []
            for key in ('Ez', 'Hx', 'Hy'):
                f = h[key][:, :, 0]
                assert f.shape == (350, 315) and f.dtype == np.float64 and np.isfinite(f).all()
                fields.append(f)
        checked.append({'file': p.name, 'sha256': row['sha256']})
        return np.array(fields)

    alpha = .5 / 34
    methods = ['forward_linear', 'backward_linear', 'quadratic']
    gates = [(150, 230), (300, 370), (330, 350)]
    centres = [200, 240, 280, 300, 320, 340]
    sums = {(method, str(radius), lo, hi): np.zeros(2) for method in methods
            for radius in (.15, .5, 1.) for lo, hi in gates}
    maps = np.zeros((len(centres), 2, 350, 315))
    map_counts = np.zeros(len(centres), dtype=int)
    traces, probes, times, steps = [], [], [], []
    previous, current = load(records[0]), load(records[1])
    for i in range(1, len(records)-1):
        following = load(records[i+1])
        t = records[i]['iteration'] * dt * 1e9
        ez = current[0]
        for method in methods:
            hx, hy = align_h(previous[1:], current[1:], following[1:], alpha, method)
            sx, sy = transport(ez, hx, hy)
            for radius, mask in aperture_masks.items():
                for lo, hi in gates:
                    if lo <= t <= hi:
                        sums[method, radius, lo, hi] += np.array([sx[mask].sum(), sy[mask].sum()])
            if method == 'forward_linear':
                probes.append(np.stack([f[np.arange(350)[None, :], iy] for f in (ez, hx, hy, sx, sy)]))
                traces.append([ez[point], hx[point], hy[point], sx[point], sy[point]])
                for j, centre in enumerate(centres):
                    if abs(t - centre) <= 4:
                        maps[j] += [sx, sy]
                        map_counts[j] += 1
        times.append(t); steps.append(records[i]['iteration'])
        previous, current = current, following
    assert len(checked) == 250 and np.all(map_counts > 0)
    maps /= map_counts[:, None, None, None]
    a.cache.mkdir(parents=True)
    cache = a.cache / 'aligned_probes_and_maps.npz'
    np.savez_compressed(cache, x_m=xx, y_m=yy, surface_m=surface, bottom_m=bottom,
                        desired_y_m=desired, actual_y_m=actual, iy=iy, time_ns=times,
                        iterations=steps, fields=np.array(probes), rx_traces=np.array(traces),
                        maps=maps, map_centres_ns=centres, map_counts=map_counts)
    angles = []
    for (method, radius, lo, hi), vector in sums.items():
        mask = aperture_masks[radius]
        ntime = sum(lo <= t <= hi for t in times)
        average = vector / (ntime * mask.sum())
        angles.append({'alignment': method, 'radius_m': float(radius), 'gate_ns': [lo, hi],
                       'spatial_samples': int(mask.sum()), 'time_samples': ntime,
                       'mean_Sx_W_m2': float(average[0]), 'mean_Sy_W_m2': float(average[1]),
                       'angle_from_positive_x_deg': float(np.degrees(np.arctan2(vector[1], vector[0]))),
                       'straight_air_backprojection_surface_x_m': intersections(xx, surface, rx, *vector)})
    # Compare the nearest snapshot point with native Rx; they are not identical locations.
    with h5py.File(a.source / 'profile.h5', 'r') as h:
        native = h['rxs/rx1/Ez'][np.array(steps)]
    observed = np.array(traces)[:, 0]
    keep = (np.array(times) >= 300) & (np.array(times) <= 370)
    correlation = float(np.dot(observed[keep], native[keep]) /
                        (np.linalg.norm(observed[keep]) * np.linalg.norm(native[keep])))
    a.out.mkdir(parents=True)
    result = {'status': 'EXPLORATORY_NATIVE_DIRECTION_LOCAL250_HASH_VERIFIED', 'solver_runs': 0,
              'script_sha256': sha(__file__), 'contract_sha256': sha(a.source / 'execution_contract.json'),
              'certificate_sha256': sha(a.source / 'completed_verification.json'), 'geometry_sha256': sha(a.geometry),
              'runtime_snapshot_source_identities': source_ids, 'locally_read_snapshots': checked,
              'cache_sha256': sha(cache), 'methods': methods, 'magnetic_alignment_fraction': alpha,
              'snapshot_spacing_ns': 34 * dt * 1e9, 'max_probe_rounding_m': float(abs(actual-desired).max()),
              'receiver_m': rx.tolist(), 'nearest_snapshot_point_m': [float(xx[point[0]]), float(yy[point[1]])],
              'late_point_vs_native_rx_signed_correlation': correlation, 'aperture_gate_directions': angles,
              'gates_status': 'Exploratory/post hoc native-time gates, not predeclared physical acceptance',
              'units': 'Total sampled E cross aligned H in W/m^2; time/area means, not separated pulse energies',
              'limits': ['Native Ricker wideband time, not source-normalized 20-170MHz SFCW time.',
                         'Coarse0.1m spatial and2.005ns temporal samples; interpolation sensitivity is not an FDTD error budget.',
                         'No conserved energy fractions, separated rays, unique bounce count or unique interface attribution.',
                         'Straight air backprojection is a geometric hypothesis from total local transport; interference may bias it.',
                         'Only190m high-loss H0, same mother model; no site, whole-line, antenna or finite3D validation.']}
    (a.out / 'analysis.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'correlation': correlation,
                      'late_directions': [r for r in angles if r['gate_ns'] == [330, 350]]}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'geometry', 'runtime-sources', 'cache', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    main(p.parse_args())
