"""Prepare four bounded cover-branch controls; no freeze, CUDA or solver launch."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists() and not a.evidence.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    parent = read(a.package / 'manifest.json'); evidence = read(a.candidate / 'analysis.json')
    audit = read(a.candidate / 'independent_audit.json'); prefix = read(a.public / 'analysis.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(a.candidate / 'analysis.json')
    assert evidence['analysis_sha256'] == sha(a.public / 'analysis.json')
    assert evidence['manifest_sha256'] == sha(a.package / 'manifest.json')
    base = next(g for g in parent['groups'] if g['id'] == 'high_x19000_H0')
    for key in ['input', 'geometry', 'material']:
        assert sha(a.package / base[key]) == base[key + '_sha256']
    receipt = next(g for g in prefix['native'] if g['id'] == base['id'])
    native = a.source / (base['id'] + '.h5'); assert sha(native) == receipt['native_sha256']
    with h5py.File(a.package / base['geometry']) as h:
        data = h['data'][:]; assert data.shape == (8400, 1700, 1)
        assert set(np.unique(data)) == {0, 1, 2}
        attributes = {k: np.asarray(v).tolist() for k, v in h.attrs.items() if k == 'dx_dy_dz'}
    card = (a.package / base['input']).read_text('utf-8'); assert '#snapshot:' not in card
    dt = parent['dt_s']; timeline = [j for j in range(0, 20352, 34) if j * dt <= 500e-9]
    suffix = ''.join(f'#snapshot: 145 10 0 180 41.5 0.025 0.1 0.1 0.025 {j} snap{j:05d}.h5\n' for j in timeline)
    groups = []; a.out.mkdir(parents=True)
    baseline = a.out / 'baseline'; baseline.mkdir()
    shutil.copyfile(native, baseline / 'native_H0.h5')
    shutil.copyfile(a.package / base['input'], baseline / 'profile.in')
    shutil.copyfile(a.package / base['geometry'], baseline / 'geometry.h5')
    shutil.copyfile(a.package / base['material'], baseline / 'materials.json')
    for sid, snapshots, roi in [('base_snapshot_H0', True, None), ('flat_snapshot_H0', True, None),
                                ('far_cover_removed_H0', False, [154., 164.]), ('near_cover_removed_H0', False, [164., 174.])]:
        folder = a.out / sid / 'geometries'; folder.mkdir(parents=True)
        geometry = folder / Path(base['geometry']).name; material = folder / Path(base['material']).name
        shutil.copyfile(a.package / base['geometry'], geometry); shutil.copyfile(a.package / base['material'], material)
        changed = data.copy()
        if sid.startswith('flat'):
            changed[:] = data[6800:6801]
        if roi is not None:
            x = np.arange(8400) * .025
            region = ((x >= roi[0]) & (x < roi[1]))[:, None, None]
            changed[(data == 2) & region] = 1
        if not np.array_equal(data, changed):
            with h5py.File(geometry, 'r+') as h: h['data'][:] = changed
        inp = a.out / sid / 'cases' / sid / 'profile.in'; inp.parent.mkdir(parents=True)
        lines = [f'#title: {sid}; 190m high-loss H0 cover-branch diagnostic' if s.startswith('#title:') else s for s in card.splitlines()]
        inp.write_text('\n'.join(lines) + '\n' + (suffix if snapshots else ''), encoding='utf-8')
        g = copy.deepcopy(base); g.update(id=sid, input=inp.relative_to(a.out).as_posix(),
            geometry=geometry.relative_to(a.out).as_posix(), material=material.relative_to(a.out).as_posix(),
            snapshots=snapshots, changed_voxels=int(np.count_nonzero(changed != data)),
            cover_removal_local_x_roi_m=roi, cover_removal_chainage_roi_m=[v + 20 for v in roi] if roi else None)
        for key in ['input', 'geometry', 'material']: g[key + '_sha256'] = sha(a.out / g[key])
        groups.append(g)
    history = 350 * 315 * len(timeline) * 6 * 8
    manifest = dict(status='PREPARED_NEXT_DIAGNOSTIC_NOT_FROZEN_OR_RUN', generator_sha256=sha(__file__),
        approval_basis='Standing user autonomous research/SSH/snapshot authorization; no permission request. Four new configurations, explicitly separate from current46attempt batch.',
        groups=groups, max_new_runs=4, no_retry=True, parent_manifest_sha256=sha(a.package / 'manifest.json'),
        candidate_analysis_sha256=sha(a.candidate / 'analysis.json'), candidate_audit_sha256=sha(a.candidate / 'independent_audit.json'),
        baseline=dict(native='baseline/native_H0.h5', native_sha256=sha(baseline / 'native_H0.h5'),
            input='baseline/profile.in', input_sha256=sha(baseline / 'profile.in'),
            geometry='baseline/geometry.h5', geometry_sha256=sha(baseline / 'geometry.h5'),
            material='baseline/materials.json', material_sha256=sha(baseline / 'materials.json')),
        source_runtime_version='4.0.1', required_native_dtype='float64', dt_s=dt, expected_samples=20352,
        snapshot_iterations=timeline, snapshot_shape=[350, 315, 1], snapshot_origin_m=[145., 10., 0.],
        snapshot_spacing_m=[.1, .1, .025], fields=['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'],
        six_field_snapshot_history_bytes_per_observed_case=history,
        payload_device_reserve_bytes=6 * 2**30, estimated_peak_device_bytes=history + 6 * 2**30,
        total_two_case_snapshot_history_bytes=2 * history,
        start_requires_current_46_attempt_batch_terminal_and_audited=True, execution_contract_frozen=False,
        invariant='All cards preserve 210x42.5m,2.5cm,1200ns,40A100MHz Ricker,source/receiver coordinates,high-loss materials,HORIPML80,averaging andFP64. Geometry-only controls; observers only in base/flat. Same exact501-tone SFCW/Hann/Blackman/full-source normalization.',
        controls='Flat copies the whole190m midpoint H0 column to everyx,including PML continuation. Far/near changes mudstone2 to cover1 throughout specifiedROI below the original cover interface;air/ground/topcover voxels preserved. Introduces vertical edges and changes lower continuation,so dependence is regional,not a pure multiple isolation.',
        planned_checks='Observer baseline receiver must match existingnative bitwise. Independently audit all geometry/keys/cards/receipts. Compare fixed300-450ns and inherited basal common345.618097-374.608117ns; retain early0-120 and late450-1100 diagnostics. Do not select a metric window from a changed result.',
        delivery='Chinese shared-scale SFCW configuration panels,raw Ricker wavefield GIF and unchanged source/timing annotations; no invented continuous Bscan from one station.',
        limits='No solver is launched by this preparer. Fresh target-runtime/resource contract and independent input audit required. Snapshots native Ricker,not SFCW;0.1m/34steps visualization decimation,not FDTD coarsening. Regional ablation includes new edges/all interactions;flat changes whole terrain/cover geometry. No unique ray,field calibration,material replacement,3D or training claim.')
    (a.out / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    a.evidence.mkdir(parents=True)
    (a.evidence / 'preparation.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=manifest['status'], snapshots_each=len(timeline), device_payload_GiB=history / 2**30,
        device_estimate_GiB=(history + 6 * 2**30) / 2**30, changes={g['id']: g['changed_voxels'] for g in groups})))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['package', 'candidate', 'source', 'public', 'out', 'evidence']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
