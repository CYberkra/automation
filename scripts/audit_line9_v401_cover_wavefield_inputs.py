"""Read-only independent voxel/card/observer audit; execution remains unfrozen."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def audit(package):
    m = json.loads((package / 'manifest.json').read_text('utf-8'))
    assert m['status'] == 'PREPARED_NEXT_DIAGNOSTIC_NOT_FROZEN_OR_RUN' and not m['execution_contract_frozen']
    assert m['max_new_runs'] == 4 and m['no_retry']
    assert not (package / 'execution_contract.json').exists() and not (package / 'execution.jsonl').exists()
    base = m['baseline']
    for key in ['native', 'input', 'geometry', 'material']:
        assert sha(package / base[key]) == base[key + '_sha256']
    with h5py.File(package / base['geometry']) as h:
        original = h['data'][:]; attributes = dict(h.attrs); keys = h['material_keys'][:]
        assert original.shape == (8400, 1700, 1) and np.unique(original).tolist() == [0, 1, 2]
    with h5py.File(package / base['native']) as h:
        assert h.attrs['gprMax'] == '4.0.1' and h['rxs/rx1/Ez'].dtype == np.float64
        assert h['rxs/rx1/Ez'].shape == (20352,) and h.attrs['dt'] == m['dt_s']
        assert np.isfinite(h['rxs/rx1/Ez'][:]).all()
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [8400, 1700, 1])
    baseline_lines = (package / base['input']).read_text('utf-8').splitlines()
    invariant_lines = [s for s in baseline_lines if not s.startswith('#title:')]
    expected_steps = list(range(0, int(np.floor(500e-9 / m['dt_s'])) + 1, 34))
    assert m['snapshot_iterations'] == expected_steps and len(expected_steps) == 250
    assert m['snapshot_shape'] == [350, 315, 1] and m['snapshot_origin_m'] == [145., 10., 0.]
    assert m['snapshot_spacing_m'] == [.1, .1, .025]
    payload = 350 * 315 * 250 * 6 * 8
    assert m['six_field_snapshot_history_bytes_per_observed_case'] == payload
    assert m['estimated_peak_device_bytes'] == payload + 6 * 2**30
    assert m['total_two_case_snapshot_history_bytes'] == 2 * payload
    assert [g['id'] for g in m['groups']] == ['base_snapshot_H0', 'flat_snapshot_H0', 'far_cover_removed_H0', 'near_cover_removed_H0']
    records = []
    for g in m['groups']:
        for key in ['input', 'geometry', 'material']:
            assert sha(package / g[key]) == g[key + '_sha256']
        assert g['material_sha256'] == base['material_sha256']
        with h5py.File(package / g['geometry']) as h:
            data = h['data'][:]; assert data.shape == original.shape and set(h.attrs) == set(attributes)
            for key, value in attributes.items(): np.testing.assert_array_equal(h.attrs[key], value)
            np.testing.assert_array_equal(h['material_keys'][:], keys)
        expected = original.copy()
        if g['id'] == 'flat_snapshot_H0':
            # Independent explicit horizontal layer realization from midpoint cuts.
            column = original[6800, :, 0]; cuts = np.flatnonzero(np.diff(column) != 0) + 1
            assert len(cuts) == 2
            expected[:, :cuts[0], 0] = 2; expected[:, cuts[0]:cuts[1], 0] = 1; expected[:, cuts[1]:, 0] = 0
        elif g['id'] == 'far_cover_removed_H0':
            slab = expected[6160:6560]; slab[slab == 2] = 1
            assert g['cover_removal_local_x_roi_m'] == [154., 164.]
        elif g['id'] == 'near_cover_removed_H0':
            slab = expected[6560:6960]; slab[slab == 2] = 1
            assert g['cover_removal_local_x_roi_m'] == [164., 174.]
        np.testing.assert_array_equal(data, expected)
        assert np.count_nonzero(data != original) == g['changed_voxels']
        if 'removed' in g['id']:
            assert np.array_equal(data[original == 0], original[original == 0])
            assert np.array_equal(data[original == 1], original[original == 1])
        card = (package / g['input']).read_text('utf-8').splitlines()
        normal = [s for s in card if not s.startswith(('#title:', '#snapshot:'))]
        assert normal == invariant_lines
        observations = [s for s in card if s.startswith('#snapshot:')]
        wanted = 250 if g['id'] in ['base_snapshot_H0', 'flat_snapshot_H0'] else 0
        assert g['snapshots'] == bool(wanted) and len(observations) == wanted
        for line, step in zip(observations, expected_steps):
            tokens = line.split()
            assert tokens[0] == '#snapshot:' and len(tokens) == 12
            np.testing.assert_array_equal(np.array(tokens[1:10], float), [145, 10, 0, 180, 41.5, .025, .1, .1, .025])
            assert int(tokens[10]) == step and tokens[11] == f'snap{step:05d}.h5'
        from gprMax.hash_cmds_file import get_user_objects
        get_user_objects(card, input_dir=(package / g['input']).parent)
        for role in ['tx_m', 'rx_m']:
            pos = g[role]; ix, iy = np.rint(np.array(pos[:2]) / .025).astype(int)
            assert data[ix, iy, 0] == 0
        records.append(dict(id=g['id'], changed_voxels=g['changed_voxels'], snapshot_count=wanted,
                            input_sha256=g['input_sha256'], geometry_sha256=g['geometry_sha256']))
    return dict(status='PASS_INDEPENDENT_PREPARED_COVER_CONTROLS_NO_EXECUTION_CONTRACT',
        auditor_sha256=sha(__file__), manifest_sha256=sha(package / 'manifest.json'),
        groups=records, calls_solver=False, calls_CUDA=False, calls_training=False,
        input_ready=True, execution_ready=False,
        pending='Wait for current46attempt run terminal and audit; then freeze a fresh target-specific runtime/resource/identity/heartbeat contract. No new permission needed under standing authorization.',
        snapshot_payload_bytes_each=payload, snapshot_device_estimate_bytes=payload + 6 * 2**30,
        limits='Parser/voxel checks, not constructed-solver import or physics. Actual baseline observer invariance and all snapshot fields require completed outputs. Regional edits add edges and change lower continuation; only regional dependence may be claimed.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args(); assert not a.out.exists()
    result = audit(a.package); a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result))
