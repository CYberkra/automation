"""One 190m H0 observer-only run, depending on audited stopped partial anchors."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import gprMax
import h5py
import numpy as np
import psutil
import hs4_station_grid_controls as supervisor
import line9_v401_version_controls as native
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz']
ROIS = dict(global_view=dict(bounds=[0, 0, 0, 210, 42.5, .025], spacing=[.25, .25, .025], shape=[840, 170, 1]),
            local_view=dict(bounds=[145, 10, 0, 180, 41.5, .025], spacing=[.1, .1, .025], shape=[350, 315, 1]))


def read(p):
    return json.loads(p.read_text('utf-8'))


def dependency(package):
    p = package / 'parent'
    c, n, v = [read(p / name) for name in ['execution_contract.json', 'continuation_execution_contract.json', 'snapshot_manifest.json']]
    analysis, audit = read(p / 'analysis.json'), read(p / 'independent_audit.json')
    assert audit['status'].startswith('PASS') and audit['analysis_sha256'] == sha(p / 'analysis.json')
    assert analysis['contract_sha256'] == v['contract_sha256'] == sha(p / 'execution_contract.json')
    assert analysis['snapshot_sha256'] == sha(p / 'snapshot_manifest.json')
    assert n['parent_contract_sha256'] == sha(p / 'execution_contract.json')
    journals = [[json.loads(s) for s in (p / name).read_text('utf-8').splitlines()] for name in ['execution.jsonl', 'continuation_execution.jsonl']]
    assert journals[0][-1] == dict(status='FAILED', error='SESSION_LEASE_EXPIRED')
    assert journals[1][-1] == dict(status='FAILED', error='USER_CANCEL_REQUESTED')
    assert journals[0][0]['contract_sha256'] == sha(p / 'execution_contract.json')
    assert journals[1][0]['contract_sha256'] == sha(p / 'continuation_execution_contract.json')
    completed = [e for events in journals for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    assert len(completed) == len({e['group'] for e in completed}) == v['completed_new'] == 37
    assert {e['group']: e['raw_sha256'] for e in completed} == {r['id']: r['native_sha256'] for r in v['records']}
    records = {r['id']: r['native_sha256'] for r in analysis['native']}
    sid = 'high_x19000_H0'
    assert records[sid] == next(e['raw_sha256'] for e in completed if e['group'] == sid)
    baseline = package / 'baseline/native_H0.h5'
    assert sha(baseline) == records[sid]
    return dict(status='PASS_TERMINAL_PARTIAL_RELEVANT_NATIVE_AND_SFCW', completed_new=37,
        baseline_native_sha256=sha(baseline), analysis_sha256=sha(p/'analysis.json'),
        audit_sha256=sha(p/'independent_audit.json'), continuation_contract_sha256=sha(p/'continuation_execution_contract.json'),
        limits='Related190m H0 dependency only; no assertion of full46 completion or permission to restart prior attempts.')


def observations(dt):
    steps = [j for j in range(0, 20352, 34) if j*dt <= 500e-9]
    assert len(steps) == 250
    return steps


def audit_input(package):
    m = read(package/'manifest.json')
    dep = dependency(package)
    g = m['groups'][0]
    assert len(m['groups']) == 1 and g['id'] == 'observer_H0' and m['no_retry']
    for key in ['input', 'geometry', 'material']:
        assert sha(package/g[key]) == g[key+'_sha256']
    card = (package/g['input']).read_text('utf-8').splitlines()
    old = (package/'baseline/profile.in').read_text('utf-8').splitlines()
    assert [s for s in card if not s.startswith('#snapshot:')] == old
    assert sha(package/g['geometry']) == sha(package/'baseline/geometry.h5')
    assert sha(package/g['material']) == sha(package/'baseline/materials.json')
    parent = read(package/'parent/execution_contract.json')
    origin = next(x for x in parent['groups'] if x['id'] == 'high_x19000_H0')
    for key in ['input', 'geometry', 'material']:
        file = package/'baseline'/dict(input='profile.in', geometry='geometry.h5', material='materials.json')[key]
        assert sha(file) == origin[key+'_sha256']
    wanted = []
    for name, roi in ROIS.items():
        for j in observations(m['dt_s']):
            wanted.append('#snapshot: '+' '.join(map(str, roi['bounds']+roi['spacing']))+f' {j} {name}_{j:05d}.h5')
    assert [s for s in card if s.startswith('#snapshot:')] == wanted
    from gprMax.hash_cmds_file import get_user_objects
    get_user_objects(card, input_dir=(package/g['input']).parent)
    with h5py.File(package/g['geometry']) as h:
        assert h['data'].shape == (8400, 1700, 1)
        assert set(np.unique(h['data'][:])) == {0, 1, 2}
    assert m['snapshot_history_bytes'] == sum(np.prod(r['shape'])*250*6*8 for r in ROIS.values()) == 3036600000
    return dict(status='PASS_OBSERVER_ONLY_INPUT_AND_PARTIAL_DEPENDENCY', manifest_sha256=sha(package/'manifest.json'),
        dependency=dep, snapshots=500, history_bytes=m['snapshot_history_bytes'], solver_called=False)


def audit(path, completed=False):
    result = native.audit(path, completed)
    c = read(path); package = Path(c['package']); m = c['study_manifest']
    assert c['max_runs'] == 1 and c['no_retry'] and c['max_lease_age_s'] == 600
    assert m == read(package/'manifest.json')
    audit_input(package)
    if completed:
        row = result['groups'][0]; raw = Path(row['native_path'])
        with h5py.File(raw) as new, h5py.File(package/'baseline/native_H0.h5') as old:
            for key in ['rxs/rx1/Ez', 'srcs/src1/excitation/samples']:
                assert new[key][:].tobytes() == old[key][:].tobytes(), 'Observer changed native trace'
        folder = raw.parent/'profile_snaps'; expected = []
        records = []
        for name, roi in ROIS.items():
            for j in observations(m['dt_s']):
                snap = folder/f'{name}_{j:05d}.h5'; expected.append(snap.name)
                with h5py.File(snap) as h:
                    assert h.attrs['gprMax'] == '4.0.1' and h.attrs['iteration'] == j
                    assert abs(h.attrs['time']-j*m['dt_s']) < 1e-20
                    assert abs(h.attrs['magnetic_time']-(j-.5)*m['dt_s']) < 1e-20
                    np.testing.assert_array_equal(h.attrs['nx_ny_nz'], roi['shape'])
                    np.testing.assert_allclose(h.attrs['origin'], roi['bounds'][:3], atol=1e-12, rtol=0)
                    np.testing.assert_array_equal(h.attrs['dx_dy_dz'], roi['spacing'])
                    assert set(h.keys()) == set(FIELDS)
                    for field in FIELDS:
                        x = h[field][:]
                        assert x.dtype == np.float64 and list(x.shape) == roi['shape'] and np.isfinite(x).all()
                records.append(dict(roi=name, iteration=j, file=str(snap), sha256=sha(snap)))
        assert sorted(p.name for p in folder.iterdir()) == sorted(expected)
        storage = read(raw.parent/'snapshot_storage.json')
        assert storage['native_gpu_to_host_streaming'] and storage['snapshot_count'] == 500
        assert storage['snapshot_history_bytes'] == m['snapshot_history_bytes']
        row.update(observer_receiver_and_source_bitwise_equal=True, snapshots=records, snapshot_count=500,
            snapshot_storage_sha256=sha(raw.parent/'snapshot_storage.json'))
        result['status'] = 'PASS_SINGLE_NATIVE_FP64_DUAL_ROI500_AND_BITWISE_OBSERVER_INVARIANCE'
    return result


def freeze(package, out):
    assert not out.exists()
    inputs = audit_input(package)
    m = read(package/'manifest.json')
    assert gprMax.__version__ == '4.0.1'
    # All predecessor-owned solver roots must be idle; do not terminate unrelated work.
    for proc in psutil.process_iter(['name', 'cmdline']):
        command = ' '.join(proc.info['cmdline'] or [])
        if 'python' in (proc.info['name'] or '').lower():
            assert not any(root in command for root in ['v401_dense_loss_5m_r1', 'v401_dense_continuation_r1', 'v401_cover_wavefield_r4'])
    nvcc, cl = shutil.which('nvcc'), shutil.which('cl'); assert nvcc and cl
    hardware = supervisor.live_resources(dict(min_available_RAM_GiB=24, min_free_VRAM_GiB=9))
    assert shutil.disk_usage(package).free > m['snapshot_history_bytes']+2*2**30
    native_source = Path(gprMax.__file__).parent
    scripts = ['line9_single_wavefield.py', 'line9_v401_version_controls.py', 'hs4_station_grid_controls.py',
        'hs_capsule_identity.py', 'gprmax_cached_cuda_entry.py', 'gprmax_snapshot_cuda_entry.py']
    c = dict(status='FROZEN_APPROVED', expected_version='4.0.1', python=str(Path(sys.executable).resolve()),
        approval_basis='User: 做吧, one190m high-loss H0 A-scan plus global/local wavefield; no other solver cases.',
        package=str(package.resolve()), study_manifest=m,
        groups=[dict(g, input=str((package/g['input']).resolve()), dt_s=m['dt_s']) for g in m['groups']],
        max_runs=1, no_retry=True, parent_dependency=inputs['dependency'],
        source_identities={p.relative_to(native_source).as_posix():sha(p) for p in native_source.rglob('*') if p.is_file() and p.suffix in ('.py','.pyd','.so','.tmpl')},
        code_identities={str(ROOT/'scripts'/name):sha(ROOT/'scripts'/name) for name in scripts},
        file_identities={str(p.resolve()):sha(p) for p in package.rglob('*') if p.is_file()},
        binary_identities={str(Path(p).resolve()):sha(p) for p in [sys.executable,nvcc,cl]},
        min_available_RAM_GiB=24, min_free_VRAM_GiB=9, max_owned_RSS_GiB=28,
        min_system_available_during_run_GiB=.5, max_group_wall_s=1800, max_batch_wall_s=2100,
        max_lease_age_s=600, gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock',
        cancel_file=str(out/'USER_STOP'), lease_file=str(out/'session_heartbeat'),
        solver_entrypoint='gprmax_snapshot_cuda_entry.py', additional_solver_args=[], hardware_at_freeze=hardware)
    out.mkdir(parents=True)
    native.save(out/'execution_contract.json', c)
    (out/'session_heartbeat').write_text('Authorized single observer task\n')
    native.save(out/'input_audit.json', inputs)
    native.save(out/'preflight_verification.json', audit(out/'execution_contract.json'))
    print(json.dumps(dict(status='FROZEN_SINGLE_OBSERVER', contract_sha256=sha(out/'execution_contract.json'))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze', 'run', 'verify', 'input-audit'])
    p.add_argument('--package', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.action == 'freeze': freeze(a.package.resolve(), a.out.resolve())
    elif a.action == 'run':
        supervisor.audit = audit
        supervisor.run(a.out.resolve()/'execution_contract.json')
    elif a.action == 'verify': print(json.dumps(audit(a.out.resolve()/'execution_contract.json', True)))
    else: native.save(a.out, audit_input(a.package.resolve()))
