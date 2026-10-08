"""Freeze/run/audit four cover controls after the dense batch completes; no retries."""
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
import line9_v401_version_controls as native_runner
from audit_line9_v401_cover_wavefield_inputs import audit as input_audit
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]
IDS = ['base_snapshot_H0', 'flat_snapshot_H0', 'far_cover_removed_H0', 'near_cover_removed_H0']


def read(path):
    return json.loads(path.read_text('utf-8'))


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def completed_parent(path, check_live=True):
    """Require actual 46 receipts and matching native audit, not a status file alone."""
    c = read(path)
    assert c['max_runs'] == len(c['groups']) == 46 and c['no_retry']
    events = [json.loads(s) for s in (path.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(path)
    assert events[-1]['status'] == 'COMPLETED' and events[-1]['traces'] == 46
    rows = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    assert len(rows) == 46 and len({e['group'] for e in rows}) == 46
    certificate = path.parent / 'completed_verification.json'
    assert events[-1]['verification_sha256'] == sha(certificate)
    audit = native_runner.audit(path, True)
    assert audit == read(certificate)
    assert {r['id']: r['native_sha256'] for r in audit['groups']} == {r['group']: r['raw_sha256'] for r in rows}
    if check_live:
        # A just-written terminal event can precede process exit. Wait for actual exit.
        roots = [str(Path(g['input']).parent).casefold() for g in c['groups']]
        for process in psutil.process_iter(['pid', 'name', 'cmdline']):
            if 'python' not in (process.info['name'] or '').lower():
                continue
            command = ' '.join(process.info['cmdline'] or []).casefold()
            assert not any(folder in command for folder in roots), 'Previous owned solver still alive'
            assert not (process.pid == events[0]['runner_pid'] and str(path.parent).casefold() in command), 'Previous runner still alive'
    return dict(contract=str(path), contract_sha256=sha(path), verification_sha256=sha(certificate),
                completed_groups=46, status='PASS_TERMINAL_RECEIPTS_NATIVE_AND_PROCESS_CHECK')


def audit_snapshot(path, iteration, manifest):
    dt = manifest['dt_s']
    with h5py.File(path) as h:
        assert h.attrs['gprMax'] == '4.0.1'
        assert int(h.attrs['iteration']) == iteration
        assert abs(float(h.attrs['time']) - iteration * dt) < 1e-20
        assert abs(float(h.attrs['magnetic_time']) - (iteration - .5) * dt) < 1e-20
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'], manifest['snapshot_shape'])
        np.testing.assert_allclose(h.attrs['origin'], manifest['snapshot_origin_m'], rtol=0, atol=1e-12)
        np.testing.assert_allclose(h.attrs['dx_dy_dz'], manifest['snapshot_spacing_m'], rtol=0, atol=1e-12)
        assert set(h.keys()) == set(manifest['fields'])
        for key in manifest['fields']:
            x = h[key][:]
            assert x.dtype == np.float64 and list(x.shape) == manifest['snapshot_shape']
            assert np.isfinite(x).all(), (path, key)
    return dict(file=str(path), iteration=iteration, sha256=sha(path))


def audit(path, completed=False):
    result = native_runner.audit(path, completed)
    c = read(path)
    m = c['study_manifest']
    assert c['max_runs'] == 4 and c['no_retry'] and [g['id'] for g in c['groups']] == IDS
    assert c['solver_entrypoint'] == 'gprmax_snapshot_cuda_entry.py'
    assert c['parent_completion']['completed_groups'] == 46
    assert sha(Path(c['package']) / 'manifest.json') == c['package_manifest_sha256']
    assert c['min_free_VRAM_GiB'] * 2**30 >= m['estimated_peak_device_bytes']
    assert c['max_lease_age_s'] == 600
    if completed:
        base = Path(c['package']) / m['baseline']['native']
        for g, row in zip(c['groups'], result['groups']):
            folder = Path(g['input']).parent
            files = sorted((folder / 'profile_snaps').glob('snap*.h5'))
            steps = m['snapshot_iterations'] if g['snapshots'] else []
            assert [p.name for p in files] == [f'snap{j:05d}.h5' for j in steps]
            row['snapshots'] = [audit_snapshot(p, j, m) for p, j in zip(files, steps)]
            row['snapshot_count'] = len(files)
            storage = read(folder / 'snapshot_storage.json')
            assert storage['snapshot_count'] == len(steps)
            assert storage['snapshot_history_bytes'] == (m['six_field_snapshot_history_bytes_per_observed_case'] if steps else 0)
            if steps:
                assert storage['native_gpu_to_host_streaming']
            row['snapshot_storage_sha256'] = sha(folder / 'snapshot_storage.json')
            if g['id'] == 'base_snapshot_H0':
                with h5py.File(base) as old, h5py.File(Path(g['input']).with_suffix('.h5')) as new:
                    for key in ['rxs/rx1/Ez', 'srcs/src1/excitation/samples']:
                        np.testing.assert_array_equal(new[key][:], old[key][:])
                        assert new[key][:].tobytes() == old[key][:].tobytes()
                row['observer_receiver_and_source_bitwise_equal'] = True
        result['status'] = 'PASS_NATIVE_FP64_ALL_SNAPSHOT_FIELDS_TIMING_AND_OBSERVER_INVARIANCE'
    return result


def collect(execution, out):
    """Copy small audited native results and receipts; keep full fields on solver host."""
    assert not out.exists()
    path = execution / 'execution_contract.json'
    verification = audit(path, True)
    assert verification == read(execution / 'completed_verification.json')
    events = [json.loads(s) for s in (execution / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[-1]['status'] == 'COMPLETED' and events[-1]['traces'] == 4
    assert events[-1]['verification_sha256'] == sha(execution / 'completed_verification.json')
    records = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    assert len(records) == 4 and {r['group']: r['raw_sha256'] for r in records} == {r['id']: r['native_sha256'] for r in verification['groups']}
    c = read(path)
    out.mkdir(parents=True)
    for g in c['groups']:
        source = Path(g['input'])
        shutil.copyfile(source.with_suffix('.h5'), out / (g['id'] + '.h5'))
        for name in ['snapshot_storage.json', 'stdout.log', 'stderr.log']:
            shutil.copyfile(source.parent / name, out / (g['id'] + '_' + name))
    for name in ['execution_contract.json', 'preflight_verification.json', 'completed_verification.json', 'execution.jsonl', 'independent_input_audit.json']:
        shutil.copyfile(execution / name, out / name)
    save(out / 'export_receipt.json', dict(status='PASS_FOUR_NATIVE_AND_FULL_REMOTE_SNAPSHOT_AUDIT_EXPORT',
        contract_sha256=sha(path), collector_sha256=sha(__file__),
        files={p.name: sha(p) for p in out.iterdir()},
        limits='Six-field snapshot files remain on solver host. This receipt is not a local independent re-read of the field histories.'))
    print(json.dumps(dict(status='PASS_EXPORT', out=str(out), groups=4)))


def freeze(package, out, parent):
    assert not out.exists(), 'Fresh isolated execution directory required'
    assert gprMax.__version__ == '4.0.1'
    inputs = input_audit(package)
    assert inputs['input_ready'] and not inputs['execution_ready']
    m = read(package / 'manifest.json')
    previous = completed_parent(parent)
    parent_contract = read(parent)
    assert m['parent_manifest_sha256'] in parent_contract['file_identities'].values(), 'Wrong predecessor manifest'
    parent_verification = read(parent.parent / 'completed_verification.json')
    baseline_receipt = next(r for r in parent_verification['groups'] if r['id'] == 'high_x19000_H0')
    assert baseline_receipt['native_sha256'] == m['baseline']['native_sha256'], 'Wrong baseline native'
    assert m['start_requires_current_46_attempt_batch_terminal_and_audited']
    nvcc, cl = shutil.which('nvcc'), shutil.which('cl')
    assert nvcc and cl, 'CUDA/MSVC environment required'
    reserve_ram = 20 + m['six_field_snapshot_history_bytes_per_observed_case'] / 2**30
    device_bytes = m['estimated_peak_device_bytes']
    hardware = supervisor.live_resources(dict(min_available_RAM_GiB=reserve_ram, min_free_VRAM_GiB=device_bytes / 2**30))
    disk_bytes = m['total_two_case_snapshot_history_bytes'] + 2 * 2**30
    assert shutil.disk_usage(package).free > disk_bytes, 'Two snapshot histories plus2GiB disk reserve required'
    groups = [dict(g, input=str((package / g['input']).resolve()), dt_s=m['dt_s']) for g in m['groups']]
    scripts = ['line9_v401_cover_wavefield_run.py', 'audit_line9_v401_cover_wavefield_inputs.py',
               'line9_v401_version_controls.py', 'hs4_station_grid_controls.py', 'hs_capsule_identity.py',
               'gprmax_cached_cuda_entry.py', 'gprmax_snapshot_cuda_entry.py']
    native = Path(gprMax.__file__).parent
    files = {str(p.resolve()): sha(p) for p in package.rglob('*') if p.is_file()}
    for p in [parent, parent.parent / 'completed_verification.json', parent.parent / 'execution.jsonl']:
        files[str(p.resolve())] = sha(p)
    contract = dict(status='FROZEN_APPROVED', expected_version='4.0.1', python=str(Path(sys.executable).resolve()),
        approval_basis=m['approval_basis'], package=str(package), package_manifest_sha256=sha(package / 'manifest.json'),
        study_manifest=m, groups=groups, max_runs=4, no_retry=True, parent_completion=previous,
        source_identities={p.relative_to(native).as_posix(): sha(p) for p in native.rglob('*')
                           if p.is_file() and p.suffix in ('.py', '.pyd', '.so', '.tmpl')},
        code_identities={str(ROOT / 'scripts' / name): sha(ROOT / 'scripts' / name) for name in scripts},
        file_identities=files, binary_identities={str(Path(p).resolve()): sha(p) for p in [sys.executable, nvcc, cl]},
        min_available_RAM_GiB=reserve_ram, min_free_VRAM_GiB=device_bytes / 2**30,
        max_owned_RSS_GiB=28, min_system_available_during_run_GiB=.5,
        max_group_wall_s=1800, max_batch_wall_s=9000, max_lease_age_s=600,
        gpu_lock=str((ROOT / 'artifacts/local_checks/hs4_gpu_exclusive.lock').resolve()),
        cancel_file=str(out / 'USER_STOP'), lease_file=str(out / 'session_heartbeat'),
        solver_entrypoint='gprmax_snapshot_cuda_entry.py', additional_solver_args=[],
        hardware_at_freeze=hardware, conservative_device_estimate_bytes=device_bytes,
        disk_budget_bytes=disk_bytes, minimum_free_disk_bytes_at_freeze=shutil.disk_usage(package).free,
        compatibility_scope='Process-local lazy snapshot size precheck; official streaming/update kernels unchanged.')
    out.mkdir(parents=True)
    (ROOT / 'artifacts/local_checks').mkdir(parents=True, exist_ok=True)
    save(out / 'execution_contract.json', contract)
    (out / 'session_heartbeat').write_text('Active authorized four-control wavefield study\n', encoding='utf-8')
    save(out / 'independent_input_audit.json', inputs)
    save(out / 'preflight_verification.json', audit(out / 'execution_contract.json'))
    print(json.dumps(dict(status=contract['status'], groups=4, contract_sha256=sha(out / 'execution_contract.json'))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze', 'run', 'verify', 'collect'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--package', type=Path)
    p.add_argument('--parent-contract', type=Path)
    p.add_argument('--export-out', type=Path)
    a = p.parse_args()
    if a.action == 'freeze':
        assert a.package is not None and a.parent_contract is not None
        freeze(a.package.resolve(), a.out.resolve(), a.parent_contract.resolve())
    elif a.action == 'run':
        supervisor.audit = audit
        supervisor.run(a.out.resolve() / 'execution_contract.json')
    elif a.action == 'collect':
        assert a.export_out is not None
        collect(a.out.resolve(), a.export_out.resolve())
    else:
        print(json.dumps(audit(a.out.resolve() / 'execution_contract.json', True)))
