"""Prepare, freeze and audit one passive full-domain Line9 wavefield replay."""
import argparse
import copy
import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json
import run_line9_2d_first as base

ROOT = Path(__file__).resolve().parents[1]
HOP = 34
SHAPE = [2000, 375, 1]
SPACING = [.2, .2, .025]
FIELDS = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz']


def prepare(package, out, evidence, reference_id='full2d_pilot_1'):
    from gprMax.hash_cmds_file import get_user_objects
    if out.exists() or evidence.exists():
        raise ValueError('Fresh snapshot package/evidence required')
    m = json.loads((package/'manifest.json').read_text('utf-8'))
    case = next(c for c in m['cases'] if c['id'] == reference_id)
    if reference_id not in m['full2d_pilot_ids']:
        raise ValueError('Use a declared full-domain pilot station')
    original = package/case['input']
    if digest(original) != case['input_sha256']:
        raise ValueError('Original station input changed')
    task = json.loads(base.TASK.read_text('utf-8'))
    for name, expected in [('full2d.h5', task['invariants']['geometry_sha256']),
                           ('line9_research_materials_v1.json', task['invariants']['material_database_sha256'])]:
        source = package/'geometries'/name
        if digest(source) != expected:
            raise ValueError('Original geometry/material changed')
    out.mkdir(parents=True)
    (out/'geometries').mkdir()
    for name in ['full2d.h5', 'line9_research_materials_v1.json']:
        shutil.copyfile(package/'geometries'/name, out/'geometries'/name)
    shutil.copyfile(package/'contacts_m.json', out/'contacts_m.json')
    dt = case['dt_CFL_initial_s']
    iterations = round(1.2e-6/dt)+1
    timeline = list(range(0, iterations, HOP))
    suffix = ''.join(f'#snapshot: 0 0 0 400 75 0.025 0.2 0.2 0.025 {j} snap{j:05d}.h5\n' for j in timeline)
    path = out/case['input']
    path.parent.mkdir(parents=True)
    original_bytes = original.read_bytes()
    if not original_bytes.endswith(b'\n'):
        raise ValueError('Expected original terminal newline')
    path.write_bytes(original_bytes+suffix.encode('ascii'))
    objects = get_user_objects(path.read_text('utf-8').splitlines(), input_dir=path.parent)
    if len(objects) != case['command_parser_object_count']+len(timeline):
        raise ValueError('Snapshot commands did not parse')
    history = int(np.prod(SHAPE))*6*8*len(timeline)
    gpu_buffer = int(np.prod(SHAPE))*6*8
    plan = dict(status='PREPARED_NOT_RUN', calls_solver=False, case=case,
        input_sha256=digest(path), source_input_sha256=digest(original),
        package_manifest_sha256=digest(package/'manifest.json'), task_sha256=digest(base.TASK),
        snapshot_iterations=timeline, snapshot_shape=SHAPE, snapshot_spacing_m=SPACING,
        dt_s=dt, native_iterations=iterations, snapshot_period_s=HOP*dt,
        six_field_history_bytes=history, six_field_GPU_single_buffer_bytes=gpu_buffer,
        # CUDA dtoh retains six host arrays per frame even if fewer are written.
        minimum_available_RAM_bytes=task['limits']['min_available_RAM_bytes']+int(np.ceil(history*1.3)),
        minimum_free_VRAM_bytes=task['limits']['min_free_VRAM_bytes']+int(np.ceil(gpu_buffer*1.3)),
        minimum_free_disk_bytes=int(np.ceil(history*1.15))+3*2**30,
        native_command_parser_count=len(objects), passive_receiver_equivalence='NOT_TESTED',
        approval_basis='User: 大域的等会再跑一个带波场快照的版本',
        scope=f'One full400x75m replay atX{case["profile_x_m"]}; native2.5cm grid unchanged; snapshot0.2m/~2ns is visualization sampling, not solver coarsening or a full-source-spectrum phase certificate.')
    save_json(out/'snapshot_plan.json', plan)
    evidence.mkdir(parents=True)
    public = {k: v for k, v in plan.items() if k not in ['case', 'snapshot_iterations']}
    public['snapshot_count'] = len(timeline)
    public['script_sha256'] = digest(Path(__file__))
    save_json(evidence/'snapshot_preparation.json', public)
    print(f'Prepared {len(timeline)} full-domain frames; {history/2**30:.2f}GiB six-field history; no solver')


def audit(path, completed=False):
    result = base.audit(path, completed)
    c = json.loads(path.read_text('utf-8'))
    g = c['groups'][0]
    p = Path(g['input'])
    original = Path(c['passive_reference_input'])
    if digest(original) != c['passive_reference_input_sha256'] or digest(Path(c['passive_reference_raw'])) != c['passive_reference_raw_sha256']:
        raise ValueError('Passive reference changed')
    expected = original.read_bytes()+''.join(f'#snapshot: 0 0 0 400 75 0.025 0.2 0.2 0.025 {j} snap{j:05d}.h5\n' for j in c['snapshot_iterations']).encode('ascii')
    if p.read_bytes() != expected:
        raise ValueError('Observer changed original input bytes beyond snapshot append')
    if completed:
        storage_path = p.parent/'snapshot_storage.json'
        storage = json.loads(storage_path.read_text('utf-8'))
        if storage['native_gpu_to_host_streaming'] is not True or storage['snapshot_count'] != len(c['snapshot_iterations']) or storage['snapshot_history_bytes'] != int(np.prod(SHAPE))*6*8*len(c['snapshot_iterations']):
            raise ValueError('Snapshot storage policy or native shape differs from frozen budget')
        with h5py.File(p.with_suffix('.h5')) as h, h5py.File(c['passive_reference_raw']) as ref:
            for name in ['rxs/rx1/Ez', 'srcs/src1/excitation/samples']:
                np.testing.assert_array_equal(h[name][:], ref[name][:])
            for name in ['dt', 'Iterations', 'nx_ny_nz', 'dx_dy_dz']:
                np.testing.assert_array_equal(h.attrs[name], ref.attrs[name])
        files = sorted(p.parent.glob('profile_snaps/snap*.h5'))
        if len(files) != len(c['snapshot_iterations']):
            raise ValueError('Incomplete snapshot outputs')
        rows = []
        for j, file in zip(c['snapshot_iterations'], files):
            with h5py.File(file) as h:
                if str(h.attrs['gprMax']) != '4.0.0' or int(h.attrs['iteration']) != j:
                    raise ValueError('Snapshot runtime/time mismatch')
                if float(h.attrs['time']) != j*c['dt_s'] or float(h.attrs['magnetic_time']) != (j-.5)*c['dt_s']:
                    raise ValueError('Snapshot E/H time staggering mismatch')
                np.testing.assert_array_equal(h.attrs['origin'], [0., 0., 0.])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'], SPACING)
                for name in FIELDS:
                    v = h[name]
                    if v.dtype != np.float64 or list(v.shape) != SHAPE or not np.isfinite(v[:]).all():
                        raise ValueError('Invalid native snapshot field')
            rows.append(dict(file=str(file.resolve()), sha256=digest(file), iteration=j))
        result.update(snapshot_count=len(files), snapshots=rows, passive_receiver_bit_identical=True,
                      snapshot_storage_sha256=digest(storage_path), snapshot_storage=storage)
    return result


def freeze(prepared, out, pilots):
    if out.exists():
        raise ValueError('Fresh execution capsule required')
    plan = json.loads((prepared/'snapshot_plan.json').read_text('utf-8'))
    rc = pilots/'execution_contract.json'
    old = json.loads(rc.read_text('utf-8'))
    base.check_files(old)
    case = plan['case']
    original_group = next(g for g in old['groups'] if g['id'] == case['id'])
    events = [json.loads(s) for s in (pilots/'execution.jsonl').read_text('utf-8').splitlines()]
    completed = [e for e in events if e.get('status') == 'COMPLETED' and e.get('group') == case['id']]
    reference_raw = Path(original_group['input']).with_suffix('.h5')
    if len(completed) != 1 or digest(reference_raw) != completed[0]['raw_sha256']:
        raise ValueError('Require one recorded completed reference solve, never an interrupted trace')
    hardware = base.resources(json.loads(base.TASK.read_text('utf-8')))
    if hardware['gpu_uuid'] != old['hardware_at_freeze']['gpu_uuid'] or hardware['driver_version'] != old['hardware_at_freeze']['driver_version']:
        raise ValueError('Reference hardware changed')
    if Path(sys.executable).resolve() != Path(old['python']).resolve():
        raise ValueError('Use the reference native runtime')
    if hardware['available_RAM_bytes'] < plan['minimum_available_RAM_bytes'] or hardware['free_VRAM_bytes'] < plan['minimum_free_VRAM_bytes']:
        raise RuntimeError('Snapshot history resource preflight rejected')
    if shutil.disk_usage(out.parent).free < plan['minimum_free_disk_bytes']:
        raise RuntimeError('Snapshot disk budget rejected')
    prepared_input = prepared/case['input']
    if digest(prepared_input) != plan['input_sha256']:
        raise ValueError('Prepared snapshot input changed')
    if original_group['input_sha256'] != plan['source_input_sha256']:
        raise ValueError('Snapshot/reference station mismatch')
    out.mkdir(parents=True)
    # Read-only view audits a completed group from an interrupted batch. It is
    # not an execution contract and cannot pass the supervisor approval check.
    view = copy.deepcopy(old)
    view.update(status='READONLY_COMPLETED_GROUP_AUDIT', groups=[original_group], max_runs=1,
                source_contract_sha256=digest(rc), purpose='Audit completed trace; no old attempt retry')
    view_path = out/'reference_audit_view.json'
    save_json(view_path, view)
    reference = base.audit(view_path, True)
    original_raw = reference['groups'][0]
    if max(plan['snapshot_iterations']) >= original_raw['iterations'] or not np.isclose(plan['dt_s'], original_raw['dt_s'], rtol=1e-12, atol=0):
        raise ValueError('Prepared snapshot timeline differs from native reference')
    save_json(out/'reference_native_verification.json', reference)
    shutil.copytree(prepared/'geometries', out/'geometries')
    shutil.copyfile(prepared/'contacts_m.json', out/'contacts_m.json')
    destination = out/case['input']
    destination.parent.mkdir(parents=True)
    shutil.copyfile(prepared_input, destination)
    c = copy.deepcopy(old)
    group = copy.deepcopy(original_group)
    group.update(input=str(destination.resolve()), input_sha256=digest(destination))
    sources = ['line9_large_domain_snapshot.py', 'run_line9_2d_first.py', 'hs4_station_grid_controls.py',
        'gprmax_cached_cuda_entry.py', 'gprmax_snapshot_cuda_entry.py', 'hs_capsule_identity.py', 'build_pdf_profile_geometry.py',
        'analyze_line9_2d_sfcw.py', 'run_line9_2d_v4.cmd', 'analyze_line9_large_wavefield.py']
    c.update(stage='snapshot', approval_basis=plan['approval_basis'], groups=[group], max_runs=1,
        reuse={}, hardware_at_freeze=hardware, additional_solver_args=[],
        solver_entrypoint='gprmax_snapshot_cuda_entry.py',
        snapshot_storage_policy='Process-local lazy nbytes accounting correction before native V4 memory check; restored before allocation; native GPU-to-host streaming and field kernels unchanged',
        code_identities={str(ROOT/'scripts'/n): digest(ROOT/'scripts'/n) for n in sources},
        file_identities={str(p.resolve()): digest(p) for p in out.rglob('*') if p.is_file()},
        passive_reference_input=original_group['input'], passive_reference_input_sha256=original_group['input_sha256'],
        passive_reference_raw=original_raw['raw_path'], passive_reference_raw_sha256=original_raw['raw_sha256'],
        snapshot_iterations=plan['snapshot_iterations'], snapshot_shape=SHAPE, dt_s=original_raw['dt_s'],
        min_available_RAM_GiB=plan['minimum_available_RAM_bytes']/2**30,
        min_free_VRAM_GiB=plan['minimum_free_VRAM_bytes']/2**30,
        max_owned_RSS_GiB=40, max_group_wall_s=3600, max_batch_wall_s=3600,
        expected_preview_total=None, snapshot_plan_sha256=digest(prepared/'snapshot_plan.json'),
        minimum_free_disk_bytes=plan['minimum_free_disk_bytes'],
        reference_source_contract_sha256=digest(rc), reference_completion_event=completed[0])
    save_json(out/'execution_contract.json', c)
    save_json(out/'preflight_verification.json', audit(out/'execution_contract.json', False))
    print('Frozen one passive large-domain replay; no solver called')


def run(out):
    import hs4_station_grid_controls as supervisor
    c = json.loads((out/'execution_contract.json').read_text('utf-8'))
    hw = base.resources(json.loads(base.TASK.read_text('utf-8')))
    if hw['gpu_uuid'] != c['hardware_at_freeze']['gpu_uuid'] or hw['driver_version'] != c['hardware_at_freeze']['driver_version']:
        raise ValueError('Frozen snapshot GPU/driver changed')
    if shutil.disk_usage(out).free < c['minimum_free_disk_bytes']:
        raise RuntimeError('Snapshot disk budget changed')
    supervisor.audit = audit
    supervisor.run(out/'execution_contract.json')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--package', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--evidence', type=Path)
    p.add_argument('--pilots', type=Path)
    p.add_argument('--reference-id', default='full2d_pilot_1')
    a = p.parse_args()
    if a.action == 'prepare':
        prepare(a.package.resolve(), a.out.resolve(), a.evidence.resolve(), a.reference_id)
    elif a.action == 'freeze':
        freeze(a.package.resolve(), a.out.resolve(), a.pilots.resolve())
    elif a.action == 'run':
        run(a.out.resolve())
    else:
        print(json.dumps(audit(a.out.resolve()/'execution_contract.json', True), indent=2))
