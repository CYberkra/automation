"""CPU-only real memory-selection and fail-closed control-audit checks; no solver."""
import argparse
import copy
from contextlib import contextmanager
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

import gprMax
import h5py
import numpy as np
from gprMax.utilities import host_info

import line9_v401_cover_wavefield_run as runner
from gprmax_snapshot_cuda_entry import with_lazy_snapshot_sizes
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_cover_wavefield import metrics


@contextmanager
def temporarily(obj, name, value):
    """Scoped fixture replacement; avoids Windows-blocked asyncio test imports."""
    old = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, old)


def reject(function):
    try:
        function()
    except (AssertionError, ValueError, KeyError, MemoryError):
        return
    raise AssertionError('Unsafe fixture accepted')


def main(out):
    assert not out.exists() and gprMax.__version__ == '4.0.1'
    checks = []
    model = SimpleNamespace(device={'dev': SimpleNamespace(total_memory=lambda: 16 * 2**30,
                                                         name=lambda: '16GiB capacity fixture'),
                                    'snapsgpu2cpu': False})
    with temporarily(host_info.config, 'sim_config', SimpleNamespace(general={'solver': 'cuda'}, dtypes={'float_or_double': 'float64'})), \
         temporarily(host_info.config, 'get_model_config', lambda: model):
        snaps = [SimpleNamespace(nbytes=0, grid_view=SimpleNamespace(size=(350, 315, 1)),
                 outputs={k: True for k in ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz']}) for _ in range(250)]
        grid = SimpleNamespace(snapshots=snaps)

        def native_capacity(grids):
            payload = sum(s.nbytes for g in grids for s in g.snapshots)
            assert payload == 1_323_000_000
            host_info.mem_check_device_snaps(6 * 2**30 + payload, payload)
            return (6 * 2**30 + payload, [])

        with_lazy_snapshot_sizes(native_capacity, [grid])
        assert model.device['snapsgpu2cpu'] and all(s.nbytes == 0 for s in snaps)
        checks.append('Installed4.0.1 official capacity selection streams250six-field frames; lazy sizes restored')
        model.device['snapsgpu2cpu'] = False

        def excessive_capacity(grids):
            payload = sum(s.nbytes for g in grids for s in g.snapshots)
            host_info.mem_check_device_snaps(17 * 2**30 + payload, payload)
            return (17 * 2**30 + payload, [])

        reject(lambda: with_lazy_snapshot_sizes(excessive_capacity, [grid]))
        assert all(s.nbytes == 0 for s in snaps)
        checks.append('Oversized non-snapshot model rejected; lazy sizes restored on failure')
        with_lazy_snapshot_sizes(lambda grids: (6 * 2**30, []), [SimpleNamespace(snapshots=[])])
        checks.append('Same wrapper accepts non-observer controls without streaming requirement')
    with tempfile.TemporaryDirectory(prefix='line9_cover_audit_') as temp:
        folder = Path(temp)
        p = folder / 'snapshot.h5'
        manifest = dict(dt_s=5.896635841874211e-11, snapshot_shape=[2, 3, 1],
                        snapshot_origin_m=[145., 10., 0.], snapshot_spacing_m=[.1, .1, .025],
                        fields=['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'])
        with h5py.File(p, 'w') as h:
            h.attrs.update(gprMax='4.0.1', iteration=34, time=34 * manifest['dt_s'], magnetic_time=33.5 * manifest['dt_s'],
                           nx_ny_nz=manifest['snapshot_shape'], origin=manifest['snapshot_origin_m'],
                           dx_dy_dz=manifest['snapshot_spacing_m'])
            for key in manifest['fields']:
                h[key] = np.zeros(manifest['snapshot_shape'], np.float64)
        runner.audit_snapshot(p, 34, manifest)
        checks.append('Six-field fixture accepted at native electric/magnetic staggered time levels')
        for key, bad in [('gprMax', '4.0.0'), ('time', 0.), ('magnetic_time', 34 * manifest['dt_s']), ('origin', [0., 0., 0.])]:
            with h5py.File(p, 'r+') as h:
                old = h.attrs[key]
                h.attrs[key] = bad
            reject(lambda: runner.audit_snapshot(p, 34, manifest))
            with h5py.File(p, 'r+') as h:
                h.attrs[key] = old
            checks.append('Snapshot rejects wrong ' + key)
        with h5py.File(p, 'r+') as h:
            h['Ez'][0, 0, 0] = np.inf
        reject(lambda: runner.audit_snapshot(p, 34, manifest))
        checks.append('Snapshot rejects non-finite field')
        with h5py.File(p, 'r+') as h:
            del h['Ez']
            h['Ez'] = np.zeros(manifest['snapshot_shape'], np.float32)
        reject(lambda: runner.audit_snapshot(p, 34, manifest))
        checks.append('Snapshot rejectsFP32 instead of nativeFP64')
        with h5py.File(p, 'r+') as h:
            del h['Ez']
        reject(lambda: runner.audit_snapshot(p, 34, manifest))
        checks.append('Snapshot rejects missing component')
        cp = folder / 'execution_contract.json'
        c = dict(max_runs=46, no_retry=True, groups=[dict(id=f'g{i}', input=str(folder / f'g{i}/profile.in')) for i in range(46)])
        runner.save(cp, c)
        cert = dict(status='MOCK_NATIVE_CERTIFICATE_NOT_REAL_EXECUTION', groups=[dict(id=f'g{i}', native_sha256=str(i)) for i in range(46)])
        vp = folder / 'completed_verification.json'
        runner.save(vp, cert)
        events = [dict(status='STARTED', contract_sha256=sha(cp), runner_pid=987654)]
        events += [dict(status='COMPLETED', group=f'g{i}', raw_sha256=str(i)) for i in range(46)]
        events += [dict(status='COMPLETED', traces=46, verification_sha256=sha(vp))]
        log = folder / 'execution.jsonl'

        def write_events(rows):
            log.write_text(''.join(json.dumps(e) + '\n' for e in rows), encoding='utf-8')

        with temporarily(runner.native_runner, 'audit', lambda *args: cert), temporarily(runner.psutil, 'process_iter', lambda *args: []):
            write_events(events)
            runner.completed_parent(cp)
            checks.append('Parent requires46unique receipts matching mocked native certificate')
            write_events(events[:-1])
            reject(lambda: runner.completed_parent(cp))
            checks.append('Running parent rejected even with46trace files/receipts')
            changed = copy.deepcopy(events)
            changed[46]['group'] = 'g0'
            write_events(changed)
            reject(lambda: runner.completed_parent(cp))
            checks.append('Duplicate/missing parent receipt rejected')
            changed = copy.deepcopy(events)
            changed[1]['raw_sha256'] = 'wrong'
            write_events(changed)
            reject(lambda: runner.completed_parent(cp))
            checks.append('Receipt/native identity mismatch rejected')
            write_events(events)
            live = SimpleNamespace(pid=987654, info=dict(name='python.exe', cmdline=['python', str(folder), 'run']))
            with temporarily(runner.psutil, 'process_iter', lambda *args: [live]):
                reject(lambda: runner.completed_parent(cp))
            checks.append('Terminal event with still-live owned runner rejected')
    time_ns = np.arange(1101.)
    pulse = np.exp(-((time_ns - 350) / 150)**2).astype(complex)
    profiles = np.column_stack([pulse, 2 * pulse, -pulse, .5j * pulse])
    for gate in metrics(profiles, time_ns).values():
        np.testing.assert_allclose([r['norm_over_base'] for r in gate['rows']], [1., 2., 1., .5], atol=1e-13, rtol=0)
        np.testing.assert_allclose([r['difference_over_base'] for r in gate['rows']], [0., 1., 2., np.sqrt(1.25)], atol=1e-13, rtol=0)
        assert gate['rows'][2]['correlation_with_base'] > 1 - 1e-13
    checks.append('Four fixed-window metrics preserve amplitude and complex differences,including opposite-phase2x residual')
    out.mkdir(parents=True)
    runner.save(out / 'checks.json', dict(status='PASS_CPU_MEMORY_AND_EXECUTION_AUDIT_GUARDS',
        native_version=gprMax.__version__, checks=checks, check_count=len(checks), calls_solver=False,
        calls_CUDA=False, auditor_sha256=sha(__file__), runner_sha256=sha(runner.__file__),
        native_memory_selection_sha256=sha(host_info.__file__),
        snapshot_compatibility_sha256=sha(Path(__file__).with_name('gprmax_snapshot_cuda_entry.py')),
        analyzer_sha256=sha(Path(__file__).with_name('analyze_line9_v401_cover_wavefield.py')),
        scope='Actual installed4.0.1 memory selection with mocked16GiB capacity; tiny HDF5 and mocked parent execution fixtures. No native allocation, FDTD, completed46claim or physical validation.'))
    print(json.dumps(dict(status='PASS', checks=len(checks), calls_solver=False)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args().out.resolve())
