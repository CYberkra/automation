"""Exercise installed V4 snapshot memory selection without allocating a GPU or running FDTD."""
import argparse
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import gprMax
from gprMax.utilities import host_info
from build_pdf_profile_geometry import digest, save_json
from gprmax_snapshot_cuda_entry import with_lazy_snapshot_sizes


def check(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('Fresh evidence and installed V4.0.0 required')
    device = SimpleNamespace(total_memory=lambda: 16*2**30, name=lambda: 'capacity fixture')
    model = SimpleNamespace(device={'dev': device, 'snapsgpu2cpu': False})
    with patch.object(host_info.config, 'sim_config', SimpleNamespace(general={'solver': 'cuda'}, dtypes={'float_or_double': 'float64'})), \
         patch.object(host_info.config, 'get_model_config', return_value=model):
        host_info.mem_check_device_snaps(30*2**30, 20*2**30)
        assert model.device['snapsgpu2cpu'] is True
        model.device['snapsgpu2cpu'] = False
        host_info.mem_check_device_snaps(40*2**30, 20*2**30)
        assert model.device['snapsgpu2cpu'] is False
        snap = SimpleNamespace(nbytes=0, grid_view=SimpleNamespace(size=(2000, 375, 1)),
                               outputs={k: True for k in ['Ex','Ey','Ez','Hx','Hy','Hz']})
        grid = SimpleNamespace(snapshots=[snap])
        def check_lazy(grids):
            assert grids[0].snapshots[0].nbytes == 36_000_000
            host_info.mem_check_device_snaps(10*2**30+snap.nbytes, snap.nbytes)
            return (10*2**30+snap.nbytes, [])
        with_lazy_snapshot_sizes(check_lazy, [grid])
        assert snap.nbytes == 0 and model.device['snapsgpu2cpu'] is True
        def fail(grids):
            raise ValueError('fixture failure')
        try:
            with_lazy_snapshot_sizes(fail, [grid])
        except ValueError:
            assert snap.nbytes == 0
        else:
            raise AssertionError('Expected check failure')
    help_result = subprocess.run([sys.executable, '-m', 'gprMax', '--help'], capture_output=True,
                                 text=True, encoding='utf-8', errors='replace', check=True)
    assert '-gpu_precision' in help_result.stdout and '-snapsgpu2cpu' not in help_result.stdout
    out.mkdir(parents=True)
    save_json(out/'checks.json', dict(status='PASS', calls_solver=False,
        native_version=gprMax.__version__, native_memory_source_sha256=digest(Path(host_info.__file__)),
        checks=['Native V4 enables streaming when non-snapshot arrays fit device',
                'Native V4 does not enable streaming for an oversized non-snapshot model',
                'Native CLI exposes GPU precision and has no legacy snapshot flag',
                'Lazy field size correction enables native streaming and restores accounting',
                'Lazy field accounting restored on failure'],
        scope='Actual installed memory-selection function and CLI only; mocked capacity, no allocation/FDTD'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    check(p.parse_args().out.resolve())
