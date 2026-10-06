"""Process-local compatibility fix for V4 lazy snapshot memory accounting.

Snapshot.nbytes is zero until initialise_snapfields(), which occurs after the
model's memory check. Supply the actual output sizes only during that check,
then restore them so native allocation accounting is unchanged. Native V4
selects GPU-to-host storage; field updates and snapshot kernels are untouched.
"""
import json
import os
from pathlib import Path
import runpy

import numpy as np
from gprMax import config, model


def with_lazy_snapshot_sizes(original, grids, log=None):
    saved = []
    sizes = []
    try:
        for grid in grids:
            for snap in grid.snapshots:
                expected = int(np.prod(snap.grid_view.size))*np.dtype(config.sim_config.dtypes['float_or_double']).itemsize*sum(bool(v) for v in snap.outputs.values())
                if snap.nbytes not in (0, expected):
                    raise ValueError('Unexpected native snapshot memory state')
                saved.append((snap, snap.nbytes))
                sizes.append(expected)
                snap.nbytes = expected
        result = original(grids)
        streaming = bool(config.get_model_config().device['snapsgpu2cpu'])
        if sizes and config.sim_config.general['solver'] == 'cuda' and not streaming:
            raise MemoryError('Native snapshot streaming capacity check failed before allocation')
        if log is not None:
            log.write_text(json.dumps(dict(snapshot_count=len(sizes), snapshot_history_bytes=sum(sizes),
                native_gpu_to_host_streaming=streaming, native_estimated_model_bytes=result[0],
                compatibility_scope='Lazy nbytes estimate only; restored before allocation; native update kernels unchanged'), indent=2)+'\n', encoding='utf-8')
        return result
    finally:
        for snap, old in saved:
            snap.nbytes = old


if __name__ == '__main__':
    original = model.mem_check_run_all
    log = Path(os.environ['HS4_CUDA_CACHE_LOG']).with_name('snapshot_storage.json')
    model.mem_check_run_all = lambda grids: with_lazy_snapshot_sizes(original, grids, log)
    runpy.run_path(str(Path(__file__).with_name('gprmax_cached_cuda_entry.py')), run_name='__main__')
