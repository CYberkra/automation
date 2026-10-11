"""Approved snapshot output selection through official Snapshot user objects.

The hash syntax has no component selector. Set the documented API outputs on
its parsed Snapshot objects, then reuse the existing capacity/cache wrapper.
No source, receiver, field update, or snapshot interpolation kernel is changed.
"""
from pathlib import Path
import runpy

from gprMax import hash_cmds_file
from gprMax.user_objects.cmds_output import Snapshot

original = hash_cmds_file.get_user_objects


def selected_snapshots(*args, **kwargs):
    objects = original(*args, **kwargs)
    snapshots = [obj for obj in objects if isinstance(obj, Snapshot)]
    if not snapshots:
        raise ValueError('Fine snapshot entry requires explicit snapshot inputs')
    for obj in snapshots:
        if obj.outputs is not None:
            raise ValueError('Unexpected prior snapshot output selection')
        obj.outputs = ['Ez', 'Hx', 'Hy']
    return objects


if __name__ == '__main__':
    hash_cmds_file.get_user_objects = selected_snapshots
    runpy.run_path(str(Path(__file__).with_name('gprmax_snapshot_cuda_entry.py')), run_name='__main__')
