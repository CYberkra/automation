"""Small synthetic snapshot passivity/rejection checks; never calls a solver."""
import argparse
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json
import line9_large_domain_snapshot as wave


def check(out):
    if out.exists():
        raise ValueError('Fresh check evidence required')
    checks = []
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        original = root/'reference.in'
        original.write_text('#title: synthetic array check\n', encoding='utf-8')
        input_path = root/'profile.in'
        valid_input = original.read_bytes()+b'#snapshot: 0 0 0 400 75 0.025 0.2 0.2 0.025 0 snap00000.h5\n'
        input_path.write_bytes(valid_input)
        raw, reference = root/'profile.h5', root/'reference.h5'
        for file in [raw, reference]:
            with h5py.File(file, 'w') as h:
                h['rxs/rx1/Ez'] = np.array([1., -.25, .125])
                h['srcs/src1/excitation/samples'] = np.array([0., 1., 0.])
                h.attrs['dt'] = 1e-9
                h.attrs['Iterations'] = 3
                h.attrs['nx_ny_nz'] = [2, 3, 1]
                h.attrs['dx_dy_dz'] = [.025]*3
        directory = root/'profile_snaps'
        directory.mkdir()
        snap = directory/'snap00000.h5'
        with h5py.File(snap, 'w') as h:
            for name in wave.FIELDS:
                h[name] = np.zeros((2, 3, 1), dtype=np.float64)
            h.attrs['gprMax'] = '4.0.0'
            h.attrs['iteration'] = 0
            h.attrs['time'] = 0.
            h.attrs['magnetic_time'] = -.5e-9
            h.attrs['origin'] = [0., 0., 0.]
            h.attrs['dx_dy_dz'] = wave.SPACING
        contract = root/'audit_contract.json'
        save_json(contract, dict(groups=[dict(input=str(input_path))],
            passive_reference_input=str(original), passive_reference_input_sha256=digest(original),
            passive_reference_raw=str(reference), passive_reference_raw_sha256=digest(reference),
            snapshot_iterations=[0], dt_s=1e-9))
        with patch.object(wave.base, 'audit', return_value={'status': 'PASS'}), patch.object(wave, 'SHAPE', [2, 3, 1]):
            result = wave.audit(contract, True)
            assert result['passive_receiver_bit_identical'] and result['snapshot_count'] == 1
            checks.append('identical native receiver and source accepted')
            with h5py.File(raw, 'r+') as h:
                h['rxs/rx1/Ez'][1] = .25
            try:
                wave.audit(contract, True)
            except AssertionError:
                checks.append('receiver polarity change rejected')
            else:
                raise AssertionError('Receiver change accepted')
            with h5py.File(raw, 'r+') as h:
                h['rxs/rx1/Ez'][1] = -.25
            with h5py.File(snap, 'r+') as h:
                del h['Ez']
                h['Ez'] = np.zeros((2, 3, 1), dtype=np.float32)
            try:
                wave.audit(contract, True)
            except ValueError:
                checks.append('native float32 snapshot rejected')
            else:
                raise AssertionError('Snapshot cast accepted')
            with h5py.File(snap, 'r+') as h:
                del h['Ez']
                h['Ez'] = np.zeros((2, 3, 1), dtype=np.float64)
                h.attrs['magnetic_time'] = 0.
            try:
                wave.audit(contract, True)
            except ValueError:
                checks.append('incorrect magnetic half-step rejected')
            else:
                raise AssertionError('Wrong staggering accepted')
            input_path.write_bytes(valid_input+b'#time_window: 1e-7\n')
            try:
                wave.audit(contract, False)
            except ValueError:
                checks.append('physical input change beyond observer append rejected')
            else:
                raise AssertionError('Changed physical input accepted')
    out.mkdir(parents=True)
    save_json(out/'checks.json', dict(status='PASS_SYNTHETIC_AUDIT_NOT_FDTD', calls_solver=False,
        checks=checks, scope='Snapshot-specific checks only; base native dtype/time tests are separate.',
        script_sha256=digest(Path(__file__)), snapshot_auditor_sha256=digest(Path(wave.__file__))))
    print(f'{len(checks)} snapshot audit checks PASS; no solver')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    check(p.parse_args().out.resolve())
