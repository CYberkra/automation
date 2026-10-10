"""Frozen, bounded CUDA FP64 smoke and exact-input version controls."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import gprMax
import h5py
import numpy as np
import psutil

import hs4_station_grid_controls as supervisor
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]


def save(p, value):
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def audit_source(excitation, group, dt):
    """Check the source declared by this batch; legacy batches remain Ricker."""
    attrs=excitation.attrs
    kind=group.get('source_type','ricker')
    assert kind in ('ricker','impulse')
    assert attrs['WaveformType']==kind
    assert attrs['WaveformFrequency']==group.get('source_frequency_Hz',100e6)
    assert attrs['WaveformAmplitude']==group.get('source_amplitude_A',40.)
    if kind=='impulse':
        samples=excitation['samples'][:]
        np.testing.assert_array_equal(np.flatnonzero(samples),[0])
        assert samples[0]==group['source_amplitude_A']
        assert attrs['SourceStartTime']==0 and attrs['TimeSampleOffset']==.5*dt


def audit(path, completed=False):
    c = json.loads(path.read_text('utf-8'))
    assert gprMax.__version__ == c['expected_version']
    assert Path(sys.executable).resolve() == Path(c['python']).resolve()
    native = Path(gprMax.__file__).parent
    for category in ['code_identities', 'file_identities', 'binary_identities']:
        for p, digest in c[category].items():
            assert sha(p) == digest, (category, p)
    for p, digest in c['source_identities'].items():
        assert sha(native/p) == digest, p
    rows = []
    for g in c['groups']:
        row = dict(id=g['id'], input_sha256=sha(g['input']))
        if completed:
            raw = Path(g['input']).with_suffix('.h5')
            with h5py.File(raw) as h:
                assert str(h.attrs['gprMax']) == c['expected_version']
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'], g['native_shape'])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'], g.get('spacing_m', [.025]*3))
                assert abs(float(h.attrs['dt']) / g['dt_s'] - 1) < 1e-13
                assert h.attrs['Iterations'] == g['expected_samples']
                x = h['rxs/rx1/Ez'][:]; s = h['srcs/src1/excitation/samples'][:]
                assert x.dtype == s.dtype == np.float64
                assert x.shape == s.shape == (g['expected_samples'],)
                assert np.isfinite(x).all() and np.isfinite(s).all() and np.any(x)
                audit_source(h['srcs/src1/excitation'],g,float(h.attrs['dt']))
                for key, position in [('srcs/src1', g['tx_m']), ('rxs/rx1', g['rx_m'])]:
                    np.testing.assert_allclose(h[key].attrs['Position'], position, atol=1e-12, rtol=0)
            row.update(native_sha256=sha(raw), native_path=str(raw), dtype=str(x.dtype),
                       samples=len(x), dt_s=g['dt_s'])
        rows.append(row)
    return dict(status='PASS_NATIVE_IDENTITY_NOT_PHYSICAL_CERTIFICATION', completed=completed,
                contract_sha256=sha(path), version=c['expected_version'], groups=rows)


def prepare(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    if a.action == 'smoke':
        p = a.out/'smoke/profile.in'; p.parent.mkdir()
        p.write_text('''#title: V4 version smoke; not geological validation
#domain_mode: TM
#domain: 2 2.5 inf
#dx_dy_dz: 0.025 0.025 0.025
#time_window: 4e-8
#omp_threads: 2
#pml_cells: 10 10 0 10 10 0
#pml_formulation: HORIPML
#dispersive_averaging: y
#material: 11 0.001 1 0 cover
#add_dispersion_debye: 1 0.5 6.4567e-9 cover
#box: 0 0 0 2 1.5 0.025 cover
#waveform: ricker 40 100000000 pulse
#hertzian_dipole: z 0.35 1.8 0.0125 pulse
#rx: 1.65 1.85 0.0125 smoke_rx Ez
''', encoding='utf-8')
        dt = .025/(299792458*np.sqrt(2))
        groups = [dict(id='smoke', input=str(p.resolve()), native_shape=[80,100,1],
                       expected_samples=int(np.ceil(4e-8/dt))+1, dt_s=float(dt),
                       tx_m=[.35,1.8,0], rx_m=[1.65,1.85,0])]
        manifest = dict(approval_basis='User: 在两个电脑上配置好，然后继续研究; bounded FP64 smoke.',
                        limits='Small model verifies runtime only; no geological inference.')
    else:
        manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
        groups = []
        for g in manifest['groups']:
            gp = a.out/'prepared'/g['id']; gp.mkdir(parents=True)
            for key in ['input','geometry','material']:
                src = a.package/g[key]; dst = a.out/'prepared'/g[key]
                assert sha(src) == g[key+'_sha256']
                dst.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(src,dst)
            row = dict(g, input=str((a.out/'prepared'/g['input']).resolve()), dt_s=manifest['dt_s'])
            groups.append(row)
        manifest.update(approval_basis='User authorizes 4.0.1 configuration on both PCs and continued research; exact three input/geometry/material files from completed 4.0.0 batch.',
                        limits='Same station configuration columns, not a new spatial scan; no AGC, amplitude/phase fit or field validation.')
    native = Path(gprMax.__file__).parent
    nvcc, cl = shutil.which('nvcc'), shutil.which('cl')
    assert nvcc and cl, 'CUDA/MSVC environment required'
    hw = supervisor.live_resources(dict(min_available_RAM_GiB=.5, min_free_VRAM_GiB=.5))
    cells = max(int(np.prod(g['native_shape'])) for g in groups)
    estimate = cells*420 + (2**30 if a.action == 'smoke' else 2*2**30)
    assert hw['free_VRAM_bytes'] > estimate
    assert psutil.virtual_memory().available > (.5 if a.action == 'smoke' else 20)*2**30
    c = dict(status='FROZEN_APPROVED', expected_version=gprMax.__version__, python=str(Path(sys.executable).resolve()),
             approval_basis=manifest['approval_basis'], study_manifest=manifest, groups=groups, max_runs=len(groups), no_retry=True,
             source_identities={p.relative_to(native).as_posix():sha(p) for p in native.rglob('*')
                                if p.is_file() and p.suffix in ('.py','.pyd','.so','.tmpl')},
             code_identities={str(ROOT/'scripts'/name):sha(ROOT/'scripts'/name) for name in
                              ['line9_v401_version_controls.py','hs4_station_grid_controls.py','hs_capsule_identity.py','gprmax_cached_cuda_entry.py']},
             file_identities={str(p.resolve()):sha(p) for p in a.out.rglob('*') if p.is_file()},
             binary_identities={str(Path(p).resolve()):sha(p) for p in [sys.executable,nvcc,cl]},
             min_available_RAM_GiB=.5 if a.action=='smoke' else 20, min_free_VRAM_GiB=estimate/2**30,
             max_owned_RSS_GiB=1.5 if a.action=='smoke' else 24, min_system_available_during_run_GiB=.2 if a.action=='smoke' else .5,
             max_group_wall_s=600 if a.action=='smoke' else 1800, max_batch_wall_s=3600,
             gpu_lock=str((ROOT/'artifacts/local_checks/hs4_gpu_exclusive.lock').resolve()),
             cancel_file=str((a.out/'USER_STOP').resolve()), lease_file=str((a.out/'session_heartbeat').resolve()),
             max_lease_age_s=600, solver_entrypoint='gprmax_cached_cuda_entry.py', additional_solver_args=[],
             hardware_at_freeze=hw, conservative_device_estimate_bytes=estimate)
    (ROOT/'artifacts/local_checks').mkdir(parents=True, exist_ok=True)
    save(a.out/'execution_contract.json',c)
    (a.out/'session_heartbeat').write_text('Active authorized version verification\n')
    save(a.out/'preflight_verification.json',audit(a.out/'execution_contract.json'))
    print(json.dumps(dict(status=c['status'], version=c['expected_version'], groups=len(groups))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['smoke','prepare','run','verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--package', type=Path)
    a = p.parse_args(); a.out = a.out.resolve()
    if a.action in ('smoke','prepare'): prepare(a)
    elif a.action == 'run':
        supervisor.audit=audit; supervisor.run(a.out/'execution_contract.json')
    else: print(json.dumps(audit(a.out/'execution_contract.json',True)))
