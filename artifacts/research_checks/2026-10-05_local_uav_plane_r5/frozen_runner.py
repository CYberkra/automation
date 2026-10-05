"""Small bounded V4 normal-incidence check; NOT a suspended-dipole benchmark."""
import argparse
import json
import os
from pathlib import Path
import sys
import gprMax
import h5py
import numpy as np
import psutil
from bounded_windows_process import supervise
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
CASES = ('air', 'halfspace', 'slab', 'conductive', 'debye')
PYTHON = Path(sys.executable).resolve()


def save(p, value):
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def recipe(case):
    lines = ['#title: Normal-incidence Fresnel/slab diagnostic; no antenna port',
             '#domain: 12 0.05 24', '#dx_dy_dz: 0.05 0.05 0.05',
             '#time_window: 300e-9', '#omp_threads: 2', '#pml_cells: 20 0 20 20 0 20',
             '#pml_formulation: HORIPML', '#waveform: ricker 1 100e6 pulse',
             '#plane_wave_axial: 2 inf 3 10 inf 22 0 z pulse',
             '#rx: 6 0 8 upper Ex Ey Ez', '#rx: 5.9 0 8 transverse Ex Ey Ez',
             '#rx: 6 0 16 lower Ex Ey Ez']
    if case != 'air':
        sigma = .001 if case in ('conductive', 'debye') else 0
        lines += [f'#material: 9 {sigma:g} 1 0 cover', '#material: 4 0 1 0 rock']
        if case == 'debye':
            lines += ['#add_dispersion_debye: 1 1 6e-9 cover']
        lines += ['#box: 0 0 12 12 0.05 24 cover n']
        if case != 'halfspace':
            lines += ['#box: 0 0 15 12 0.05 24 rock n']
    return lines


def freeze(out):
    if out.exists() or gprMax.__version__ != '4.0.0':
        raise ValueError('new capsule and actual V4 required')
    out.mkdir(parents=True)
    groups = []
    for case in CASES:
        p = out/(case+'.in')
        p.write_text('\n'.join(recipe(case))+'\n', encoding='utf-8')
        groups.append(dict(id=case, input=str(p.resolve()), input_sha256=sha256(p)))
    package = Path(gprMax.__file__).parent
    identities = {str(p.relative_to(package)): sha256(p) for p in package.rglob('*.py')}
    identities.update({str(p.relative_to(package)): sha256(p) for p in package.rglob('*.pyd')})
    c = dict(status='FROZEN_APPROVED', approval_basis='2026-10-05 user requests quick local checks',
             scope='Normal-incidence material/source/phase checks only; no finite-height point-source certification',
             python=str(PYTHON), groups=groups, source_identities=identities,
             code_identities={str(Path(__file__).resolve()): sha256(__file__),
                              str(ROOT/'scripts/bounded_windows_process.py'): sha256(ROOT/'scripts/bounded_windows_process.py')},
             backend='CPU', precision='double', max_cases=5, no_retry=True,
             min_available_RAM_GiB=.95, max_job_commit_GiB=.8, max_case_wall_s=180,
             child_thread_environment=dict(OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2', MKL_NUM_THREADS='2'),
             max_batch_wall_s=900, max_case_output_bytes=20*2**20,
             diagnostic_tolerances=dict(complex_relative_L2=.05, maximum_phase_error_deg=5),
             tolerance_scope='Predeclared finite-grid engineering diagnostic, not field or convergence certification',
             reference='Continuous normal-incidence Fresnel/slab with exp(+jwt), epsilon*=epsilon_inf+delta/(1+jwtau)-j*sigma/(omega*epsilon0)',
             source_normalisation='Matched air incident receiver, NOT a circuit port or S21',
             physical_interfaces_m=[12,15], receiver_z_m=[8,8,16], fit_parameters=False)
    save(out/'execution_contract.json', c)


def run(out):
    cp = out/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    if (out/'execution.jsonl').exists():
        raise ValueError('attempt already consumed; no retries')
    for p, digest in c['code_identities'].items():
        if sha256(p) != digest:
            raise ValueError('code changed')
    package = Path(gprMax.__file__).parent
    for p, digest in c['source_identities'].items():
        if sha256(package/p) != digest:
            raise ValueError('runtime changed')
    available = psutil.virtual_memory().available
    if available < c['min_available_RAM_GiB']*2**30:
        save(out/'capacity_rejection.json', dict(status='PREPARED_NOT_RUN_RESOURCE_PREFLIGHT_REJECTED',
             available_RAM_bytes=available, execution_attempt_consumed=False))
        return
    record = out/'execution.jsonl'
    os.environ.update(c['child_thread_environment'])
    rows = []
    for g in c['groups']:
        p = Path(g['input'])
        if sha256(p) != g['input_sha256'] or p.read_text('utf-8').splitlines() != recipe(g['id']):
            raise ValueError('input changed')
        with record.open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(dict(event='STARTED', case=g['id'], contract_sha256=sha256(cp)))+'\n')
        result = supervise([str(PYTHON), '-m', 'gprMax', str(p), '-cpu_precision', 'double'],
                           out/g['id'], wall_s=c['max_case_wall_s'],
                           memory_bytes=round(c['max_job_commit_GiB']*2**30),
                           output_bytes=c['max_case_output_bytes'])
        save(out/(g['id']+'_supervision.json'), result)
        if result['reason'] != 'completed':
            raise RuntimeError(f"{g['id']} stopped: {result['reason']}; preserve failed attempt")
        raw = p.with_suffix('.h5')
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax']) != '4.0.0' or not np.array_equal(h.attrs['nx_ny_nz'], [240,1,480]):
                raise ValueError('native identity differs')
            if not np.isclose(h.attrs['dt'], .05/(299792458*np.sqrt(2)), rtol=1e-12):
                raise ValueError('native time differs')
            for i, xyz in enumerate(([6,0,8],[5.9,0,8],[6,0,16]), 1):
                r = h[f'rxs/rx{i}']
                if not np.array_equal(r.attrs['GridPosition'], np.rint(np.array(xyz)/.05).astype(int)):
                    raise ValueError('native receiver differs')
                for k in ('Ex','Ey','Ez'):
                    v = r[k][:]
                    if v.dtype != np.float64 or not np.isfinite(v).all():
                        raise ValueError('native finite double required')
        row = dict(case=g['id'], raw_sha256=sha256(raw), supervision=result,
                   dtype='float64', input_sha256=sha256(p))
        rows.append(row)
        with record.open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(dict(event='COMPLETED', case=g['id'], raw_sha256=sha256(raw)))+'\n')
        print(g['id'], 'completed', round(result['wall_s'],2), 's', flush=True)
    save(out/'completed_verification.json', dict(status='PASS', contract_sha256=sha256(cp), groups=rows))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('freeze','run'))
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    (freeze if args.action == 'freeze' else run)(args.out.resolve())
