"""Independent direct-DFT check of local wavefield spectral diagnostics.

Uses archived gathers rather than solver calls or the analysis FFT helpers.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT/'artifacts/research_checks/2026-10-04_hs4_local_'
ROWS = ('cover_z11.975', 'air_z16.025', 'air_z19.925')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new evidence directory required')
    capsule = Path(str(PREFIX)+'wavefield_validation_r1')
    cp = capsule/'execution_contract.json'
    contract = json.loads(cp.read_text('utf-8'))
    verification = json.loads((capsule/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (capsule/'execution.jsonl').read_text('utf-8').splitlines()]
    if not (verification['status'] == 'PASS' and verification['completed']
            and verification['contract_sha256'] == sha256(cp)
            and events[-1]['status'] == 'COMPLETED'
            and events[-1]['verification_sha256'] == sha256(capsule/'completed_verification.json')):
        raise ValueError('completed execution and audit identities required')
    completed = [e for e in events if e['status'] == 'COMPLETED' and 'group' in e]
    if len(completed) != 3 or [e['group'] for e in completed] != [g['id'] for g in contract['groups']]:
        raise ValueError('three completed cases required')
    for event, group in zip(completed, verification['groups']):
        if event['raw_sha256'] != group['raw_sha256'] or sha256(capsule/group['id']/'profile.h5') != group['raw_sha256']:
            raise ValueError('native output identity mismatch')
    windows = []
    for suffix in ('aircone_r1', 'aircone_w80_180', 'aircone_w100_180'):
        folder = Path(str(PREFIX)+suffix)
        summary = json.loads((folder/'summary.json').read_text('utf-8'))
        if summary['contract_sha256'] != sha256(cp) or summary['physical_attribution_certified']:
            raise ValueError('contract mismatch or unsupported physical promotion')
        with np.load(folder/'spectral_arrays.npz', allow_pickle=False) as arrays:
            ts = arrays['time_ns']*1e-9
            frequency, kx = arrays['frequency_hz'], arrays['kx_rad_per_m']
            nt, nx = len(ts), len(kx)
            dt, dx = ts[1]-ts[0], contract['snapshot_spacing_m'][0]
            # Direct matrix DFT at saved physical frequency/wavenumber bins.
            temporal = np.exp(-2j*np.pi*frequency[:, None]*dt*np.arange(nt)[None, :])
            spatial = np.exp(-1j*dx*np.arange(nx)[:, None]*kx[None, :])
            taper = np.hanning(nt)[:, None]*np.hanning(nx)[None, :]
            cone = abs(kx[None, :]) <= 2*np.pi*frequency[:, None]/299792458.
            np.testing.assert_array_equal(cone, arrays['air_cone'])
            old_mask = abs(kx) > 2*np.pi*95e6/299792458.
            current, fixed, errors = {}, {}, {}
            for row in ROWS:
                square = abs(temporal @ (arrays[row+'_gather']*taper) @ spatial)**2
                saved = arrays[row+'_field_spectral_square']
                error = float(np.linalg.norm(square-saved)/np.linalg.norm(saved))
                if error > 1e-12:
                    raise ValueError(f'direct DFT and FFT differ: {row}, {error}')
                errors[row] = error
                current[row] = float(square[~cone].sum()/square.sum())
                fixed[row] = float(square[:, old_mask].sum()/square.sum())
                np.testing.assert_allclose(current[row], summary['field_spectral_square_fraction_outside_air_cone'][row], rtol=1e-12)
            windows.append({'analysis': str(folder.relative_to(ROOT)), 'summary_sha256': sha256(folder/'summary.json'),
                            'arrays_sha256': sha256(folder/'spectral_arrays.npz'), 'window_ns': summary['window_ns'],
                            'direct_DFT_vs_FFT_relative_L2': errors,
                            'per_frequency_air_cone_fraction': current, 'fixed95MHz_fraction': fixed})
    script = ROOT/'scripts/analyze_hs4_height8m_surface_filter.py'
    guards = []
    for window in (['170', '90'], ['nan', '170'], ['-1', '170']):
        pending = ROOT/'artifacts/local_checks/hs4_invalid_window_must_not_exist'
        result = subprocess.run([sys.executable, str(script), '--capsule', str(capsule), '--out', str(pending),
                                 '--window-ns', *window], capture_output=True, text=True)
        if result.returncode == 0 or pending.exists() or 'finite increasing nonnegative time window required' not in result.stderr:
            raise ValueError('invalid-window guard failed')
        guards.append(window)
    result = subprocess.run([sys.executable, str(script), '--capsule', str(capsule), '--out', str(Path(str(PREFIX)+'aircone_r1'))],
                            capture_output=True, text=True)
    if result.returncode == 0 or 'historical results must not be overwritten' not in result.stderr:
        raise ValueError('historical overwrite guard failed')
    report = {'status': 'PASS', 'code_sha256': sha256(__file__), 'solver_called_by_this_check': False,
              'actual_completed_solver_cases': 3, 'case_runtime': completed, 'window_results': windows,
              'invalid_window_rejections': guards, 'historical_overwrite_rejected': True,
              'physical_attribution_certified': False}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
