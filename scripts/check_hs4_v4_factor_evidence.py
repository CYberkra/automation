"""Independent native-clock direct DFT and fixed-window result checks (no solve)."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]


def direct_response(path, frequency):
    with h5py.File(path) as h:
        excitation = h['srcs/src1/excitation']
        s = excitation['samples'][:]
        dt = float(excitation.attrs['SampleInterval'])
        source_offset = float(excitation.attrs['TimeSampleOffset'])
        dataset = h['rxs/rx1/Ey']
        y = dataset[:]
        receiver_offset = float(dataset.attrs['TimeSampleOffset'])
        if s.dtype != np.float64 or y.dtype != np.float64 or dt != dataset.attrs['SampleInterval']:
            raise ValueError('native double and matching clocks required')
        tail = round(200e-9/dt)
        y[-tail:] *= .5*(1+np.cos(np.linspace(0, np.pi, tail)))
        # Independent sum of actual staggered sample histories; no CZT/loader.
        source_matrix = np.exp(-2j*np.pi*frequency[:, None]*(source_offset+dt*np.arange(len(s))[None, :]))
        receiver_matrix = np.exp(-2j*np.pi*frequency[:, None]*(receiver_offset+dt*np.arange(len(y))[None, :]))
        return (dt*(receiver_matrix@y))/(dt*(source_matrix@s))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new independent evidence directory required')
    cases = [('factor_controls', 11), ('material_controls', 4)]
    reports = []
    for stem, expected_count in cases:
        capsule = ROOT/f'artifacts/research_checks/2026-10-04_hs4_v4_{stem}_r1'
        result_folder = ROOT/f'artifacts/research_checks/2026-10-04_hs4_v4_{stem}_results_r1'
        cp = capsule/'execution_contract.json'
        c = json.loads(cp.read_text('utf-8'))
        v = json.loads((capsule/'completed_verification.json').read_text('utf-8'))
        summary = json.loads((result_folder/'summary.json').read_text('utf-8'))
        events = [json.loads(s) for s in (capsule/'execution.jsonl').read_text('utf-8').splitlines()]
        if (v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp)
                or summary['contract_sha256'] != sha256(cp) or summary['physical_attribution_certified']
                or events[-1]['status'] != 'COMPLETED' or events[-1]['traces'] != expected_count
                or events[-1]['verification_sha256'] != sha256(capsule/'completed_verification.json')):
            raise ValueError('completed numerical evidence identities required')
        rows = {r['id']: r for r in v['groups']}
        with np.load(result_folder/'comparison_arrays.npz', allow_pickle=False) as arrays:
            take = np.array([0, 83, 167, 250, 333, 417, 500])
            frequency = arrays['frequency_hz'][take]
            dft_errors = {}
            for g in c['groups']:
                path = capsule/g['id']/'profile.h5'
                if sha256(path) != rows[g['id']]['raw_sha256']:
                    raise ValueError('raw identity differs')
                key = g['id']+('_frequency_response' if stem == 'factor_controls' else '_native_frequency_response')
                saved = arrays[key][take]
                exact = direct_response(path, frequency)
                error = float(np.linalg.norm(exact-saved)/np.linalg.norm(exact))
                if error > 1e-9:
                    raise ValueError(f'phase-preserving direct DFT mismatch: {g["id"]} {error}')
                dft_errors[g['id']] = error
            time = arrays['time_ns']
            metric_errors = []
            for name, windows in summary['fixed_window_comparisons'].items():
                if stem == 'factor_controls':
                    station = name.rsplit('_', 1)[-1]
                    reference, candidate = 'baseline_'+station, name
                else:
                    reference, candidate = {'debye_vs_matched95': ('debye_n', 'matched95'),
                                            'matched95_vs_low_sigma': ('matched95', 'low_sigma'),
                                            'debye_vs_low_sigma': ('debye_n', 'low_sigma')}[name]
                for window, limits in {'underground': (160, 220), 'early': (160, 180), 'late': (180, 220)}.items():
                    mask = (time >= limits[0])&(time <= limits[1])
                    a = arrays[reference+'_signed_waveform'][mask]
                    b = arrays[candidate+'_signed_waveform'][mask]
                    # Scalar sums are independent of the analysis norm helper.
                    value = float(np.sqrt(sum(float(q)*float(q) for q in b-a)/sum(float(q)*float(q) for q in a)))
                    error = abs(value-windows[window]['signed_waveform_relative_L2'])
                    if error > 1e-12:
                        raise ValueError('fixed-window waveform metric mismatch')
                    metric_errors.append(error)
            flat_translation = {}
            if stem == 'factor_controls':
                mask = (time >= 160)&(time <= 220)
                centre = arrays['flat_centre_signed_waveform'][mask]
                for station in ('left', 'right'):
                    candidate = arrays['flat_'+station+'_signed_waveform'][mask]
                    flat_translation[station] = float(np.linalg.norm(candidate-centre)/np.linalg.norm(centre))
        reports.append({'study': stem, 'completed_solver_cases': expected_count,
                        'contract_sha256': sha256(cp), 'summary_sha256': sha256(result_folder/'summary.json'),
                        'native_staggered_direct_DFT_relative_L2_by_case': dft_errors,
                        'flat_three_station_translation_signed_relative_L2_to_centre': flat_translation,
                        'maximum_independent_signed_metric_absolute_difference': max(metric_errors)})
    output = {'status': 'PASS', 'code_sha256': sha256(__file__), 'solver_called': False,
              'studies': reports, 'physical_attribution_certified': False,
              'scope': 'Independent processing/identity checks, not FDTD continuum or field truth.'}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
