"""Read-only HS review; replay/mutation writes occur only in fresh local copies.

Run with the existing gprMax Python environment. No solver invocation.
Default report goes to ignored local_checks; --output must name a NEW file.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np
from scipy.signal import hilbert

from sfcw_official_loader_v0_2 import signed_official, time_ns, trace_time_response

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'
ACCEPTANCE = ROOT / 'scripts/run_hs_acceptance_v0_1.py'
NAMES = ('hs1_flat_halfspace', 'hs2_coveronly_halfspace', 'hs3_domain16_halfspace')
WINDOWS = {'direct': (0, 10), 'ground': (90, 115), 'interface': (160, 220)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def replay(work, name, scale=None):
    dest = work / name
    shutil.copytree(CAPSULE, dest)
    if scale is not None:
        with h5py.File(dest / (NAMES[2] + '.h5'), 'r+') as h:
            data = h['rxs/rx1/Ex']
            data[:] = data[:] * scale
    outputs = ['hs_acceptance_metrics.json', 'hs_sfcw_official_traces.npz']
    before = {p: sha(dest / p) for p in outputs}
    result = subprocess.run(
        [sys.executable, str(ACCEPTANCE), '--dir', str(dest),
         '--fig', str(dest / 'review_replay.png')],
        cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
    (dest / 'stdout.txt').write_text(result.stdout, encoding='utf-8')
    (dest / 'stderr.txt').write_text(result.stderr, encoding='utf-8')
    require(result.returncode == 0, 'Replay failed: ' + result.stderr)
    return {
        'exit_code': result.returncode,
        'metrics': read_json(dest / outputs[0]),
        'existing_outputs_changed': {p: before[p] != sha(dest / p) for p in outputs},
        'input_manifest_matches': sha(dest / (NAMES[2] + '.h5')) ==
            read_json(dest / 'manifest.json')[NAMES[2] + '.h5'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output is not None:
        require(not args.output.exists(), 'Review output already exists')
    local = ROOT / 'artifacts/local_checks'
    local.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='hs_review_20261002_', dir=local))
    manifest = read_json(CAPSULE / 'manifest.json')
    initial = {p.name: sha(p) for p in CAPSULE.iterdir() if p.is_file()}
    require(all(initial.get(n) == digest for n, digest in manifest.items()),
            'Original capsule has manifest mismatches')
    traces, metadata = [], {}
    for name in NAMES:
        path = CAPSULE / (name + '.h5')
        with h5py.File(path, 'r') as h:
            metadata[name] = {
                'version': str(h.attrs['gprMax']),
                'dt': float(h.attrs['dt']),
                'iterations': int(h.attrs['Iterations']),
                'grid': np.asarray(h.attrs['nx_ny_nz']).tolist(),
                'receiver_position': np.asarray(h['rxs/rx1'].attrs['Position']).tolist(),
                'dtype': str(h['rxs/rx1/Ex'].dtype),
            }
        traces.append(trace_time_response(path))
    t = time_ns(traces[0])
    require(all(np.array_equal(t, time_ns(tr)) for tr in traces), 'Axes differ')
    signed = [signed_official(tr) for tr in traces]
    with np.load(CAPSULE / 'hs_sfcw_official_traces.npz') as archived:
        require(np.array_equal(t, archived['t']), 'Archived axis differs')
        require(all(np.array_equal(s, archived[f's{i}'])
                    for i, s in enumerate(signed, 1)), 'Archived signed traces differ')
    official_envelope = 2 * np.abs(traces[0].complex_envelope)
    official_difference = 2 * np.abs(
        traces[0].complex_envelope - traces[1].complex_envelope)
    hilbert_envelope = np.abs(hilbert(signed[0]))
    hilbert_difference = np.abs(hilbert(signed[0] - signed[1]))
    comparison = {}
    for name, (lo, hi) in WINDOWS.items():
        mask = (t >= lo) & (t <= hi)
        analytic, real_hilbert = ((official_difference, hilbert_difference)
                                 if name == 'interface' else
                                 (official_envelope, hilbert_envelope))
        comparison[name] = {
            'official_peak': float(analytic[mask].max()),
            'official_peak_ns': float(t[mask][np.argmax(analytic[mask])]),
            'hilbert_peak': float(real_hilbert[mask].max()),
            'hilbert_peak_ns': float(t[mask][np.argmax(real_hilbert[mask])]),
            'max_envelope_difference': float(np.abs(analytic[mask] - real_hilbert[mask]).max()),
        }
    ratio = comparison['interface']['official_peak'] / comparison['ground']['official_peak']
    comparison['official_interface_over_ground'] = ratio
    comparison['official_interface_over_ground_dB'] = float(20 * np.log10(ratio))
    normal = replay(work, 'normal')
    changed = replay(work, 'hs3_amplitude_x2', scale=2)
    require(not any(normal['existing_outputs_changed'].values()), 'Normal replay drifted')
    require(max(changed['metrics']['hs3_vs_hs1_max_rel_diff_by_window'].values()) > 0.8,
            'Mutation did not produce the intended disagreement')
    require(changed['metrics']['hs3_vs_hs1_verdict'] == normal['metrics']['hs3_vs_hs1_verdict'],
            'Acceptance behavior changed; revisit this review finding')
    ladder_path = ROOT / 'artifacts/research_checks/2026-10-02_t3_ladder_carrierfix_v0_2_r1.json'
    ladder = read_json(ladder_path)['records']
    noise = []
    for row in ladder:
        if row['representation'] != 'legacy95' or row['damage'] != 'nc_amp':
            continue
        counterpart = next(r for r in ladder if r['representation'] == 'official20'
                           and all(r[k] == row[k] for k in ('family', 'damage', 'level')))
        noise.append({'family': row['family'], 'level': row['level'],
                      'legacy_Nb_ratio': row['Nb_ratio'],
                      'official_Nb_ratio': counterpart['Nb_ratio']})
    ref_path = ROOT / 'artifacts/research_checks/2026-10-02_s1s3_refwindow_carrierfix_v0_2_r1.json'
    rows = read_json(ref_path)['families']['S1X']['rows']
    rms = {key: float(np.mean([r[key] for r in rows]))
           for key in ('nc_rms_legacy', 'nc_rms_official', 'floor_rms_legacy', 'floor_rms_official')}
    require(all(sha(CAPSULE / n) == digest for n, digest in initial.items()),
            'Original capsule changed during review')
    result = {
        'schema': 'hs_storage_review/1',
        'reviewed_commit': 'c921881',
        'scope': 'Archived HS postprocessing only; no solver, training or C5/C8 reads',
        'capsule_hashes': initial,
        'manifest_matches': len(manifest),
        'original_capsule_unchanged': True,
        'h5_metadata': metadata,
        'archived_signed_arrays_exact': True,
        'envelope_comparison': comparison,
        'normal_replay': normal,
        'scaled_hs3_replay': changed,
        'legacy_archive_counterexamples': {'noise_Nb': noise, 'S1X_mean_RMS': rms},
        'source_hashes': {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p)
                          for p in (Path(__file__), ACCEPTANCE,
                                    ROOT / 'scripts/sfcw_official_loader_v0_2.py',
                                    ladder_path, ref_path)},
    }
    output = args.output if args.output is not None else work / 'review.json'
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({'report': str(output.resolve()), 'replicas': str(work),
                      'manifest_matches': len(manifest),
                      'mutation_exit_code': changed['exit_code'],
                      'mutation_max_relative_difference': max(changed['metrics'][
                          'hs3_vs_hs1_max_rel_diff_by_window'].values())}, indent=2))


if __name__ == '__main__':
    main()
