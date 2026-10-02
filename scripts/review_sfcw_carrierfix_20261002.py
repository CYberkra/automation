"""Independent, CPU-only review of 5af6826; no solver or test-family input.

Compares submitted JSON with frozen archives and exercises the installed
official loader against a tone and a direct-sum reconstruction of one locally
available development control. Full t3/S1S3/B2 H5 reruns are not claimed.
Output defaults to ignored local_checks; exclusive creation preserves evidence.
"""

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
CHECKS = ROOT / 'artifacts/research_checks'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def audit_archives():
    provenance = {}
    submitted = {}
    names = {'t3': 't3_ladder', 's1s3': 's1s3_refwindow', 'b2': 'b2_reward'}
    for key, name in names.items():
        p1 = CHECKS / f'2026-10-02_{name}_carrierfix_v0_2_r1.json'
        p2 = CHECKS / f'2026-10-02_{name}_carrierfix_v0_2_r2.json'
        assert p1.read_bytes() == p2.read_bytes(), key
        provenance[key] = {'r1_r2_identical': True, 'sha256': digest(p1)}
        submitted[key] = read(p1)

    old_t3 = CHECKS / '2026-09-28_t3_damage_ladder_r1.json'
    old_b2 = CHECKS / '2026-10-01_reward_protocol_b2_pilot_r1.json'
    ref_path = ROOT / 'configs/research/s1s3_reference_window_contract_v0.1.json'
    tol_path = ROOT / 'configs/research/reward_tolerance_contract_v0.2.json'
    assert submitted['t3']['archived_r1_sha256'] == digest(old_t3)
    assert submitted['b2']['archived_pilot_r1_sha256'] == digest(old_b2)
    for p in (old_t3, old_b2, ref_path, tol_path):
        provenance[p.relative_to(ROOT).as_posix()] = digest(p)
    tau_d = read(tol_path)['tolerances']['tau_D']['value']

    def key(r):
        return r['family'], r['damage'], r['level']

    keys = ('D', 'a', 'A', 'H', 'rho', 'Nb_ratio',
            'arrival_drift_ns_median', 'arrival_drift_ns_max_abs')
    old = {key(r): r for r in read(old_t3)['records']}
    legacy = {key(r): r for r in submitted['t3']['records']
              if r['representation'] == 'legacy95'}
    official = {key(r): r for r in submitted['t3']['records']
                if r['representation'] == 'official20'}
    assert set(old) == set(legacy) == set(official)
    assert all(old[k][m] == legacy[k][m] for k in old for m in keys)
    shift2 = [{'family': k[0], 'D': r['D'], 'tau_D': tau_d,
               'passes_D_gate': r['D'] <= tau_d,
               'margin_above_tau_D': round(r['D'] - tau_d, 6)}
              for k, r in sorted(official.items()) if k[1:] == ('shift_smp', 2)]
    nonshift_changes = []
    for k in sorted(old):
        if k[1] == 'shift_smp':
            continue
        changed = {m: {'legacy': legacy[k][m], 'official': official[k][m]}
                   for m in keys if legacy[k][m] != official[k][m]}
        if changed:
            nonshift_changes.append({'family': k[0], 'damage': k[1],
                                     'level': k[2], 'changes': changed})
    for delta in submitted['t3']['paired_delta']:
        k = key(delta)
        for m in keys:
            a, b = legacy[k][m], official[k][m]
            expected = None if a is None or b is None else round(b - a, 6)
            assert delta[f'{m}_legacy'] == a
            assert delta[f'{m}_official'] == b
            assert delta[f'{m}_delta'] == expected

    old_rewards = {r['family']: r for r in read(old_b2)['results']}
    rank_changes, winner_changes, feasible_changes = [], [], []
    for fam in submitted['b2']['results']:
        name = fam['family']
        a, b = (fam['representations'][r] for r in ('legacy95', 'official20'))
        old_rows = {r['config']: r for r in old_rewards[name]['rows']}
        assert len(a['rows']) == len(old_rows)
        for r in a['rows']:
            for m, value in r.items():
                assert old_rows[r['config']][m] == value, (name, r['config'], m)
        for m in ('selection', 'selection_R_bg_db', 'ranking'):
            assert a[m] == old_rewards[name][m]
        if a['ranking'] != b['ranking']:
            rank_changes.append({'family': name, 'legacy': a['ranking'],
                                 'official': b['ranking']})
        if a['selection'] != b['selection']:
            winner_changes.append(name)
        assert fam['selection_flip'] == (name in winner_changes)
        by_a = {r['config']: r for r in a['rows'] if r.get('available')}
        by_b = {r['config']: r for r in b['rows'] if r.get('available')}
        assert set(by_a) == set(by_b)
        for cid in by_a:
            if by_a[cid]['feasible'] != by_b[cid]['feasible']:
                feasible_changes.append([name, cid])

    ref = read(ref_path)
    windows, moved = {}, 0
    for fam, result in submitted['s1s3']['families'].items():
        rows = result['rows']
        frozen = ref['families'][fam]
        assert len(rows) == len(frozen['t_measured_ns']) == 33
        assert all(r['t_measured_legacy_ns'] == frozen['t_measured_ns'][i]
                   == r['t_measured_frozen_ns'] for i, r in enumerate(rows))
        assert result['legacy_matches_frozen_contract'] is True
        tf = np.asarray(frozen['t_fermat_er18_ns'])
        peaks = np.asarray([r['t_measured_official_ns'] for r in rows])
        max_delta = float(np.max(np.abs(peaks - tf)))
        clearance = float(ref['constants']['nc_window_ns'][0]
                          - np.max(peaks) - ref['constants']['ev_half_ns'])
        assert max_delta < ref['constants']['ev_half_ns']
        assert clearance > 0
        count = sum(r['t_measured_official_ns'] != r['t_measured_legacy_ns']
                    for r in rows)
        assert count == result['n_traces_peak_moved']
        moved += count
        ratios = {m: [r[f'{m}_rms_official'] / r[f'{m}_rms_legacy']
                      for r in rows] for m in ('nc', 'floor')}
        windows[fam] = {'moved': count,
                        'max_abs_official_minus_fermat_ns': max_delta,
                        'measured_window_nc_clearance_ns': clearance,
                        'rms_ratio_ranges': {m: [min(v), max(v)]
                                             for m, v in ratios.items()}}
    return {'provenance': provenance, 't3_legacy_metrics_equal': len(old) * len(keys),
            't3_shift_2samples': shift2, 't3_nonshift_metric_changes': nonshift_changes,
            'b2_legacy_reported_fields_equal': True, 'b2_ranking_changes': rank_changes,
            'b2_winner_changes': winner_changes, 'b2_feasibility_changes': feasible_changes,
            's1s3_peak_moved_total': moved, 's1s3_window_recheck': windows}


def audit_actual_loader():
    import sfcw_official_loader_v0_2 as loader
    import gprMax.toolboxes.SFCW.processing as processing

    processing_path = Path(processing.__file__)
    expected = 'adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b'
    assert digest(processing_path) == expected, 'official module differs from reviewed version'
    freq = np.linspace(20e6, 170e6, 501)
    response = np.zeros(501, dtype=complex)
    response[200] = 1
    tone = processing.reconstruct_time_response(
        SimpleNamespace(frequency=freq, response=response),
        window='rectangular', zero_pad_factor=8)
    off, leg = loader.signed_official(tone), loader.signed_legacy(tone)
    tone_error = float(np.max(np.abs(off - 2 / 501 * np.cos(2 * np.pi * 80e6 * tone.time))))
    legacy_error = float(np.max(np.abs(leg - 1 / 501 * np.cos(2 * np.pi * 155e6 * tone.time))))
    assert max(tone_error, legacy_error) < 1e-14

    path = CHECKS / '2026-10-02_benchmark3d_r2_acceptance/ctl3dA_t07_asis.h5'
    # Official H5 parsers only; DFT, taper and inverse sum below are independent.
    src = processing.load_source(path)
    rx = processing.load_receiver(path, component='Ex')
    n = min(len(rx.samples), int(np.floor(1200e-9 / rx.dt)) + 1)
    raw = rx.samples[:n].copy()
    ntail = max(2, round(200e-9 / rx.dt))
    raw[-ntail:] *= (1 + np.cos(np.linspace(0, np.pi, ntail))) / 2
    def direct(samples, dt, offset):
        result = []
        t = offset + np.arange(len(samples)) * dt
        for f in freq:
            result.append(dt * np.sum(samples * np.exp(-2j * np.pi * f * t)))
        return np.asarray(result)
    h = direct(raw, rx.dt, rx.time_offset) / direct(src.samples, src.dt, src.time_offset)
    window = np.hanning(501)
    window /= window.mean()
    tr = loader.trace_time_response(path)
    indices = np.unique(np.r_[np.arange(0, 4008, 173), 0, 4007])
    expected_samples = np.asarray([
        2 / 501 * np.real(np.sum(h * window * np.exp(2j * np.pi * freq * tr.time[i])))
        for i in indices])
    actual = loader.signed_official(tr)[indices]
    err = float(np.max(np.abs(actual - expected_samples)) / np.max(np.abs(expected_samples)))
    assert err < 1e-9, err
    return {'official_processing_sha256': digest(processing_path),
            'loader_sha256': digest(ROOT / 'scripts/sfcw_official_loader_v0_2.py'),
            'python': platform.python_version(), 'numpy': np.__version__,
            'tone_official_peak_mhz': loader.dominant_freq_mhz(off, tone.time[1]),
            'tone_legacy_peak_mhz': loader.dominant_freq_mhz(leg, tone.time[1]),
            'tone_peak_amplitude_ratio': float(np.max(np.abs(off)) / np.max(np.abs(leg))),
            'tone_analytic_max_error': tone_error, 'legacy_analytic_max_error': legacy_error,
            'control_path': path.relative_to(ROOT).as_posix(),
            'control_h5_sha256': digest(path), 'control_raw_dtype': str(rx.samples.dtype),
            'control_source_time_offset_s': src.time_offset,
            'control_receiver_time_offset_s': rx.time_offset,
            'control_sample_count': len(tr.time), 'control_independent_samples': len(indices),
            'control_direct_sum_max_error_relative_to_selected_peak': err}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'artifacts/local_checks/carrierfix_rereview_20261002.json')
    args = parser.parse_args()
    doc = {'schema': 'carrierfix-independent-review/1', 'reviewed_commit': '5af6826',
           'scope': 'JSON/archive crosscheck + analytic tone + one archived dev H5; '
                    'no solver/training/test-family inputs; full three-line raw rerun unavailable',
           'archives': audit_archives(), 'actual_loader': audit_actual_loader(),
           'review_script_sha256': digest(__file__)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(doc, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print('Independent archive comparisons, tone and archived control checks passed.')
    print('2-sample D gate:', doc['archives']['t3_shift_2samples'])
    print('B2 ranking changed:', [r['family'] for r in doc['archives']['b2_ranking_changes']])
    print('Output:', args.output)


if __name__ == '__main__':
    main()
