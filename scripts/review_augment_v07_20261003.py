"""CPU-only review probes for commit 18422a6, not an augmentation acceptance suite.

No field CSV, full simulation batch, solver, training or test-family access.
The main-path probe supplies synthetic arrays and redirects every output into
a temporary directory. Historical conclusions apply to the reviewed source.
"""

import argparse
import contextlib
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


def record(path):
    path = Path(path)
    return {'file': path.relative_to(ROOT).as_posix(),
            'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def synthetic_main_probe():
    import augment_sim2real_v0_7 as module
    # A known pulse: replacing it with its negative must fail correlation >= .95.
    time = np.linspace(0, 600, 721)
    pulse = np.exp(-((time - 10) / 4) ** 2) * np.cos(2 * np.pi * .08 * (time - 10))
    signals = np.tile(pulse, (33, 1))
    fake_field = np.ones((33, 501))  # No real field data are loaded.
    stdout = io.StringIO()
    figures = []
    with tempfile.TemporaryDirectory(prefix='augment_v07_review_',
                                     dir=ROOT / 'artifacts/local_checks') as temp:
        destination = Path(temp) / 'metrics.json'
        destination.write_text('historical-sentinel', encoding='utf-8')
        def redirected_open(name, mode='r', **kwargs):
            require(str(name) == r'E:\automation_djh\augment_v07_metrics.json',
                    'unexpected output path from reviewed main')
            return destination.open(mode, **kwargs)
        def fake_save(_figure, name, *args, **kwargs):
            require(str(name) in (r'E:\automation_djh\fig_augment_v07_bscan.png',
                                 r'E:\automation_djh\fig_augment_v07_profile.png'),
                    'unexpected figure output path')
            figures.append(str(name))
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(module, 'load_bscan_both',
                return_value=(time, signals, signals / 2)))
            stack.enter_context(patch.object(module, 'load_line9', return_value=fake_field))
            stack.enter_context(patch.object(module, 'calibrate_q',
                return_value={'p25': .05, 'p75': .09}))
            stack.enter_context(patch.object(module, 'op2_surface_edit',
                side_effect=lambda values, *_: values.copy()))
            stack.enter_context(patch.object(module, 'op1_ringing_v06',
                side_effect=lambda values, *_: (-values, 0)))
            stack.enter_context(patch.object(module, 'env_profile', return_value=np.ones(5)))
            stack.enter_context(patch.object(module, 'lat_env_diffrms', return_value=.18))
            stack.enter_context(patch.object(module, 'open', redirected_open, create=True))
            from matplotlib.figure import Figure
            stack.enter_context(patch.object(Figure, 'savefig', fake_save))
            stack.enter_context(contextlib.redirect_stdout(stdout))
            stack.enter_context(warnings.catch_warnings())
            warnings.simplefilter('ignore')
            # The real direct_fidelity, JSON writer and main control flow execute.
            module.main()
        result = json.loads(destination.read_text(encoding='utf-8'))
        require(result['direct_fidelity']['corr_min'] < -.99,
                'synthetic polarity reversal was not established')
        return {'kind': 'synthetic main-path control-flow probe; not physical validation',
                'actual_direct_fidelity': result['direct_fidelity'],
                'correlation_gate_failed': True,
                'main_returned_normally': True,
                'existing_output_overwritten': True,
                'figure_save_calls_after_failed_gate': len(figures),
                'stdout_finished_with_done': stdout.getvalue().rstrip().endswith('done'),
                'mocked': ['input arrays, field calibration, augmentation operators',
                           'envelope/texture diagnostics, output destinations'],
                'not_mocked': ['direct_fidelity', 'main control flow', 'JSON serialization']}


def flat_response_probe():
    import sfcw_official_loader_v0_2 as loader
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    loader.verify_official_runtime()
    flat = SimpleNamespace(frequency=np.linspace(20e6, 170e6, 501),
                           response=np.ones(501, dtype=complex))
    rows = []
    for window in ('rectangular', 'hann'):
        for padding in (1, 8):
            response = reconstruct_time_response(flat, window=window,
                                                  zero_pad_factor=padding)
            envelope = float(np.max(np.abs(response.complex_envelope)))
            signed = float(np.max(np.abs(response.real_bandpass)))
            require(abs(envelope - 1) < 1e-15 and abs(signed - 2) < 1e-15,
                    'flat-response amplitude convention differs')
            rows.append({'window': window, 'zero_pad_factor': padding,
                         'complex_envelope_abs_peak': envelope,
                         'official_real_bandpass_abs_peak': signed})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'review output already exists')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    (ROOT / 'artifacts/local_checks').mkdir(parents=True, exist_ok=True)
    archive = json.loads((ROOT / 'docs/research/augment_v07_metrics.json').read_text('utf-8'))
    profile = archive['env_profile']
    ratios = np.asarray(profile['v07']) / np.asarray(profile['line9'])
    require(np.allclose(ratios, profile['v07_over_line9'], rtol=1e-14, atol=0),
            'archived ratios do not reproduce')
    tracked = subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
    changed = subprocess.check_output(['git', 'diff', '--name-only', '9daeeae', '18422a6'],
                                     cwd=ROOT, text=True).splitlines()
    claimed = {}
    for name in ('plot_sfcw_bscan.py', 'render_sfcw_compare.py',
                 'render_v6_ringing.py', 'analyze_direct_fidelity.py'):
        paths = [p for p in tracked if Path(p).name == name]
        claimed[name] = {'tracked_paths': paths,
                         'changed_in_reviewed_commit': [p for p in paths if p in changed]}
    document = {
        'schema': 'augment-v07-review/1',
        'reviewed_commit': '18422a6686987086c3477466c1ea48999ffaeec6',
        'scope': 'source/committed metrics + synthetic control flow + pinned official flat H; '
                 'no solver, training, field data or full augmentation rerun',
        'archive_ratios_independently_recomputed': ratios.tolist(),
        'windows_outside_declared_0_8_to_1_25': [
            profile['wins'][i] for i, ratio in enumerate(ratios) if not .8 <= ratio <= 1.25],
        'claimed_repairs_in_git': claimed,
        'local_missing_full_rerun_inputs': {
            'first_b2_h5_present': (ROOT / 'artifacts/simulations/'
                '2026-10-01_B2D-C1mX-BG-CO33-t01/B2D-C1mX-BG-CO33-t01.h5').exists(),
            'hardcoded_field_csv_present': Path(
                r'E:\automation_djh\temporary product\Line9origin(36).csv').exists()},
        'synthetic_main_probe': synthetic_main_probe(),
        'official_flat_response': flat_response_probe(),
        'provenance': [record(ROOT / path) for path in (
            'scripts/augment_sim2real_v0_3.py', 'scripts/augment_sim2real_v0_4.py',
            'scripts/augment_sim2real_v0_5.py', 'scripts/augment_sim2real_v0_6.py',
            'scripts/augment_sim2real_v0_7.py', 'scripts/sfcw_official_loader_v0_2.py',
            'scripts/review_augment_v07_20261003.py',
            'docs/research/augment_v07_metrics.json',
            'docs/research/2026-10-03_sfcw_chain_review_and_repairs.md',
            'artifacts/research_checks/2026-10-02_benchmark3d_r2_acceptance/plot_sfcw_bscan.py')]
    }
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(document, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print('Review evidence saved:', args.output)


if __name__ == '__main__':
    main()
