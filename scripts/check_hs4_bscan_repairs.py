"""CPU regression: real HS inputs, phase/identity failures, and ns plots.

Only temporary copies are mutated. No solver, no historical output writes.
Also runnable with python -O: rejection checks use unittest, not assert.
"""
import argparse
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
import json
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from hs_capsule_identity import read_manifest, sha256, verify_file
import run_hs_acceptance_v0_2 as acceptance
from plot_hs4_bscan_v0_2 import diagnostic_arrays

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'
T2 = ROOT / 'artifacts/research_checks/2026-10-02_hs4t2d_transect'
WORK = None


def snapshot(folder):
    return {p.name: sha256(p) for p in folder.iterdir() if p.is_file()}


class Repairs(unittest.TestCase):
    def setUp(self):
        self.folder = Path(tempfile.mkdtemp(dir=WORK))
        self.base = self.folder / 'input'
        self.base.mkdir()
        for tag in acceptance.H5S:
            shutil.copy2(CAPSULE / f'{tag}.h5', self.base)
        self.records = [{'file': p.name, 'sha256': sha256(p), 'bytes': p.stat().st_size}
                        for p in sorted(self.base.iterdir())]
        self.manifest(self.records)

    def manifest(self, records):
        (self.base / 'manifest.json').write_text(json.dumps(records), encoding='utf-8')

    def invoke(self, expected, base=None):
        dest = self.folder / 'output'
        argv = ['acceptance', '--dir', str(base or self.base), '--out', str(dest),
                '--fig', str(dest / 'acceptance.png')]
        log = StringIO()
        with patch.object(sys, 'argv', argv), redirect_stdout(log), redirect_stderr(log):
            if expected:
                with self.assertRaises(SystemExit) as raised:
                    acceptance.main()
                self.assertEqual(raised.exception.code, expected)
            else:
                acceptance.main()
        (self.folder / 'log.txt').write_text(log.getvalue(), encoding='utf-8')
        return dest

    def test_actual_current_capsule_pass(self):
        dest = self.invoke(0, CAPSULE)
        metrics = json.loads((dest / 'hs_acceptance_metrics_v0_3.json').read_text(encoding='utf-8'))
        self.assertEqual(metrics['verdict'], 'PASS')
        self.assertAlmostEqual(metrics['ground_env_peak_ns'], 99.8003992015968, places=10)
        self.assertEqual(metrics['domain_control_threshold_over_direct_peak'], 1e-5)

    def test_historical_dict_pass(self):
        self.manifest({r['file']: r['sha256'] for r in self.records})
        self.invoke(0)

    def test_list_pass(self):
        self.invoke(0)

    def mutate_receiver(self, function):
        p = self.base / f'{acceptance.H5S[2]}.h5'
        with h5py.File(p, 'r+') as h:
            h['rxs/rx1/Ex'][:] = function(h['rxs/rx1/Ex'][:])
        for r in self.records:
            if r['file'] == p.name:
                r['sha256'] = sha256(p)
        self.manifest(self.records)

    def test_polarity_reversal_rejected_with_matching_identity(self):
        self.mutate_receiver(lambda v: -v)
        dest = self.invoke(1)
        metrics = json.loads((dest / 'hs_acceptance_metrics_v0_3.json').read_text(encoding='utf-8'))
        self.assertEqual(metrics['verdict'], 'FAIL')
        self.assertGreater(metrics['hs3_vs_hs1_rel_diff_over_direct_peak']['direct_ns'], 1.9)
        self.assertLess(metrics['hs3_vs_hs1_amplitude_only_over_direct_DIAGNOSTIC_ONLY']['direct_ns'], 1e-5)

    def test_amplitude_change_rejected(self):
        self.mutate_receiver(lambda v: v * 2)
        self.invoke(1)

    def test_temporal_shift_rejected(self):
        self.mutate_receiver(lambda v: np.roll(v, 8))
        self.invoke(1)

    def test_hash_mismatch_no_output(self):
        self.records[0]['sha256'] = '0' * 64
        self.manifest(self.records)
        self.invoke(2)
        self.assertFalse((self.folder / 'output').exists())

    def test_byte_count_mismatch(self):
        self.records[0]['bytes'] += 1
        self.manifest(self.records)
        self.invoke(2)

    def test_missing_record(self):
        self.manifest(self.records[1:])
        self.invoke(2)

    def test_missing_file(self):
        (self.base / self.records[0]['file']).unlink()
        self.invoke(2)

    def test_duplicate_list_record(self):
        self.manifest(self.records + [self.records[0]])
        self.invoke(2)

    def test_duplicate_dict_key(self):
        r = self.records[0]
        (self.base / 'manifest.json').write_text('{"a":"' + r['sha256'] + '","a":"' + r['sha256'] + '"}')
        self.invoke(2)

    def test_malformed_schema(self):
        for value in [42, None, ['bad'], [{'file': 'a', 'sha256': 'abc'}]]:
            with self.subTest(value=value):
                self.manifest(value)
                self.invoke(2)

    def test_unsafe_paths(self):
        for name in ['../x', '/x', 'D:/x', 'a\\b', 'a/../b', 'a//b']:
            with self.subTest(name=name):
                self.manifest([{'file': name, 'sha256': '0' * 64}])
                self.invoke(2)

    def test_existing_output_unchanged(self):
        dest = self.folder / 'output'
        dest.mkdir()
        (dest / 'old.txt').write_text('history')
        before = snapshot(dest)
        self.invoke(1)
        self.assertEqual(snapshot(dest), before)

    def test_existing_figure_unchanged(self):
        fig = self.folder / 'old.png'
        fig.write_bytes(b'historical figure')
        with patch.object(sys, 'argv', ['acceptance', '--dir', str(self.base), '--out',
                                       str(self.folder / 'output'), '--fig', str(fig)]):
            with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                acceptance.main()
        self.assertEqual(fig.read_bytes(), b'historical figure')
        self.assertFalse((self.folder / 'output').exists())

    def test_output_inside_input_rejected(self):
        before = snapshot(self.base)
        with patch.object(sys, 'argv', ['acceptance', '--dir', str(self.base), '--out',
                                       str(self.base / 'new'), '--fig', str(self.folder / 'fig.png')]):
            with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                acceptance.main()
        self.assertEqual(snapshot(self.base), before)
        self.assertFalse((self.base / 'new').exists())

    def test_figure_inside_input_rejected(self):
        before = snapshot(self.base)
        with patch.object(sys, 'argv', ['acceptance', '--dir', str(self.base), '--out',
                                       str(self.folder / 'output'), '--fig', str(self.base / 'new.png')]):
            with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                acceptance.main()
        self.assertEqual(snapshot(self.base), before)

    def test_changed_during_loading(self):
        original = acceptance.trace_time_response
        def concurrent_change(path):
            tr = original(path)
            if path.name == f'{acceptance.H5S[2]}.h5':
                with h5py.File(path, 'r+') as h:
                    h['rxs/rx1/Ex'][0] = 1
            return tr
        with patch.object(acceptance, 'trace_time_response', concurrent_change):
            self.invoke(2)

    def test_reconstructed_axes_and_envelopes(self):
        t = np.arange(301) * 1e-9
        good = {k: SimpleNamespace(time=t.copy(), complex_envelope=np.ones(301, dtype=complex))
                for k in acceptance.H5S}
        acceptance.checked_envelopes(good)
        for field, replacement in [('time', t + 1e-12), ('complex_envelope', np.ones(300)),
                                   ('complex_envelope', np.full(301, complex(1, np.nan)))]:
            with self.subTest(field=field):
                bad = dict(good)
                bad[acceptance.H5S[2]] = SimpleNamespace(time=t.copy(), complex_envelope=np.ones(301, dtype=complex))
                setattr(bad[acceptance.H5S[2]], field, replacement)
                with self.assertRaises(ValueError):
                    acceptance.checked_envelopes(bad)

    def test_plot_ns_and_signed_identity(self):
        t = np.arange(401, dtype=float)
        data = np.tile(np.sin(t[:, None] / 10), (1, 3))
        tm, raw, centered, residual, gain, agc, samples, floor = diagnostic_arrays(t, data)
        np.testing.assert_array_equal(tm, np.arange(251))
        np.testing.assert_array_equal(raw, data[:251])
        np.testing.assert_allclose(residual, 0, atol=1e-14)
        np.testing.assert_array_equal(agc, residual * gain)
        self.assertEqual(samples, 20)

    def test_plot_invalid_axes(self):
        data = np.zeros((301, 3))
        for t, v in [(np.arange(301)[::-1], data), (np.arange(301), data.T),
                     (np.arange(301), data.astype(complex)), (np.arange(301) ** 2, data),
                     (np.arange(301), data + np.nan)]:
            with self.subTest():
                with self.assertRaises(ValueError):
                    diagnostic_arrays(t, v)


def main():
    global WORK
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', help='NEW directory; defaults to ignored unique local_checks')
    args = ap.parse_args()
    WORK = Path(args.out) if args.out else Path(tempfile.mkdtemp(prefix='hs4_repairs_', dir=ROOT / 'artifacts/local_checks'))
    if args.out:
        WORK.mkdir(parents=True)
    before = {str(p): snapshot(p) for p in (CAPSULE, T2)}
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Repairs))
    unchanged = before == {str(p): snapshot(p) for p in (CAPSULE, T2)}
    summary = {'tests': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
               'input_capsules_unchanged': unchanged, 'optimized_python': not __debug__,
               'script_sha256': sha256(Path(__file__)), 'solver_executed': False}
    (WORK / 'checks.json').write_text(json.dumps(summary, indent=1), encoding='utf-8')
    print(json.dumps(summary, indent=1))
    if not result.wasSuccessful() or not unchanged:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
