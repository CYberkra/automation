"""Acceptance regressions: reject corrupt archives, inputs and window claims.

Only development JSON, temporary byte fixtures and one archived control H5.
No solver, training, test-family or original field-data access.
"""

import argparse
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import sfcw_carrierfix_acceptance as audit

CHECKS = ROOT / 'artifacts/research_checks'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


class AcceptanceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t3 = read(CHECKS / '2026-10-02_t3_ladder_carrierfix_v0_2_r1.json')['records']
        cls.old_t3 = read(CHECKS / '2026-09-28_t3_damage_ladder_r1.json')
        cls.keys = ('D', 'a', 'A', 'H', 'rho', 'Nb_ratio',
                    'arrival_drift_ns_median', 'arrival_drift_ns_max_abs')
        cls.b2 = read(CHECKS / '2026-10-02_b2_reward_carrierfix_v0_2_r1.json')['results']
        for r in cls.b2:
            r['ranking_changed'] = (r['representations']['legacy95']['ranking']
                                     != r['representations']['official20']['ranking'])
        cls.old_b2 = read(CHECKS / '2026-10-01_reward_protocol_b2_pilot_r1.json')
        cls.s1 = read(CHECKS / '2026-10-02_s1s3_refwindow_carrierfix_v0_2_r1.json')
        cls.contract = read(ROOT / 'configs/research/s1s3_reference_window_contract_v0.1.json')

    def test_t3_archive_matches(self):
        self.assertEqual(audit.accept_t3(self.t3, self.old_t3, self.keys)['metrics_compared'], 264)

    def test_t3_altered_legacy_metric_rejected(self):
        rows = copy.deepcopy(self.t3)
        next(r for r in rows if r['representation'] == 'legacy95')['D'] += 0.001
        with self.assertRaisesRegex(ValueError, 'legacy/archive mismatch'):
            audit.accept_t3(rows, self.old_t3, self.keys)

    def test_t3_missing_or_duplicate_rows_rejected(self):
        for rows in (self.t3[:-1], self.t3 + [self.t3[0]]):
            with self.assertRaises(ValueError):
                audit.accept_t3(rows, self.old_t3, self.keys)

    def test_b2_archive_matches_and_rank_changes_recorded(self):
        result = audit.accept_b2(self.b2, self.old_b2)
        self.assertEqual(result['ranking_changed_families'], ['S2X', 'S2TZX'])

    def test_b2_changed_selection_or_score_rejected(self):
        for field, value in [('selection', 'wrong_config'), ('selection_R_bg_db', -99)]:
            rows = copy.deepcopy(self.b2)
            rows[0]['representations']['legacy95'][field] = value
            with self.assertRaisesRegex(ValueError, 'legacy/archive mismatch'):
                audit.accept_b2(rows, self.old_b2)

    def test_b2_missing_metric_or_false_rank_flag_rejected(self):
        rows = copy.deepcopy(self.b2)
        del rows[0]['representations']['legacy95']['rows'][-1]['D_e']
        with self.assertRaisesRegex(ValueError, 'field coverage'):
            audit.accept_b2(rows, self.old_b2)
        rows = copy.deepcopy(self.b2)
        rows[1]['ranking_changed'] = False
        with self.assertRaisesRegex(ValueError, 'ranking flag'):
            audit.accept_b2(rows, self.old_b2)

    def test_s1s3_frozen_registration_and_windows(self):
        for fam, result in self.s1['families'].items():
            self.assertTrue(audit.accept_s1s3(result['rows'], self.contract['families'][fam],
                            self.contract['constants'])['peak_window_and_nc_pass'])

    def test_s1s3_changed_legacy_peak_rejected(self):
        rows = copy.deepcopy(self.s1['families']['S1X']['rows'])
        rows[0]['t_measured_legacy_ns'] += .1
        with self.assertRaisesRegex(ValueError, 'legacy/frozen'):
            audit.accept_s1s3(rows, self.contract['families']['S1X'], self.contract['constants'])

    def test_s1s3_outside_event_window_rejected(self):
        rows = copy.deepcopy(self.s1['families']['S1X']['rows'])
        rows[0]['t_measured_official_ns'] += 20
        with self.assertRaisesRegex(ValueError, 'outside Fermat window'):
            audit.accept_s1s3(rows, self.contract['families']['S1X'], self.contract['constants'])

    def test_s1s3_nc_overlap_rejected(self):
        constants = dict(self.contract['constants'], nc_window_ns=[180, 400])
        with self.assertRaisesRegex(ValueError, 'overlaps NC'):
            audit.accept_s1s3(self.s1['families']['S1X']['rows'],
                            self.contract['families']['S1X'], constants)

    def test_input_identity_and_mutation(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / 'dev_run_att2'
            source.mkdir()
            path = source / 'dev_run.h5'
            path.write_bytes(b'byte fixture; not scientific H5')
            expected = {'run_id': 'dev_run', 'source_dir': source.name,
                        'h5_sha256': audit.sha256(path)}
            self.assertEqual(audit.verify_input(path, expected)['sha256'], expected['h5_sha256'])
            with self.assertRaisesRegex(ValueError, 'directory mismatch'):
                audit.verify_input(path, dict(expected, source_dir='dev_run_att3'))
            with self.assertRaisesRegex(ValueError, 'run_id/path mismatch'):
                audit.verify_input(path, dict(expected, run_id='wrong_run'))
            path.write_bytes(b'changed bytes')
            with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
                audit.verify_input(path, expected)

    def test_manifest_missing_duplicate_and_scope(self):
        with self.assertRaisesRegex(ValueError, 'required'):
            audit.input_manifest(None, ['dev_run'])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'manifest.json'
            row = {'run_id': 'dev_run', 'source_dir': 'dev_run_att2', 'h5_sha256': 'a' * 64}
            for rows, ids in [([row, row], ['dev_run']), ([row], ['wrong_run']),
                              ([dict(row, source_dir='../outside')], ['dev_run']),
                              ([dict(row, h5_sha256='')], ['dev_run'])]:
                path.write_text(json.dumps({'inputs': rows}), encoding='utf-8')
                with self.assertRaises(ValueError):
                    audit.input_manifest(path, ids)
            path.write_text(json.dumps({'inputs': [row]}), encoding='utf-8')
            self.assertEqual(audit.input_manifest(path, ['dev_run'])['dev_run'], row)

    def test_official_code_drift_rejected_before_h5_read(self):
        import sfcw_official_loader_v0_2 as loader
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'changed_processing.py'
            path.write_text('# changed official source', encoding='utf-8')
            with patch.object(loader._processing, '__file__', str(path)):
                with self.assertRaisesRegex(ValueError, 'differs from reviewed version'):
                    loader.trace_time_response('not_a_real_input.h5')

    def test_writer_rejects_mismatch_and_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            p1, p2 = Path(temp) / 'r1.json', Path(temp) / 'r2.json'
            with self.assertRaisesRegex(ValueError, 'r1/r2 mismatch'):
                audit.write_new_pair(p1, p2, {'x': 1}, {'x': 2})
            self.assertFalse(p1.exists() or p2.exists())
            audit.write_new_pair(p1, p2, {'x': 1}, {'x': 1})
            before = p1.read_bytes()
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                audit.write_new_pair(p1, p2, {'x': 2}, {'x': 2})
            self.assertEqual(p1.read_bytes(), before)

    def test_t3_main_stops_before_writing_bad_legacy(self):
        import study_t3_ladder_carrierfix_v0_2 as module
        rows = copy.deepcopy(self.t3)
        rows[0]['D'] += .01  # legacy row deliberately corrupt
        with tempfile.TemporaryDirectory() as temp:
            p1, p2 = Path(temp) / 'r1.json', Path(temp) / 'r2.json'
            with patch.object(sys, 'argv', ['check', '--input-manifest', 'unused.json']), \
                 patch.object(module, 'input_manifest', return_value={}), \
                 patch.object(module, 'build_records', return_value=rows), \
                 patch.object(module, 'OUT_R1', p1), patch.object(module, 'OUT_R2', p2):
                with self.assertRaisesRegex(ValueError, 'legacy/archive mismatch'):
                    module.main()
            self.assertFalse(p1.exists() or p2.exists())

    def test_b2_main_stops_before_writing_bad_legacy(self):
        import study_b2_reward_carrierfix_v0_2 as module
        rows = copy.deepcopy(self.b2)
        rows[0]['representations']['legacy95']['selection'] = 'wrong_config'
        with tempfile.TemporaryDirectory() as temp:
            p1, p2 = Path(temp) / 'r1.json', Path(temp) / 'r2.json'
            with patch.object(sys, 'argv', ['check', '--input-manifest', 'unused.json']), \
                 patch.object(module, 'input_manifest', return_value={}), \
                 patch.object(module, 'build', return_value=rows), \
                 patch.object(module, 'OUT_R1', p1), patch.object(module, 'OUT_R2', p2):
                with self.assertRaisesRegex(ValueError, 'legacy/archive mismatch'):
                    module.main()
            self.assertFalse(p1.exists() or p2.exists())

    def test_s1s3_main_stops_before_writing_bad_legacy(self):
        import numpy as np
        import study_s1s3_refwindow_carrierfix_v0_2 as module
        with tempfile.TemporaryDirectory() as temp:
            p1, p2 = Path(temp) / 'r1.json', Path(temp) / 'r2.json'
            # Fake zero waveforms intentionally cannot reproduce frozen peaks.
            # Exercise the real builder/acceptance path without reading raw H5.
            with patch.object(module, 'TIERS', ['S1X']), \
                 patch.object(module, 'fermat_times', return_value=np.asarray(
                     self.contract['families']['S1X']['t_fermat_er18_ns'])), \
                 patch.object(module, 'verify_input', return_value={}), \
                 patch.object(module.L, 'trace_time_response', return_value=None), \
                 patch.object(module.L, 'time_ns', return_value=np.arange(4008) / 1.2024), \
                 patch.object(module.L, 'signed_official', return_value=np.zeros(4008)), \
                 patch.object(module.L, 'signed_legacy', return_value=np.zeros(4008)), \
                 patch.object(module, 'OUT_R1', p1), patch.object(module, 'OUT_R2', p2):
                with self.assertRaisesRegex(ValueError, 'legacy/frozen peak mismatch'):
                    module.main()
            self.assertFalse(p1.exists() or p2.exists())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'artifacts/local_checks/carrierfix_acceptance_checks.json')
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(AcceptanceChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    from review_sfcw_carrierfix_20261002 import audit_actual_loader, audit_archives
    doc = {'schema': 'carrierfix-acceptance-checks/1', 'checks_passed': result.testsRun,
           'full_raw_rerun': 'not performed: original full dev H5 sets unavailable',
           'archives': audit_archives(), 'actual_loader': audit_actual_loader(),
           'code': [audit.file_record(ROOT / 'scripts' / name) for name in (
               'sfcw_carrierfix_acceptance.py', 'sfcw_official_loader_v0_2.py',
               'study_t3_ladder_carrierfix_v0_2.py', 'study_b2_reward_carrierfix_v0_2.py',
               'study_s1s3_refwindow_carrierfix_v0_2.py', 'check_sfcw_carrierfix_acceptance.py')]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(doc, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print('Acceptance tests, archive crosschecks and actual-loader checks passed:', args.output)


if __name__ == '__main__':
    main()
