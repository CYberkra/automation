"""Regression counterexamples for HS4 frequency masks and matched shape metrics."""
import json
import unittest
import numpy as np
from hs4_analysis_metrics import air_propagation_mask, matched_profile_metrics


class MetricChecks(unittest.TestCase):
    def setUp(self):
        self.x = np.arange(128) / 128 * 8
        self.low = np.cos(2*np.pi*self.x/8)
        self.high = .1*np.sin(2*np.pi*self.x/.5)
        self.truth = 9+self.low+self.high

    def metrics(self, ridge):
        return matched_profile_metrics(ridge, self.truth, 8/128, 1.6)

    def test_perfect_details(self):
        m = self.metrics(self.truth)
        for band in ('lowpass_matched', 'highpass_matched'):
            self.assertAlmostEqual(m[band]['corr'], 1)
            self.assertEqual(m[band]['relative_L2_error'], 0)

    def test_erased_details(self):
        m = self.metrics(9+self.low)['highpass_matched']
        self.assertLess(m['amplitude_norm_ratio'], 1e-12)
        self.assertAlmostEqual(m['relative_L2_error'], 1)

    def test_reversed_details(self):
        m = self.metrics(9+self.low-self.high)['highpass_matched']
        self.assertAlmostEqual(m['corr'], -1)
        self.assertAlmostEqual(m['signed_projection_gain'], -1)
        self.assertAlmostEqual(m['relative_L2_error'], 2)

    def test_overamplified_details(self):
        m = self.metrics(9+self.low+5*self.high)['highpass_matched']
        self.assertAlmostEqual(m['corr'], 1)
        self.assertAlmostEqual(m['amplitude_norm_ratio'], 5)
        self.assertAlmostEqual(m['relative_L2_error'], 4)

    def test_bias_preserved(self):
        m = self.metrics(self.truth+.2)
        self.assertAlmostEqual(m['mean_bias_m'], .2)
        self.assertAlmostEqual(m['rmse_m'], .2)

    def test_constant_missing_direction(self):
        m = self.metrics(np.full_like(self.truth, 9))
        self.assertIsNone(m['highpass_matched']['corr'])

    def test_frequency_not_carrier(self):
        # Fixed 95 MHz gives the wrong classification in both directions.
        m = air_propagation_mask([70e6, 120e6], [1.7, 2.3])
        np.testing.assert_array_equal(m, [[False, False], [True, True]])

    def test_mask_symmetry_and_units(self):
        m = air_propagation_mask([100e6], [-2.2, 0, 2.2])
        np.testing.assert_array_equal(m, [[False, True, False]])

    def test_invalid_frequency(self):
        for f in ([0], [-1], [np.nan]):
            with self.assertRaises(ValueError):
                air_propagation_mask(f, [0])

    def test_nonfinite_profile(self):
        with self.assertRaises(ValueError):
            self.metrics(np.full_like(self.truth, np.nan))

    def test_axis_mismatch(self):
        with self.assertRaises(ValueError):
            self.metrics(self.truth[:-1])

    def test_invalid_scale(self):
        with self.assertRaises(ValueError):
            matched_profile_metrics(self.truth, self.truth, .1, 0)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MetricChecks)
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    print(json.dumps({'tests': result.testsRun, 'passed': result.wasSuccessful(),
                      'scope': 'Array counterexamples, not FDTD or physical validation'}))
    raise SystemExit(0 if result.wasSuccessful() else 1)
