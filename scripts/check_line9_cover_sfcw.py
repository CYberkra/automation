"""CPU regression: physical sample origins, trace axes and inverse convention."""
import unittest

import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from analyze_line9_cover_planes import field_response
from analyze_line9_v401_version_controls import FREQ, inverse


class CoverSFCWChecks(unittest.TestCase):
    def test_known_delay_with_half_step_source(self):
        dt=1e-10
        samples=np.zeros(256);samples[0]=1.
        source=sf.SampledSignal('synthetic/impulse',samples,dt,.5*dt)
        fields=np.zeros((2,3,256))
        amplitudes=np.array([[1.,-2.,.25],[3.,-.5,2.]])
        fields[:,:,16]=amplitudes
        result=field_response(source,fields,dt)
        # Analytic delta pair, independent of any source/receiver DFT code.
        expected=np.exp(-2j*np.pi*FREQ*(16-.5)*dt)[:,None]*amplitudes.reshape(1,-1)/.025
        np.testing.assert_allclose(result.response,expected,rtol=2e-12,atol=1e-12)
        wrong=np.exp(-2j*np.pi*FREQ*16*dt)[:,None]*amplitudes.reshape(1,-1)/.025
        self.assertGreater(np.linalg.norm(result.response-wrong)/np.linalg.norm(expected),.01)

    def test_official_inverse_preserves_existing_complex_convention(self):
        samples=np.zeros(128);samples[0]=1.
        source=sf.SampledSignal('synthetic/impulse',samples,1e-10,.5e-10)
        fields=np.zeros((1,2,128));fields[0,:,30]=[1.,-2.]
        result=field_response(source,fields,1e-10)
        for window in ['hann','blackman']:
            official=sf.reconstruct_time_response(result,window=window,zero_pad_factor=8)
            legacy,t=inverse(result.response,sf.spectral_window(window,501))
            np.testing.assert_array_equal(t,official.time)
            np.testing.assert_allclose(legacy,official.complex_bandpass,rtol=2e-13,atol=1e-13)
            np.testing.assert_array_equal(official.real_bandpass,2*official.complex_bandpass.real)
            # Explicit finite sinusoidal sum also checks first-frequency carrier.
            pick=np.array([0,30,100,1000,4007])
            expected=np.exp(2j*np.pi*t[pick,None]*FREQ)@(result.response*official.weights[:,None])/501
            np.testing.assert_allclose(official.complex_bandpass[pick],expected,rtol=1e-9,atol=1e-11)

    def test_inconsistent_sampling_is_rejected(self):
        source=sf.SampledSignal('synthetic/impulse',np.r_[1.,np.zeros(63)],1e-10,.5e-10)
        with self.assertRaises(ValueError):
            field_response(source,np.ones((1,2,64)),2e-10)


if __name__=='__main__':
    unittest.main(verbosity=2)
