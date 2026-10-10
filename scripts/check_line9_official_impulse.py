"""CPU checks of source-only edits, built-in impulse and native-source rejection."""
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from gprMax.waveforms import Waveform
from gprMax.toolboxes.SFCW import processing as sf
from prepare_line9_official_impulse import replace_source
from line9_v401_version_controls import audit_source


BASE='''#domain_mode: TM
#domain: 210 42.5 inf
#waveform: ricker 40 100000000 pulse
#hertzian_dipole: z 169.35 38.95 0.0125 pulse
#rx: 170.65 39.025 0.0125 rx1 Ez
'''


class ImpulseChecks(unittest.TestCase):
    def test_only_source_commands_change(self):
        for ending in ['\n','\r\n']:
            old=BASE.replace('\n',ending);new=replace_source(old)
            before=old.splitlines(keepends=True);after=new.splitlines(keepends=True)
            self.assertEqual([i for i,(a,b) in enumerate(zip(before,after)) if a!=b],[2,3])
            self.assertEqual(after[2],'#waveform: impulse 1 1 impulse'+ending)
            self.assertEqual(after[3],'#hertzian_dipole: z 169.35 38.95 0.0125 impulse'+ending)

    def test_ambiguous_or_delayed_source_rejected(self):
        for bad in [BASE+'#waveform: ricker 40 1 other\n',
                    BASE.replace('0.0125 pulse','0.0125 pulse 1e-9'),
                    BASE.replace('0.0125 pulse','0.0125 missing'),
                    BASE.replace('#hertzian_dipole: z','#hertzian_dipole: x')]:
            with self.assertRaises(ValueError):replace_source(bad)

    def test_official_builtin_is_single_sample_at_half_step(self):
        dt=5.896635841874211e-11
        wave=Waveform();wave.type='impulse';wave.amp=1.;wave.freq=1.
        samples=np.array([wave.calculate_value((i+.5)*dt,dt) for i in range(256)])
        np.testing.assert_array_equal(samples,np.r_[1.,np.zeros(255)])
        frequencies=20e6+np.arange(501)*300000.
        spectrum=sf.engineering_dft(samples,dt,frequencies,time_offset=.5*dt)
        np.testing.assert_allclose(spectrum,dt*np.exp(-2j*np.pi*frequencies*.5*dt),rtol=1e-13,atol=1e-24)

    def test_native_impulse_shape_and_metadata_checked(self):
        dt=1e-10;spec={'source_type':'impulse','source_frequency_Hz':1.,'source_amplitude_A':1.}
        with tempfile.TemporaryDirectory(prefix='line9-impulse-check-') as temp:
            with h5py.File(Path(temp)/'synthetic.h5','w') as h:
                g=h.create_group('excitation');g.attrs.update(WaveformType='impulse',WaveformFrequency=1.,
                    WaveformAmplitude=1.,SourceStartTime=0.,TimeSampleOffset=.5*dt)
                data=g.create_dataset('samples',data=np.r_[1.,np.zeros(15)])
                audit_source(g,spec,dt)
                data[1]=1.
                with self.assertRaises(AssertionError):audit_source(g,spec,dt)
                data[1]=0.;g.attrs['TimeSampleOffset']=0.
                with self.assertRaises(AssertionError):audit_source(g,spec,dt)
                g.attrs['TimeSampleOffset']=.5*dt;g.attrs['WaveformType']='ricker'
                with self.assertRaises(AssertionError):audit_source(g,spec,dt)

    def test_legacy_ricker_still_passes_source_audit(self):
        with tempfile.TemporaryDirectory(prefix='line9-ricker-check-') as temp:
            with h5py.File(Path(temp)/'synthetic.h5','w') as h:
                g=h.create_group('excitation');g.attrs.update(WaveformType='ricker',
                    WaveformFrequency=100e6,WaveformAmplitude=40.)
                audit_source(g,{},1e-10)


if __name__=='__main__':
    unittest.main(verbosity=2)
