"""Audit official SFCW output; compare the M00 transverse dipole to continuum theory.

This does not synthesise SFCW, invoke FDTD, fit phase/amplitude, or certify convergence.
Theory: Cornell ECE303 Lecture 28, pp.10-11, exp(+j omega t) convention.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c, epsilon_0, mu_0


def transverse_dipole(frequency, distance, length):
    """Ex/I for x current and receiver displaced along +y, units V/(m A)."""
    k = 2 * np.pi * np.asarray(frequency) / c
    kr = k * distance
    eta = np.sqrt(mu_0 / epsilon_0)
    # theta points opposite to the source axis in its equatorial plane.
    return (-1j * eta * k * length / (4 * np.pi * distance)
            * (1 + 1 / (1j * kr) + 1 / (1j * kr)**2) * np.exp(-1j * kr))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--raw-audit', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    raw = json.loads(args.raw_audit.read_text(encoding='utf-8'))
    if not raw['passed']:
        raise ValueError('Raw metadata audit must pass first')
    with h5py.File(args.input, 'r') as f:
        freq, response = f['frequency'][:], f['response'][:]
        src, rec = f['source_spectrum'][:], f['receiver_spectrum'][:]
        valid = f['source_valid'][:].astype(bool)
        t = f['time_response/time'][:]
        env = f['time_response/complex_envelope'][:]
        a = dict(f.attrs)
        checks = {
            'format': a['Format'] == 'gprMax SFCW toolbox',
            'method_direct': a['Method'] == 'direct',
            '501_frequency_points': freq.shape == (501,),
            'exact_requested_grid': bool(np.array_equal(freq, np.linspace(20e6, 170e6, 501))),
            'finite_response': bool(np.all(np.isfinite(response))),
            'nonzero_response': bool(np.any(response != 0)),
            'valid_sources': bool(valid.shape == (501,) and np.all(valid)),
            'source_normalisation': bool(np.allclose(response * src, rec, rtol=1e-12, atol=0)),
            'IQ_consistent': bool(np.array_equal(f['I'][:] + 1j * f['Q'][:], response)),
            'source_identity': a['SourcePath'] == raw['source_path'],
            'receiver_identity': a['ReceiverPath'] == raw['receiver_path'] + '/Ex',
            'source_units_A': a['SourceUnits'] == 'A',
            'receiver_units_V_m': a['ReceiverUnits'] == 'V/m',
            'source_time_offset': bool(np.isclose(a['SourceTimeSampleOffset'], raw['root_attrs']['dt']/2, rtol=1e-12, atol=0)),
            'receiver_time_offset': float(a['ReceiverTimeSampleOffset']) == 0,
            'rectangular_window': f['time_response'].attrs['Window'] == 'rectangular',
            'no_zero_padding': int(f['time_response'].attrs['ZeroPadFactor']) == 1,
            'no_time_shift': float(f['time_response'].attrs['TimeShift']) == 0,
            'no_tail_taper': float(a['TailTaperFraction']) == 0,
            'periodic_time_grid': bool(len(t)==501 and np.allclose(t, np.arange(501)/(501*300e3), rtol=1e-12, atol=1e-18)),
        }
    distance = float(np.linalg.norm(np.array(raw['receiver_attrs']['Position']) - raw['source_attrs']['Position']))
    length = float(a['SourceSpatialScale'])
    theory = transverse_dipole(freq, distance, length)
    relative = np.abs(response-theory)/np.abs(theory)
    amplitude_db = 20*np.log10(np.abs(response)/np.abs(theory))
    phase_deg = np.angle(response/theory, deg=True)
    with (args.output_dir/'frequency_comparison.csv').open('w', newline='', encoding='utf-8') as out:
        writer=csv.writer(out)
        writer.writerow(['frequency_Hz','official_real_V_per_m_A','official_imag_V_per_m_A','theory_real_V_per_m_A','theory_imag_V_per_m_A','relative_complex_error','amplitude_difference_dB','wrapped_phase_difference_deg'])
        writer.writerows(zip(freq,response.real,response.imag,theory.real,theory.imag,relative,amplitude_db,phase_deg))
    result = {
        'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
        'raw_audit_sha256': hashlib.sha256(args.raw_audit.read_bytes()).hexdigest(),
        'checks': checks, 'passed': all(checks.values()),
        'theory_url': 'https://courses.cit.cornell.edu/ece303/Lectures/lecture28.pdf',
        'theory_pages_one_based': [10,11],
        'theory_assumptions': 'Infinite homogeneous vacuum; ideal point current element; transverse geometry; exp(+j omega t). No fitting. Finite-record FDTD need not equal infinite-time harmonic theory.',
        'distance_m': distance, 'dipole_length_m': length,
        'response_units': 'V/(m A); not hardware S21',
        'max_relative_complex_error': float(np.max(relative)),
        'median_relative_complex_error': float(np.median(relative)),
        'max_absolute_amplitude_difference_dB': float(np.max(np.abs(amplitude_db))),
        'max_absolute_wrapped_phase_difference_deg': float(np.max(np.abs(phase_deg))),
        'receiver_tail_relative_dB': float(a['ReceiverTailRelativeDB']),
        'period_s': 1/300e3, 'physical_FDTD_record_s': raw['root_attrs']['dt']*raw['root_attrs']['Iterations'],
        'physical_accuracy_threshold': None, 'convergence_certified': False,
        'solver_invoked': False,
    }
    (args.output_dir/'audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2,figsize=(11,7),layout='constrained')
    mhz = freq/1e6
    axes[0,0].plot(mhz, np.abs(response), label='Official direct; finite record')
    axes[0,0].plot(mhz, np.abs(theory), '--', label='Infinite-medium harmonic theory')
    axes[0,0].set(ylabel='|Ex/I| [V/(m A)]',xlabel='Frequency [MHz]')
    axes[0,0].legend(fontsize=8)
    axes[0,1].plot(mhz, amplitude_db)
    axes[0,1].set(ylabel='Amplitude difference [dB]',xlabel='Frequency [MHz]')
    axes[1,0].plot(mhz, phase_deg)
    axes[1,0].set(ylabel='Wrapped phase difference [degree]',xlabel='Frequency [MHz]')
    axes[1,1].plot(t*1e6,np.abs(env))
    axes[1,1].set(xlabel='Periodic reconstructed time [us]',ylabel='|Complex envelope| [V/(m A)]')
    for ax in axes.flat: ax.grid(alpha=.25)
    fig.suptitle('M00: 1.3 m transverse baseline, 5 cm dipole; no fitted correction')
    fig.savefig(args.output_dir/'sfcw_diagnostics.png',dpi=160)
    plt.close(fig)
    print(json.dumps(result))
    if not result['passed']: raise SystemExit(1)


if __name__ == '__main__':
    main()
