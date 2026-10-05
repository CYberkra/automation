"""Independent continuum dipole sanity checks, not an FDTD or device validation.

MIT 6.013 notes section 10.2: exact Hertzian E includes 1/r, 1/r^2,
and 1/r^3 terms; constant-current radiated power scales as frequency^2.
The source is an ideal fixed-length current element in homogeneous air.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def axial_to_broadside(kr):
    kr = np.asarray(kr, dtype=float)
    return 2*np.sqrt(1+kr**2)/np.sqrt((1-kr**2)**2+kr**2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('fresh output directory required')
    # Static near field has a 2:1 axial/broadside amplitude ratio. At large
    # kr, the axial residual falls as 1/r^2 while broadside radiation is 1/r.
    assert abs(float(axial_to_broadside(1e-4))-2) < 1e-6
    assert abs(float(axial_to_broadside(1e4))/(2/1e4)-1) < 1e-6
    c = 299792458.
    distance = 1.3
    freq = np.array([20., 50., 95., 170.])*1e6
    kr = 2*np.pi*freq*distance/c
    rows = [{'frequency_MHz': float(f/1e6), 'r_over_wavelength': float(q/(2*np.pi)),
             'kr': float(q), 'axial_over_broadside_dB': float(20*np.log10(v))}
            for f, q, v in zip(freq, kr, axial_to_broadside(kr))]
    result = {'status': 'PASS_ANALYTIC_SANITY_CHECK_NOT_FDTD_VALIDATION',
              'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'baseline_m': distance, 'frequency_rows': rows,
              'fixed_current_radiated_power_170_over20_dB': float(20*np.log10(170/20)),
              'equal_radiated_power_current_170_over20_dB': float(20*np.log10(20/170)),
              'source': 'https://ocw.mit.edu/courses/6-013-electromagnetics-and-applications-spring-2009/d3be4ea78b036a6362230fb41780cf54_MIT6_013S09_notes.pdf',
              'source_section': '10.2.1 exact electric field and 10.2.13 radiated power',
              'asymptotic_checks': {'static_ratio_2': True, 'far_ratio_2_over_kr': True},
              'scope': 'Homogeneous air, ideal 3D infinitesimal dipole, same co-polar field component; no ground/antenna/port.',
              'limits': ['No source strength or spectral weighting correction inferred for the actual UAV.',
                         'These 3D power relations are not 2D line-source power calibration.',
                         'Axial far-field radiation null does not mean zero near-field coupling.',
                         'No confirmation that layout is the main B-scan cause.'],
              'solver_called': False, 'field_data_read': False, 'training_called': False}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
