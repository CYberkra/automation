"""Independent Yee mask identity and saved official source-history audits.

No solver runs. Checks signs, electric/magnetic moments and half-step indexing.
This does not certify a finite-aperture plane or terrain/return propagation.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def mask_identity():
    r = np.random.default_rng(20261010)
    nx, ny, p = 19, 23, 9
    dx, dy, dz = .05, .05, .05
    a, b = .731, .123
    e = r.normal(size=(nx, ny+1))
    hx = r.normal(size=(nx, ny))
    hy = r.normal(size=(nx+1, ny+1))
    me = (np.arange(ny+1) <= p).astype(float)
    mhx = (np.arange(ny) < p).astype(float)
    # Advance H on full and masked incident fields, then inject M_x.
    hx_true = hx-b*np.diff(e, axis=1)/dy
    hy_true = hy.copy()
    hy_true[1:-1] += b*np.diff(e, axis=0)/dx
    hx_raw = hx*mhx-b*np.diff(e*me, axis=1)/dy
    hy_raw = hy*me
    hy_raw[1:-1] += b*np.diff(e*me, axis=0)/dx
    moment = dx*dz*e[:,p]
    hx_raw[:,p] -= b*moment/(dx*dy*dz)
    h_error = max(np.max(abs(hx_raw-hx_true*mhx)), np.max(abs(hy_raw-hy_true*me)))
    # Advance E at n+1 from corrected H(n+1/2), then inject J_z.
    true_new = e[:,1:-1] + a*(np.diff(hy_true, axis=0)[:,1:-1]/dx-np.diff(hx_true, axis=1)/dy)
    raw_new = e[:,1:-1]*me[1:-1] + a*(np.diff(hy_raw, axis=0)[:,1:-1]/dx-np.diff(hx_raw, axis=1)/dy)
    current = dx*hx_true[:,p]
    raw_new[:,p-1] -= a*current*dz/(dx*dy*dz)
    e_error = float(np.max(abs(raw_new-true_new*me[1:-1])))
    assert h_error < 1e-13 and e_error < 1e-13
    return dict(mask_magnetic_max_abs_error=float(h_error), mask_electric_max_abs_error=e_error)


def histories(directory):
    result = []
    with h5py.File(directory/'full_air/model.h5') as reference:
        dt = float(reference.attrs['dt']); dl = float(reference.attrs['dx_dy_dz'][0])
        incident = {str(g.attrs['Name']):g for g in reference['rxs'].values()}
        for name in ['compact_air', 'compact_cover']:
            filename = directory/name/'model.h5'
            if not filename.exists():
                continue
            with h5py.File(filename) as compact:
                count = int(compact.attrs['Iterations'])
                worst = 0.; seen = 0
                for g in compact['srcs'].values():
                    waveform = str(g['excitation'].attrs['WaveformID'])
                    index = int(waveform[1:])
                    samples = g['excitation/samples'][:]
                    assert samples.dtype == np.float64 and samples.size == count
                    if waveform.startswith('J'):
                        target = dl*incident[f'inc{index}']['Hx'][1:count+1]
                        offset = .5*dt
                    elif waveform.startswith('M'):
                        target = dl*dl*incident[f'inc{index}']['Ez'][:count]
                        offset = 0.
                    else:
                        raise AssertionError('Unexpected virtual source')
                    assert float(g['excitation'].attrs['TimeSampleOffset']) == offset
                    peak = np.max(abs(target))
                    err = np.max(abs(samples-target))/peak if peak else np.max(abs(samples-target))
                    worst = max(worst, float(err)); seen += 1
                assert worst < 1e-11
                result.append(dict(case=name, sources_checked=seen, max_normalized_sample_error=worst))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path)
    args = parser.parse_args()
    results = mask_identity()
    if args.directory:
        results['saved_source_histories'] = histories(args.directory)
    print(json.dumps(results, ensure_ascii=False, indent=2))
