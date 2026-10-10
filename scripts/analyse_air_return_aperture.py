"""Separate padding/quadrature checks for existing flat-ground return spectra.

No solver, new source, time window, SFCW tones, taper or amplitude fit.
The angular spectrum is a project diagnostic, not an official gprMax tool.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def analyse(root):
    start = time.perf_counter()
    with np.load(root/'spectra.npz') as d:
        f = d['frequency_hz']; xs = d['x_m']
        full = d['full_cover_return']-d['full_air_return']
        compact = d['compact_cover_return']-d['compact_air_return']
        truth = d['full_cover_original']-d['full_air_original']
    c = json.loads((root/'contract.json').read_text(encoding='utf-8'))
    dl = c['mesh_m']
    dt = 1/(299792458.*np.sqrt(2)/dl)
    distance = c['original_rx_logical_m'][1]-c['return_plane_y_m']
    ix = int(np.argmin(abs(xs-c['original_rx_logical_m'][0])))
    wt = 2*np.sin(np.pi*f*dt)/dt
    rows = []; values = {}
    for factor in [1,2,4,8,16,32]:
        count = factor*len(xs)
        q = 2*np.pi*np.fft.fftfreq(count,d=dl)
        qt = 2*np.sin(q*dl/2)/dl
        rootky = np.sqrt((wt[:,None]/299792458.)**2-qt[None,:]**2+0j)
        ky = 2/dl*np.arcsin(dl*rootky/2)
        ky = ky.real-1j*np.abs(ky.imag)
        kernel = np.exp(-1j*ky*distance)
        for label, field in [('full',full),('compact',compact)]:
            predicted = np.fft.ifft(np.fft.fft(field,n=count,axis=1)*kernel,axis=1)[:,ix]
            err = np.linalg.norm(predicted-truth)/np.linalg.norm(truth)
            low=f<=40e6
            lowerr=np.linalg.norm((predicted-truth)[low])/np.linalg.norm(truth[low])
            rows.append(dict(padding_factor=factor,plane=label,fullband_relative_l2=float(err),low20_40_relative_l2=float(lowerr)))
            values[f'{label}_pad{factor}'] = predicted
    target=root/'return_aperture_checks.json'
    target.write_text(json.dumps(dict(rows=rows,analysis_wall_s=time.perf_counter()-start,
                                     interpretation='Padding tests periodic-image sensitivity, not missing physical aperture or time convergence'),indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(root/'return_aperture_checks.npz',frequency_hz=f,truth_scattered=truth,**values)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    analyse(p.parse_args().root)
