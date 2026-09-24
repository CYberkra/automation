"""Sensitivity of one completed M00 record using only official SFCW transforms.

Not FDTD, not measured hardware, not a universal taper-selection rule.
Predeclared windows 200/300/full ns and taper fractions 0/.1/.25/.5.
"""
import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import epsilon_0
from gprMax.toolboxes.SFCW.processing import (load_source, load_receiver,
    direct_frequency_response, reconstruct_time_response, write_sfcw_output)
from audit_official_sfcw import transverse_dipole


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw',type=Path)
    p.add_argument('--raw-audit',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    audit=json.loads(a.raw_audit.read_text(encoding='utf-8'))
    if not audit['passed']: raise SystemExit('Require passing raw audit')
    if hashlib.sha256(a.raw.read_bytes()).hexdigest()!=audit['input_sha256']:
        raise SystemExit('Raw file changed since audit')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    source=load_source(a.raw,audit['source_path'])
    receiver=load_receiver(a.raw,receiver_path='name:measurement',component='Ex')
    frequency=np.linspace(20e6,170e6,501)
    theory=transverse_dipole(frequency,1.3,source.spatial_scale)
    rows=[]
    frequency_rows=[]
    for requested_ns in (200.,300.,None):
        n=len(receiver.samples) if requested_ns is None else min(len(receiver.samples),int(np.floor(requested_ns*1e-9/receiver.dt))+1)
        selected=replace(receiver,samples=receiver.samples[:n])
        for taper in (0.,.1,.25,.5):
            response=direct_frequency_response(source,selected,frequency,tail_taper_fraction=taper)
            if not np.all(response.source_valid) or not np.all(np.isfinite(response.response)):
                raise ValueError('Invalid frequency output')
            relative=np.abs(response.response-theory)/abs(theory)
            amp=20*np.log10(abs(response.response)/abs(theory))
            phase=np.angle(response.response/theory,deg=True)
            label=f'n{n}_taper{taper:g}'
            time=reconstruct_time_response(response,window='rectangular',zero_pad_factor=1,time_shift=0)
            write_sfcw_output(a.output_dir/f'{label}.h5',response,time)
            tail=selected.samples[-max(8,int(np.ceil(.05*n))):]
            rows.append(dict(label=label,n=n,last_sample_ns=float(selected.times[-1]*1e9),taper_fraction=taper,
                median_tail_V_m=float(np.median(tail)),tail_peak_to_peak_V_m=float(np.ptp(tail)),
                max_relative_complex_error=float(max(relative)),median_relative_complex_error=float(np.median(relative)),
                max_absolute_amplitude_difference_dB=float(max(abs(amp))),max_absolute_phase_difference_deg=float(max(abs(phase)))))
            frequency_rows.extend(zip([label]*501,frequency,relative,amp,phase))
    with (a.output_dir/'frequency_sensitivity.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['case','frequency_Hz','relative_complex_error','amplitude_difference_dB','wrapped_phase_difference_deg']);writer.writerows(frequency_rows)
    q=float(np.sum(source.samples)*source.dt)
    result=dict(raw_sha256=audit['input_sha256'],solver_invoked=False,official_transform=True,
        source_charge_integral_A_s=q,continuous_static_tail_prediction_V_m=-source.spatial_scale*q/(4*np.pi*epsilon_0*1.3**3),
        rows=rows,physical_accuracy_threshold=None,convergence_certified=False,
        caution='Truncation/taper sensitivity on one air record only; no fitted correction and no field-data performance claim.')
    (a.output_dir/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':main()
