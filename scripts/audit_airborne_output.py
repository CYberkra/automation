"""Read-only raw V4 HDF5 audit for the agreed single M00_x_3d geometry.

No solver or SFCW reimplementation. Does not certify grid/boundary convergence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np


def scalar(value):
    if isinstance(value,bytes):return value.decode('utf-8')
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    return value


def attrs(group):return {k:scalar(v) for k,v in group.attrs.items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise SystemExit('Refusing overwrite')
    with h5py.File(a.input,'r') as f:
        root=attrs(f)
        receivers=[g for g in f['rxs'].values() if scalar(g.attrs.get('Name',''))=='measurement']
        if len(receivers)!=1:raise ValueError('Expected unique named receiver')
        rx=receivers[0]; field=rx['Ex']; y=field[:]
        sources=[g for g in f['srcs'].values() if 'excitation' in g and 'samples' in g['excitation']]
        if len(sources)!=1:raise ValueError('Expected unique sampled source')
        src=sources[0]; excitation=src['excitation']; x=excitation['samples'][:]
        dt=float(root['dt'])
        checks={'version':scalar(root['gprMax'])=='4.0.0',
                'receiver_position':bool(np.allclose(rx.attrs['Position'],[12,12.6,19],rtol=0,atol=1e-10)),
                'source_position':bool(np.allclose(src.attrs['Position'],[12,11.3,19],rtol=0,atol=1e-10)),
                'receiver_float64':y.dtype==np.dtype('float64'),
                'receiver_finite_nonzero':bool(np.all(np.isfinite(y)) and np.any(y!=0)),
                'source_finite_impulse':bool(np.all(np.isfinite(x)) and np.count_nonzero(x)==1 and x[0]!=0),
                'sample_count':len(y)==int(root['Iterations'])==len(x),
                'source_half_step':bool(np.isclose(float(excitation.attrs['TimeSampleOffset']),dt/2,rtol=1e-10,atol=0)),
                'source_interval':bool(np.isclose(float(excitation.attrs['SampleInterval']),dt,rtol=1e-10,atol=0)),
                'dipole_length':bool(np.isclose(float(excitation.attrs['SpatialScale']),.05,rtol=1e-10,atol=0))}
        peak=float(np.max(np.abs(y))); tail=float(np.max(np.abs(y[-max(1,int(np.ceil(.05*len(y)))):])))
        result={'input':str(a.input),'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),
            'root_attrs':root,'receiver_path':rx.name,'receiver_attrs':attrs(rx),
            'receiver_dataset_attrs':attrs(field),'source_path':src.name,'source_attrs':attrs(src),
            'excitation_attrs':attrs(excitation),'checks':checks,'passed':all(checks.values()),
            'max_abs_Ex_V_m':peak,'peak_sample_time_ns':int(np.argmax(np.abs(y)))*dt*1e9,
            'last_5pct_peak_relative_dB':20*np.log10(tail/peak) if tail>0 and peak>0 else None,
            'solver_called_by_audit':False,'convergence_certified':False}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':result['passed'],'checks':checks}))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
