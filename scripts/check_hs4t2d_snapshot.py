"""Independent portable audit of passive snapshot outputs and source identity."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new audit file required')
    c=json.loads((args.run/'execution_contract.json').read_text('utf-8'))
    e=json.loads((args.run/'execution.json').read_text('utf-8'))
    if e['status']!='COMPLETED' or e['contract_sha256']!=sha256(args.run/'execution_contract.json'):
        raise ValueError('completed contract identity differs')
    if sha256(args.source)!=c['source_input_sha256'] or sha256(args.reference)!=c['reference_raw_sha256']:
        raise ValueError('reference identities differ')
    original=args.source.read_bytes()
    modified=(args.run/'profile.in').read_bytes()
    if not modified.startswith(original) or sha256(args.run/'profile.in')!=c['input_sha256']:
        raise ValueError('model directives changed')
    expected='\n'.join(f'#snapshot: 0 0 0 12 0.05 33 0.05 0.05 0.05 {n}e-9 wave_{n:03d}.h5'
                       for n in c['snapshot_times_ns'])+'\n'
    if modified[len(original):]!=expected.encode('ascii'):
        raise ValueError('appended directives differ')
    with h5py.File(args.run/'profile.h5','r') as h,h5py.File(args.reference,'r') as old:
        if not np.array_equal(h['rxs/rx1/Ey'][:],old['rxs/rx1/Ey'][:]):
            raise ValueError('observer changed receiver')
        dt=float(h.attrs['dt'])
        for path in ['srcs/src1','rxs/rx1']:
            if not np.array_equal(h[path].attrs['GridPosition'],old[path].attrs['GridPosition']):
                raise ValueError('source/receiver grid position changed')
    if len(e['snapshots'])!=9:
        raise ValueError('incomplete snapshots')
    for row in e['snapshots']:
        p=args.run/row['file']
        if sha256(p)!=row['sha256']:
            raise ValueError('snapshot identity differs')
        with h5py.File(p,'r') as h:
            if str(h.attrs['gprMax'])!='4.0.0' or h.attrs['time']!=h.attrs['iteration']*dt:
                raise ValueError('snapshot version/electric time differs')
            if h.attrs['magnetic_time']!=(h.attrs['iteration']-.5)*dt:
                raise ValueError('magnetic time staggering differs')
            if not np.array_equal(h.attrs['nx_ny_nz'],[240,1,660]):
                raise ValueError('snapshot dimensions differ')
            for component in ['Ex','Ey','Ez','Hx','Hy','Hz']:
                value=h[component][:]
                if value.shape!=(240,1,660) or value.dtype!=np.float64 or not np.isfinite(value).all():
                    raise ValueError('snapshot component invalid')
    report={'status':'PASS','snapshots':9,'float64_finite_components_verified':54,
            'unchanged_source_bytes_except_snapshot_append':True,'raw_receiver_bit_identical':True,
            'electric_time_n_dt_and_magnetic_half_step_verified':True,
            'checker_sha256':sha256(__file__),'contract_sha256':sha256(args.run/'execution_contract.json'),
            'scope':'output provenance, passive observer, dtype and time staggering; not propagation acceptance'}
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
