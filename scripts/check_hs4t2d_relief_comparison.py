"""Independently verify completed V4 raw data and actual cell-material exports."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new verification file required')
    c=json.loads(args.contract.read_text('utf-8'))
    root=Path(__file__).resolve().parents[1]
    run=Path(c['run_directory'])
    result={}
    for g in c['groups']:
        folder=Path(g['input']).parent
        if sha256(g['input'])!=g['sha256']:
            raise ValueError('series input changed')
        expected=np.full((660,1,240),2,dtype=np.uint32)
        expected[:240]=3
        for line in Path(g['input']).read_text('utf-8').splitlines():
            words=line.split()
            if line.startswith('#box:') and words[-1]=='cover':
                x0,y0,z0,x1,y1,z1=map(float,words[1:7])
                expected[round(z0/.05):round(z1/.05),:,round(x0/.05):round(x1/.05)]=4
        audit=json.loads((folder/'audit.json').read_text('utf-8'))
        if len(audit)!=121:
            raise ValueError('incomplete audit')
        raw_difference=[]
        for k,row in enumerate(audit,1):
            path=folder/f'profile{k}.h5'
            if sha256(path)!=row['sha256']:
                raise ValueError('output hash differs')
            with h5py.File(path,'r') as h:
                raw=h['rxs/rx1/Ey'][:]
                if raw.shape!=(5089,) or raw.dtype!=np.float64 or not np.isfinite(raw).all():
                    raise ValueError('raw validity failed')
                if str(h.attrs['gprMax'])!='4.0.0':
                    raise ValueError('not V4.0.0')
                with h5py.File(root/f'artifacts/research_checks/2026-10-02_hs4t2d_transect/hs4t2d_t{k:02d}.h5','r') as old:
                    for key in ['srcs/src1','rxs/rx1']:
                        if not np.array_equal(h[key].attrs['GridPosition'],old[key].attrs['GridPosition']):
                            raise ValueError('grid station differs')
                    if g['id']=='baseline':
                        raw_difference.append(float(np.max(np.abs(raw-old['rxs/rx1/Ey'][:]))))
            with h5py.File(folder/f'hs4t2d_geom{k}.vtkhdf','r') as h:
                actual=h['VTKHDF/CellData/Material'][:]
                if not np.array_equal(actual,expected):
                    raise ValueError('actual cell material geometry differs')
                if not np.array_equal(h['VTKHDF'].attrs['Spacing'],[.05,.05,.05]):
                    raise ValueError('geometry spacing differs')
        result[g['id']]={'raw_h5_verified':121,'actual_cell_material_grids_verified':121,
                         'cells_per_grid':int(expected.size),'material_grid_mismatches':0,
                         'version':'4.0.0','raw_dtype':'float64',
                         'archived_baseline_raw_max_abs_difference':max(raw_difference) if raw_difference else None}
    reference=Path(c['cache_verification_reference'])
    for k in range(1,8):
        with h5py.File(reference/f'profile{k}.h5','r') as a,h5py.File(run/'baseline'/f'profile{k}.h5','r') as b:
            if not np.array_equal(a['rxs/rx1/Ey'][:],b['rxs/rx1/Ey'][:]):
                raise ValueError('compile cache changed raw output')
    report={'status':'PASS','groups':result,'cache_uncached_raw_bit_identical_traces':7,
            'checker_sha256':sha256(__file__),'contract_sha256':sha256(args.contract),
            'scope':'execution, actual cell material geometry and numerical reproducibility; not physical waveform validation'}
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
