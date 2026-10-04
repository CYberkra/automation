"""Verify portable completed-case catalog using tracked native files only."""
import argparse
from pathlib import Path
import h5py
import numpy as np
from hs4_cross_pc import ROOT, read, save
from hs_capsule_identity import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--catalog',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.out.exists(): raise ValueError('new verification output required')
    catalog=read(a.catalog); records=[]
    for g in catalog['completed_cases']:
        for item in ('input','raw'):
            if sha256(ROOT/g[item])!=g[item+'_sha256']: raise ValueError('completed identity changed')
        with h5py.File(ROOT/g['raw']) as h:
            v=h['rxs/rx1/Ey'][:]
            if str(h.attrs['gprMax'])!='4.0.0' or v.dtype!=np.float64 or not np.isfinite(v).all():
                raise ValueError('raw version/float64/finite failed')
            for key,expected in [('srcs/src1',g['tx_grid']),('rxs/rx1',g['rx_grid'])]:
                if not np.array_equal(h[key].attrs['GridPosition'],expected): raise ValueError('actual pose differs')
            if not np.array_equal(h.attrs['dx_dy_dz'],g['spacing_m']): raise ValueError('actual grid differs')
        records.append(g['case_id'])
    for g in catalog['reconstruction_files']:
        if sha256(ROOT/g['file'])!=g['sha256']: raise ValueError('reconstruction changed')
    save(a.out,{'status':'PASS','count':len(records),'case_ids':records,'catalog_sha256':sha256(a.catalog),
                'scope':'Tracked native identities, actual version/precision/grid/pose; historical material-map acceptance remains archived. No new solve.'})
    print('PASS',len(records),'completed cases')


if __name__=='__main__': main()
