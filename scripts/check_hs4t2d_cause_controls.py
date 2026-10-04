"""Independently check cause-control factors and actual material grids."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--completed',action='store_true')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new verification required')
    c=json.loads(args.contract.read_text('utf-8'))
    root=Path(__file__).resolve().parents[1]
    z=np.genfromtxt(root/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',delimiter=',',names=True)['grid_z_m']
    fixed=['#title:','#time_window:','#omp_threads:','#pml_formulation:','#waveform:','#material:','#add_dispersion_debye:']
    groups={}
    for g in c['groups']:
        path=args.contract.parent/g['id']/'profile.in'
        old=root/g['source_input']
        if sha256(path)!=g['sha256'] or sha256(old)!=g['source_sha256']:
            raise ValueError('input identity differs')
        lines=path.read_text('utf-8').splitlines()
        expected_source=old.read_text('utf-8')
        if g.get('source_condition')=='ricker95':
            expected_source=expected_source.replace('#waveform: impulse 1 1 impulse','#waveform: ricker 1 95e6 impulse')
        elif g.get('source_condition','impulse')!='impulse':
            raise ValueError('undeclared source condition')
        if [s for s in lines if any(s.startswith(p) for p in fixed)]!=[s for s in expected_source.splitlines() if any(s.startswith(p) for p in fixed)]:
            raise ValueError('source/material/time invariants differ')
        nx,nz=round(g['width_m']/.05),round(33/g['dz_m'])
        expected=np.full((nz,1,nx),2,dtype=np.uint32)
        expected[:round(12/g['dz_m'])]=3
        x=(np.arange(nx)+.5)*.05-g['x_translation_m']
        lower=np.zeros(nx) if g['halfspace'] else z[np.clip(np.floor(x/.25).astype(int),0,47)]
        for j,bottom in enumerate(lower):
            expected[round(bottom/g['dz_m']):round(12/g['dz_m']),0,j]=4
        parsed=np.full_like(expected,2)
        for s in lines:
            if s.startswith('#box:'):
                w=s.split()
                x0,y0,z0,x1,y1,z1=map(float,w[1:7])
                parsed[round(z0/g['dz_m']):round(z1/g['dz_m']),:,round(x0/.05):round(x1/.05)]={'rock':3,'cover':4}[w[-1]]
        if not np.array_equal(parsed,expected):
            raise ValueError('input material map is not the independently declared control')
        if f"#pml_cells: 20 0 {round(1/g['dz_m'])} 20 0 {round(1/g['dz_m'])}" not in lines:
            raise ValueError('physical PML thickness differs')
        material_hash=__import__('hashlib').sha256(expected.tobytes()).hexdigest()
        grids=[]
        if args.completed:
            audit=json.loads((path.parent/'audit.json').read_text('utf-8'))
            if len(audit)!=g['traces']:
                raise ValueError('incomplete raw audit')
            for k,row in enumerate(audit,1):
                p=path.parent/row['file']
                if sha256(p)!=row['sha256']:
                    raise ValueError('raw identity changed')
                with h5py.File(p,'r') as h:
                    a=h['rxs/rx1/Ey'][:]
                    if str(h.attrs['gprMax'])!='4.0.0' or a.dtype!=np.float64 or not np.isfinite(a).all():
                        raise ValueError('raw V4/FP64 validity differs')
                geom=path.parent/(f'hs4t2d_geom{k}.vtkhdf' if g['traces']>1 else 'hs4t2d_geom.vtkhdf')
                with h5py.File(geom,'r') as h:
                    if not np.array_equal(h['VTKHDF/CellData/Material'][:],expected):
                        raise ValueError('actual solver grid differs from independent declaration')
                grids.append({'file':geom.name,'sha256':sha256(geom)})
        groups[g['id']]={'input_factor_and_invariants':'PASS','material_array_sha256':material_hash,
                         'actual_material_grids_verified':len(grids),'geometry_exports':grids}
    if args.completed:
        events=[json.loads(s) for s in (args.contract.parent/'execution.jsonl').read_text('utf-8').splitlines()]
        if events[-1].get('status')!='COMPLETED' or events[-1].get('traces')!=c['max_runs'] or events[0]['contract_sha256']!=sha256(args.contract):
            raise ValueError('incomplete or changed frozen run')
    report={'status':'PASS','groups':groups,'completed_outputs_checked':args.completed,
            'code_sha256':sha256(__file__),'contract_sha256':sha256(args.contract),
            'scope':'control factors, execution and actual cell-material grid; not physical attribution'}
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('PASS',len(groups),'groups;',sum(g['actual_material_grids_verified'] for g in groups.values()),'actual grids')


if __name__=='__main__':
    main()
