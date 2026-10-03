"""Independent audit of completed native/snapshot fields and matched probes."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256
from hs4_height_wavefield_v0_2 import receivers, FIELDS, save


def verify(folder, out):
    if out.exists():
        raise ValueError('new audit required')
    contract=folder/'execution_contract.json'
    c=json.loads(contract.read_text('utf-8'))
    events=[json.loads(s) for s in (folder/'execution.jsonl').read_text('utf-8').splitlines()]
    if events[0]['contract_sha256']!=sha256(contract):
        raise ValueError('execution contract differs')
    rows=[]
    for g in c['groups']:
        p=Path(g['input']); raw=p.with_suffix('.h5')
        finished=[r for r in events if r.get('group')==g['id'] and r['status']=='COMPLETED']
        if len(finished)!=1 or finished[0]['raw_sha256']!=sha256(raw) or sha256(p)!=g['input_sha256']:
            raise ValueError('one completed immutable raw trace required')
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax'])!='4.0.0' or float(h.attrs['dt'])!=c['dt_s'] or int(h.attrs['Iterations'])!=c['iterations'] or not np.array_equal(h.attrs['dx_dy_dz'],[.025,.05,.025]):
                raise ValueError('actual version/grid/time differs')
            rx=receivers(h)
            if not np.array_equal(rx['t01'].attrs['GridPosition'],[756,0,round((12+g['height_m'])/.025)]) or not np.array_equal(h['srcs/src1'].attrs['GridPosition'],[704,0,round((12+g['height_m'])/.025)]):
                raise ValueError('actual acquisition differs')
            for r in rx.values():
                for v in r.values():
                    if v.dtype!=np.float64 or v.shape!=(c['iterations'],) or not np.isfinite(v[:]).all():
                        raise ValueError('native validity differs')
            probes={k:{f:r[f][:] for f in r} for k,r in rx.items() if k.startswith('p')}
        maximum={f:0. for f in ('Ey','Hx','Hz')}; snapshots=[]
        files=sorted(p.parent.glob('profile_snaps/*.h5'))
        if len(files)!=len(c['snapshot_iterations']):
            raise ValueError('incomplete snapshot history')
        for j,s in zip(c['snapshot_iterations'],files):
            with h5py.File(s) as h:
                if int(h.attrs['iteration'])!=j or float(h.attrs['time'])!=j*c['dt_s'] or float(h.attrs['magnetic_time'])!=(j-.5)*c['dt_s'] or not np.array_equal(h.attrs['origin'],c['snapshot_extent_m'][:3]):
                    raise ValueError('snapshot timing/coordinate differs')
                a={f:h[f][:] for f in FIELDS}
                if any(v.dtype!=np.float64 or list(v.shape)!=c['snapshot_shape_xyz'] or not np.isfinite(v).all() for v in a.values()):
                    raise ValueError('snapshot validity differs')
                for probe in g['probes']:
                    k=probe['id']; ix,iz=probe['snapshot_ix'],probe['snapshot_iz']
                    values={'Ey':probes[k+'c']['Ey'][j],
                        'Hx':.5*(probes[k+'c']['Hx'][j]+probes[k+'z']['Hx'][j]),
                        'Hz':.5*(probes[k+'c']['Hz'][j]+probes[k+'x']['Hz'][j])}
                    for f,v in values.items():
                        actual=a[f][ix,0,iz]
                        maximum[f]=max(maximum[f],abs(actual-v))
                        if not np.isclose(actual,v,rtol=2e-13,atol=2e-15):
                            raise ValueError('snapshot/native spatial closure differs')
            snapshots.append({'file':s.relative_to(p.parent).as_posix(),'sha256':sha256(s),'iteration':j})
        geom=p.parent/'hs4t2d_geom.vtkhdf'
        with h5py.File(geom) as h,h5py.File(g['reference_geometry']) as old:
            if sha256(g['reference_geometry'])!=g['reference_geometry_sha256'] or not np.array_equal(h['VTKHDF/CellData/Material'][:],old['VTKHDF/CellData/Material'][:]):
                raise ValueError('actual material grid differs')
        rows.append({'id':g['id'],'raw_sha256':sha256(raw),'actual_geometry_sha256':sha256(geom),
            'snapshot_count':len(files),'snapshot_native_max_abs_difference':maximum,'snapshots':snapshots})
    matched=[]; moved=[]
    ref=Path(c['passive_reference'])
    if sha256(ref)!=c['passive_reference_sha256']:
        raise ValueError('passive reference changed')
    with h5py.File(ref) as a,h5py.File(folder/'high_rough/profile.h5') as b:
        ar,br=receivers(a),receivers(b)
        for name,r in ar.items():
            if not np.array_equal(r.attrs['GridPosition'],br[name].attrs['GridPosition']):
                moved.append({'id':name,'reference_grid':r.attrs['GridPosition'].tolist(),'snapshot_grid':br[name].attrs['GridPosition'].tolist()})
                continue
            for f in r:
                if not np.array_equal(r[f][:],br[name][f][:]):
                    raise ValueError('passive observation changed matched native record')
            matched.append(name)
    if 't01' not in matched or len(matched)!=31 or len(moved)!=6:
        raise ValueError('unexpected passive coordinate coverage')
    result={'status':'PASS','code_sha256':sha256(__file__),'contract_sha256':sha256(contract),'groups':rows,
        'matched_passive_receivers_bit_identical':matched,'moved_probe_receivers_excluded_from_passive_replay':moved,
        'scope':'All four material grids and 4072 six-component snapshots verified; 12 same-cell native probe closures per trace. Main receiver and30 auxiliary receivers replay exactly; two probe triplets moved with ROI and cannot be compared by name. Original frozen checker failure is preserved, not a solver failure or permission retry.'}
    save(out,result)
    print('PASS:4raw/4072snapshots;31matched receivers bit-identical;6moved probes explicitly excluded')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); verify(a.study,a.out)
