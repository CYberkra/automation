"""Bounded one-station H0 boundary-distance controls; preserve native grid and exact source."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np

import line9_basal_pair as pair
import hs4_station_grid_controls as supervisor

ROOT=Path(__file__).resolve().parents[1]
CASES=[('base',0,0),('top40',0,40),('sides40',40,0),('joint40',40,40)]


def prepare(source,reference,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh outputs required')
    m=json.loads((source/'manifest.json').read_text('utf-8'))
    old=source/'H0/geometries/full2d_compact.h5';db=source/'H0/geometries/line9_research_materials_v1_smoothed.json'
    assert pair.sha(old)==m['groups'][1]['geometry_sha256'] and pair.sha(db)==m['material_sha256']
    with h5py.File(old) as h:data=h['data'][:]
    assert data.shape==(4400,1600,1) and np.all(data[:,-1,0]==0)
    original=(source/'H0/cases/full2d_c0008/profile.in').read_text('utf-8')
    base='\n'.join(line for line in original.splitlines() if not line.startswith('#snapshot:'))+'\n'
    out.mkdir(parents=True);shutil.copyfile(reference,out/'reference_H0.h5');pair.save(out/'parent_manifest.json',m)
    groups=[]
    for name,padding,top in CASES:
        folder=out/name;gp=folder/'geometries';gp.mkdir(parents=True)
        p=folder/'cases/full2d_c0008/profile.in';p.parent.mkdir(parents=True)
        shift=round(padding/.025);height=round(top/.025)
        v=np.pad(data,((shift,shift),(0,0),(0,0)),mode='edge')
        v=np.pad(v,((0,0),(0,height),(0,0)),mode='constant',constant_values=0)
        np.testing.assert_array_equal(v[shift:shift+4400,:1600],data)
        if height:assert np.all(v[:,1600:]==0)
        shutil.copyfile(old,gp/old.name);shutil.copyfile(db,gp/db.name)
        with h5py.File(gp/old.name,'r+') as h:
            del h['data'];h.create_dataset('data',data=v,compression='gzip',compression_opts=4)
            h.attrs['shape_nxyz']=v.shape
            h.attrs['BoundaryControl']='Original interior exact; side edge-column continuation and/or top free space extension'
        card=[]
        for line in base.splitlines():
            words=line.split()
            if line.startswith('#domain:'):line=f'#domain: {110+2*padding} {40+top} inf'
            elif line.startswith('#hertzian_dipole:'):words[2]=str(78.6+padding);line=' '.join(words)
            elif line.startswith('#rx:'):words[1]=str(79.9+padding);line=' '.join(words)
            card.append(line)
        p.write_text('\n'.join(card)+'\n',encoding='utf-8')
        groups.append(dict(id=name,input=str(p.relative_to(out)),input_sha256=pair.sha(p),geometry=str((gp/old.name).relative_to(out)),geometry_sha256=pair.sha(gp/old.name),
            material=str((gp/db.name).relative_to(out)),material_sha256=pair.sha(db),padding_x_m=padding,top_extension_m=top,
            native_shape=list(v.shape),domain_m=[110+2*padding,40+top,.025],tx_m=[78.6+padding,35.,0.],rx_m=[79.9+padding,35.,0.]))
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        path=out/g['input'];get_user_objects(path.read_text('utf-8').splitlines(),input_dir=path.parent)
    manifest=dict(status='PREPARED_APPROVED_BOUNDARY_CONTROLS_NOT_RUN',approval_basis='User persistent objective: 研究明白; existing autonomous research and SSH simulation authorization.',
        parent_manifest_sha256=pair.sha(source/'manifest.json'),parent_H0_native_sha256=pair.sha(reference),groups=groups,
        invariants='Exact original interior voxel map, 2.5cm native grid, dt/800ns,40A impulse,Yee timing,material database,Tx/Rx relative to ground,FP64,HORIPML80cells; snapshots omitted in all4cases.',
        factors='Top extension: only add40m free space; side extension: add40m both sides by constant edge-column geological continuation, translate source/rx40m; joint tests interaction.',
        processing='Exact501tones20-170MHz/0.3MHz,source-spatial/time normalization,no tail taper,Hann+Blackman,complex comparison; no AGC.',
        predefined_gates_ns={'bottom':[332.31137724550894,356.31137724550894],'upper':[217.54091816367264,267.32268795741845]},
        snapshot_count=0,limits='One selected H0 station; side extension changes exterior continuation and boundary distance together. Top is an air-only distance control. Not whole-line or real-data validation.')
    pair.save(out/'manifest.json',manifest);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',manifest)


def audit(path,completed=False):
    c=json.loads(path.read_text('utf-8'));pair.check_files(c);m=c['study_manifest'];rows=[]
    with h5py.File(Path(c['package'])/'reference_H0.h5') as ref:
        for g in c['groups']:
            raw=Path(g['input']).with_suffix('.h5');r=dict(id=g['id'],input_sha256=pair.sha(g['input']))
            if completed:
                with h5py.File(raw) as h:
                    assert str(h.attrs['gprMax'])=='4.0.0'
                    for key in ['Iterations','dt','dx_dy_dz']:np.testing.assert_array_equal(h.attrs[key],ref.attrs[key])
                    np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape'])
                    for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
                    x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
                    assert x.dtype==s.dtype==np.float64 and np.isfinite(x).all() and x.shape==ref['rxs/rx1/Ez'].shape
                    np.testing.assert_array_equal(s,ref['srcs/src1/excitation/samples'][:])
                    r.update(native_sha256=pair.sha(raw),native_path=str(raw),dtype=str(x.dtype),dt_s=float(h.attrs['dt']),samples=len(x))
            rows.append(r)
    return dict(status='PASS_NATIVE_IDENTITY_NOT_CAUSAL_CLASSIFICATION',completed=completed,contract_sha256=pair.sha(path),groups=rows)


def freeze(package,out,parent):
    if out.exists():raise ValueError('Fresh execution capsule required')
    c=json.loads(parent.read_text('utf-8'));pair.check_files(c);hardware=pair.resources()
    m=json.loads((package/'manifest.json').read_text('utf-8'))
    groups=[dict(g,input=str((package/g['input']).resolve())) for g in m['groups']]
    cells=max(np.prod(g['native_shape']) for g in groups)
    estimate=int(cells)*420+2*2**30
    if hardware['free_VRAM_bytes']<estimate:raise RuntimeError('Conservative device budget rejected')
    for g in groups:
        for k in ['input','geometry','material']:
            p=Path(g['input']) if k=='input' else package/g[k]
            assert pair.sha(p)==g[k+'_sha256']
    c.update(package=str(package),groups=groups,max_runs=len(groups),study_manifest=m,hardware_at_freeze=hardware,
        parent_execution_contract_sha256=pair.sha(parent),approval_basis=m['approval_basis'],max_batch_wall_s=3600,
        cancel_file=str(out/'USER_STOP'),lease_file=str(out/'session_heartbeat'),max_lease_age_s=600,
        file_identities={str(p.resolve()):pair.sha(p) for p in package.rglob('*') if p.is_file()},
        code_identities={str(ROOT/'scripts'/name):pair.sha(ROOT/'scripts'/name) for name in ['line9_boundary_controls.py',*pair.CODE]},
        conservative_device_estimate_bytes=estimate,limits=m['limits'])
    c.pop('prepared_manifest',None)
    out.mkdir(parents=True);pair.save(out/'execution_contract.json',c);(out/'session_heartbeat').write_text('Active authorized boundary controls\n')
    pair.save(out/'preflight_verification.json',audit(out/'execution_contract.json'))
    print(json.dumps(dict(contract_sha256=pair.sha(out/'execution_contract.json'),device_estimate_GiB=estimate/2**30,groups=len(groups))))


def run(out):
    supervisor.audit=audit;supervisor.run(out/'execution_contract.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for name in ['source','reference','out','evidence','package','parent']:p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.source,a.reference,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':run(a.out.resolve())
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
