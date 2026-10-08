"""Single-station geological/direct-field contrasts for the persistent H0 event."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np

import line9_basal_pair as pair
import line9_boundary_controls as boundary
import hs4_station_grid_controls as supervisor


def prepare(source,reference,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh package required')
    m=json.loads((source/'manifest.json').read_text('utf-8'))
    old=source/'H0/geometries/full2d_compact.h5';db=source/'H0/geometries/line9_research_materials_v1_smoothed.json'
    assert pair.sha(old)==m['groups'][1]['geometry_sha256']
    with h5py.File(old) as h:data=h['data'][:]
    assert data.shape==(4400,1600,1)
    no_interbed=data.copy();no_interbed[no_interbed==3]=2
    homogeneous=data.copy();homogeneous[homogeneous!=0]=1
    flat=np.repeat(data[3170:3171],4400,axis=0)
    bottom=np.pad(data,((0,0),(1600,0),(0,0)),mode='edge')
    np.testing.assert_array_equal(bottom[:,1600:],data)
    cases=[('no_interbed',no_interbed,'上覆砂岩夹层→泥岩'),('homogeneous_cover',homogeneous,'保留地表，地下全为覆盖层'),('free_space',np.zeros_like(data),'全空气：仅直耦/数值响应'),('flat_layers',flat,'各层横向拉平，取原中点体素列'),('bottom40',bottom,'只移远底边40m，原地质不变')]
    original=(source/'H0/cases/full2d_c0008/profile.in').read_text('utf-8')
    card='\n'.join(line for line in original.splitlines() if not line.startswith('#snapshot:'))+'\n'
    out.mkdir(parents=True);shutil.copyfile(reference,out/'reference_H0.h5')
    groups=[]
    for name,volume,label in cases:
        gp=out/name/'geometries';gp.mkdir(parents=True);p=out/name/'cases/full2d_c0008/profile.in';p.parent.mkdir(parents=True)
        shutil.copyfile(old,gp/old.name);shutil.copyfile(db,gp/db.name)
        with h5py.File(gp/old.name,'r+') as h:
            del h['data'];h.create_dataset('data',data=volume,compression='gzip',compression_opts=4)
            h.attrs['shape_nxyz']=volume.shape
        content=card
        if name=='bottom40':
            content=content.replace('#domain: 110 40 inf','#domain: 110 80 inf').replace('#hertzian_dipole: z 78.6 35 ', '#hertzian_dipole: z 78.6 75 ').replace('#rx: 79.9 35 ','#rx: 79.9 75 ')
        p.write_text(content,encoding='utf-8')
        if name=='no_interbed':np.testing.assert_array_equal(volume[data!=3],data[data!=3]);assert np.all(volume[data==3]==2)
        if name=='homogeneous_cover':np.testing.assert_array_equal(volume==0,data==0)
        groups.append(dict(id=name,input=str(p.relative_to(out)),input_sha256=pair.sha(p),geometry=str((gp/old.name).relative_to(out)),geometry_sha256=pair.sha(gp/old.name),
            material=str((gp/db.name).relative_to(out)),material_sha256=pair.sha(db),native_shape=list(volume.shape),domain_m=[110,80 if name=='bottom40' else 40,.025],tx_m=[78.6,75. if name=='bottom40' else 35.,0.],rx_m=[79.9,75. if name=='bottom40' else 35.,0.],
            description=label,changed_voxel_count=int(np.count_nonzero(volume!=data)) if name!='bottom40' else 0,added_lower_m=40 if name=='bottom40' else 0))
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        p=out/g['input'];get_user_objects(p.read_text('utf-8').splitlines(),input_dir=p.parent)
    manifest=dict(status='PREPARED_APPROVED_GEOLOGICAL_CONTRASTS_NOT_RUN',approval_basis='User persistent objective: 研究明白; autonomous research and SSH simulation authorization.',
        parent_manifest_sha256=pair.sha(source/'manifest.json'),parent_H0_native_sha256=pair.sha(reference),groups=groups,
        case_labels=['原H0：保留上覆砂岩夹层']+[c[2] for c in cases],reference_is_prior_H0=True,
        invariants='2.5cm,800ns,HORIPML80cells,FP64,relative source/rx/current/Yee timing/material database. First4cases110x40m change declared geological voxels; bottom40 adds40m lower material and translates original interior/source/rx40m up. No snapshots.',
        factors='Remove upper sandstone interbed; remove all subsurface contrast while preserving surface; remove ground entirely; flatten layer stack using midpoint column3170; independently extend bottom40m.',
        predefined_gates_ns={'bottom':[332.31137724550894,356.31137724550894],'upper':[217.54091816367264,267.32268795741845]},
        snapshot_count=0,limits='Causal responses to declared material/geometry contrasts, including changed PML material continuation and interactions; not pure isolated reflections or measured clean truth. Flat case changes whole-domain geological shape, not only slope at source.')
    pair.save(out/'manifest.json',manifest);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',manifest)


def freeze(package,out,parent):
    boundary.freeze(package,out,parent)
    path=out/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    c['code_identities'][str(Path(__file__).resolve())]=pair.sha(__file__)
    pair.save(path,c);pair.save(out/'preflight_verification.json',boundary.audit(path))
    print(json.dumps(dict(final_contract_sha256=pair.sha(path),groups=len(c['groups']))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for key in ['source','reference','out','evidence','package','parent']:p.add_argument('--'+key,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.source,a.reference,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':boundary.run(a.out.resolve())
    else:print(json.dumps(boundary.audit(a.out/'execution_contract.json',True)))
