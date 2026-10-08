"""Frozen localized upper-interbed contrasts and same-domain low-AGL basal pair."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import line9_basal_pair as pair
import line9_boundary_controls as boundary


def prepare(source,reference,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh outputs required')
    m=json.loads((source/'manifest.json').read_text('utf-8'))
    db=source/'H0/geometries/line9_research_materials_v1_smoothed.json'
    old=source/'H0/geometries/full2d_compact.h5'
    assert pair.sha(old)==m['groups'][1]['geometry_sha256'] and pair.sha(db)==m['material_sha256']
    with h5py.File(old) as h:h0=h['data'][:]
    with h5py.File(source/'H1/geometries/full2d_compact.h5') as h:h1=h['data'][:]
    x=np.arange(h0.shape[0])*.025
    far=h0.copy();fm=(h0==3)&((x>=45)&(x<60))[:,None,None];far[fm]=2
    near=h0.copy();nm=(h0==3)&((x>=65)&(x<90))[:,None,None];near[nm]=2
    cases=[('far_interbed_removed',far,35.,'高航高H0：仅x45–60m夹层→泥岩'),
           ('near_interbed_removed',near,35.,'高航高H0：仅x65–90m夹层→泥岩'),
           ('low_H1',h1,25.725,'AGL1m：完整原地质H1'),('low_H0',h0,25.725,'AGL1m：仅去底砂H0')]
    card='\n'.join(l for l in (source/'H0/cases/full2d_c0008/profile.in').read_text('utf-8').splitlines() if not l.startswith('#snapshot:'))+'\n'
    out.mkdir(parents=True);shutil.copyfile(reference,out/'reference_H0.h5');groups=[]
    for name,volume,y,label in cases:
        gp=out/name/'geometries';gp.mkdir(parents=True);p=out/name/'cases/full2d_c0008/profile.in';p.parent.mkdir(parents=True)
        shutil.copyfile(old,gp/old.name);shutil.copyfile(db,gp/db.name)
        with h5py.File(gp/old.name,'r+') as h:h['data'][:]=volume
        content=card.replace('#hertzian_dipole: z 78.6 35 ',f'#hertzian_dipole: z 78.6 {y} ').replace('#rx: 79.9 35 ',f'#rx: 79.9 {y} ')
        p.write_text(content,encoding='utf-8')
        groups.append(dict(id=name,input=p.relative_to(out).as_posix(),input_sha256=pair.sha(p),geometry=(gp/old.name).relative_to(out).as_posix(),geometry_sha256=pair.sha(gp/old.name),
            material=(gp/db.name).relative_to(out).as_posix(),material_sha256=pair.sha(db),native_shape=list(volume.shape),domain_m=[110,40,.025],tx_m=[78.6,y,0.],rx_m=[79.9,y,0.],
            description=label,midpoint_agl_m=y-24.725,changed_voxels_vs_H0=int(np.count_nonzero(volume!=h0))))
    np.testing.assert_array_equal(far[~fm],h0[~fm]);np.testing.assert_array_equal(near[~nm],h0[~nm])
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        p=out/g['input'];get_user_objects(p.read_text('utf-8').splitlines(),input_dir=p.parent)
    manifest=dict(status='PREPARED_APPROVED_LOCALIZED_AND_HEIGHT_CONTROLS_NOT_RUN',approval_basis='User objective 研究明白 and explicit SSH simulation/autonomous research authorization.',
        parent_manifest_sha256=pair.sha(source/'manifest.json'),parent_H0_native_sha256=pair.sha(reference),groups=groups,reference_is_prior_H0=True,
        case_labels=['高航高原H0']+[v[3] for v in cases],snapshot_count=0,
        invariants='110x40m,2.5cm,800ns,FP64,HORIPML80cells,40A impulse,material database,relative horizontal Tx/Rx; no snapshots. Local edits affect sandstoneID3 only. Low pair same geometry/domain as original high pair, only y35->25.725m.',
        predefined_gates_ns={'high_basal':[332.31137724550894,356.31137724550894], 'low_basal_air_shift_proxy':[270.433,294.433]},
        limits='Far/near removals introduce ROI edges and alter interactions; they locate sensitivity regions, not isolated pinchout reflection. Low gate fixed before solve by approximate air delay; evaluate full complex traces and model-derived template as well, no adaptive production correction.')
    pair.save(out/'manifest.json',manifest);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',manifest)


def freeze(package,out,parent):
    boundary.freeze(package,out,parent)
    path=out/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    c['code_identities'][str(Path(__file__).resolve())]=pair.sha(__file__)
    pair.save(path,c);pair.save(out/'preflight_verification.json',boundary.audit(path))
    print(json.dumps(dict(final_contract_sha256=pair.sha(path),groups=len(c['groups']))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for name in ['source','reference','out','evidence','package','parent']:p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.source,a.reference,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':boundary.run(a.out.resolve())
    else:print(json.dumps(boundary.audit(a.out/'execution_contract.json',True)))
