"""Dense Ricker wavefield pair restricted to the confirmed interbed/receiver region."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import line9_basal_pair as pair
import line9_boundary_controls as boundary
import line9_source_station_controls as source_controls
import hs4_station_grid_controls as supervisor


def prepare(source,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh package required')
    old=json.loads((source/'manifest.json').read_text('utf-8'));parent=next(g for g in old['groups'] if g['id']=='ricker_800')
    n=parent['expected_samples'];timeline=list(range(0,n,17));origin=[40.,10.,0.];shape=[450,270,1];spacing=[.1,.1,.025]
    suffix=''.join(f'#snapshot: 40 10 0 85 37 0.025 0.1 0.1 0.025 {j} snap{j:05d}.h5\n' for j in timeline)
    base=source/parent['input'];card=base.read_text('utf-8');assert '#snapshot:' not in card
    out.mkdir(parents=True);shutil.copyfile(source/'reference_H0.h5',out/'reference_H0.h5');groups=[]
    for name in ['H0','far_removed']:
        gp=out/name/'geometries';gp.mkdir(parents=True);p=out/name/'cases/full2d_c0008/profile.in';p.parent.mkdir(parents=True)
        geo=source/parent['geometry'];db=source/parent['material'];shutil.copyfile(geo,gp/geo.name);shutil.copyfile(db,gp/db.name)
        changed=0
        if name=='far_removed':
            with h5py.File(gp/geo.name,'r+') as h:
                data=h['data'][:];x=np.arange(data.shape[0])*.025;mask=(data==3)&((x>=45)&(x<60))[:,None,None];changed=int(mask.sum());data[mask]=2;h['data'][:]=data
        p.write_text(card+suffix,encoding='utf-8');g=dict(parent,id=name,input=p.relative_to(out).as_posix(),input_sha256=pair.sha(p),geometry=(gp/geo.name).relative_to(out).as_posix(),geometry_sha256=pair.sha(gp/geo.name),material=(gp/db.name).relative_to(out).as_posix(),material_sha256=pair.sha(db),changed_voxels=changed)
        groups.append(g)
    assert pair.sha(out/groups[0]['input'])==pair.sha(out/groups[1]['input'])
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        p=out/g['input'];get_user_objects(p.read_text('utf-8').splitlines(),input_dir=p.parent)
    history=int(np.prod(shape))*len(timeline)*6*8
    manifest=dict(status='PREPARED_APPROVED_NOT_RUN',approval_basis='User objective 研究明白 and standing autonomous wavefield/SSH simulation authorization.',
        parent_manifest_sha256=pair.sha(source/'manifest.json'),parent_H0_native_sha256=old['parent_H0_native_sha256'],groups=groups,
        snapshot_iterations=timeline,snapshot_count=len(timeline),snapshot_shape=shape,snapshot_origin_m=origin,snapshot_spacing_m=spacing,snapshot_six_field_history_bytes=history,
        fields=pair.FIELDS,invariants='Same c0008/H0,Ricker100MHz40A,800ns,110x40m,2.5cm,FP64,HORIPML80cells,material/source. Only upper sandstoneID3 insideROI45<=x<60 replaced by mudstoneID2; exact same observers in both.',
        limits='Native Ricker movie is pulse-weighted, not a501-tone SFCW movie. Spatial0.1m and time17native steps(~1ns) are visualization sampling, not FDTD coarsening. Window40-85m/10-37m excludes exterior wavefield. ROI deletion adds edges/changes interactions; no claim of one isolated diffraction or full-line/3D/field certification.')
    pair.save(out/'manifest.json',manifest);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',manifest)
    print(json.dumps(dict(snapshot_count=len(timeline),history_GiB=history/2**30)))


def audit(path,completed=False):
    result=source_controls.audit(path,completed);c=json.loads(path.read_text('utf-8'));m=c['study_manifest']
    if completed:
        for g,row in zip(c['groups'],result['groups']):
            files=sorted(Path(g['input']).parent.glob('profile_snaps/snap*.h5'));assert len(files)==m['snapshot_count'];rows=[]
            for file,j in zip(files,m['snapshot_iterations']):
                with h5py.File(file) as h:
                    assert h.attrs['iteration']==j and abs(h.attrs['time']-j*row['dt_s'])<1e-20
                    np.testing.assert_array_equal(h.attrs['origin'],m['snapshot_origin_m']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],m['snapshot_spacing_m'])
                    for key in m['fields']:
                        values=h[key][:];assert values.dtype==np.float64 and list(values.shape)==m['snapshot_shape'] and np.isfinite(values).all()
                rows.append(dict(file=str(file),iteration=j,sha256=pair.sha(file)))
            row.update(snapshot_count=len(files),snapshots=rows)
    return result


def freeze(package,out,parent):
    hardware=pair.resources();m=json.loads((package/'manifest.json').read_text('utf-8'));estimate=m['snapshot_six_field_history_bytes']+6*2**30
    if hardware['free_VRAM_bytes']<estimate:raise RuntimeError('Snapshot history plus6GiB native reserve does not fit')
    boundary.audit=audit;boundary.freeze(package,out,parent)
    path=out/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    c['code_identities'].update({str(Path(__file__).resolve()):pair.sha(__file__),str(Path(source_controls.__file__).resolve()):pair.sha(source_controls.__file__)})
    c.update(conservative_device_estimate_bytes=estimate,snapshot_budget_basis='Six-field observer history plus6GiB native/dispersive reserve, no field-kernel changes/streaming.')
    pair.save(path,c);pair.save(out/'preflight_verification.json',audit(path));print(json.dumps(dict(final_contract_sha256=pair.sha(path),history_GiB=m['snapshot_six_field_history_bytes']/2**30)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for key in ['source','out','evidence','package','parent']:p.add_argument('--'+key,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.source,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':supervisor.audit=audit;supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
