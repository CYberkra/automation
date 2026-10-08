"""Bounded low-loss v5 side/bottom distance pairs; reuse completed baseline."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner


def prepare(a):
    assert not a.out.exists()
    parent=json.loads((a.parent/'manifest.json').read_text('utf-8'))
    v=json.loads((a.reference/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256']==sha(a.reference/'execution_contract.json')
    groups={g['id']:g for g in parent['groups']};records={r['id']:r for r in v['groups']}
    a.out.mkdir(parents=True);new=[];reused=[]
    for role in ['H0','H1']:
        base=groups['span03_dc0003_'+role]
        for key in ['input','geometry','material']:assert sha(a.parent/base[key])==base[key+'_sha256']
        raw=a.reference/(base['id']+'.h5');assert sha(raw)==records[base['id']]['native_sha256']
        dest=a.out/'reused'/('base_'+role+'.h5');dest.parent.mkdir(exist_ok=True);shutil.copyfile(raw,dest)
        reused.append(dict(id='base_'+role,path=dest.relative_to(a.out).as_posix(),native_sha256=sha(dest),parent_group=base))
        original=(a.parent/base['input']).read_text('utf-8').splitlines()
        with h5py.File(a.parent/base['geometry']) as h:data=h['data'][:]
        assert data.shape==(8400,1700,1)
        for factor,sx,sy in [('sides40',1600,0),('bottom20',0,800)]:
            name=factor+'_'+role;folder=a.out/name;gp=folder/'geometries';gp.mkdir(parents=True)
            outdata=np.pad(data,((sx,sx),(sy,0),(0,0)),mode='edge')
            np.testing.assert_array_equal(outdata[sx:sx+8400,sy:sy+1700],data)
            geometry=gp/'full2d_compact.h5'
            with h5py.File(a.parent/base['geometry']) as src,h5py.File(geometry,'x') as dst:
                for key,value in src.attrs.items():dst.attrs[key]=value
                for key in src:
                    if key!='data':src.copy(key,dst)
                ds=dst.create_dataset('data',data=outdata,compression='gzip',compression_opts=4)
                for key,value in src['data'].attrs.items():ds.attrs[key]=value
                dst.attrs['shape_nxyz']=outdata.shape
                dst.attrs['profile_x_offset_m']=float(src.attrs['profile_x_offset_m'])-sx*.025
                dst.attrs['elevation_offset_m']=float(src.attrs['elevation_offset_m'])-sy*.025
                dst.attrs['BoundaryControl']='Original interior exact; edge-column side continuation or basal-row bottom continuation; source translated with original interior'
            material=gp/Path(base['material']).name;shutil.copyfile(a.parent/base['material'],material)
            card=[]
            for line in original:
                words=line.split()
                if line.startswith('#title:'):line='#title: '+name+'; low-loss boundary sensitivity; no snapshots'
                elif line.startswith('#domain:'):line=f'#domain: {outdata.shape[0]*.025:g} {outdata.shape[1]*.025:g} inf'
                elif line.startswith('#hertzian_dipole:'):
                    words[2]=f'{float(words[2])+sx*.025:.12g}';words[3]=f'{float(words[3])+sy*.025:.12g}';line=' '.join(words)
                elif line.startswith('#rx:'):
                    words[1]=f'{float(words[1])+sx*.025:.12g}';words[2]=f'{float(words[2])+sy*.025:.12g}';line=' '.join(words)
                card.append(line)
            inp=folder/'cases/full2d_r0098/profile.in';inp.parent.mkdir(parents=True);inp.write_text('\n'.join(card)+'\n',encoding='utf-8')
            g=copy.deepcopy(base);g.update(id=name,input=inp.relative_to(a.out).as_posix(),geometry=geometry.relative_to(a.out).as_posix(),material=material.relative_to(a.out).as_posix(),native_shape=list(outdata.shape),translation_m=[sx*.025,sy*.025,0],domain_height_m=outdata.shape[1]*.025)
            for key in ['tx_m','rx_m']:g[key]=(np.array(base[key])+[sx*.025,sy*.025,0]).tolist()
            for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
            new.append(g)
    new.sort(key=lambda g:(0 if g['id'].startswith('sides') else 1,g['id']))
    m=copy.deepcopy(parent);m.pop('reused',None)
    m.update(groups=new,reused=reused,generator_sha256=sha(__file__),parent_manifest_sha256=sha(a.parent/'manifest.json'),
        reference_verification_sha256=sha(a.reference/'completed_verification.json'),reference_contract_sha256=sha(a.reference/'execution_contract.json'),
        approval_basis='Persistent user goal research fully; standing autonomous SSH simulation authorization. Four low-loss v5 boundary controls frozen before solving, reuse completed baseline pair.',
        diagnostic='Same span0.3/DC0.0003 low-loss spectrum; independent sides+40m and bottom+20m controls; no joint extension in this batch.',
        invariants='Exact original interior voxels,material JSON,2.5cm,1200ns,dt,FP64,40A/100MHz Ricker,AGL8m,Tx/Rx relative to geology,HORIPML80,averaging and exact501-tone SFCW.',
        limits='One low-loss counterfactual station. Extensions change outside continuation together with PML distance; finite sensitivity only, not PML-only mechanism or global convergence certificate. No joint control, whole-line,3D,field or production material validation.')
    (a.out/'manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PREPARED_FOUR_NEW_TWO_REUSED',groups=[g['id'] for g in new])))


def freeze(a):
    assert runner.gprMax.__version__=='4.0.1'
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=json.loads(path.read_text('utf-8'));m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    c['study_manifest']=m;c['approval_basis']=m['approval_basis']
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    c['file_identities'][str((a.package/'manifest.json').resolve())]=sha(a.package/'manifest.json')
    for r in m['reused']:
        p=(a.package/r['path']).resolve();assert sha(p)==r['native_sha256'];c['file_identities'][str(p)]=r['native_sha256']
    runner.save(path,c);runner.save(a.out/'preflight_verification.json',runner.audit(path))
    print(json.dumps(dict(status='FROZEN_FOUR_NEW_TWO_REUSED',contract_sha256=sha(path))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze'])
    for key in ['parent','reference','package','out']:p.add_argument('--'+key,type=Path,required=key=='out')
    a=p.parse_args();a.out=a.out.resolve();prepare(a) if a.action=='prepare' else freeze(a)
