"""Bounded nonflat material sensitivity: reuse completed H0, solve only three new cases."""
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
    old = json.loads((a.parent/'manifest.json').read_text('utf-8'))
    pair = json.loads((a.pair/'manifest.json').read_text('utf-8'))
    h0 = next(g for g in old['groups'] if g['id']=='nonflat_ricker_H0')
    h1 = next(g for g in pair['groups'] if g['id']=='base_H1')
    ver = json.loads((a.reference/'completed_verification.json').read_text('utf-8'))
    contract = json.loads((a.reference/'execution_contract.json').read_text('utf-8'))
    assert ver['completed'] and ver['contract_sha256']==sha(a.reference/'execution_contract.json')
    record = next(g for g in ver['groups'] if g['id']==h0['id'])
    reused = a.reference/(h0['id']+'.h5')
    assert sha(reused)==record['native_sha256']
    reference_group = next(g for g in contract['groups'] if g['id']==h0['id'])
    for key in ['input','geometry','material']:
        assert sha(a.parent/h0[key])==h0[key+'_sha256']==reference_group[key+'_sha256']
        assert sha(a.pair/h1[key])==h1[key+'_sha256']
    assert h0['material_sha256']==h1['material_sha256']
    with h5py.File(a.parent/h0['geometry']) as b, h5py.File(a.pair/h1['geometry']) as t:
        x=b['data'][:];y=t['data'][:]
        changed=x!=y
        assert np.any(changed) and np.all(x[changed]==2) and np.all(y[changed]==3)
        assert x.shape==y.shape==(8400,1700,1)
    a.out.mkdir(parents=True)
    shutil.copyfile(reused,a.out/'reused_nonflat_H0.h5')
    original=(a.parent/h0['input']).read_text('utf-8')
    groups=[]
    for name, low, source, package in [('baseline_H1',False,h1,a.pair),('low_H0',True,h0,a.parent),('low_H1',True,h1,a.pair)]:
        gp=a.out/name/'geometries';gp.mkdir(parents=True)
        geometry=gp/Path(source['geometry']).name
        shutil.copyfile(package/source['geometry'],geometry)
        material=gp/Path(source['material']).name
        db=json.loads((package/source['material']).read_text('utf-8'))
        if low:
            mud=db['materials']['material_002_mudstone']
            assert mud['base']['electric_conductivity_s_per_m']==.003
            mud['base']['electric_conductivity_s_per_m']=.0003
            mud['metadata']['parent_parameter_design_sha256']=mud['metadata'].pop('parameter_design_sha256')
            mud['metadata']['parameter_status']='DC_CONDUCTIVITY_COUNTERFACTUAL_NOT_SITE_CALIBRATION'
            db['database']['description']='Nonflat diagnostic only: mudstone sigmaDC .003 -> .0003 S/m'
            material.write_text(json.dumps(db,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        else:shutil.copyfile(package/source['material'],material)
        p=a.out/name/'cases/full2d_r0098/profile.in';p.parent.mkdir(parents=True)
        lines=original.splitlines();lines[0]='#title: '+name+'; nonflat v5 mudstone DC diagnostic; no snapshots'
        p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        g=copy.deepcopy(h0)
        g.update(id=name,input=p.relative_to(a.out).as_posix(),geometry=geometry.relative_to(a.out).as_posix(),material=material.relative_to(a.out).as_posix())
        for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
        groups.append(g)
    m=copy.deepcopy(old)
    m.update(groups=groups,generator_sha256=sha(__file__),parent_manifest_sha256=sha(a.parent/'manifest.json'),pair_manifest_sha256=sha(a.pair/'manifest.json'),
             reused_H0=dict(native_sha256=sha(reused),verification_sha256=sha(a.reference/'completed_verification.json'),contract_sha256=sha(a.reference/'execution_contract.json'),group=reference_group),
             approval_basis='User continues authorized autonomous research after two-host V4.0.1 setup; three bounded nonflat DC-loss controls, reusing completed H0.',
             diagnostic='Only mudstone DC conductivity differs between baseline and low pairs; same original nonflat geometry, poles, source, heights, domain and numerical parameters.',
             limits='Single nonflat v5 station; not whole-line, finite-3D or field certification. Low sigma is counterfactual, never a selected/calibrated production material.')
    (a.out/'manifest.json').write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(new_groups=3,reused_native_sha256=sha(reused))))


def freeze(a):
    assert runner.gprMax.__version__=='4.0.1'
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert sha(a.package/'reused_nonflat_H0.h5')==m['reused_H0']['native_sha256']
    c['study_manifest']=m;c['approval_basis']=m['approval_basis']
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    c['file_identities'][str((a.package/'manifest.json').resolve())]=sha(a.package/'manifest.json')
    c['file_identities'][str((a.package/'reused_nonflat_H0.h5').resolve())]=m['reused_H0']['native_sha256']
    runner.save(path,c);runner.save(a.out/'preflight_verification.json',runner.audit(path))
    print(json.dumps(dict(status='FROZEN_THREE_NEW_REUSE_ONE',contract_sha256=sha(path))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['prepare','freeze'])
    for name in ['parent','pair','reference','out','package']:p.add_argument('--'+name,type=Path,required=name=='out')
    a=p.parse_args();a.out=a.out.resolve()
    prepare(a) if a.action=='prepare' else freeze(a)
