"""Freeze a 2x2 mudstone polarization/DC diagnostic; reuse four, solve four."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import numpy as np
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner


def prepare(a):
    assert not a.out.exists()
    m=json.loads((a.parent/'manifest.json').read_text('utf-8'))
    v=json.loads((a.reference/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256']==sha(a.reference/'execution_contract.json')
    groups={g['id']:g for g in m['groups']}
    records={g['id']:g for g in v['groups']}
    a.out.mkdir(parents=True)
    reused=[]
    for name,source,record in [
        ('span3_dc003_H0',a.parent/'reused_nonflat_H0.h5',m['reused_H0']),
        ('span3_dc003_H1',a.reference/'baseline_H1.h5',records['baseline_H1']),
        ('span3_dc0003_H0',a.reference/'low_H0.h5',records['low_H0']),
        ('span3_dc0003_H1',a.reference/'low_H1.h5',records['low_H1'])]:
        assert sha(source)==record['native_sha256']
        dst=a.out/'reused'/(name+'.h5');dst.parent.mkdir(exist_ok=True)
        shutil.copyfile(source,dst)
        reused.append(dict(id=name,path=dst.relative_to(a.out).as_posix(),native_sha256=sha(dst)))
    new=[]
    for dc,label in [(.003,'003'),(.0003,'0003')]:
        for role in ['H0','H1']:
            base=groups['low_H0' if role=='H0' else 'baseline_H1']
            name=f'span03_dc{label}_{role}'
            target=a.out/name
            geometry=target/'geometries/full2d_compact.h5';geometry.parent.mkdir(parents=True)
            assert sha(a.parent/base['geometry'])==base['geometry_sha256']
            shutil.copyfile(a.parent/base['geometry'],geometry)
            db=json.loads((a.parent/groups['baseline_H1']['material']).read_text('utf-8'))
            mud=db['materials']['material_002_mudstone']
            pole=mud['poles'][0];tau=pole['relaxation_time_s']
            assert tau==8e-9 and len(mud['poles'])==1
            old_delta=pole['relative_permittivity_difference']
            assert abs(mud['base']['relative_permittivity']+old_delta/(1+(2*np.pi*95e6*tau)**2)-12)<1e-12
            pole['relative_permittivity_difference']=old_delta*.1
            mud['base']['relative_permittivity']=12-pole['relative_permittivity_difference']/(1+(2*np.pi*95e6*tau)**2)
            mud['base']['electric_conductivity_s_per_m']=dc
            mud['metadata']['parent_parameter_design_sha256']=mud['metadata'].pop('parameter_design_sha256')
            mud['metadata']['parameter_status']='POLARIZATION_DC_FACTORIAL_COUNTERFACTUAL_NOT_SITE_CALIBRATION'
            db['database']['description']='Diagnostic span0.3 at epsilonReal95=12,tau8ns; DC factorial, not field fit'
            material=target/'geometries/line9_research_materials_v1_smoothed.json'
            material.write_text(json.dumps(db,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
            original=(a.parent/base['input']).read_text('utf-8').splitlines()
            original[0]='#title: '+name+'; mudstone polarization/DC diagnostic; no snapshots'
            inp=target/'cases/full2d_r0098/profile.in';inp.parent.mkdir(parents=True)
            inp.write_text('\n'.join(original)+'\n',encoding='utf-8')
            g=copy.deepcopy(base);g.update(id=name,input=inp.relative_to(a.out).as_posix(),geometry=geometry.relative_to(a.out).as_posix(),material=material.relative_to(a.out).as_posix())
            for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
            new.append(g)
    out=copy.deepcopy(m)
    out.pop('reused_H0',None)
    out.update(groups=new,reused=reused,parent_manifest_sha256=sha(a.parent/'manifest.json'),
        reference_verification_sha256=sha(a.reference/'completed_verification.json'),generator_sha256=sha(__file__),
        approval_basis='Standing autonomous research authorization; bounded four missing controls after completed original nonflat DC pair. Freeze before solver, no retries.',
        diagnostic='2x2 real endpoint span {3,.3} x DC {.003,.0003}; epsilonReal95=12,tau8ns; same original v5 nonflat station and all other physics/numerics.',
        invariants='Original 210x42.5m,2.5cm,1200ns,FP64,AGL8m,source/rx,PML80,other materials and exact H1/H0 voxels; epsilonReal95=12,tau8ns. Span change necessarily changes off-centre real dispersion and polarization loss.',
        limits='Mechanism diagnostic, not site confidence interval or fitted material. H0 is geological response, delta only defined bottom replacement. One station; no whole-line, 3D or field validation. Formal materials unchanged.')
    (a.out/'manifest.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PREPARED_FOUR_NEW_FOUR_REUSED',groups=[g['id'] for g in new])))


def freeze(a):
    assert runner.gprMax.__version__=='4.0.1'
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=json.loads(path.read_text('utf-8'))
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    c['study_manifest']=m;c['approval_basis']=m['approval_basis']
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    c['file_identities'][str((a.package/'manifest.json').resolve())]=sha(a.package/'manifest.json')
    for r in m['reused']:
        p=(a.package/r['path']).resolve();assert sha(p)==r['native_sha256']
        c['file_identities'][str(p)]=r['native_sha256']
    runner.save(path,c);runner.save(a.out/'preflight_verification.json',runner.audit(path))
    print(json.dumps(dict(status='FROZEN_FOUR_NEW_FOUR_REUSED',contract_sha256=sha(path))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['prepare','freeze'])
    for name in ['parent','reference','out','package']:p.add_argument('--'+name,type=Path,required=name=='out')
    a=p.parse_args();a.out=a.out.resolve()
    prepare(a) if a.action=='prepare' else freeze(a)
