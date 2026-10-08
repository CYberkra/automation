"""Six bounded original-v5 low-loss station solves; reuse the completed anchor."""
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
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(a.parent/'manifest.json');d=read(a.design);v=read(a.reference/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(a.reference/'execution_contract.json')
    assert d['parent_manifest_sha256']==sha(a.parent/'manifest.json') and d['new_solves_started']==0
    assert [s['chainage_m'] for s in d['stations']]==[190.,180.,146.4]
    groups={g['id']:g for g in m['groups']};records={g['id']:g for g in v['groups']}
    a.out.mkdir(parents=True);new=[];reused=[]
    for role in ['H0','H1']:
        base=groups['span03_dc0003_'+role]
        for key in ['input','geometry','material']:assert sha(a.parent/base[key])==base[key+'_sha256']
        raw=a.reference/(base['id']+'.h5');assert sha(raw)==records[base['id']]['native_sha256']
        dst=a.out/'reused'/('anchor_'+role+'.h5');dst.parent.mkdir(exist_ok=True);shutil.copyfile(raw,dst)
        reused.append(dict(id='anchor_'+role,chainage_m=198.6,path=dst.relative_to(a.out).as_posix(),native_sha256=sha(dst),parent_group=base))
        original=(a.parent/base['input']).read_text('utf-8').splitlines()
        for s in d['stations']:
            name=s['id']+'_'+role;gp=a.out/name/'geometries';gp.mkdir(parents=True)
            geometry=gp/'full2d_compact.h5';material=gp/Path(base['material']).name
            shutil.copyfile(a.parent/base['geometry'],geometry);shutil.copyfile(a.parent/base['material'],material)
            lines=[]
            for line in original:
                if line.startswith('#title:'):line=f'#title: {name}; X{s["chainage_m"]:g}; low-loss spatial control; no snapshots'
                elif line.startswith('#hertzian_dipole:'):line='#hertzian_dipole: z '+' '.join(f'{q:.12g}' for q in s['tx_m'])+' pulse'
                elif line.startswith('#rx:'):line='#rx: '+' '.join(f'{q:.12g}' for q in s['rx_m'])+' '+s['id']+'_rx1 Ez'
                lines.append(line)
            inp=a.out/name/'cases'/s['id']/'profile.in';inp.parent.mkdir(parents=True);inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            g=copy.deepcopy(base);g.update(id=name,station_id=s['id'],chainage_m=s['chainage_m'],
                input=inp.relative_to(a.out).as_posix(),geometry=geometry.relative_to(a.out).as_posix(),material=material.relative_to(a.out).as_posix(),
                tx_m=s['tx_m'][:2]+[0.],rx_m=s['rx_m'][:2]+[0.],input_tx_m=s['tx_m'],input_rx_m=s['rx_m'],
                basal_gate_ns=s['proposed_basal_gate_ns'],wide_gate_ns=s['proposed_wide_gate_ns'])
            for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
            new.append(g)
    new.sort(key=lambda g:(-g['chainage_m'],g['id']))
    out=copy.deepcopy(m);out.update(groups=new,reused=reused,station_design=d,station_design_sha256=sha(a.design),
        parent_manifest_sha256=sha(a.parent/'manifest.json'),reference_verification_sha256=sha(a.reference/'completed_verification.json'),
        reference_contract_sha256=sha(a.reference/'execution_contract.json'),generator_sha256=sha(__file__),
        approval_basis='User explicitly says do after three-station 190/180/146.4m plan; standing SSH simulation authorization. Six new single-attempt solves, existing 198.6m pair reused.',
        diagnostic='Four sparse physical stations including anchor; unchanged original full-domain geology and low-loss material, only Tx/Rx station coordinates change.',
        invariants='Original 210x42.5m/2.5cm/1200ns/dt/FP64/40A100MHz Ricker/HORIPML80/averaging and exact501-tone SFCW. Same H0 and H1 geometry bytes at every site. AGL8m grid-quantized original planned coordinates.',
        limits='Low-loss counterfactual, not calibrated field materials. Four sparse stations are not a dense line; H1-H0 includes bottom-replacement interactions. Geometry-guided gates/templates are not blind event truth. No new mesh/boundary/top convergence, finite3D or field validation; production material unchanged.')
    runner.save(a.out/'manifest.json',out);print(json.dumps(dict(status='PREPARED_SIX_NEW_TWO_REUSED',groups=[g['id'] for g in new])))


def freeze(a):
    assert runner.gprMax.__version__=='4.0.1'
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=json.loads(path.read_text('utf-8'));m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    c['study_manifest']=m;c['approval_basis']=m['approval_basis'];c['max_batch_wall_s']=3600
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    c['file_identities'][str((a.package/'manifest.json').resolve())]=sha(a.package/'manifest.json')
    for r in m['reused']:
        p=(a.package/r['path']).resolve();assert sha(p)==r['native_sha256'];c['file_identities'][str(p)]=sha(p)
    audit=a.package/'independent_input_audit.json'
    assert json.loads(audit.read_text('utf-8'))['manifest_sha256']==sha(a.package/'manifest.json')
    c['file_identities'][str(audit.resolve())]=sha(audit)
    runner.save(path,c);runner.save(a.out/'preflight_verification.json',runner.audit(path))
    print(json.dumps(dict(status='FROZEN_SIX_NEW_TWO_REUSED',contract_sha256=sha(path))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze'])
    for key in ['parent','reference','design','package','out']:p.add_argument('--'+key,type=Path,required=key=='out')
    a=p.parse_args();a.out=a.out.resolve();prepare(a) if a.action=='prepare' else freeze(a)
