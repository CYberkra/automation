"""One190m high-loss H0 right-edge extension; all prior attempts stay stopped."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import psutil

from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner


def read(p):
    return json.loads(p.read_text('utf-8'))


def prepare(a):
    assert not a.out.exists()
    c=read(a.parent/'execution_contract.json');v=read(a.parent/'completed_verification.json')
    x=read(a.parent/'analysis.json')
    assert v['completed'] and v['contract_sha256']==x['contract_sha256']==sha(a.parent/'execution_contract.json')
    assert x['verification_sha256']==sha(a.parent/'completed_verification.json')
    assert x['observer_native_and_spectrum_bitwise_equal'] and c['max_runs']==1
    old=c['study_manifest']['groups'][0]
    assert old['chainage_m']==190 and old['role']=='H0' and old['variant']=='high'
    assert sha(a.baseline/'native_H0.h5')==v['groups'][0]['native_sha256']==x['new_native_sha256']
    a.out.mkdir(parents=True);base=a.out/'baseline';base.mkdir();parent=a.out/'parent';parent.mkdir()
    for name in ['profile.in','geometry.h5','materials.json','native_H0.h5']:
        shutil.copyfile(a.baseline/name,base/name)
    for name in ['execution_contract.json','completed_verification.json','analysis.json']:
        shutil.copyfile(a.parent/name,parent/name)
    shutil.copyfile(a.observer_input,parent/'observer_profile.in')
    group='right40_H0';gp=a.out/group/'geometries';gp.mkdir(parents=True)
    geom=gp/'full2d_compact.h5'
    with h5py.File(base/'geometry.h5') as src,h5py.File(geom,'x') as dst:
        olddata=src['data'][:];assert olddata.shape==(8400,1700,1)
        after=np.pad(olddata,((0,1600),(0,0),(0,0)),mode='edge')
        for key,value in src.attrs.items():dst.attrs[key]=value
        for key in src:
            if key!='data':src.copy(key,dst)
        ds=dst.create_dataset('data',data=after,compression='gzip',compression_opts=4)
        for key,value in src['data'].attrs.items():ds.attrs[key]=value
        dst.attrs['shape_nxyz']=list(after.shape)
        dst.attrs['BoundaryControl']='Right+40m edge-column continuation; original8400x1700 interior exact; no coordinate translation'
    mat=gp/'line9_research_materials_v1_smoothed.json';shutil.copyfile(base/'materials.json',mat)
    lines=(base/'profile.in').read_text('utf-8').splitlines()
    assert sum(line.startswith('#domain:') for line in lines)==1
    lines=['#domain: 250 42.5 inf' if line.startswith('#domain:') else line for line in lines]
    inp=a.out/group/'cases'/group/'profile.in';inp.parent.mkdir(parents=True)
    inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    g=copy.deepcopy(old);g.update(id=group,factor='right40',input=inp.relative_to(a.out).as_posix(),
        geometry=geom.relative_to(a.out).as_posix(),material=mat.relative_to(a.out).as_posix(),
        native_shape=[10000,1700,1],domain_height_m=42.5)
    for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
    m=dict(status='PREPARED_APPROVED_ONE_RIGHT40_H0',groups=[g],dt_s=c['study_manifest']['dt_s'],
        no_retry=True,baseline_sha256={p.name:sha(p) for p in base.iterdir()},
        parent_sha256={p.name:sha(p) for p in parent.iterdir()},generator_sha256=sha(__file__),
        approval_basis='User active goal: 持续推进直到研究明白; next one190m high-loss H0 right-edge+40m A-scan control; prior stopped line/four-case pipeline not resumed; no new snapshots, H1 or retries.',
        gates_ns=dict(early=[0,120],wide=[300,450],basal=[345.6180971390552,369.6180971390552],late=[450,1100]),
        limits='One H0 finite right-extension sensitivity, including added exterior geology; not a PML-only reflection coefficient; no matched new H1/basal difference, full boundary convergence or field calibration.')
    runner.save(a.out/'manifest.json',m)
    print(json.dumps(dict(status=m['status'],cases=1,native_shape=g['native_shape'])))


def freeze(a):
    from audit_line9_high_h0_right40_inputs import check
    assert runner.gprMax.__version__=='4.0.1'
    m=read(a.package/'manifest.json');audit=check(a.package)
    frozen_audit=read(a.package/'independent_input_audit.json')
    assert frozen_audit==audit
    previous=read(a.package/'parent/execution_contract.json')
    native=Path(runner.gprMax.__file__).parent
    for path,h in previous['source_identities'].items():assert sha(native/path)==h
    roots=['v401_dense_loss_5m_r1','v401_dense_continuation_r1','v401_cover_wavefield_r4','v401_single_wavefield_r1','v401_high_h0_top20_r1']
    for proc in psutil.process_iter(['name','cmdline']):
        try:
            cmd=' '.join(proc.info['cmdline'] or []).lower()
            if 'python' in (proc.info['name'] or '').lower():assert not any(r in cmd for r in roots),'Prior owned task still active'
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    p=a.out/'execution_contract.json';c=read(p)
    c.update(study_manifest=m,approval_basis=m['approval_basis'],max_batch_wall_s=2100,
        gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock')
    for name in ['line9_high_h0_right40.py','audit_line9_high_h0_right40_inputs.py']:
        path=Path(__file__).parent/name;c['code_identities'][str(path.resolve())]=sha(path)
    for path in a.package.rglob('*'):
        if path.is_file():c['file_identities'][str(path.resolve())]=sha(path)
    runner.save(p,c);runner.save(a.out/'preflight_verification.json',runner.audit(p))
    print(json.dumps(dict(status='FROZEN_ONE_RIGHT40_H0',contract_sha256=sha(p))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['prepare','freeze']);p.add_argument('--out',type=Path,required=True)
    for key in ['parent','baseline','observer-input','package']:p.add_argument('--'+key,type=Path)
    a=p.parse_args();a.out=a.out.resolve()
    prepare(a) if a.action=='prepare' else freeze(a)
