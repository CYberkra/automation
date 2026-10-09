"""Two bounded190m H0 interface controls, independent of the stopped line batch."""
import argparse,copy,json,shutil
from pathlib import Path
import h5py,numpy as np,psutil
import line9_v401_version_controls as runner
from hs_capsule_identity import sha256 as sha

read=lambda p:json.loads(p.read_text('utf-8'))


def prepare(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    for folder in ['baseline','parent']:shutil.copytree(a.baseline_package/folder,a.out/folder)
    dep=a.out/'right_dependency';dep.mkdir()
    for name in ['execution_contract.json','completed_verification.json','analysis.json','independent_audit.json','peak_diagnostic.json']:
        shutil.copyfile(a.right_result/name,dep/name)
    b=read(a.out/'parent/execution_contract.json')['study_manifest']['groups'][0]
    with h5py.File(a.out/'baseline/geometry.h5') as h:data=h['data'][:]
    assert data.shape==(8400,1700,1)
    groups=[]
    for ident in ['no_cover_contrast_H0','cover_bottom_down050_H0']:
        gp=a.out/ident/'geometries';gp.mkdir(parents=True)
        geom=gp/'full2d_compact.h5';after=data.copy()
        if ident=='no_cover_contrast_H0':after[data==2]=1
        else:
            for i in range(8400):
                mud=np.flatnonzero(data[i,:,0]==2)
                assert len(mud)>100 and np.array_equal(mud,np.arange(mud[-1]+1))
                after[i,mud[-20:],0]=1
        with h5py.File(a.out/'baseline/geometry.h5') as src,h5py.File(geom,'x') as dst:
            for k,val in src.attrs.items():dst.attrs[k]=val
            for k in src:
                if k!='data':src.copy(k,dst)
            ds=dst.create_dataset('data',data=after,compression='gzip',compression_opts=4)
            for k,val in src['data'].attrs.items():ds.attrs[k]=val
            dst.attrs['InterfaceControl']=ident
        mat=gp/'line9_research_materials_v1_smoothed.json';shutil.copyfile(a.out/'baseline/materials.json',mat)
        inp=a.out/ident/'cases'/ident/'profile.in';inp.parent.mkdir(parents=True)
        shutil.copyfile(a.out/'baseline/profile.in',inp)
        g=copy.deepcopy(b);g.update(id=ident,factor=ident,input=inp.relative_to(a.out).as_posix(),geometry=geom.relative_to(a.out).as_posix(),material=mat.relative_to(a.out).as_posix(),native_shape=[8400,1700,1],domain_height_m=42.5,changed_voxels=int(np.count_nonzero(after!=data)))
        for k in ['input','geometry','material']:g[k+'_sha256']=sha(a.out/g[k])
        groups.append(g)
    m=dict(status='PREPARED_TWO_COVER_INTERFACE_CONTROLS',groups=groups,no_retry=True,dt_s=read(a.out/'parent/execution_contract.json')['study_manifest']['dt_s'],
        generator_sha256=sha(__file__),approval_basis='User active goal 持续推进直到研究明白 and 做吧; two new bounded single190m high-loss H0 interface controls, no line resumption or retry.',
        baseline_sha256={p.name:sha(p) for p in (a.out/'baseline').iterdir()},parent_sha256={p.name:sha(p) for p in (a.out/'parent').iterdir()},right_dependency_sha256={p.name:sha(p) for p in dep.iterdir()},
        gates_ns=dict(early=[0,120],wide=[300,450],original_peak=[318.1729873586161,342.1729873586161],basal=[345.6180971390552,369.6180971390552],late=[450,1100]),
        limits='No-contrast replacement changes bulk mud propagation and bottom PML medium; thickness shift tests interface participation/path length, not uniquely reflection count. H0 controls only, no new H1, snapshots, whole line or field calibration.')
    runner.save(a.out/'manifest.json',m)
    print(json.dumps(dict(status=m['status'],cases=2,changes=[g['changed_voxels'] for g in groups])))


def freeze(a):
    from audit_line9_cover_interface_inputs import check
    m=read(a.package/'manifest.json');assert read(a.package/'independent_input_audit.json')==check(a.package)
    assert runner.gprMax.__version__=='4.0.1'
    old=read(a.package/'right_dependency/execution_contract.json');native=Path(runner.gprMax.__file__).parent
    for name,h in old['source_identities'].items():assert sha(native/name)==h
    roots=['v401_dense_loss_5m_r1','v401_dense_continuation_r1','v401_cover_wavefield_r4','v401_single_wavefield_r1','v401_high_h0_top20_r1','v401_high_h0_right40_r1']
    for proc in psutil.process_iter(['name','cmdline']):
        try:
            if 'python' in (proc.info['name'] or '').lower():assert not any(x in ' '.join(proc.info['cmdline'] or []).lower() for x in roots)
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    p=a.out/'execution_contract.json';c=read(p)
    c.update(study_manifest=m,approval_basis=m['approval_basis'],max_batch_wall_s=3900,gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock')
    for name in ['line9_cover_interface_controls.py','audit_line9_cover_interface_inputs.py']:c['code_identities'][str(Path(__file__).parent/name)]=sha(Path(__file__).parent/name)
    for f in a.package.rglob('*'):
        if f.is_file():c['file_identities'][str(f.resolve())]=sha(f)
    assert c['max_runs']==2 and c['no_retry'];runner.save(p,c);runner.save(a.out/'preflight_verification.json',runner.audit(p))
    print(json.dumps(dict(status='FROZEN_TWO_COVER_CONTROLS',contract_sha256=sha(p))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze']);p.add_argument('--out',type=Path,required=True)
    for k in ['baseline-package','right-result','package']:p.add_argument('--'+k,type=Path)
    a=p.parse_args();a.out=a.out.resolve();prepare(a) if a.action=='prepare' else freeze(a)
