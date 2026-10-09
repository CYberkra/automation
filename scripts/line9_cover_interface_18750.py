"""Two187.5m H0 controls replicate completed190m interface experiments."""
import argparse,copy,json,shutil
from pathlib import Path
import h5py,numpy as np,psutil
import line9_v401_version_controls as runner
from hs_capsule_identity import sha256 as sha
read=lambda p:json.loads(p.read_text('utf-8'))


def prepare(a):
    assert not a.out.exists();a.out.mkdir(parents=True)
    baseline=a.out/'baseline';baseline.mkdir();dep=a.out/'dependency190';dep.mkdir();parent=a.out/'parent';parent.mkdir()
    original=read(a.original_package/'manifest.json');g=next(x for x in original['groups'] if x['id']=='high_x18750_H0')
    for key,name in [('input','profile.in'),('geometry','geometry.h5'),('material','materials.json')]:
        assert sha(a.original_package/g[key])==g[key+'_sha256'];shutil.copyfile(a.original_package/g[key],baseline/name)
    shutil.copyfile(a.original_source/(g['id']+'.h5'),baseline/'native_H0.h5')
    for name in ['execution_contract.json','execution.jsonl','stopped_certificate.json']:shutil.copyfile(a.original_source/name,parent/name)
    for name in ['analysis.json','independent_audit.json','stopped_export_audit.json']:shutil.copyfile(a.original_result/name,parent/name)
    for name in ['execution_contract.json','completed_verification.json','analysis.json','independent_audit.json']:shutil.copyfile(a.result190/name,dep/name)
    old=read(a.package190/'manifest.json');groups=[]
    assert g['chainage_m']==187.5 and g['geometry_sha256']==old['baseline_sha256']['geometry.h5'] and g['material_sha256']==old['baseline_sha256']['materials.json']
    for og in old['groups']:
        gp=a.out/og['id']/'geometries';gp.mkdir(parents=True)
        geom=gp/'full2d_compact.h5';mat=gp/'line9_research_materials_v1_smoothed.json'
        shutil.copyfile(a.package190/og['geometry'],geom);shutil.copyfile(baseline/'materials.json',mat)
        inp=a.out/og['id']/'cases'/og['id']/'profile.in';inp.parent.mkdir(parents=True);shutil.copyfile(baseline/'profile.in',inp)
        row=copy.deepcopy(g);row.update(id=og['id'],factor=og['factor'],input=inp.relative_to(a.out).as_posix(),geometry=geom.relative_to(a.out).as_posix(),material=mat.relative_to(a.out).as_posix(),changed_voxels=og['changed_voxels'])
        for key in ['input','geometry','material']:row[key+'_sha256']=sha(a.out/row[key])
        assert row['geometry_sha256']==og['geometry_sha256'];groups.append(row)
    # Old baseline only, before either new solve; fixed501/Hann peak defines its replication gate.
    from analyze_line9_v401_version_controls import response,inverse
    z,error=response(baseline/'native_H0.h5',.025);w=np.hanning(501);w/=w.mean();x,t=inverse(z[:,None],w)
    k=(t*1e9>=300)&(t*1e9<=450);peak=float(t[k][abs(x[k,0]).argmax()]*1e9)
    m=dict(status='PREPARED_TWO18750_INTERFACE_REPLICATION_CONTROLS',groups=groups,original_group=g,no_retry=True,dt_s=g['dt_s'],generator_sha256=sha(__file__),
        approval_basis='User继续 after190m result: next planned187.5m two new interface-control A-scans; no stopped line resumption or retry.',
        gates_ns=dict(early=[0,120],wide=[300,450],original_peak=[peak-12,peak+12],basal=g['basal_gate_ns'],late=[450,1100]),
        baseline_peak_reference=dict(hann_peak_ns=peak,native_sha256=sha(baseline/'native_H0.h5'),direct_DFT_relative_L2=error),
        first_difference_replication_gate_ns=[120,250],limits='Two187.5m controls replicate interface dependence; same global geometries as190m. Replacement changes bulk propagation/bottom PML medium; differences are not pure echoes or clean labels; no full line,3D or field validation.')
    for folder in ['baseline','parent','dependency190']:m[folder+'_sha256']={p.name:sha(p) for p in (a.out/folder).iterdir()}
    runner.save(a.out/'manifest.json',m);print(json.dumps(dict(status=m['status'],cases=2,original_peak_ns=peak)))


def freeze(a):
    from audit_line9_cover_interface_18750_inputs import check
    m=read(a.package/'manifest.json');assert read(a.package/'independent_input_audit.json')==check(a.package)
    assert runner.gprMax.__version__=='4.0.1'
    old=read(a.package/'dependency190/execution_contract.json');native=Path(runner.gprMax.__file__).parent
    for name,h in old['source_identities'].items():assert sha(native/name)==h
    roots=['v401_dense_loss_5m_r1','v401_dense_continuation_r1','v401_cover_wavefield_r4','v401_single_wavefield_r1','v401_high_h0_top20_r1','v401_high_h0_right40_r1','v401_cover_interface_r1']
    for proc in psutil.process_iter(['name','cmdline']):
        try:
            if 'python' in (proc.info['name'] or '').lower():assert not any(x in ' '.join(proc.info['cmdline'] or []).lower() for x in roots)
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package));p=a.out/'execution_contract.json';c=read(p)
    c.update(study_manifest=m,approval_basis=m['approval_basis'],max_batch_wall_s=3900,gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock')
    for name in ['line9_cover_interface_18750.py','audit_line9_cover_interface_18750_inputs.py']:c['code_identities'][str(Path(__file__).parent/name)]=sha(Path(__file__).parent/name)
    for f in a.package.rglob('*'):
        if f.is_file():c['file_identities'][str(f.resolve())]=sha(f)
    assert c['max_runs']==2 and c['no_retry'];runner.save(p,c);runner.save(a.out/'preflight_verification.json',runner.audit(p));print(json.dumps(dict(status='FROZEN_TWO18750_COVER_CONTROLS',contract_sha256=sha(p))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze']);p.add_argument('--out',type=Path,required=True)
    for k in ['original-package','original-source','original-result','package190','result190','package']:p.add_argument('--'+k,type=Path)
    a=p.parse_args();a.out=a.out.resolve();prepare(a) if a.action=='prepare' else freeze(a)
