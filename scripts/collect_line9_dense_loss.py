"""Export only receipt-completed, native-audited dense-line traces; no solver start."""
import argparse,json,shutil,time,zipfile
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner

def main(root):
    read=lambda p:json.loads(p.read_text('utf-8'))
    execution=root/'execution';c=read(execution/'execution_contract.json');m=c['study_manifest']
    assert runner.audit(execution/'execution_contract.json')==read(execution/'preflight_verification.json')
    # Ignore an in-progress final line instead of guessing its terminal status.
    text=(execution/'execution.jsonl').read_text('utf-8');lines=text[:text.rfind('\n')+1].splitlines();events=[json.loads(s) for s in lines]
    assert events[0]['contract_sha256']==sha(execution/'execution_contract.json')
    receipts={e['group']:e for e in events if e['status']=='COMPLETED' and 'group' in e}
    assert len(receipts)==sum(e['status']=='COMPLETED' and 'group' in e for e in events)
    out=root/'exports'/str(time.time_ns());out.mkdir(parents=True);records=[]
    for g in c['groups']:
        sid=g['id']
        if sid not in receipts:continue
        path=Path(g['input']).with_suffix('.h5');digest=sha(path);assert digest==receipts[sid]['raw_sha256']
        with h5py.File(path) as h:
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==g['expected_samples']==20352
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
            for node,role in [('srcs/src1','tx'),('rxs/rx1','rx')]:np.testing.assert_allclose(h[node].attrs['Position'],g[role+'_m'],rtol=0,atol=1e-11)
            for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
                x=h[key][:];assert x.dtype==np.float64 and x.shape==(20352,) and np.isfinite(x).all() and np.any(x)
        shutil.copyfile(path,out/(sid+'.h5'))
        for name in ['stdout.log','stderr.log']:shutil.copyfile(path.parent/name,out/(sid+'_'+name))
        records.append(dict(id=sid,native_sha256=digest,elapsed_s=receipts[sid]['elapsed_s']))
    for r in m['reused']:
        p=root/'package'/r['path'];assert sha(p)==r['native_sha256'];shutil.copyfile(p,out/(r['id']+'.h5'))
    for name in ['execution_contract.json','preflight_verification.json']:shutil.copyfile(execution/name,out/name)
    (out/'execution.jsonl').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    completed=len(records)==46 and (execution/'completed_verification.json').exists()
    if completed:
        v=read(execution/'completed_verification.json');assert v==runner.audit(execution/'execution_contract.json',True)
        shutil.copyfile(execution/'completed_verification.json',out/'completed_verification.json')
    manifest=dict(status='COMPLETE' if completed else 'PARTIAL',contract_sha256=sha(execution/'execution_contract.json'),exporter_sha256=sha(__file__),completed_new=len(records),requested_new=46,reused=2,records=records,terminal_event_observed=events[-1]['status'],limits='Only completed receipts exported. Missing columns remain missing; no rerun,interpolation or full-batch completion claim from file presence.')
    (out/'snapshot_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    zip_path=out.with_suffix('.zip')
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.iterdir()):z.write(p,p.name)
    print(json.dumps(dict(status=manifest['status'],new=len(records),path=str(zip_path),sha256=sha(zip_path))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root.resolve())
