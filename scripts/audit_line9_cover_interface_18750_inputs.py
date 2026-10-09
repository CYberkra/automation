"""Independent full-voxel audit of the two cover-bottom experiments; no solve."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha


def check(p):
    read=lambda f:json.loads(f.read_text('utf-8'))
    m=read(p/'manifest.json')
    for folder in ['baseline','parent','dependency190']:
        for name,h in m[folder+'_sha256'].items():assert sha(p/folder/name)==h
    oc=read(p/'parent/execution_contract.json');ov=read(p/'parent/stopped_certificate.json');oa=read(p/'parent/analysis.json');oi=read(p/'parent/independent_audit.json');oe=read(p/'parent/stopped_export_audit.json')
    assert ov['owned_solver_processes']==0 and ov['original_contract_sha256']==sha(p/'parent/execution_contract.json')
    assert oe['certificate_sha256']==sha(p/'parent/stopped_certificate.json') and oi['analysis_sha256']==sha(p/'parent/analysis.json')
    old=next(g for g in oc['groups'] if g['id']=='high_x18750_H0')
    assert old['role']=='H0' and old['variant']=='high' and old['chainage_m']==187.5
    assert m['original_group']==next(g for g in oc['study_manifest']['groups'] if g['id']==old['id'])
    receipts=[json.loads(x) for x in (p/'parent/execution.jsonl').read_text('utf-8').splitlines()]
    completed=[r for r in receipts if r.get('group')==old['id'] and r['status']=='COMPLETED'];assert len(completed)==1
    native_hash=next(g['native_sha256'] for g in ov['groups'] if g['id']==old['id'])
    assert sha(p/'baseline/native_H0.h5')==native_hash==completed[0]['raw_sha256']
    assert sha(p/'baseline/profile.in')==old['input_sha256'] and sha(p/'baseline/geometry.h5')==old['geometry_sha256'] and sha(p/'baseline/materials.json')==old['material_sha256']
    dc=read(p/'dependency190/execution_contract.json');dv=read(p/'dependency190/completed_verification.json');da=read(p/'dependency190/analysis.json');di=read(p/'dependency190/independent_audit.json')
    assert dv['completed'] and dv['contract_sha256']==da['contract_sha256']==di['contract_sha256']==sha(p/'dependency190/execution_contract.json')
    assert da['verification_sha256']==sha(p/'dependency190/completed_verification.json') and di['analysis_sha256']==sha(p/'dependency190/analysis.json') and dc['max_runs']==2
    assert old['geometry_sha256']==dc['study_manifest']['baseline_sha256']['geometry.h5'] and old['material_sha256']==dc['study_manifest']['baseline_sha256']['materials.json']
    assert m['first_difference_replication_gate_ns']==[120,250]
    assert m['no_retry'] and len(m['groups'])==2
    assert [g['id'] for g in m['groups']]==['no_cover_contrast_H0','cover_bottom_down050_H0']
    before=(p/'baseline/profile.in').read_bytes()
    assert not any(s.startswith('#snapshot:') for s in before.decode().splitlines())
    from gprMax.hash_cmds_file import get_user_objects
    rows=[]
    with h5py.File(p/'baseline/geometry.h5') as src:
        x=src['data'][:];assert x.shape==(8400,1700,1) and set(np.unique(x))=={0,1,2}
        for g in m['groups']:
            for k in ['input','geometry','material']:assert sha(p/g[k])==g[k+'_sha256']
            assert (p/g['input']).read_bytes()==before and g['material_sha256']==sha(p/'baseline/materials.json')
            assert g['native_shape']==[8400,1700,1] and g['domain_height_m']==42.5 and g['expected_samples']==20352
            assert g['tx_m']==old['tx_m'] and g['rx_m']==old['rx_m'] and g['chainage_m']==187.5
            assert g['geometry_sha256']==next(d['geometry_sha256'] for d in dc['study_manifest']['groups'] if d['id']==g['id'])
            with h5py.File(p/g['geometry']) as dst:
                y=dst['data'][:];assert y.shape==x.shape and set(src)==set(dst)
                np.testing.assert_array_equal(y[x==0],x[x==0]);np.testing.assert_array_equal(y[x==1],x[x==1])
                for name in src:
                    if name!='data':np.testing.assert_array_equal(src[name][:],dst[name][:])
                for name,val in src.attrs.items():np.testing.assert_array_equal(val,dst.attrs[name])
                changed=y!=x
                if g['id']=='no_cover_contrast_H0':
                    np.testing.assert_array_equal(changed,x==2);assert np.all(y[x==2]==1)
                else:
                    for i in range(8400):
                        top=np.flatnonzero(x[i,:,0]==2)[-1]
                        np.testing.assert_array_equal(np.flatnonzero(changed[i,:,0]),np.arange(top-19,top+1))
                        assert np.all(y[i,top-19:top+1,0]==1) and np.all(y[i,:top-19,0]==x[i,:top-19,0])
                    assert np.count_nonzero(changed)==8400*20
                assert int(np.count_nonzero(changed))==g['changed_voxels']
                for q in [g['tx_m'],g['rx_m']]:assert y[round(q[0]/.025),round(q[1]/.025),0]==0
            get_user_objects(before.decode().splitlines(),input_dir=(p/g['input']).parent)
            rows.append(dict(id=g['id'],changed_voxels=g['changed_voxels'],air_surface_source_material_cards_preserved=True))
    assert m['dt_s']==oc['study_manifest']['dt_s']
    return dict(status='PASS_TWO18750_INDEPENDENT_VOXELS_CARD_RECEIPT_REPLICATION',script_sha256=sha(__file__),manifest_sha256=sha(p/'manifest.json'),rows=rows,solver_runs=0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['package','out']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();assert not a.out.exists();r=check(a.package);a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))
