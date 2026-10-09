"""Independent voxel/card invariants for a single190m high-loss H0 right control."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def check(package):
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(package/'manifest.json');base=package/'baseline';parent=package/'parent'
    c=read(parent/'execution_contract.json');v=read(parent/'completed_verification.json');a=read(parent/'analysis.json')
    for folder,key in [(base,'baseline_sha256'),(parent,'parent_sha256')]:
        for name,h in m[key].items():assert sha(folder/name)==h
    assert v['completed'] and c['max_runs']==1 and v['groups'][0]['snapshot_count']==500
    assert v['contract_sha256']==a['contract_sha256']==sha(parent/'execution_contract.json')
    assert a['verification_sha256']==sha(parent/'completed_verification.json') and a['observer_native_and_spectrum_bitwise_equal']
    b=c['study_manifest']['groups'][0]
    assert sha(base/'native_H0.h5')==v['groups'][0]['native_sha256']==a['new_native_sha256']
    assert sha(parent/'observer_profile.in')==b['input_sha256']
    before=(base/'profile.in').read_text('utf-8').splitlines()
    observer=(parent/'observer_profile.in').read_text('utf-8').splitlines()
    assert [line for line in observer if not line.startswith('#snapshot:')]==before
    assert sha(base/'geometry.h5')==b['geometry_sha256'] and sha(base/'materials.json')==b['material_sha256']
    assert len(m['groups'])==1 and m['no_retry'];g=m['groups'][0]
    assert g['id']=='right40_H0' and g['role']=='H0' and g['variant']=='high' and g['chainage_m']==190
    for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
    assert g['native_shape']==[10000,1700,1] and g['domain_height_m']==42.5
    assert g['material_sha256']==sha(base/'materials.json')
    after=(package/g['input']).read_text('utf-8').splitlines()
    strip=lambda lines:[line for line in lines if not line.startswith('#domain:')]
    assert strip(after)==strip(before) and '#domain: 250 42.5 inf' in after
    assert not any(line.startswith('#snapshot:') for line in after)
    with h5py.File(base/'geometry.h5') as src,h5py.File(package/g['geometry']) as dst:
        x=src['data'][:];y=dst['data'][:]
        assert y.shape==(10000,1700,1);np.testing.assert_array_equal(y[:8400],x)
        np.testing.assert_array_equal(y[8400:],np.repeat(x[-1:],1600,axis=0))
        assert set(np.unique(x))=={0,1,2}
        assert set(src)==set(dst)
        for key in src:
            if key!='data':np.testing.assert_array_equal(src[key][:],dst[key][:])
        for key in ['profile_x_offset_m','elevation_offset_m','dx_dy_dz']:
            np.testing.assert_array_equal(src.attrs[key],dst.attrs[key])
        for q in [g['tx_m'],g['rx_m']]:assert y[round(q[0]/.025),round(q[1]/.025),0]==0
    assert g['tx_m']==b['tx_m'] and g['rx_m']==b['rx_m'] and g['expected_samples']==20352
    assert m['dt_s']==c['study_manifest']['dt_s']
    from gprMax.hash_cmds_file import get_user_objects
    get_user_objects(after,input_dir=(package/g['input']).parent)
    return dict(status='PASS_INDEPENDENT_ONE_RIGHT40_INTERIOR_EDGE_CARD_AND_PARENT',script_sha256=sha(__file__),
        manifest_sha256=sha(package/'manifest.json'),original_interior_bitwise_equal=True,
        added_edge_continuation_cells=1600*1700,no_coordinate_translation=True,material_exact=True,
        original_H0_native_sha256=sha(base/'native_H0.h5'),solver_called=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();assert not a.out.exists();r=check(a.package)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))
