"""Independent immutable prefix and homogeneous flat-cover stencil audit."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from audit_line9_native_probe_inputs import check as check_parent


def check(package):
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(package/'manifest.json');old=read(package/'parent/manifest.json');check_parent(package/'parent')
    assert m['no_retry'] and len(m['groups'])==1 and m['receiver_count']==1231 and len(m['probes'])==410
    assert m['parent_manifest_sha256']==sha(package/'parent/manifest.json')
    assert m['parent_contract_sha256']==sha(package/'parent_execution_contract.json')
    assert m['parent_completed_sha256']==sha(package/'parent_completed_verification.json')
    c=read(package/'parent_execution_contract.json');v=read(package/'parent_completed_verification.json')
    assert c['study_manifest']==old and v['completed'] and v['contract_sha256']==m['parent_contract_sha256']
    assert m['reference_native_sha256']==v['groups'][0]['native_sha256']=='46288bae706a5534940baf3b87748c0f68eea02ffc6e1ba6ae1ae17c6efbe68e'
    assert m['geometry_proposal_sha256']==sha(package/'geometry_proposal.json')
    g=m['groups'][0];og=old['groups'][0];assert g['id']=='native_cover_planes_H0'
    for key in ['native_shape','tx_m','rx_m','expected_samples','dt_s','time_window_ns','source_type','source_frequency_Hz','source_amplitude_A']:
        assert g[key]==og[key],key
    for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
    for key in ['geometry','material']:assert g[key+'_sha256']==og[key+'_sha256']
    original=(package/'parent'/og['input']).read_bytes();card=(package/g['input']).read_bytes();assert card.startswith(original)
    lines=card[len(original):].decode('utf-8').splitlines();assert len(lines)==369
    assert m['probes'][:287]==old['probes'] and m['previous_receiver_count']==862 and m['previous_logical_points']==287
    assert m['dt_s']==old['dt_s'] and m['collocation']==old['collocation'] and m['SFCW']==old['SFCW']
    assert m['cover_plane_y_m']==[27.025,28.525,30.025]
    names={q['name'] for p in old['probes'] for q in p['anchors']};counter=0
    with h5py.File(package/g['geometry']) as h:geometry=h['data'][:,:,0]
    for yi,y in enumerate(m['cover_plane_y_m']):
        for xi,x in enumerate(np.arange(154,174.001,.5)):
            p=m['probes'][287+yi*41+xi];i=round(x/.025);j=round(y/.025)
            assert p['x_m']==x and p['y_m']==j*.025 and p['material_id']==1 and p['wanted_y_m']==p['y_m']
            assert p['band']=='flat_cover_'+str(j)
            assert p['surface_m']==(np.flatnonzero(geometry[i]!=0)[-1]+1)*.025
            assert p['bottom_m']==(np.flatnonzero(geometry[i]==2)[-1]+1)*.025
            for k,(di,dj,outputs) in enumerate([(0,0,['Ez','Hx','Hy']),(0,-1,['Hx']),(-1,0,['Hy'])]):
                q=p['anchors'][k];ii=i+di;jj=j+dj
                assert q['coord']==[ii,jj,0] and q['outputs']==outputs and q['receiver_index']==863+counter
                np.testing.assert_allclose(q['position_m'],[ii*.025,jj*.025,0],rtol=0,atol=1e-12)
                assert q['name'] not in names;names.add(q['name'])
                assert 80<ii<8320 and 80<jj<1620 and (geometry[ii-1:ii+2,jj-1:jj+2]==1).all()
                line=lines[counter].split();assert line[0]=='#rx:' and line[4]==q['name'] and line[5:]==outputs
                np.testing.assert_allclose([float(v) for v in line[1:4]],[ii*.025,jj*.025,.0125],rtol=0,atol=1e-12)
                counter+=1
    assert counter==369 and len(names)==1230
    assert m['device_receiver_history_bytes']==1231*20352*6*8 and m['saved_receiver_history_bytes']==20352*8*(1+5*410)
    return {'status':'PASS_IMMUTABLE1231_RX_AND_THREE_UNIFORM_COVER_PLANES','manifest_sha256':sha(package/'manifest.json'),
            'geometry_sha256':g['geometry_sha256'],'original_points_unchanged':287,'new_points':123,'raw_receivers':1231,'new_solver_runs':0}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not a.out.exists();v=check(a.package);a.out.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');print(json.dumps(v))
