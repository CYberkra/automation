"""Independent prefix, geometry and native collocation-stencil audit; no solver."""
import argparse
import json
from pathlib import Path, PureWindowsPath
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def check(package):
    m=json.loads((package/'manifest.json').read_text('utf-8'))
    assert m['no_retry'] and len(m['groups'])==1
    group=m['groups'][0]
    assert group['id']=='native_probes_H0' and group['native_shape']==[8400,1700,1]
    assert group['expected_samples']==20352 and group['tx_m']==[169.35,38.95,0.0]
    np.testing.assert_allclose(group['rx_m'],[170.65,39.025,0],rtol=0,atol=1e-12)
    for k in ['input','geometry','material']:assert sha(package/group[k])==group[k+'_sha256']
    for k,digest in m['baseline_sha256'].items():assert sha(package/'baseline'/k)==digest
    for k,digest in m['dependency_sha256'].items():assert sha(package/'dependency'/k)==digest
    assert sha(package/group['geometry'])==sha(package/'baseline/geometry.h5')
    assert sha(package/group['material'])==sha(package/'baseline/materials.json')
    with h5py.File(package/'baseline/native_H0.h5','r') as h:
        assert float(h.attrs['dt'])==m['dt_s'] and h['rxs/rx1/Ez'].shape==(20352,)
    parent=json.loads((package/'dependency/execution_contract.json').read_text('utf-8'))
    cert=json.loads((package/'dependency/completed_verification.json').read_text('utf-8'))
    assert cert['completed'] and cert['contract_sha256']==sha(package/'dependency/execution_contract.json')
    assert sha(package/'baseline/native_H0.h5')==cert['groups'][0]['native_sha256']
    for name in ['profile.in','geometry.h5','materials.json','native_H0.h5']:
        identities=[digest for file,digest in parent['file_identities'].items()
                    if 'baseline' in PureWindowsPath(file).parts and PureWindowsPath(file).name==name]
        assert len(identities)==1 and sha(package/'baseline'/name)==identities[0]
    assert group['geometry_sha256']==parent['groups'][0]['geometry_sha256']
    assert group['material_sha256']==parent['groups'][0]['material_sha256']
    original=(package/'baseline/profile.in').read_bytes();card=(package/group['input']).read_bytes()
    assert card.startswith(original)
    lines=card[len(original):].decode('utf-8').splitlines()
    assert len(lines)==861 and len(m['probes'])==287 and m['receiver_count']==862
    expected_names=[];counter=0
    with h5py.File(package/group['geometry'],'r') as h:g=h['data'][:,:,0]
    for probe in m['probes']:
        anchors=probe['anchors'];assert len(anchors)==3
        center=np.array(anchors[0]['coord']);i,j,k=center;assert k==0
        np.testing.assert_array_equal(anchors[1]['coord'],center+[0,-1,0])
        np.testing.assert_array_equal(anchors[2]['coord'],center+[-1,0,0])
        assert [x['outputs'] for x in anchors]==[['Ez','Hx','Hy'],['Hx'],['Hy']]
        np.testing.assert_allclose([probe['x_m'],probe['y_m']],center[:2]*.025,atol=1e-12,rtol=0)
        assert abs(probe['wanted_y_m']-probe['y_m'])<=.0125+1e-12
        surface=(np.flatnonzero(g[i]!=0)[-1]+1)*.025
        bottom=(np.flatnonzero(g[i]==2)[-1]+1)*.025
        assert surface==probe['surface_m'] and bottom==probe['bottom_m']
        if probe['band'].startswith('air'):assert probe['y_m']>=surface+.2
        elif probe['band']=='mud_below':assert probe['y_m']<=bottom-.2
        else:assert bottom+.2<=probe['y_m']<=surface-.2
        for anchor in anchors:
            ii,jj,kk=anchor['coord'];assert 80<ii<8320 and 80<jj<1620 and kk==0
            support=g[ii-1:ii+2,jj-1:jj+2];assert np.all(support==probe['material_id'])
            line=lines[counter].split();assert line[0]=='#rx:'
            np.testing.assert_allclose([float(x) for x in line[1:4]],[ii*.025,jj*.025,.0125],rtol=0,atol=1e-12)
            assert line[4]==anchor['name'] and line[5:]==anchor['outputs']
            assert anchor['receiver_index']==2+counter
            expected_names.append(anchor['name']);counter+=1
    assert len(set(expected_names))==861 and counter==861
    assert sorted({r['x_m'] for r in m['probes']})==list(np.arange(154,174.001,.5))
    assert m['device_receiver_history_bytes']==862*20352*6*8
    assert m['saved_receiver_history_bytes']==20352*8*(1+287*5)
    assert m['collocation']=={'Ez':'raw anchor(i,j)','Hx':'mean raw anchors(i,j) and(i,j-1)',
                             'Hy':'mean raw anchors(i,j) and(i-1,j)',
                             'H_at_E_time':'mean H[n] and H[n+1], n=0..N-2',
                             'TM_invariant_axis':'z collapsed to0; x/y Yee offsets checked, not finite3D'}
    assert m['SFCW']['tones']==501 and m['SFCW']['min_Hz']==20e6 and m['SFCW']['max_Hz']==170e6 and m['SFCW']['step_Hz']==300000
    return {'status':'PASS_ONE_OBSERVER862_RX_STENCILS_AND_UNCHANGED_SOURCE_GEOMETRY_MATERIAL',
            'manifest_sha256':sha(package/'manifest.json'),'baseline_native_sha256':sha(package/'baseline/native_H0.h5'),
            'collocated_points':287,'raw_receivers':862,'new_solver_runs':0}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();assert not a.out.exists();j=check(a.package)
    a.out.write_text(json.dumps(j,indent=2)+'\n',encoding='utf-8');print(json.dumps(j))
