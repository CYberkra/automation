"""Independent complete-card, voxel and reuse audit before diagnostic line solves."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha

def audit(package):
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(package/'manifest.json');bases=m['baseline_sources'];maps={};cards={};rows=[]
    def fields(path):
        result={}
        for line in path.read_text('utf-8').splitlines():
            key,*values=line.split();assert key not in result;result[key]=values
        return result
    for name,b in bases.items():
        parent=Path(b['parent']);g=b['group']
        for key in ['geometry','material','input']:assert sha(parent/g[key])==g[key+'_sha256']
        cards[name]=fields(parent/g['input'])
        with h5py.File(parent/g['geometry']) as h:maps[name]=h['data'][:,:,0];assert maps[name].shape==(8400,1700)
    np.testing.assert_array_equal(maps['low_H1'],maps['high_H1']);np.testing.assert_array_equal(maps['low_H0'],maps['high_H0'])
    replaced=maps['low_H1']!=maps['low_H0'];assert replaced.sum()==5606675
    assert np.all(maps['low_H1'][replaced]==3) and np.all(maps['low_H0'][replaced]==2)
    expected=sorted(set([220.-2*i for i in range(41)]+[198.6,146.4]),reverse=True)
    assert [s['chainage_m'] for s in m['stations']]==expected
    assert m['regular_line_chainage_m']==[220.-2*i for i in range(41)]
    station={s['chainage_m']:s for s in m['stations']};positions=[]
    for s in m['stations']:
        tx=np.array(s['tx_m']);rx=np.array(s['rx_m']);mid=(tx+rx)/2
        assert abs(mid[0]+20-s['chainage_m'])<1e-10 and abs(rx[0]-tx[0]-1.3)<1e-10
        for q in [tx,rx]:
            np.testing.assert_allclose(q[:2]/.025,np.rint(q[:2]/.025),atol=1e-9,rtol=0)
            assert q[2]==.0125 and 2<q[0]<208 and 2<q[1]<40.5
            assert maps['low_H1'][round(q[0]/.025),round(q[1]/.025)]==0
        col=maps['low_H1'][round(mid[0]/.025)];surface=(np.flatnonzero(col!=0)[-1]+1)*.025
        coverbottom=(np.flatnonzero(col==2)[-1]+1)*.025;sandtop=(np.flatnonzero(col==3)[-1]+1)*.025
        assert abs(mid[1]-surface-8)<=.025
        np.testing.assert_allclose([surface, surface-coverbottom, surface-sandtop],[s['geometry']['surface_y_m'],s['geometry']['cover_base']['depth_m'],s['geometry']['basal_sand']['depth_m']],atol=1e-10,rtol=0)
        positions.append(dict(chainage_m=s['chainage_m'],AGL_m=float(mid[1]-surface),top_PML_entry_clearance_m=float(40.5-max(tx[1],rx[1])),side_PML_entry_clearance_m=float(min(tx[0]-2,208-rx[0]))))
    for g in m['groups']:
        name=g['variant']+'_'+g['role'];base=bases[name]['group'];s=station[g['chainage_m']]
        for key in ['geometry','material','input']:assert sha(package/g[key])==g[key+'_sha256']
        assert g['geometry_sha256']==base['geometry_sha256'] and g['material_sha256']==base['material_sha256']
        c=fields(package/g['input']);old=cards[name];allowed={'#title:','#hertzian_dipole:','#rx:'}
        assert {k:v for k,v in c.items() if k not in allowed}=={k:v for k,v in old.items() if k not in allowed}
        assert c['#hertzian_dipole:'][0]=='z' and c['#hertzian_dipole:'][-1]=='pulse' and c['#rx:'][-1]=='Ez'
        for key,start,role in [('#hertzian_dipole:',1,'tx'),('#rx:',0,'rx')]:
            p=np.array(c[key][start:start+3],float);np.testing.assert_allclose(p,s[role+'_m'],atol=1e-11,rtol=0)
            np.testing.assert_allclose(p[:2],g[role+'_m'][:2],atol=1e-11,rtol=0);assert g[role+'_m'][2]==0
        rows.append(dict(id=g['id'],input_sha256=g['input_sha256']))
    keys=[]
    for r in m['reused']:
        p=package/r['path'];assert sha(p)==r['native_sha256']
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64
            assert h['rxs/rx1/Ez'].shape==(20352,)
            for node,role in [('srcs/src1','tx'),('rxs/rx1','rx')]:np.testing.assert_allclose(h[node].attrs['Position'],r[role+'_m'],atol=1e-11,rtol=0)
        keys.append((r['variant'],r['chainage_m'],r['role']))
    keys += [(g['variant'],g['chainage_m'],g['role']) for g in m['groups']]
    assert len(keys)==len(set(keys))==58 and len(rows)==48 and len(m['reused'])==10
    for x in expected:assert ('low',x,'H1') in keys
    return dict(status='PASS_INPUT_IDENTITIES_NOT_BOUNDARY_CONVERGENCE',script_sha256=sha(__file__),manifest_sha256=sha(package/'manifest.json'),groups=rows,positions=positions,
        new=48,reused=10,stations=43,bulk_execution_ready=False,blocking_check='Worst-clearance220m low-loss top/right extension controls before bulk freeze.',limits=m['limits'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert not a.out.exists()
    result=audit(a.package);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','new','reused','stations','bulk_execution_ready','blocking_check']}))
