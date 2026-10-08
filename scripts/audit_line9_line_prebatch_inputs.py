"""Independent voxel, card and coordinate checks for220m boundary controls."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha

def main(a):
    assert not a.out.exists();m=json.loads((a.package/'manifest.json').read_text('utf-8'));rows=[]
    fields=lambda p:{line.split()[0]:line.split()[1:] for line in p.read_text('utf-8').splitlines()}
    for g in m['groups']:
        ref=m['baseline_sources']['low_'+g['role']];p=Path(ref['parent']);b=ref['group']
        for key in ['input','material','geometry']:assert sha(a.package/g[key])==g[key+'_sha256'] and sha(p/b[key])==b[key+'_sha256']
        assert g['material_sha256']==b['material_sha256']
        with h5py.File(p/b['geometry']) as src,h5py.File(a.package/g['geometry']) as dst:
            before=src['data'][:];after=dst['data'][:];np.testing.assert_array_equal(after[:8400,:1700],before)
            shape={'base':(8400,1700,1),'top20':(8400,2500,1),'right40':(10000,1700,1)}[g['factor']];assert after.shape==shape==tuple(g['native_shape'])
            if g['factor']=='top20':assert np.all(after[:,1700:]==0)
            if g['factor']=='right40':np.testing.assert_array_equal(after[8400:],np.repeat(before[-1:],1600,axis=0))
            for key in ['material_keys']:np.testing.assert_array_equal(src[key][:],dst[key][:])
            for key in ['profile_x_offset_m','elevation_offset_m','dx_dy_dz']:np.testing.assert_array_equal(src.attrs[key],dst.attrs[key])
            for q in [g['tx_m'],g['rx_m']]:assert after[round(q[0]/.025),round(q[1]/.025),0]==0
        c=fields(a.package/g['input']);old=fields(p/b['input']);mutable={'#title:','#domain:','#hertzian_dipole:','#rx:'}
        assert {k:v for k,v in c.items() if k not in mutable}=={k:v for k,v in old.items() if k not in mutable}
        np.testing.assert_allclose(np.array(c['#domain:'][:2],float),np.array(shape[:2])*.025,rtol=0,atol=1e-11)
        for key,start,role in [('#hertzian_dipole:',1,'tx'),('#rx:',0,'rx')]:
            q=np.array(c[key][start:start+3],float);np.testing.assert_allclose(q,m['station'][role+'_m'],atol=1e-11,rtol=0)
            np.testing.assert_allclose(q[:2],g[role+'_m'][:2],atol=1e-11,rtol=0);assert q[2]==.0125 and g[role+'_m'][2]==0
        rows.append(dict(id=g['id'],shape=shape,original_voxels_exact=True,material_exact=True))
    assert len(rows)==6 and len({g['id'] for g in m['groups']})==6
    result=dict(status='PASS_INDEPENDENT_SIX_PREBATCH_INTERIOR_EXTERIOR_MATERIAL_COORDINATES',script_sha256=sha(__file__),manifest_sha256=sha(a.package/'manifest.json'),groups=rows,limits=m['limits'])
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],groups=len(rows))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
