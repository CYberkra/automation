"""Independent whole-voxel/coordinate audit for low-loss boundary controls."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def audit(package,parent):
    m=json.loads((package/'manifest.json').read_text('utf-8'))
    old=json.loads((parent/'manifest.json').read_text('utf-8'))
    assert sha(parent/'manifest.json')==m['parent_manifest_sha256']
    byid={g['id']:g for g in old['groups']};rows=[]
    for g in m['groups']:
        role=g['id'].split('_')[-1];base=byid['span03_dc0003_'+role]
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        assert g['material_sha256']==base['material_sha256']
        sx,sy=(1600,0) if g['id'].startswith('sides') else (0,800)
        np.testing.assert_array_equal(g['translation_m'],[sx*.025,sy*.025,0])
        with h5py.File(parent/base['geometry']) as src,h5py.File(package/g['geometry']) as dst:
            before=src['data'][:];after=dst['data'][:]
            assert after.shape==tuple(g['native_shape'])==(8400+2*sx,1700+sy,1)
            np.testing.assert_array_equal(after[sx:sx+8400,sy:sy+1700],before)
            if sx:
                np.testing.assert_array_equal(after[:sx],np.repeat(before[:1],sx,axis=0))
                np.testing.assert_array_equal(after[-sx:],np.repeat(before[-1:],sx,axis=0))
            if sy:np.testing.assert_array_equal(after[:,:sy],np.repeat(before[:,:1],sy,axis=1))
            np.testing.assert_array_equal(src['material_keys'][:],dst['material_keys'][:])
            np.testing.assert_array_equal(dst.attrs['shape_nxyz'],g['native_shape'])
            np.testing.assert_array_equal(src.attrs['dx_dy_dz'],dst.attrs['dx_dy_dz'])
            assert dst.attrs['profile_x_offset_m']+sx*.025==src.attrs['profile_x_offset_m']
            assert dst.attrs['elevation_offset_m']+sy*.025==src.attrs['elevation_offset_m']
        card=(package/g['input']).read_text('utf-8').splitlines()
        original=(parent/base['input']).read_text('utf-8').splitlines()
        fields=lambda lines:{x.split(':',1)[0]:x.split(':',1)[1].split() for x in lines if ':' in x}
        x=fields(card);y=fields(original)
        for key in x:
            if key not in ['#title','#domain','#hertzian_dipole','#rx']:assert x[key]==y[key]
        assert x['#domain']==[f'{after.shape[0]*.025:g}',f'{after.shape[1]*.025:g}','inf']
        for key,offset in [('#hertzian_dipole',1),('#rx',0)]:
            a=np.array(x[key][offset:offset+3],float);b=np.array(y[key][offset:offset+3],float)
            np.testing.assert_allclose(a-b,[sx*.025,sy*.025,0],atol=1e-12,rtol=0)
            assert x[key][:offset]==y[key][:offset] and x[key][offset+3:]==y[key][offset+3:]
        rows.append(dict(id=g['id'],interior_exact=True,exterior_edge_continuation_exact=True,native_shape=g['native_shape'],material_sha256=g['material_sha256']))
    assert len(rows)==4 and len(m['reused'])==2
    for r in m['reused']:assert sha(package/r['path'])==r['native_sha256']
    return dict(status='PASS_INDEPENDENT_ORIGINAL_INTERIOR_EXTERIOR_MATERIAL_COORDINATES',script_sha256=sha(__file__),manifest_sha256=sha(package/'manifest.json'),groups=rows,limits=m['limits'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','parent','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();assert not a.out.exists();r=audit(a.package,a.parent)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))
