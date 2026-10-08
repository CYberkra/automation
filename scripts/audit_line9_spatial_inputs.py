"""Independently audit all six cards against the original maps and planned sites."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def audit(package,parent):
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(package/'manifest.json');old=read(parent/'manifest.json');d=m['station_design']
    assert m['parent_manifest_sha256']==d['parent_manifest_sha256']==sha(parent/'manifest.json')
    base={g['id'].rsplit('_',1)[1]:g for g in old['groups'] if g['id'].startswith('span03_dc0003_')}
    assert len(m['groups'])==6 and len(m['reused'])==2
    cards=lambda p:{x.split()[0]:x.split()[1:] for x in p.read_text('utf-8').splitlines()}
    rows=[];maps={}
    for role,b in base.items():
        for key in ['geometry','material','input']:assert sha(parent/b[key])==b[key+'_sha256']
        with h5py.File(parent/b['geometry']) as h:maps[role]=h['data'][:]
    changed=maps['H0']!=maps['H1'];assert changed.any() and np.all(maps['H0'][changed]==2) and np.all(maps['H1'][changed]==3)
    for g in m['groups']:
        role=g['id'].rsplit('_',1)[1];b=base[role];s=next(s for s in d['stations'] if s['id']==g['station_id'])
        assert g['chainage_m']==s['chainage_m']
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        assert g['geometry_sha256']==b['geometry_sha256'] and g['material_sha256']==b['material_sha256']
        card=cards(package/g['input']);original=cards(parent/b['input'])
        mutable={'#title:','#hertzian_dipole:','#rx:'}
        assert {k:v for k,v in card.items() if k not in mutable}=={k:v for k,v in original.items() if k not in mutable}
        tx=np.array(list(map(float,card['#hertzian_dipole:'][1:4])));rx=np.array(list(map(float,card['#rx:'][:3])))
        np.testing.assert_allclose(tx,s['tx_m'],atol=1e-11,rtol=0);np.testing.assert_allclose(rx,s['rx_m'],atol=1e-11,rtol=0)
        assert abs((tx[0]+rx[0])/2+20-s['chainage_m'])<1e-10 and abs(rx[0]-tx[0]-1.3)<1e-10
        assert tx[2]==rx[2]==.0125 and card['#hertzian_dipole:'][0]=='z'
        np.testing.assert_allclose(np.array(g['tx_m'])[:2],tx[:2],atol=1e-11,rtol=0)
        np.testing.assert_allclose(np.array(g['rx_m'])[:2],rx[:2],atol=1e-11,rtol=0)
        assert g['tx_m'][2]==g['rx_m'][2]==0
        midpoint=(tx+rx)/2;column=maps['H1'][round(midpoint[0]/.025),:,0]
        # Independently count the top of each occupied material from voxel IDs.
        surface=(np.flatnonzero(column!=0)[-1]+1)*.025
        cover_base=(np.flatnonzero(column!=1)[np.flatnonzero(column!=1)<round(surface/.025)][-1]+1)*.025
        basal_top=(np.flatnonzero(column==3)[-1]+1)*.025
        np.testing.assert_allclose([surface-cover_base,cover_base-basal_top,surface-basal_top],
            [s['cover_thickness_m'],s['mud_thickness_m'],s['basal_depth_m']],atol=1e-10,rtol=0)
        agl=midpoint[1]-surface;assert abs(agl-8)<=.025
        for q in [tx,rx]:
            np.testing.assert_allclose(q[:2]/.025,np.rint(q[:2]/.025),atol=1e-9,rtol=0)
            assert maps[role][round(q[0]/.025),round(q[1]/.025),0]==0
            assert 2<q[0]<208 and 2<q[1]<40.5
        rows.append(dict(id=g['id'],chainage_m=g['chainage_m'],midpoint_AGL_m=float(agl),geometry_byte_identical=True,material_byte_identical=True,cover_m=float(surface-cover_base),mud_m=float(cover_base-basal_top)))
    return dict(status='PASS_INDEPENDENT_SIX_STATION_CARDS_MAPS_MATERIALS_VOXEL_DEPTHS',script_sha256=sha(__file__),manifest_sha256=sha(package/'manifest.json'),bottom_replacement_voxels=int(changed.sum()),groups=rows,limits='Input integrity and voxel geometry audit, not FDTD or field validation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','parent','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();assert not a.out.exists();r=audit(a.package,a.parent)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))
