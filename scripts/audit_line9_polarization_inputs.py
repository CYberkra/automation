"""Independent constitutive and voxel invariants before/after factorial solves."""
import argparse
import copy
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def audit(package,parent):
    m=json.loads((package/'manifest.json').read_text('utf-8'))
    old=json.loads((parent/'manifest.json').read_text('utf-8'))
    assert sha(parent/'manifest.json')==m['parent_manifest_sha256']
    base=next(g for g in old['groups'] if g['id']=='baseline_H1')
    original=json.loads((parent/base['material']).read_text('utf-8'))
    source_lines=[x for x in (parent/base['input']).read_text('utf-8').splitlines() if not x.startswith('#title:')]
    rows=[];geo={}
    for g in m['groups']:
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        lines=[x for x in (package/g['input']).read_text('utf-8').splitlines() if not x.startswith('#title:')]
        assert lines==source_lines
        db=json.loads((package/g['material']).read_text('utf-8'))
        mud=db['materials']['material_002_mudstone'];ei=mud['base']['relative_permittivity'];dc=mud['base']['electric_conductivity_s_per_m']
        de=mud['poles'][0]['relative_permittivity_difference'];tau=mud['poles'][0]['relaxation_time_s']
        f=np.array([20,95,170])*1e6
        # Real-valued formulas independent of generator and gprMax API.
        real=ei+de/(1+(2*np.pi*f*tau)**2)
        pol=de*2*np.pi*f*tau/(1+(2*np.pi*f*tau)**2)
        assert abs(real[1]-12)<1e-12 and abs(real[0]-real[2]-.3)<1e-12
        assert tau==8e-9 and dc==(.0003 if 'dc0003' in g['id'] else .003)
        restored=copy.deepcopy(db);restored['database']=copy.deepcopy(original['database'])
        restored['materials']['material_002_mudstone']=copy.deepcopy(original['materials']['material_002_mudstone'])
        assert restored==original
        assert mud['averagable']==original['materials']['material_002_mudstone']['averagable']
        assert mud['base']['relative_permeability']==1 and mud['base']['magnetic_conductivity_s_per_m']==0
        assert mud['name']=='mudstone' and mud['model']=='debye'
        role=g['id'].split('_')[-1];geo.setdefault(role,package/g['geometry'])
        assert sha(geo[role])==g['geometry_sha256']
        rows.append(dict(id=g['id'],epsilon_real_20_95_170=real.tolist(),epsilon_imag_polarization_20_95_170=pol.tolist(),sigma_DC_S_m=dc,epsilon_infinity=ei,delta_epsilon=de,tau_s=tau))
    changed=0
    with h5py.File(geo['H0']) as h0,h5py.File(geo['H1']) as h1:
        x=h1['data'][:,:,0];y=h0['data'][:,:,0]
        assert x.shape==y.shape==(8400,1700)
        for col,bg in zip(x,y):
            n=0
            while n<len(col) and col[n]==3:n+=1
            expect=col.copy();expect[:n]=2
            np.testing.assert_array_equal(bg,expect);changed+=n
    assert changed==5606675
    expected={r['id']:r for r in m['reused']}
    assert len(expected)==4
    for r in expected.values():assert sha(package/r['path'])==r['native_sha256']
    return dict(status='PASS_INDEPENDENT_FACTORIAL_INPUT_MATERIAL_VOXEL_AUDIT',script_sha256=sha(__file__),manifest_sha256=sha(package/'manifest.json'),new_groups=rows,reused_native_count=4,bottom_connected_voxels=changed,limits=m['limits'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','parent','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();assert not a.out.exists()
    result=audit(a.package,a.parent)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))
