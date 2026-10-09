"""Six tamper/rejection fixtures for the two interface controls, CPU only."""
import argparse,json,shutil,tempfile
from pathlib import Path
import h5py,numpy as np
from audit_line9_cover_interface_inputs import check
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists();check(a.package);passed=[]
    for name in ['changed_tx','changed_material','changed_air','wrong_thickness','unfinished_parent','extra_case']:
        with tempfile.TemporaryDirectory(prefix='line9_cover_guard_') as tmp:
            p=Path(tmp)/'package';shutil.copytree(a.package,p);m=json.loads((p/'manifest.json').read_text('utf-8'));g=m['groups'][0]
            if name=='changed_tx':
                f=p/g['input'];s=f.read_text('utf-8');assert 'z 169.35 38.95' in s;f.write_text(s.replace('z 169.35 38.95','z 169.375 38.95'),encoding='utf-8');g['input_sha256']=sha(f)
            elif name=='changed_material':
                f=p/g['material'];f.write_bytes(f.read_bytes()+b'\n');g['material_sha256']=sha(f)
            elif name in ['changed_air','wrong_thickness']:
                if name=='wrong_thickness':g=m['groups'][1]
                f=p/g['geometry']
                with h5py.File(p/'baseline/geometry.h5') as h:x=h['data'][:]
                with h5py.File(f,'r+') as h:
                    if name=='changed_air':i,j,k=np.argwhere(x==0)[0];h['data'][i,j,k]=1
                    else:top=np.flatnonzero(x[4000,:,0]==2)[-1];h['data'][4000,top-20,0]=1
                    g['changed_voxels']=int(np.count_nonzero(h['data'][:]!=x))
                g['geometry_sha256']=sha(f)
            elif name=='unfinished_parent':
                f=p/'parent/completed_verification.json';v=json.loads(f.read_text('utf-8'));v['completed']=False;f.write_text(json.dumps(v),encoding='utf-8');m['parent_sha256'][f.name]=sha(f)
                f=p/'parent/analysis.json';v=json.loads(f.read_text('utf-8'));v['verification_sha256']=sha(p/'parent/completed_verification.json');f.write_text(json.dumps(v),encoding='utf-8');m['parent_sha256'][f.name]=sha(f)
            else:m['groups'].append(dict(g))
            (p/'manifest.json').write_text(json.dumps(m),encoding='utf-8')
            try:check(p)
            except AssertionError:passed.append(name)
            else:raise AssertionError('Accepted '+name)
    result=dict(status='PASS_SIX_COVER_INTERFACE_REJECTION_FIXTURES_NO_SOLVER',script_sha256=sha(__file__),package_manifest_sha256=sha(a.package/'manifest.json'),checks=passed,solver_runs=0)
    a.out.parent.mkdir(parents=True);a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['package','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
