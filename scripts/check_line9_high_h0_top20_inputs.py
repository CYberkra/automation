"""Six rejection fixtures for the single top-boundary input; no solver execution."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import h5py
from hs_capsule_identity import sha256 as sha
from audit_line9_high_h0_top20_inputs import check


def main(a):
    assert not a.out.exists();check(a.package);results=[]
    for name in ['changed_tx','changed_material','changed_original_voxel','nonair_padding','unfinished_parent','extra_case']:
        with tempfile.TemporaryDirectory(prefix='line9_top20_guard_') as tmp:
            p=Path(tmp)/'package';shutil.copytree(a.package,p)
            m=json.loads((p/'manifest.json').read_text('utf-8'));g=m['groups'][0]
            if name=='changed_tx':
                f=p/g['input'];s=f.read_text('utf-8').replace('z 169.35 38.95','z 169.375 38.95');f.write_text(s,encoding='utf-8');g['input_sha256']=sha(f)
            elif name=='changed_material':
                f=p/g['material'];f.write_bytes(f.read_bytes()+b'\n');g['material_sha256']=sha(f)
            elif name in ['changed_original_voxel','nonair_padding']:
                f=p/g['geometry']
                with h5py.File(f,'r+') as h:
                    if name=='changed_original_voxel':h['data'][4000,1000,0]=(int(h['data'][4000,1000,0])+1)%3
                    else:h['data'][4000,1800,0]=1
                g['geometry_sha256']=sha(f)
            elif name=='unfinished_parent':
                f=p/'parent/completed_verification.json';v=json.loads(f.read_text('utf-8'));v['completed']=False
                f.write_text(json.dumps(v),encoding='utf-8');m['parent_sha256'][f.name]=sha(f)
                f=p/'parent/analysis.json';v=json.loads(f.read_text('utf-8'));v['verification_sha256']=sha(p/'parent/completed_verification.json')
                f.write_text(json.dumps(v),encoding='utf-8');m['parent_sha256'][f.name]=sha(f)
            else:m['groups'].append(dict(g))
            (p/'manifest.json').write_text(json.dumps(m),encoding='utf-8')
            try:check(p)
            except AssertionError:results.append(name)
            else:raise AssertionError('Guard accepted '+name)
    r=dict(status='PASS_SIX_TOP20_REJECTION_FIXTURES_NO_SOLVER',script_sha256=sha(__file__),
        package_manifest_sha256=sha(a.package/'manifest.json'),checks=results,solver_runs=0)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
