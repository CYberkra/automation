"""Reject physically invalid yet rehashed probe packages before consuming one attempt."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
from audit_line9_native_probe_inputs import check
from hs_capsule_identity import sha256 as sha


def main(a):
    original=check(a.package);cases=[]
    a.work.mkdir(parents=True,exist_ok=True)
    for name in ['source_changed_rehashed','Hx_wrong_neighbor','Hy_wrong_neighbor','main_rx_changed',
                 'retry_enabled','wrong_time_collocation']:
        with tempfile.TemporaryDirectory(prefix=name+'_',dir=a.work) as folder:
            package=Path(folder)/'package';shutil.copytree(a.package,package)
            path=package/'manifest.json';m=json.loads(path.read_text('utf-8'));g=m['groups'][0]
            card=package/g['input']
            if name=='source_changed_rehashed':
                p=package/'baseline/profile.in';p.write_bytes(p.read_bytes().replace(b'ricker 40 ',b'ricker 41 '))
                card.write_bytes(card.read_bytes().replace(b'ricker 40 ',b'ricker 41 '))
                m['baseline_sha256']['profile.in']=sha(p)
            elif name in ['Hx_wrong_neighbor','Hy_wrong_neighbor']:
                anchor=m['probes'][0]['anchors'][1 if name.startswith('Hx') else 2]
                anchor['coord'][1 if name.startswith('Hx') else 0]+=2
            elif name=='main_rx_changed':
                card.write_bytes(card.read_bytes().replace(b'#rx: 170.65 ',b'#rx: 170.675 ',1))
            elif name=='retry_enabled':m['no_retry']=False
            else:m['collocation']['H_at_E_time']='use H[n] without Yee time alignment'
            g['input_sha256']=sha(card);path.write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
            try:check(package)
            except (AssertionError,ValueError):cases.append(name)
            else:raise AssertionError('Invalid package accepted: '+name)
    result={'status':'PASS_SIX_REHASHED_PROBE_PACKAGE_REJECTIONS','positive':original,'rejected':cases,
            'script_sha256':sha(__file__),'new_solver_runs':0}
    assert not a.out.exists();a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['package','work','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
