"""Fresh H0-only recovery after a preserved session-lease interruption."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import line9_basal_pair as pair
import hs4_station_grid_controls as supervisor


def prepare(root):
    old=root/'execution/execution_contract.json'
    c=json.loads(old.read_text('utf-8'))
    events=[json.loads(s) for s in (old.parent/'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[-1]==dict(status='FAILED',error='SESSION_LEASE_EXPIRED')
    completed=[r for r in events if r.get('group')=='H1' and r['status']=='COMPLETED']
    assert len(completed)==1
    assert pair.sha(Path(c['groups'][0]['input']).with_suffix('.h5'))==completed[0]['raw_sha256']
    dest=root/'recovery_prepared/H0'
    if dest.exists():raise ValueError('Fresh recovery required')
    source=Path(c['package'])/'H0'
    for relative in ['cases/full2d_c0008/profile.in','geometries/full2d_compact.h5','geometries/line9_research_materials_v1_smoothed.json']:
        p=dest/relative;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/relative,p)
        assert pair.sha(p)==pair.sha(source/relative)
    c['groups'][1]['input']=str(dest/'cases/full2d_c0008/profile.in')
    c['file_identities'].update({str(p):pair.sha(p) for p in dest.rglob('*') if p.is_file()})
    c['code_identities'][str(Path(__file__).resolve())]=pair.sha(__file__)
    c['recovery_basis']='Prior H1 complete and hash verified; interrupted H0 preserved. One fresh H0 attempt, identical scientific inputs; no H1 rerun.'
    c['prior_execution_contract_sha256']=pair.sha(old)
    out=root/'recovery_execution';out.mkdir()
    c.update(cancel_file=str(out/'USER_STOP'),lease_file=str(out/'session_heartbeat'),max_lease_age_s=600)
    pair.save(out/'paired_verification_contract.json',c)
    c['groups']=c['groups'][1:];c['max_runs']=1
    pair.save(out/'execution_contract.json',c)
    (out/'session_heartbeat').write_text('Active authorized recovery\n')
    pair.save(out/'preflight_verification.json',pair.audit(out/'execution_contract.json',False))


def run(root):
    out=root/'recovery_execution'
    def audit(path,completed):
        return pair.audit(out/'paired_verification_contract.json' if completed else path,completed)
    pair.resources();supervisor.audit=audit
    supervisor.run(out/'execution_contract.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','run']);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();globals()[a.action](a.root.resolve())
