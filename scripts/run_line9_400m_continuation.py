"""New, interruptible 400m coarse-scan attempt; reuse certified X220, never retry consumed runs."""
import argparse
import json
from pathlib import Path

from build_pdf_profile_geometry import digest, save_json
import run_line9_2d_first as runner


def freeze(package,out,reference_study):
    reference_contract=reference_study/'execution_contract.json'
    verified=runner.audit(reference_contract,True,group_ids={'full2d_pilot_0'})
    row=verified['groups'][0]
    events=[json.loads(line) for line in (reference_study/'execution.jsonl').read_text('utf-8').splitlines()]
    matches=[e for e in events if e.get('status')=='COMPLETED' and e.get('group')=='full2d_pilot_0']
    if len(matches)!=1 or matches[0]['raw_sha256']!=row['raw_sha256']:
        raise ValueError('Reference must have exactly one matching completed event')
    task=json.loads(runner.TASK.read_text('utf-8'))
    if row['profile_x_m']!=220 or row['acquisition_s_m']!=0:
        raise ValueError('Only certified X220 reference can be reused')
    ids=[name for name in task['acquisition']['preview_ids'] if name!='full2d_s0000']
    runner.freeze(package,out,'continuation',station_ids=ids)
    contract=out/'execution_contract.json'
    c=json.loads(contract.read_text('utf-8'))
    c.update(approval_basis='User: run current400m large domain until interrupted; compute as many as possible. Existing99-station first-scan scope retained.',
        expected_preview_total=99,reuse={'full2d_s0000':dict(raw_path=row['raw_path'],raw_sha256=row['raw_sha256'],pilot_id='full2d_pilot_0')},
        cancel_file=str((out/'USER_STOP').resolve()),lease_file=str((out/'session_heartbeat').resolve()),
        max_lease_age_s=120,analysis_every_new_traces=5,
        reference_contract_sha256=digest(reference_contract))
    c['file_identities'][str(reference_contract.resolve())]=digest(reference_contract)
    c['file_identities'][row['raw_path']]=row['raw_sha256']
    c['code_identities'][str(Path(__file__).resolve())]=digest(Path(__file__))
    save_json(contract,c)
    save_json(out/'reference_verification.json',verified)
    save_json(out/'preflight_verification.json',runner.audit(contract,False))
    (out/'group_audits').mkdir()
    print('Frozen98 new full-domain solves +1 verified reuse;120s session lease; no solver called',flush=True)


def run(out,report_root):
    import hs4_station_grid_controls as supervisor
    from analyze_line9_2d_sfcw import main as analyze
    c=json.loads((out/'execution_contract.json').read_text('utf-8'))
    if c['stage']!='continuation': raise ValueError('Dedicated continuation contract required')
    for item in c['reuse'].values():
        if digest(Path(item['raw_path']))!=item['raw_sha256']: raise ValueError('Reference changed')
    count=0
    def completed(group):
        nonlocal count
        verified=runner.audit(out/'execution_contract.json',True,group_ids={group['id']})
        save_json(out/'group_audits'/f"{group['id']}.json",verified)
        count+=1
        if count==1 or count%c['analysis_every_new_traces']==0:
            analyze(out,report_root/f'completed_{count+1:03d}',partial=True)
    supervisor.group_completed_hook=completed
    try:
        runner.run(out)
        analyze(out,report_root/'complete_099')
    finally:
        supervisor.group_completed_hook=None


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['freeze','run'])
    p.add_argument('--out',type=Path,required=True);p.add_argument('--package',type=Path)
    p.add_argument('--reference-study',type=Path);p.add_argument('--report-root',type=Path)
    a=p.parse_args()
    if a.action=='freeze':
        if a.package is None or a.reference_study is None: p.error('freeze requires package/reference-study')
        freeze(a.package.resolve(),a.out.resolve(),a.reference_study.resolve())
    else:
        if a.report_root is None: p.error('run requires report-root')
        run(a.out.resolve(),a.report_root.resolve())
