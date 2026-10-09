"""Read-only native audit/export of a user-stopped continuation; never solve/retry."""
import argparse
import json
from pathlib import Path
import shutil
import zipfile


def stopped_prefix(contract, events):
    """Require single attempts in frozen order and a deliberate terminal stop."""
    assert events[-1] == dict(status='FAILED', error='USER_CANCEL_REQUESTED')
    expected = [g['id'] for g in contract['groups']]
    started = [e['group'] for e in events if e['status'] == 'STARTED' and e.get('group')]
    done = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    ids = [e['group'] for e in done]
    assert len(expected) == len(set(expected))
    assert len(started) == len(set(started)) and started == expected[:len(started)]
    assert len(ids) == len(set(ids)) and ids == started[:len(ids)]
    assert 0 <= len(started) - len(ids) <= 1
    # A completed receipt must follow its STARTED record and precede the next attempt.
    active = None
    for event in events:
        if event['status'] == 'STARTED' and event.get('group'):
            assert active is None
            active = event['group']
        elif event['status'] == 'COMPLETED' and event.get('group'):
            assert active == event['group']
            active = None
    return done, started[len(ids):], expected[len(started):]


def collect(execution, out):
    from hs_capsule_identity import sha256 as sha
    import line9_v401_dense_continuation as frozen

    assert not out.exists() and not out.with_suffix('.zip').exists()
    cp = execution / 'execution_contract.json'
    audit = frozen.audit(cp, False)
    n = frozen.read(cp)
    events = [json.loads(s) for s in (execution / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(cp)
    receipts, interrupted, unstarted = stopped_prefix(n, events)
    frozen.no_live_owned(n, execution)
    parent = Path(n['parent_contract'])
    old, old_events, _, old_rows = frozen.predecessor(parent)
    groups = {g['id']: g for g in n['groups']}
    new_rows = [frozen.native_record(groups[e['group']], Path(groups[e['group']]['input']).with_suffix('.h5'), e['raw_sha256']) for e in receipts]
    rows = old_rows + new_rows
    assert len(rows) == len({r['id'] for r in rows})
    consumed = n['consumed_incomplete'] + interrupted
    expected = {g['id'] for g in old['groups']}
    assert {r['id'] for r in rows}.isdisjoint(consumed + unstarted)
    assert {r['id'] for r in rows} | set(consumed + unstarted) == expected
    reused = []
    for r in old['study_manifest']['reused']:
        raw = Path(n['parent_package']) / r['path']
        assert sha(raw) == r['native_sha256']
        identity = dict(r, expected_samples=20352, dt_s=old['study_manifest']['dt_s'])
        reused.append(frozen.native_record(identity, raw, r['native_sha256']))
    out.mkdir(parents=True)
    for r in rows + reused:
        shutil.copyfile(r['native_path'], out / (r['id'] + '.h5'))
    for directory, prefix in [(parent.parent, ''), (execution, 'continuation_')]:
        for name in ['execution_contract.json', 'execution.jsonl', 'preflight_verification.json']:
            shutil.copyfile(directory / name, out / (prefix + name))
    certificate = dict(status='PASS_USER_STOPPED_NATIVE_PARTIAL_NO_RETRY',
        original_contract_sha256=sha(parent), continuation_contract_sha256=sha(cp),
        original_completed=len(old_rows), continuation_completed=len(new_rows), completed_new=len(rows),
        reused=len(reused), native_total=len(rows + reused), consumed_incomplete=consumed,
        never_started=unstarted, owned_solver_processes=0,
        terminal_event=events[-1], original_terminal_event=old_events[-1],
        continuation_preflight=audit, groups=rows, reused_groups=reused,
        limits='Completed native files only. User stop, not full batch completion. No wavefield/physical-path/site validation.')
    frozen.save(out / 'stopped_certificate.json', certificate)
    frozen.save(out / 'snapshot_manifest.json', dict(status='PARTIAL', contract_sha256=sha(parent),
        completed_new=len(rows), requested_new=len(old['groups']), reused=len(reused),
        records=[dict(id=r['id'], native_sha256=r['native_sha256']) for r in rows],
        consumed_incomplete=consumed, never_started=unstarted,
        terminal_event_observed='OLD_FAILED_LEASE_NEW_FAILED_USER_CANCEL',
        certificate_sha256=sha(out / 'stopped_certificate.json')))
    frozen.save(out / 'delivery_hashes.json', {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    with zipfile.ZipFile(out.with_suffix('.zip'), 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.iterdir()):
            z.write(p, p.name)
    print(json.dumps(dict(path=str(out.with_suffix('.zip')), sha256=sha(out.with_suffix('.zip')),
        completed_new=len(rows), reused=len(reused), continuation_completed=len(new_rows),
        consumed_incomplete=consumed, never_started=unstarted)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    collect(a.execution.resolve(), a.out.resolve())
