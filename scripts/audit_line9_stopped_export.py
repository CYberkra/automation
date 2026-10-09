"""Independent local receipt/transport audit and stopped-prefix rejection fixtures."""
import argparse
import hashlib
import json
from pathlib import Path

from collect_line9_stopped_continuation import stopped_prefix


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads(Path(p).read_text('utf-8'))


def fixtures():
    c = {'groups': [{'id': v} for v in ['a', 'b', 'c']]}
    s = lambda x: dict(status='STARTED', group=x)
    d = lambda x: dict(status='COMPLETED', group=x, raw_sha256='test')
    f = dict(status='FAILED', error='USER_CANCEL_REQUESTED')
    assert stopped_prefix(c, [s('a'), d('a'), s('b'), f])[1:] == (['b'], ['c'])
    assert stopped_prefix(c, [f])[1:] == ([], ['a', 'b', 'c'])
    assert stopped_prefix(c, [s('a'), d('a'), f])[1:] == ([], ['b', 'c'])
    bad = [[s('a'), d('a'), s('a'), f], [s('a'), d('b'), f],
           [d('a'), s('a'), f], [s('b'), f],
           [s('a'), dict(status='FAILED', error='SESSION_LEASE_EXPIRED')]]
    for events in bad:
        try:
            stopped_prefix(c, events)
        except AssertionError:
            pass
        else:
            raise AssertionError(events)
    return dict(valid_cases=3, rejection_cases=5)


def audit(source, archive, out):
    assert not out.exists()
    checks = fixtures()
    hashes = read(source / 'delivery_hashes.json')
    assert set(p.name for p in source.iterdir()) == set(hashes) | {'delivery_hashes.json'}
    assert all(sha(source / name) == digest for name, digest in hashes.items())
    c = read(source / 'execution_contract.json')
    n = read(source / 'continuation_execution_contract.json')
    cert = read(source / 'stopped_certificate.json')
    v = read(source / 'snapshot_manifest.json')
    old = [json.loads(x) for x in (source / 'execution.jsonl').read_text('utf-8').splitlines()]
    events = [json.loads(x) for x in (source / 'continuation_execution.jsonl').read_text('utf-8').splitlines()]
    assert old[0]['contract_sha256'] == cert['original_contract_sha256'] == sha(source / 'execution_contract.json')
    assert events[0]['contract_sha256'] == cert['continuation_contract_sha256'] == sha(source / 'continuation_execution_contract.json')
    assert n['parent_contract_sha256'] == sha(source / 'execution_contract.json')
    assert n['parent_execution_sha256'] == sha(source / 'execution.jsonl')
    assert old[-1] == dict(status='FAILED', error='SESSION_LEASE_EXPIRED')
    receipts, interrupted, unstarted = stopped_prefix(n, events)
    old_done = [e for e in old if e['status'] == 'COMPLETED' and e.get('group')]
    old_started = [e['group'] for e in old if e['status'] == 'STARTED' and e.get('group')]
    assert len(old_done) == len({e['group'] for e in old_done}) == 28
    assert len(old_started) == len(set(old_started)) == 29
    assert set(old_started) - {e['group'] for e in old_done} == {'low_x18675_H1'}
    assert [g['id'] for g in n['groups']] == [g['id'] for g in c['groups'] if g['id'] not in old_started]
    complete = {e['group']: e['raw_sha256'] for e in old_done + receipts}
    assert len(complete) == len(old_done) + len(receipts)
    assert complete == {r['id']: r['native_sha256'] for r in v['records']} == {r['id']: r['native_sha256'] for r in cert['groups']}
    assert set(complete) | set(v['consumed_incomplete']) | set(v['never_started']) == {g['id'] for g in c['groups']}
    assert v['consumed_incomplete'] == cert['consumed_incomplete'] == ['low_x18675_H1'] + interrupted
    assert v['never_started'] == cert['never_started'] == unstarted
    assert v['status'] == 'PARTIAL' and v['completed_new'] == cert['completed_new'] == len(complete)
    assert cert['continuation_completed'] == len(receipts) and cert['original_completed'] == len(old_done)
    assert cert['terminal_event'] == events[-1] and cert['owned_solver_processes'] == 0
    assert v['certificate_sha256'] == sha(source / 'stopped_certificate.json')
    reused = {r['id']: r['native_sha256'] for r in c['study_manifest']['reused']}
    assert reused == {r['id']: r['native_sha256'] for r in cert['reused_groups']}
    all_records = dict(complete, **reused)
    assert set(p.stem for p in source.glob('*.h5')) == set(all_records)
    assert all(sha(source / (sid + '.h5')) == digest for sid, digest in all_records.items())
    result = dict(status='PASS_TRANSPORT_RECEIPTS_PARTIAL_INVENTORY_NO_RETRY',
        auditor_sha256=sha(__file__), collector_sha256=sha(Path(__file__).with_name('collect_line9_stopped_continuation.py')),
        archive_sha256=sha(archive), certificate_sha256=sha(source / 'stopped_certificate.json'),
        completed_new=len(complete), continuation_completed=len(receipts), reused=len(reused), native_total=len(all_records),
        consumed_incomplete=v['consumed_incomplete'], never_started=unstarted, fixtures=checks,
        limits='Local immutable metadata and transport audit. Native waveform/DFT is audited separately; no physical attribution.')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'archive', 'out']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    audit(a.source, a.archive, a.out)
