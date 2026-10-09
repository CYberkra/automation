"""CPU combined-receipt negative checks; no fabricated physical evidence."""
import argparse
import copy
import json
from pathlib import Path

from hs_capsule_identity import sha256 as sha
from line9_v401_dense_continuation import read, inventory
from line9_v401_dense_join import validate_join


def main(a):
    assert not a.out.exists()
    old = read(a.source / 'execution_contract.json')
    oe = [json.loads(s) for s in (a.source / 'execution.jsonl').read_text('utf-8').splitlines()]
    _, groups = inventory(old, oe)
    new = dict(max_runs=17, no_retry=True, groups=groups)
    events = [dict(status='STARTED', group=g['id']) for g in groups]
    events += [dict(status='COMPLETED', group=g['id'], raw_sha256='fixture_' + g['id']) for g in groups]
    events += [dict(status='COMPLETED', traces=17)]
    cert = dict(completed=True, groups=[dict(id=g['id'], native_sha256='fixture_' + g['id']) for g in groups])
    assert len(validate_join(old, oe, new, events, cert)) == 45
    checks = ['Fixture join accepts45unique receipts and exactly the actual consumed missing low18675']
    def reject(n, e, v):
        try:
            validate_join(old, oe, n, e, v)
        except AssertionError:
            return
        raise AssertionError('Invalid combined lineage accepted')
    q = copy.deepcopy(new); q['groups'][0] = old['groups'][0]; reject(q, events, cert)
    checks.append('Retry of completed or consumed group rejected')
    reject(new, events[:-1], cert); checks.append('Nonterminal continuation rejected')
    q = copy.deepcopy(events); q[-1]['traces'] = 46; reject(new, q, cert)
    checks.append('False46completion rejected')
    q = copy.deepcopy(events); q.insert(-1, q[17]); reject(new, q, cert)
    checks.append('Duplicate completion cannot hide missing coverage')
    q = copy.deepcopy(events); q[0]['group'] = 'low_x18675_H1'; reject(new, q, cert)
    checks.append('Consumed incomplete trace cannot be silently retried')
    q = copy.deepcopy(cert); q['groups'][0]['native_sha256'] = 'wrong'; reject(new, events, q)
    checks.append('Certificate/native receipt mismatch rejected')
    q = copy.deepcopy(cert); q['groups'].append(q['groups'][0]); reject(new, events, q)
    checks.append('Duplicate native certificate row rejected')
    q = copy.deepcopy(cert); q['completed'] = False; reject(new, events, q)
    checks.append('Incomplete certificate rejected')
    a.out.mkdir(parents=True)
    (a.out / 'checks.json').write_text(json.dumps(dict(status='PASS_CPU_JOIN_RECEIPT_GUARDS',
        checks=checks, fixture=True, calls_solver=False, calls_CUDA=False,
        checker_sha256=sha(__file__), join_checker_sha256=sha(Path(__file__).with_name('line9_v401_dense_join.py')),
        source_execution_sha256=sha(a.source / 'execution.jsonl'),
        limits='Receipt fixtures only; actual17completion,combined45native and independent SFCW remain required.'), indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', checks=len(checks))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
