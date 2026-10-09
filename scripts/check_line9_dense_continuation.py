"""CPU rejection checks for the actual interrupted lineage; no solver or CUDA."""
import argparse
import copy
import json
from pathlib import Path

from hs_capsule_identity import sha256 as sha
from line9_v401_dense_continuation import inventory, read
from run_line9_ssh_lease_controller import probe_script


def main(a):
    assert not a.out.exists()
    c = read(a.source / 'execution_contract.json')
    e = [json.loads(s) for s in (a.source / 'execution.jsonl').read_text('utf-8').splitlines()]
    done, remaining = inventory(c, e)
    assert len(done) == 28 and len(remaining) == 17
    assert set(g['id'] for g in remaining).isdisjoint({r['group'] for r in e if r.get('group')})
    checks = ['Actual28complete/29consumed/17never-started lineage partitions without overlap']
    def reject(events):
        try:
            inventory(c, events)
        except AssertionError:
            return
        raise AssertionError('Invalid lineage accepted')
    reject(e[:-1]); checks.append('Nonterminal predecessor rejected')
    q = copy.deepcopy(e); q.append(q[-2]); reject(q); checks.append('Terminal status cannot be inferred from earlier failure')
    q = copy.deepcopy(e); q.insert(-1, next(r for r in q if r.get('group') and r['status'] == 'COMPLETED')); reject(q); checks.append('Duplicate completed receipt rejected')
    q = copy.deepcopy(e); q.insert(-1, next(r for r in q if r.get('group') and r['status'] == 'STARTED')); reject(q); checks.append('Repeated consumed attempt rejected')
    q = copy.deepcopy(e); q[-1]['error'] = 'other_failure'; reject(q); checks.append('Different failure requires separate continuation review')
    q = copy.deepcopy(e)
    for r in q:
        if r.get('group') == 'low_x18675_H1':
            r['group'] = 'high_x18500_H0'
    reject(q); checks.append('Missing anchor or different consumed trace rejected')
    quoted = probe_script("E:/case'quoted/execution", 'a' * 64)
    assert "$r='E:/case''quoted/execution'" in quoted and 'USER_STOP' not in quoted
    assert 'USER_STOP' in probe_script('E:/case/execution', 'a' * 64, True)
    checks.append('Lease controller quotes literal PS paths and separates renew/owned-stop actions')
    a.out.mkdir(parents=True)
    (a.out / 'checks.json').write_text(json.dumps(dict(status='PASS_INTERRUPTED_LINEAGE_AND_CONTROLLER_CPU_GUARDS',
        checks=checks, calls_solver=False, calls_CUDA=False,
        script_sha256=sha(__file__), continuation_sha256=sha(Path(__file__).with_name('line9_v401_dense_continuation.py')),
        controller_sha256=sha(Path(__file__).with_name('run_line9_ssh_lease_controller.py')),
        source_contract_sha256=sha(a.source / 'execution_contract.json'), source_execution_sha256=sha(a.source / 'execution.jsonl')), indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', checks=len(checks))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
