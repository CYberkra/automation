"""Verify combined immutable dense-line executions with one declared consumed gap."""
import json
from pathlib import Path
from hs_capsule_identity import sha256 as sha
from line9_v401_dense_continuation import read, inventory, predecessor, audit, no_live_owned

def validate_join(old, old_events, new, events, verification):
    done, remaining = inventory(old, old_events)
    assert new['max_runs'] == len(new['groups']) == 17 and new['no_retry']
    assert [g['id'] for g in new['groups']] == [g['id'] for g in remaining]
    started = [e['group'] for e in events if e['status'] == 'STARTED' and e.get('group')]
    receipts = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    expected_new = {g['id'] for g in remaining}
    assert len(started) == len(set(started)) == 17 and set(started) == expected_new
    assert len(receipts) == len({e['group'] for e in receipts}) == 17
    assert events[-1].get('status') == 'COMPLETED' and events[-1].get('traces') == 17
    assert verification['completed']
    assert {e['group']: e['raw_sha256'] for e in receipts} == {g['id']: g['native_sha256'] for g in verification['groups']}
    assert len(verification['groups']) == 17
    joined = done + receipts
    assert len(joined) == len({e['group'] for e in joined}) == 45
    assert {e['group'] for e in joined} == {g['id'] for g in old['groups']} - {'low_x18675_H1'}
    return joined


def verify_combined(path, check_live=True):
    """Re-read both immutable lineages and native files; no fabricated completion."""
    folder = path.parent
    certificate = read(path)
    assert certificate['status'] == 'PASS_COMBINED45_NEW_PLUS2_REUSED_ONE_CONSUMED_MISSING'
    assert certificate['completed_groups'] == 45 and certificate['consumed_incomplete'] == ['low_x18675_H1']
    for name, digest in certificate['files'].items():
        assert Path(name).name == name and sha(folder / name) == digest
    cp = folder / 'continuation_execution_contract.json'
    n = read(cp)
    parent = Path(n['parent_contract'])
    assert sha(cp) == certificate['continuation_contract_sha256']
    assert sha(parent) == sha(folder / 'execution_contract.json') == certificate['original_contract_sha256']
    assert sha(parent.parent / 'execution.jsonl') == sha(folder / 'execution.jsonl')
    original, old_events, _, old_rows = predecessor(parent)
    assert sha(Path(n['parent_package']) / 'manifest.json') == certificate['original_manifest_sha256']
    execution = Path(n['groups'][0]['input']).parents[4]
    assert sha(execution / 'execution_contract.json') == sha(cp)
    assert sha(execution / 'execution.jsonl') == sha(folder / 'continuation_execution.jsonl')
    new_audit = audit(cp, True)
    assert new_audit == read(folder / 'continuation_completed_verification.json')
    events = [json.loads(s) for s in (folder / 'continuation_execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(cp)
    assert events[-1]['verification_sha256'] == sha(folder / 'continuation_completed_verification.json')
    receipts = validate_join(original, old_events, n, events, new_audit)
    receipt_map = {e['group']: e['raw_sha256'] for e in receipts}
    assert len(certificate['groups']) == 45
    assert {g['id']: g['native_sha256'] for g in certificate['groups']} == receipt_map
    for g in certificate['groups']:
        assert sha(Path(g['native_path'])) == sha(folder / (g['id'] + '.h5')) == g['native_sha256']
    for r in original['study_manifest']['reused']:
        assert sha(folder / (r['id'] + '.h5')) == r['native_sha256']
    actual = {p.stem for p in folder.glob('*.h5')}
    assert actual == set(receipt_map) | {r['id'] for r in original['study_manifest']['reused']}
    assert certificate['all_three_H0_anchors_present']
    assert all(f'{v}_x{int(x*100):05d}_H0' in actual for x in [190,187.5,185] for v in ['high','low'])
    if check_live:
        no_live_owned(n, execution)
    # The local independent exact-tone analysis must also be complete and bound.
    source_analysis = read(folder / 'analysis.json')
    source_audit = read(folder / 'sfcw_independent_audit.json')
    assert source_audit['status'] == 'PASS_INDEPENDENT_DENSE_NATIVE_DFT_INVERSE_MASKS'
    assert source_audit['completed_new'] == source_analysis['completed_new'] == 45
    assert source_audit['analysis_sha256'] == sha(folder / 'analysis.json')
    assert source_audit['contract_sha256'] == certificate['original_contract_sha256']
    assert source_audit['manifest_sha256'] == certificate['original_manifest_sha256']
    assert {r['id']: r['native_sha256'] for r in source_analysis['native']} == {name: sha(folder / (name + '.h5')) for name in actual}
    return dict(status='PASS_TERMINAL_COMBINED45_NATIVE_AND_SFCW_ONE_DECLARED_MISSING',
        contract=str(parent), contract_sha256=sha(parent), verification_sha256=sha(path),
        completed_groups=45, consumed_incomplete=['low_x18675_H1'],
        manifest_sha256=certificate['original_manifest_sha256'], baseline_native_sha256=receipt_map['high_x19000_H0'],
        gpu_lock=n['gpu_lock'],
        source_analysis_sha256=sha(folder / 'analysis.json'), source_audit_sha256=sha(folder / 'sfcw_independent_audit.json'))


