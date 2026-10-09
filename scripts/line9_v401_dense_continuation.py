"""Single-attempt continuation: reuse28, exclude the consumed incomplete trace."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import time
import zipfile

import h5py
import numpy as np
import psutil

from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as native
import hs4_station_grid_controls as supervisor


def read(p):
    return json.loads(Path(p).read_text('utf-8'))


def save(p, value):
    native.save(p, value)


def inventory(contract, events):
    assert contract['max_runs'] == len(contract['groups']) == 46 and contract['no_retry']
    assert events[-1] == dict(status='FAILED', error='SESSION_LEASE_EXPIRED')
    expected = [g['id'] for g in contract['groups']]
    started = [e['group'] for e in events if e['status'] == 'STARTED' and e.get('group')]
    done = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    assert len(expected) == len(set(expected)) and len(started) == len(set(started)) == 29
    assert len(done) == len({e['group'] for e in done}) == 28
    assert set(e['group'] for e in done) <= set(started) <= set(expected)
    missing = set(started) - {e['group'] for e in done}
    assert missing == {'low_x18675_H1'}
    remaining = [g for g in contract['groups'] if g['id'] not in started]
    assert len(remaining) == 17
    return done, remaining


def no_live_owned(contract, directory):
    roots = [str(Path(g['input']).parent).casefold() for g in contract['groups']]
    for process in psutil.process_iter(['name', 'cmdline']):
        if 'python' not in (process.info['name'] or '').lower():
            continue
        command = ' '.join(process.info['cmdline'] or []).casefold()
        assert not any(p in command for p in roots), 'Owned predecessor solver is alive'
        assert not (str(directory).casefold() in command and ' run ' in ' ' + command + ' '), 'Owned predecessor runner alive'


def native_record(g, raw, digest):
    assert sha(raw) == digest
    with h5py.File(raw) as h:
        assert h.attrs['gprMax'] == '4.0.1' and h.attrs['Iterations'] == g['expected_samples'] == 20352
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [8400, 1700, 1])
        np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025] * 3)
        assert abs(h.attrs['dt'] / g['dt_s'] - 1) < 1e-14
        for node, role in [('srcs/src1', 'tx'), ('rxs/rx1', 'rx')]:
            np.testing.assert_allclose(h[node].attrs['Position'], g[role + '_m'], rtol=0, atol=1e-11)
            np.testing.assert_array_equal(h[node].attrs['GridPosition'], np.rint(np.array(g[role + '_m']) / .025))
        for key in ['rxs/rx1/Ez', 'srcs/src1/excitation/samples']:
            x = h[key][:]
            assert x.dtype == np.float64 and x.shape == (20352,) and np.isfinite(x).all() and np.any(x)
        s = h['srcs/src1/excitation']
        assert s.attrs['WaveformType'] == 'ricker' and s.attrs['WaveformFrequency'] == 100e6
    return dict(id=g['id'], native_path=str(raw.resolve()), native_sha256=digest)


def predecessor(path):
    c = read(path)
    native.audit(path)
    events = [json.loads(s) for s in (path.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(path)
    done, remaining = inventory(c, events)
    no_live_owned(c, path.parent)
    groups = {g['id']: g for g in c['groups']}
    rows = [native_record(groups[e['group']], Path(groups[e['group']]['input']).with_suffix('.h5'), e['raw_sha256']) for e in done]
    return c, events, remaining, rows


def freeze(parent, parent_package, out):
    from audit_line9_dense_loss_inputs import audit as input_audit
    assert not out.exists()
    c, events, remaining, rows = predecessor(parent)
    m = read(parent_package / 'manifest.json')
    assert m == c['study_manifest']
    original_groups = {g['id']: g for g in m['groups']}
    remaining = [original_groups[g['id']] for g in remaining]
    inputs = input_audit(parent_package)
    assert inputs['status'].startswith('PASS') and inputs['bulk_execution_ready']
    package = out.parent / 'package'
    assert not package.exists()
    package.mkdir(parents=True)
    for g in remaining:
        for key in ['input', 'geometry', 'material']:
            src, dst = parent_package / g[key], package / g[key]
            assert not Path(g[key]).is_absolute() and not dst.exists()
            assert sha(src) == g[key + '_sha256']
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
    continuation = copy.deepcopy(m)
    continuation.update(groups=remaining, approval_basis='User: 继续;17 never-started traces only,retain28 complete and1 consumed incomplete,no retry.',
        continuation_of_contract_sha256=sha(parent), excluded_consumed=['low_x18675_H1'])
    save(package / 'manifest.json', continuation)
    native.prepare(argparse.Namespace(action='prepare', package=package, out=out))
    cp = out / 'execution_contract.json'
    n = read(cp)
    n.update(study_manifest=continuation, approval_basis=continuation['approval_basis'], max_batch_wall_s=7200,
        parent_contract=str(parent), parent_package=str(parent_package), parent_contract_sha256=sha(parent),
        parent_execution_sha256=sha(parent.parent / 'execution.jsonl'), predecessor_native=rows,
        consumed_incomplete=['low_x18675_H1'], lease_authority='Bounded local SSH controller;600s remote expiry,7200s controller and solver batch ceilings;USER_STOP retained.',
        gpu_lock=str(Path('E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock')))
    for p in [Path(__file__), Path(__file__).with_name('audit_line9_dense_loss_inputs.py')]:
        n['code_identities'][str(p.resolve())] = sha(p)
    for p in [parent, parent.parent / 'execution.jsonl', parent_package / 'manifest.json']:
        n['file_identities'][str(p.resolve())] = sha(p)
    save(cp, n)
    save(out / 'predecessor_partial_audit.json', dict(status='PASS28_RECEIPTS_NATIVE_AND17_NEVER_STARTED',
        parent_contract_sha256=sha(parent), execution_sha256=sha(parent.parent / 'execution.jsonl'), groups=rows,
        never_started=[g['id'] for g in remaining], consumed_incomplete=['low_x18675_H1'], input_audit=inputs))
    save(out / 'preflight_verification.json', audit(cp))
    print(json.dumps(dict(status='FROZEN17_NEVER_STARTED_NO_RETRY', contract_sha256=sha(cp), out=str(out))))


def audit(path, completed=False):
    result = native.audit(path, completed)
    c = read(path)
    assert c['max_runs'] == len(c['groups']) == 17 and c['no_retry'] and c['max_lease_age_s'] == 600
    assert c['consumed_incomplete'] == ['low_x18675_H1']
    parent = Path(c['parent_contract'])
    assert sha(parent) == c['parent_contract_sha256'] and sha(parent.parent / 'execution.jsonl') == c['parent_execution_sha256']
    old = read(parent)
    events = [json.loads(s) for s in (parent.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
    _, remaining = inventory(old, events)
    assert [g['id'] for g in remaining] == [g['id'] for g in c['groups']]
    return result


def combine(execution, out):
    assert not out.exists()
    cp = execution / 'execution_contract.json'
    audit_new = audit(cp, True)
    assert audit_new == read(execution / 'completed_verification.json')
    n = read(cp)
    parent = Path(n['parent_contract'])
    old, old_events, _, old_rows = predecessor(parent)
    events = [json.loads(s) for s in (execution / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(cp)
    assert events[-1] == dict(status='COMPLETED', traces=17, verification_sha256=sha(execution / 'completed_verification.json'))
    receipts = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    assert len(receipts) == len({e['group'] for e in receipts}) == 17
    assert {e['group']: e['raw_sha256'] for e in receipts} == {g['id']: g['native_sha256'] for g in audit_new['groups']}
    no_live_owned(n, execution)
    rows = old_rows + [dict(id=g['id'], native_path=g['native_path'], native_sha256=g['native_sha256']) for g in audit_new['groups']]
    expected = {g['id'] for g in old['groups']} - {'low_x18675_H1'}
    assert len(rows) == len({r['id'] for r in rows}) == 45 and {r['id'] for r in rows} == expected
    out.mkdir(parents=True)
    for r in rows:
        shutil.copyfile(r['native_path'], out / (r['id'] + '.h5'))
    for r in old['study_manifest']['reused']:
        raw = Path(n['parent_package']) / r['path']
        assert sha(raw) == r['native_sha256']
        shutil.copyfile(raw, out / (r['id'] + '.h5'))
    for name in ['execution_contract.json', 'execution.jsonl', 'preflight_verification.json']:
        shutil.copyfile(parent.parent / name, out / name)
    for name in ['execution_contract.json', 'execution.jsonl', 'preflight_verification.json', 'completed_verification.json', 'predecessor_partial_audit.json']:
        shutil.copyfile(execution / name, out / ('continuation_' + name))
    save(out / 'snapshot_manifest.json', dict(status='PARTIAL', contract_sha256=sha(parent), completed_new=45,
        requested_new=46, reused=2, records=[dict(id=r['id'], native_sha256=r['native_sha256']) for r in rows],
        consumed_incomplete=['low_x18675_H1'], terminal_event_observed='OLD_FAILED_NEW_COMPLETED17',
        limits='Combined audited output,not completed46. Missing consumed low18675 remains missing. No synthesized receipt or retry.'))
    certificate = dict(status='PASS_COMBINED45_NEW_PLUS2_REUSED_ONE_CONSUMED_MISSING', original_contract_sha256=sha(parent),
        original_manifest_sha256=sha(Path(n['parent_package']) / 'manifest.json'), continuation_contract_sha256=sha(cp),
        consumed_incomplete=['low_x18675_H1'], groups=rows, completed_groups=45,
        files={p.name: sha(p) for p in out.iterdir() if p.is_file()},
        all_three_H0_anchors_present=all(f'{v}_x{int(x*100):05d}_H0' in {r['id'] for r in rows + old['study_manifest']['reused']} for x in [190,187.5,185] for v in ['high','low']),
        limits='Native receipt consistency only,not physical/path/site certification. Full46 completion never asserted.')
    assert certificate['all_three_H0_anchors_present']
    save(out / 'combined_certificate.json', certificate)
    zip_path = out.with_suffix('.zip')
    with zipfile.ZipFile(zip_path, 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.iterdir()):
            z.write(p, p.name)
    print(json.dumps(dict(path=str(zip_path), sha256=sha(zip_path), groups=45, reused=2, missing=1)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze', 'run', 'verify', 'combine'])
    for key in ['out', 'parent', 'parent-package', 'execution']:
        p.add_argument('--' + key, type=Path, required=key == 'out')
    a = p.parse_args()
    a.out = a.out.resolve()
    if a.action == 'freeze':
        freeze(a.parent.resolve(), a.parent_package.resolve(), a.out)
    elif a.action == 'run':
        supervisor.audit = audit
        supervisor.run(a.out / 'execution_contract.json')
    elif a.action == 'verify':
        print(json.dumps(audit(a.out / 'execution_contract.json', True)))
    else:
        combine(a.execution.resolve(), a.out)
