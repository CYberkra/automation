"""Single bounded continuation of the x190 AGL5 material ablation batch.

Parent attempt: 8 groups frozen under contract 3d174683..., 1 completed
(abl_base_x190_H0), 1 consumed with zero output (abl_base_x190_H1, killed by
SESSION_LEASE_EXPIRED before any .h5 bytes), 6 never started.

Policy recorded in the continuation contract: the six never-started groups run
unchanged; abl_base_x190_H1 is resumed once because its consumed attempt
produced zero output bytes (no result selection is possible) and the approved
deliverable requires the base H0/H1 reference pair. This is disclosed as a
deviation from the strict never-rerun-a-started-attempt precedent; the full
parent event log is preserved and hash-pinned.
"""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import psutil

from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as native
import hs4_station_grid_controls as supervisor

SOURCE = dict(source_type='impulse', source_frequency_Hz=1., source_amplitude_A=1.)
DONE = ['abl_base_x190_H0']
RESUMED_ZERO_OUTPUT = ['abl_base_x190_H1']
APPROVAL = ('2026-10-11 continuation of the approved 8-group material ablation after '
            'SESSION_LEASE_EXPIRED: 6 never-started groups unchanged; abl_base_x190_H1 '
            'resumed once (consumed attempt wrote zero output bytes; no result selection '
            'possible; base H0/H1 pair is the approved reference). No other change.')


def read(p):
    return json.loads(Path(p).read_text('utf-8'))


def inventory(contract, events):
    assert contract['max_runs'] == len(contract['groups']) == 8 and contract['no_retry']
    assert events[-1] == dict(status='FAILED', error='SESSION_LEASE_EXPIRED')
    expected = [g['id'] for g in contract['groups']]
    started = [e['group'] for e in events if e['status'] == 'STARTED' and e.get('group')]
    done = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    assert len(expected) == len(set(expected)) == 8
    assert set(started) == set(DONE) | set(RESUMED_ZERO_OUTPUT) and len(started) == len(set(started))
    assert [e['group'] for e in done] == DONE
    missing = set(started) - {e['group'] for e in done}
    assert missing == set(RESUMED_ZERO_OUTPUT)
    remaining = [g for g in contract['groups'] if g['id'] not in DONE]
    assert len(remaining) == 7
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
        assert h.attrs['gprMax'] == '4.0.1' and h.attrs['Iterations'] == g['expected_samples'] == 40703
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [8400, 1700, 1])
        np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025] * 3)
        assert abs(h.attrs['dt'] / g['dt_s'] - 1) < 1e-14
        for node, role in [('srcs/src1', 'tx'), ('rxs/rx1', 'rx')]:
            np.testing.assert_allclose(h[node].attrs['Position'], g[role + '_m'], rtol=0, atol=1e-11)
        x = h['rxs/rx1/Ez'][:]
        assert x.dtype == np.float64 and x.shape == (40703,) and np.isfinite(x).all() and np.any(x)
        exc = h['srcs/src1/excitation']
        merged = dict(g, **SOURCE)
        native.audit_source(exc, merged, float(h.attrs['dt']))
    return dict(id=g['id'], native_path=str(raw.resolve()), native_sha256=digest)


def predecessor(path):
    c = read(path)
    native.audit(path)
    events = [json.loads(s) for s in (path.parent / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(path)
    done, remaining = inventory(c, events)
    no_live_owned(c, path.parent)
    groups = {g['id']: g for g in c['groups']}
    for g in remaining:
        assert not Path(g['input']).with_suffix('.h5').exists(), 'Residual output blocks continuation'
    rows = [native_record(groups[e['group']], Path(groups[e['group']]['input']).with_suffix('.h5'), e['raw_sha256']) for e in done]
    return c, events, remaining, rows


def freeze(parent, parent_package, out):
    assert not out.exists()
    c, events, remaining, rows = predecessor(parent)
    m = read(parent_package / 'manifest.json')
    assert m == c['study_manifest']
    original_groups = {g['id']: g for g in m['groups']}
    remaining = [dict(original_groups[g['id']], **SOURCE) for g in remaining]
    package = out.parent / 'package_cont'
    assert not package.exists()
    package.mkdir(parents=True)
    for g in remaining:
        for key in ['input', 'geometry', 'material']:
            src, dst = parent_package / g[key], package / g[key]
            assert not Path(g[key]).is_absolute() and not dst.exists()
            assert sha(src) == g[key + '_sha256']
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
        line = next(s for s in (package / g['input']).read_text('utf-8').splitlines() if s.startswith('#waveform:'))
        assert line == '#waveform: impulse 1 1 impulse'
    continuation = copy.deepcopy(m)
    continuation.update(groups=remaining, approval_basis=APPROVAL,
        continuation_of_contract_sha256=sha(parent),
        resumed_zero_output=RESUMED_ZERO_OUTPUT,
        resumed_zero_output_justification='Consumed attempt killed by SESSION_LEASE_EXPIRED before any output bytes; no result selection possible; approved deliverable requires the base H0/H1 reference pair. Deviation from the strict started-attempt exclusion precedent is disclosed.',
        never_started=[g['id'] for g in remaining if g['id'] not in RESUMED_ZERO_OUTPUT])
    native.save(package / 'manifest.json', continuation)
    native.prepare(argparse.Namespace(action='prepare', package=package, out=out))
    cp = out / 'execution_contract.json'
    n = read(cp)
    n.update(study_manifest=continuation, approval_basis=APPROVAL, limits=continuation['limits'],
        max_batch_wall_s=7200,
        parent_contract=str(parent), parent_package=str(parent_package),
        parent_contract_sha256=sha(parent), parent_execution_sha256=sha(parent.parent / 'execution.jsonl'),
        predecessor_native=rows, resumed_zero_output=RESUMED_ZERO_OUTPUT,
        lease_authority='Bounded local SSH controller; 600s remote expiry, 7200s controller and solver batch ceilings; USER_STOP retained.')
    for p in [Path(__file__), Path(__file__).with_name('line9_material_ablation.py')]:
        n['code_identities'][str(p.resolve())] = sha(p)
    for p in [parent, parent.parent / 'execution.jsonl', parent_package / 'manifest.json']:
        n['file_identities'][str(p.resolve())] = sha(p)
    native.save(cp, n)
    native.save(out / 'predecessor_partial_audit.json', dict(
        status='PASS1_RECEIPT_NATIVE_AND7_CONTINUATION', parent_contract_sha256=sha(parent),
        execution_sha256=sha(parent.parent / 'execution.jsonl'), groups=rows,
        resumed_zero_output=RESUMED_ZERO_OUTPUT,
        never_started=[g['id'] for g in remaining if g['id'] not in RESUMED_ZERO_OUTPUT]))
    native.save(out / 'preflight_verification.json', audit(cp))
    print(json.dumps(dict(status='FROZEN_CONTINUATION7_NO_RETRY', contract_sha256=sha(cp), out=str(out))))


def audit(path, completed=False):
    result = native.audit(path, completed)
    c = read(path)
    assert c['max_runs'] == len(c['groups']) == 7 and c['no_retry'] and c['max_lease_age_s'] == 600
    assert c['resumed_zero_output'] == RESUMED_ZERO_OUTPUT
    assert all(g['source_type'] == 'impulse' and g['source_frequency_Hz'] == 1. and g['source_amplitude_A'] == 1. for g in c['groups'])
    parent = Path(c['parent_contract'])
    assert sha(parent) == c['parent_contract_sha256']
    assert sha(parent.parent / 'execution.jsonl') == c['parent_execution_sha256']
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
    old, _, _, old_rows = predecessor(parent)
    events = [json.loads(s) for s in (execution / 'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[0]['contract_sha256'] == sha(cp)
    assert events[-1]['status'] == 'COMPLETED' and events[-1]['traces'] == 7
    receipts = [e for e in events if e['status'] == 'COMPLETED' and e.get('group')]
    assert len(receipts) == len({e['group'] for e in receipts}) == 7
    assert {e['group']: e['raw_sha256'] for e in receipts} == {g['id']: g['native_sha256'] for g in audit_new['groups']}
    no_live_owned(n, execution)
    rows = old_rows + [dict(id=g['id'], native_path=g['native_path'], native_sha256=g['native_sha256']) for g in audit_new['groups']]
    expected = [g['id'] for g in old['groups']]
    assert len(rows) == len({r['id'] for r in rows}) == 8 and {r['id'] for r in rows} == set(expected)
    out.mkdir(parents=True)
    order = {g['id']: g for g in old['study_manifest']['groups']}
    merged_groups = []
    merged_rows = []
    for gid in expected:
        r = next(r for r in rows if r['id'] == gid)
        dst = out / (gid + '.h5')
        shutil.copyfile(r['native_path'], dst)
        assert sha(dst) == r['native_sha256']
        merged_groups.append(dict(order[gid], **SOURCE))
        merged_rows.append(dict(id=gid, native_path=str(dst.resolve()), native_sha256=r['native_sha256'],
                                dtype='float64', samples=40703, dt_s=old['study_manifest']['dt_s']))
    merged_contract = dict(status='COMBINED_AUDITED_VIEW', expected_version='4.0.1',
        approval_basis=old['approval_basis'] + ' | ' + APPROVAL,
        study_manifest=old['study_manifest'], groups=merged_groups, max_runs=8, no_retry=True,
        parent_contract_sha256=n['parent_contract_sha256'], continuation_contract_sha256=sha(cp),
        limits=old['limits'])
    native.save(out / 'execution_contract.json', merged_contract)
    verification = dict(status='PASS_COMBINED8_NATIVE_IDENTITY_NOT_PHYSICAL_CERTIFICATION', completed=True,
        contract_sha256=sha(out / 'execution_contract.json'), version='4.0.1', groups=merged_rows)
    native.save(out / 'completed_verification.json', verification)
    for name in ['execution_contract.json', 'execution.jsonl', 'preflight_verification.json']:
        shutil.copyfile(parent.parent / name, out / ('parent_' + name))
    for name in ['execution_contract.json', 'execution.jsonl', 'preflight_verification.json',
                 'completed_verification.json', 'predecessor_partial_audit.json']:
        shutil.copyfile(execution / name, out / ('continuation_' + name))
    native.save(out / 'combined_certificate.json', dict(
        status='PASS_COMBINED8_ONE_RESUMED_ZERO_OUTPUT', original_contract_sha256=n['parent_contract_sha256'],
        continuation_contract_sha256=sha(cp), groups=merged_rows, completed_groups=8,
        resumed_zero_output=RESUMED_ZERO_OUTPUT,
        files={p.name: sha(p) for p in out.iterdir() if p.is_file()},
        limits='Native receipt consistency only, not physical/path/site certification.'))
    print(json.dumps(dict(status='COMBINED8', out=str(out), groups=8,
                          certificate_sha256=sha(out / 'combined_certificate.json'))))


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
