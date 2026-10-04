"""Bounded continuation contract for a host-interrupted HS4 cross-PC attempt.

Per docs/research/2026-10-04_hs4_cross_pc_task.md: once any group has STARTED,
preserve the original capsule, analyse the cause, and freeze a new bounded
contract for only the remaining groups; never delete logs/raws to reopen an
attempt and never rerun completed groups. This script implements exactly that:
freeze-continuation creates the bounded contract for the unfinished groups of
an interrupted capsule; assemble-analysis merges the completed prior groups
with the continuation groups into one byte-verified analysis capsule and
re-audits every raw file independently.
"""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psutil

from hs_capsule_identity import sha256
from hs4_cross_pc import ROOT, audit, read, save

CHECKS = ROOT / 'artifacts/research_checks'
TASK_DOC = 'docs/research/2026-10-04_hs4_cross_pc_task.md'


def events(capsule):
    log = capsule / 'execution.jsonl'
    return [json.loads(line) for line in log.read_text('utf-8').splitlines()] if log.exists() else []


def interruption(capsule):
    contract = read(capsule / 'execution_contract.json')
    if contract['status'] != 'FROZEN_APPROVED':
        raise ValueError('prior capsule is not a frozen approved contract')
    log = events(capsule)
    if not log or log[0]['status'] != 'STARTED':
        raise ValueError('prior capsule has no started attempt')
    if any(e['status'] == 'FAILED' for e in log) or any(
            e['status'] == 'COMPLETED' and 'group' not in e for e in log):
        raise ValueError('prior attempt finished or failed; no continuation needed/allowed')
    started = {e['group'] for e in log if e['status'] == 'STARTED' and 'group' in e}
    completed = {e['group']: e for e in log if e['status'] == 'COMPLETED' and 'group' in e}
    if not started - set(completed):
        raise ValueError('no dangling group found; attempt was not interrupted mid-group')
    for gid, event in completed.items():
        group = next(g for g in contract['groups'] if g['id'] == gid)
        raw = capsule / group['input']
        if sha256(raw.with_suffix('.h5')) != event['raw_sha256']:
            raise ValueError('completed raw differs from its execution log record')
    remaining = [g for g in contract['groups'] if g['id'] not in completed]
    return contract, log, completed, remaining


def freeze_continuation(prior, out):
    if out.exists():
        raise ValueError('new continuation capsule required')
    contract, log, completed, remaining = interruption(prior)
    for process in psutil.process_iter(['cmdline']):
        cmd = ' '.join(process.info['cmdline'] or [])
        if 'gprmax_cached_cuda_entry' in cmd:
            raise ValueError('a solver from the interrupted attempt is still alive: ' + cmd)
    out.mkdir(parents=True)
    groups = []
    for group in remaining:
        target = out / group['input']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(prior / group['input'], target)
        if sha256(target) != group['input_sha256']:
            raise ValueError('copied input differs from frozen hash: ' + group['id'])
        groups.append(group)
    new = copy.deepcopy(contract)
    new['stage'] = contract['stage'] + '-continuation'
    new['groups'] = groups
    new['max_runs'] = len(groups)
    new['max_batch_wall_s'] = contract['max_group_wall_s'] * len(groups)
    new['continuation'] = {
        'reason': 'Host application exited; solver child was terminated mid-group. '
                  'Not a solver or physics failure.',
        'basis': TASK_DOC + ': interrupted attempt is preserved; a bounded contract '
                 'covers only the remaining groups; completed groups are never rerun.',
        'prior_capsule': prior.name,
        'prior_contract_sha256': sha256(prior / 'execution_contract.json'),
        'prior_execution_log_sha256': sha256(prior / 'execution.jsonl'),
        'completed_in_prior_capsule': [
            {'id': gid, 'raw_sha256': e['raw_sha256'], 'elapsed_s': e['elapsed_s'],
             'peak_owned_RSS_bytes': e['peak_owned_RSS_bytes']}
            for gid, e in sorted(completed.items())],
        'remaining_group_ids': [g['id'] for g in groups]}
    save(out / 'execution_contract.json', new)
    print('Frozen continuation:', out.name, len(groups), 'traces:',
          ', '.join(g['id'] for g in groups))


def assemble_analysis(prior, continuation, out):
    if out.exists():
        raise ValueError('new analysis capsule required')
    prior_contract, log, completed, _ = interruption(prior)
    cont_contract = read(continuation / 'execution_contract.json')
    cont_check = read(continuation / 'completed_verification.json')
    if cont_check['status'] != 'PASS' or \
            cont_check['contract_sha256'] != sha256(continuation / 'execution_contract.json'):
        raise ValueError('continuation capsule is not a verified complete capsule')
    cont_records = {g['id']: g for g in cont_check['groups']}
    out.mkdir(parents=True)
    sources = {}
    for group in prior_contract['groups']:
        gid = group['id']
        origin, record = (prior, completed[gid]) if gid in completed else (continuation, cont_records.get(gid))
        if record is None:
            raise ValueError('group has no verified raw in either capsule: ' + gid)
        target = out / group['input']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(origin / group['input'], target)
        raw = target.with_suffix('.h5')
        shutil.copyfile((origin / group['input']).with_suffix('.h5'), raw)
        if sha256(raw) != record['raw_sha256']:
            raise ValueError('merged raw differs from its execution record: ' + gid)
        geom_name = 'hs4_geom.vtkhdf' if group['component'] == 'Ex' else 'hs4t2d_geom.vtkhdf'
        shutil.copyfile((origin / group['input']).parent / geom_name, target.parent / geom_name)
        sources[gid] = {'capsule': origin.name, 'raw_sha256': record['raw_sha256']}
    merged = copy.deepcopy(prior_contract)
    merged['assembly'] = {
        'reason': 'Analysis view over byte-identical copies from the interrupted capsule '
                  'and its bounded continuation; every raw re-verified by audit().',
        'group_sources': sources,
        'prior_capsule': prior.name,
        'continuation_capsule': continuation.name,
        'prior_execution_log_sha256': sha256(prior / 'execution.jsonl'),
        'continuation_execution_log_sha256': sha256(continuation / 'execution.jsonl'),
        'continuation_verification_sha256': sha256(continuation / 'completed_verification.json')}
    save(out / 'execution_contract.json', merged)
    verification = audit(out)
    save(out / 'completed_verification.json', verification)
    print('Assembled analysis capsule:', out.name, 'audit', verification['status'],
          len(verification['groups']), 'groups')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze-continuation', 'assemble-analysis'))
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--continuation', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'freeze-continuation':
        freeze_continuation(args.prior, args.out)
    else:
        if args.continuation is None:
            raise ValueError('--continuation required')
        assemble_analysis(args.prior, args.continuation, args.out)


if __name__ == '__main__':
    main()
