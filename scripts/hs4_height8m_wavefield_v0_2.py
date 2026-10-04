"""Continuation of the 8 m height wavefield study after a launcher-only failure.

The first capsule failed inside the VRAM probe before any solver process
started: the stripped launch environment lacked SYSTEMROOT/base PATH, so
nvidia-smi could not initialise NVML. No group was started and no solver
output exists. This continuation freezes byte-identical physics inputs under
a new contract whose only change is the launcher identity, and retains the
failed capsule untouched.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hs4_height8m_wavefield_v0_1 as v01

ROOT = v01.ROOT


def prepare_continuation(folder, previous):
    records = [json.loads(s) for s in (previous/'execution.jsonl').read_text('utf-8').splitlines()]
    if (len(records)!=2 or records[0].get('status')!='STARTED' or 'group' in records[0]
            or records[1].get('status')!='FAILED' or 'nvidia-smi' not in records[1].get('error','')):
        raise ValueError('previous capsule must show a launcher-only failure with no started group')
    if any(previous.glob('*/profile.h5')):
        raise ValueError('solver output exists in previous capsule')
    old_inputs = {p.parent.name: v01.sha256(p) for p in previous.glob('*/profile.in')}
    v01.prepare(folder)
    path = folder/'execution_contract.json'
    c = json.loads(path.read_text('utf-8'))
    for g in c['groups']:
        if old_inputs[g['id']] != g['input_sha256']:
            raise ValueError('continuation input differs from failed attempt')
    c['code_identities'] = {str(ROOT/name): v01.sha256(ROOT/name) for name in
        ['scripts/hs4_height8m_wavefield_v0_2.py','scripts/hs4_height8m_wavefield_v0_1.py',
         'scripts/run_hs4_height8m_wavefield_v0_2.cmd','scripts/gprmax_cached_cuda_entry.py']}
    c['previous_failed_capsule'] = str(previous.resolve())
    c['previous_contract_sha256'] = v01.sha256(previous/'execution_contract.json')
    c['launcher_fix'] = ('Launcher-only fix: runner cmd now seeds SYSTEMROOT and a base PATH '
        '(C:\\Windows\\System32;C:\\Windows) before vcvars, so nvidia-smi/NVML works under a '
        'stripped environment. Physics inputs byte-identical to the failed capsule; no solver ran there.')
    v01.save(path, c)
    v01.check(path, False)
    print('FrozenContinuation3runs;', c['snapshot_six_field_payload_GiB'], 'GiB snapshot payload per snapshot trace')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare','run','check'])
    parser.add_argument('--out', type=Path)
    parser.add_argument('--contract', type=Path)
    parser.add_argument('--previous', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.action=='prepare':
        prepare_continuation(args.out, args.previous)
    elif args.action=='run':
        if not args.execute: raise ValueError('--execute required')
        v01.run(args.contract)
    else:
        print(json.dumps(v01.check(args.contract, True), ensure_ascii=False)[:500])
