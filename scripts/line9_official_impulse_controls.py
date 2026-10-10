"""Freeze/run the user-authorized, two-case official-impulse comparison."""
import argparse
import json
from pathlib import Path

import line9_v401_version_controls as native


def freeze(a):
    manifest_path = a.package / 'manifest.json'
    m = json.loads(manifest_path.read_text('utf-8'))
    assert [g['id'] for g in m['groups']] == ['high_x19000_H0', 'high_x19000_H1']
    for g in m['groups']:
        assert g['source_type'] == 'impulse' and g['source_amplitude_A'] == 1
        assert g['source_start_s'] == 0 and g['source_frequency_slot_unused']
        text = (a.package / g['input']).read_text('utf-8')
        assert '#waveform: impulse 1 1 impulse' in text
    a.action = 'prepare'
    native.prepare(a)
    p = a.out / 'execution_contract.json'
    c = json.loads(p.read_text('utf-8'))
    c['approval_basis'] = '用户：重新用新的去开始正演；官方 impulse，190 m 同站位 H0/H1 两项。'
    c['study_manifest']['preparation_approval_basis'] = m['approval_basis']
    c['study_manifest'].update(approval_basis=c['approval_basis'], execution_authorized=True,
        status='FROZEN_TWO_OFFICIAL_IMPULSE_CASES', preparation_manifest_sha256=native.sha(manifest_path))
    c['limits'] = 'Two cases only; unchanged geometry/material/height/grid; FP64; no retry; no snapshots or extra observers. Official direct SFCW 20-170 MHz, 0.3 MHz, 501 tones; Hann/Blackman, pad8, no tail taper or time shift.'
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    native.save(p, c)
    native.save(a.out / 'preflight_verification.json', native.audit(p))
    print(json.dumps({'status': c['status'], 'groups': len(c['groups']), 'contract_sha256': native.sha(p)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'run', 'verify'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--package', type=Path)
    args = parser.parse_args()
    args.out = args.out.resolve()
    if args.action == 'freeze':
        args.package = args.package.resolve()
        freeze(args)
    elif args.action == 'run':
        native.supervisor.audit = native.audit
        native.supervisor.run(args.out / 'execution_contract.json')
    else:
        print(json.dumps(native.audit(args.out / 'execution_contract.json', True)))
