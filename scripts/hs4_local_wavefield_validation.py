"""Fresh local replay of the archived 8 m pair; no historical attempt reuse.

Use the existing owned-process supervisor and native/snapshot/passive audit.
The local replay supplies missing frames and measures cross-build differences.
"""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys

import gprMax
import hs4_height8m_wavefield_v0_1 as wavefield
import hs4_station_grid_controls as supervisor
from hs_capsule_identity import sha256

ROOT = wavefield.ROOT
OLD = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
REFERENCE = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre'


def freeze(out):
    if out.exists():
        raise ValueError('new local capsule required')
    if str(gprMax.__version__) != '4.0.0':
        raise ValueError('V4.0.0 required')
    old = json.loads((OLD/'execution_contract.json').read_text('utf-8'))
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    identities = {name: sha256(package/name) for name in old['source_identities']}
    changed = [name for name in identities if identities[name] != old['source_identities'][name]]
    if any(name.endswith('.py') for name in changed):
        raise ValueError('reviewed Python source differs')
    c = copy.deepcopy(old)
    c.pop('wheel', None); c.pop('wheel_sha256', None)
    c['python'] = str(Path(sys.executable).resolve())
    c['source_identities'] = identities
    c['source_identity_note'] = 'Actual local package hashes frozen; Python source identical to reference.'
    c['changed_binary_files_from_reference'] = changed
    c['reference_contract_sha256'] = sha256(OLD/'execution_contract.json')
    c['reference_capsule'] = str(OLD)
    c['approval_basis'] = 'User: 你本地跑个仿真验证一下. Fresh 3-case local V4/double replay; same physical inputs, new attempt and local snapshot archive.'
    c['min_available_RAM_GiB'] = 1.25
    c['max_batch_wall_s'] = 1800
    c['min_system_available_during_run_GiB'] = .2
    c['code_identities'] = {str(ROOT/name): sha256(ROOT/name) for name in (
        'scripts/hs4_local_wavefield_validation.py', 'scripts/run_hs4_local_wavefield_validation.cmd',
        'scripts/hs4_station_grid_controls.py', 'scripts/hs4_height8m_wavefield_v0_1.py',
        'scripts/gprmax_cached_cuda_entry.py')}
    c['processing'] = 'Ricker95MHz raw wavefield; per-frequency60-140MHz spatial diagnostic. Receivers separately use reviewed20-170MHz501-tone chain. No power attribution from Ey squared alone.'
    out.mkdir(parents=True)
    for g in c['groups']:
        folder = out/g['id']; folder.mkdir()
        p = folder/'profile.in'
        p.write_bytes((OLD/g['id']/'profile.in').read_bytes())
        if sha256(p) != g['input_sha256']:
            raise ValueError('replay input not byte-identical')
        g['input'] = str(p.resolve())
        ref = REFERENCE/('centre_'+g['role'])
        g['reference_input'] = str(ref/'profile.in')
        g['reference_geometry'] = str(ref/'hs4t2d_geom.vtkhdf')
    wavefield.save(out/'execution_contract.json', c)
    wavefield.check(out/'execution_contract.json', False)
    print('Frozen3localcases; changed binaries', len(changed), '; snapshots', len(c['snapshot_iterations']))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['freeze', 'run', 'verify'])
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    path = args.out/'execution_contract.json'
    if args.action == 'freeze':
        freeze(args.out)
    elif args.action == 'verify':
        print(json.dumps(wavefield.check(path, True), ensure_ascii=False)[:200])
    else:
        if not args.execute:
            raise ValueError('--execute required')
        # This study has wavefield probes/snapshots. Keep the common supervisor,
        # substituting only the study-specific audit, not its resource guards.
        supervisor.audit = wavefield.check
        supervisor.run(path)


if __name__ == '__main__':
    main()
