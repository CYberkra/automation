"""Line9 时窗收敛验证批：夜间41道批尾部最差3道的4800ns延长窗对照。

仅改 #time_window: 2.4e-06 -> 4.8e-06，其余与夜间批逐字节一致。
parent 包为夜间批准备包（含 manifest.json/cases/geometries）。
"""
import argparse
import copy
import json
from pathlib import Path
import shutil

import line9_v401_version_controls as native

APPROVAL = ('2026-10-12 用户明确指示「之前说的时窗收敛那个你也解决了」。'
            '对应方案：夜间41道批（契约ed9791f1…）尾部末5%峰超过官方−60dB建议的11道H1中，'
            '取最差3道（X170.0/170.5/171.0，−57.25/−57.30/−57.30dB）各跑一发4800ns延长窗H1，'
            '其余输入逐字节不变；用既有crop-vs-full加性检验量化2400ns截断对501音调复谱的影响。')

STATIONS = [170.0, 170.5, 171.0]
NEW_WINDOW_S = 4.8e-06
CROP_SAMPLES = 40703
_DT = 5.896635841874211e-11
NEW_SAMPLES = int(__import__('math').ceil(NEW_WINDOW_S/_DT))+1


def prepare(a):
    assert not a.out.exists()
    old = json.loads((a.parent/'manifest.json').read_text('utf-8'))
    assert old['status'] == 'PREPARED_APPROVED_OVERNIGHT_LINE'
    by_id = {g['id']: g for g in old['groups']}
    groups = []
    for x in STATIONS:
        src_id = f'line_x{round(x*100):05d}_H1'
        g = copy.deepcopy(by_id[src_id])
        name = f'window_x{round(x*100):05d}_H1'
        g.update(id=name, input=f'cases/{name}/profile.in', time_window_ns=4800.0,
                 expected_samples=NEW_SAMPLES, comparison_crop_samples=CROP_SAMPLES,
                 phase='window_convergence', parent_group_id=src_id)
        lines = []
        for line in (a.parent/by_id[src_id]['input']).read_text('utf-8').splitlines():
            if line.startswith('#title:'):
                line = f'#title: {name}; window convergence 4800ns; parent {src_id}'
            elif line.startswith('#time_window:'):
                line = '#time_window: 4.8e-06'
            lines.append(line)
        assert any(l == '#time_window: 4.8e-06' for l in lines)
        p = a.out/g['input']; p.parent.mkdir(parents=True)
        p.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        for key in ['geometry', 'material']:
            src = a.parent/g[key]; dst = a.out/g[key]
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists():
                shutil.copyfile(src, dst)
            assert native.sha(dst) == g[key+'_sha256']
        g['input_sha256'] = native.sha(p)
        groups.append(g)
    import math
    dt = old['dt_s']
    manifest = dict(status='PREPARED_APPROVED_WINDOW_CONVERGENCE', approval_basis=APPROVAL,
                    execution_authorized=True, groups=groups, no_retry=True,
                    parent_manifest_sha256=native.sha(a.parent/'manifest.json'),
                    generator_sha256=native.sha(__file__), dt_s=dt, SFCW=old['SFCW'],
                    stations=[s for s in old['stations'] if s['profile_x_m'] in STATIONS],
                    limits='3 stations, 4800ns window only; same geometry/materials/source/mesh as parent batch. '
                           'Truncation effect on the 501-tone complex spectrum is quantified by crop-vs-full; '
                           'not a physical convergence certification for other stations or field data.')
    native.save(a.out/'manifest.json', manifest)
    print(json.dumps(dict(groups=len(groups), expected_samples=NEW_SAMPLES)))


def freeze(a):
    manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert manifest['approval_basis'] == APPROVAL and manifest['execution_authorized']
    assert len(manifest['groups']) == 3
    a.action = 'prepare'; native.prepare(a)
    p = a.out/'execution_contract.json'; c = json.loads(p.read_text('utf-8'))
    c.update(approval_basis=APPROVAL, study_manifest=manifest, limits=manifest['limits'],
             max_batch_wall_s=10800, max_group_wall_s=2400)
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    native.save(p, c); native.save(a.out/'preflight_verification.json', native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p), groups=3, max_batch_wall_s=10800)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--parent', type=Path)
    p.add_argument('--package', type=Path)
    a = p.parse_args(); a.out = a.out.resolve()
    if a.action == 'prepare':
        a.parent = a.parent.resolve(); prepare(a)
    elif a.action == 'freeze':
        freeze(a)
    elif a.action == 'run':
        native.supervisor.audit = native.audit
        native.supervisor.run(a.out/'execution_contract.json')
    else:
        print(json.dumps(native.audit(a.out/'execution_contract.json', True)))
