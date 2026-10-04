"""CPU checks of portable geometry and resource/identity guards; no solver."""
import argparse
import json
import os
from pathlib import Path
import tempfile
from unittest.mock import patch
import numpy as np

import hs4_cross_pc as task
from hs_capsule_identity import sha256


def rejected(action):
    try:
        action()
    except (ValueError, RuntimeError):
        return
    raise AssertionError('invalid request accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--design',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.out.exists(): raise ValueError('new test output required')
    design=task.verify_design(a.design); runtime=task.version_runtime(); results=[]
    base=task.ROOT/'artifacts/local_checks'; base.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='hs4_portable_cpu_',dir=base) as name:
        tmp=Path(name); regenerated=tmp/'design'; task.prepare(regenerated)
        other=task.verify_design(regenerated)
        assert design['stages']==other['stages']; results.append('relative design regenerates exactly')
        for stage in task.STAGES:
            for g in design['stages'][stage]: assert not Path(g['input']).is_absolute()
        results.append('all20 input locations portable')
        for g in design['stages']['3d-centre']+design['stages']['3d-ends']:
            assert g['tx_m'][0]==1.6 and g['rx_m'][0]==1.6
        results.append('3D actual grid X retained rather than nominal1.625')
        normal=task.commands(a.design/design['stages']['3d-centre'][0]['input'])['#box:'][1:]
        wide=task.commands(a.design/design['stages']['3d-xwide-centre'][0]['input'])['#box:'][1:]
        assert len(normal)==2304 and len(wide)==2400
        for old,new in zip(normal,wide[:2304]):
            old=np.asarray(list(map(float,old[:6]))); new=np.asarray(list(map(float,new[:6])))
            assert np.array_equal(old+np.array([2,0,0,2,0,0]),new)
        results.append('Xwide full original interior preserved by translation')
        rejected(lambda:task.prepare(regenerated)); results.append('existing design refuses overwrite')
        original=a.design/design['stages']['replay'][0]['input']
        text=original.read_text('utf-8').replace('17.6','17.601')
        rejected(lambda:task.case_record(tmp,'bad','offgrid',text,task.FINE/'base/profile.in'))
        results.append('off-grid active source coordinate rejected')
        fixture=tmp/'fixture'; fixture.mkdir(); groups=[]
        for g in design['stages']['replay']:
            g=dict(g); directory=fixture/g['id']; directory.mkdir()
            (directory/'profile.in').write_bytes((a.design/g['input']).read_bytes())
            g['input']=g['id']+'/profile.in'; groups.append(g)
            for filename in ('profile.h5','hs4t2d_geom.vtkhdf'):
                os.link(task.FINE/g['id']/filename,directory/filename)
        task.save(fixture/'execution_contract.json',{'stage':'replay','groups':groups,'runtime':runtime})
        verified=task.audit(fixture); assert len(verified['groups'])==4
        results.append('all4 archived float64 fields, poses, CFL and entire material maps independently audited')
        clean=tmp/'resource_guard'; clean.mkdir(); g=dict(groups[0]); (clean/g['id']).mkdir()
        (clean/g['input']).write_bytes((fixture/g['input']).read_bytes())
        c={'status':'FROZEN_APPROVED','stage':'replay','groups':[g],'max_runs':1,'runtime':runtime,
           'code_identities':{'scripts/hs4_cross_pc.py':sha256(task.__file__)},
           'min_available_RAM_GiB':2.25,'min_free_VRAM_GiB':3.5}
        task.save(clean/'execution_contract.json',c)
        with patch.object(task,'version_runtime',return_value=runtime),patch.object(task,'resources',return_value={'available_RAM_bytes':0,'free_VRAM_bytes':0}),patch.object(task.subprocess,'Popen') as launch:
            rejected(lambda:task.run(clean)); assert not launch.called and not (clean/'execution.jsonl').exists()
        results.append('capacity rejection starts no solver and consumes no attempt')
        (clean/g['input']).write_text('changed',encoding='utf-8')
        with patch.object(task,'version_runtime',return_value=runtime),patch.object(task.subprocess,'Popen') as launch:
            rejected(lambda:task.run(clean)); assert not launch.called
        results.append('frozen input tampering rejected before launch')
        with patch.object(task,'version_runtime',return_value=runtime):
            rejected(lambda:task.freeze(a.design,'centre1cm',tmp/'unreplayed',None))
        assert not (tmp/'unreplayed').exists(); results.append('target replay mandatory before new stage')
    task.save(a.out,{'status':'PASS','solver_executed':False,'checks':results,'count':len(results),
                    'code_sha256':sha256(__file__),'runner_sha256':sha256(task.__file__),
                    'design_sha256':sha256(a.design/'task.json'),
                    'scope':'CPU portable-design and execution-guard checks using archived4 raw controls; not a new replay, simulation or physical validation.'})
    print(json.dumps({'status':'PASS','count':len(results)},indent=2))


if __name__=='__main__': main()
