"""Independent physical invariants and sampling audit; no solver invocation."""
import argparse
import json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha


def check(package):
    m=json.loads((package/'manifest.json').read_text('utf-8'));rows=[]
    dt=m['dt_s'];assert m['execution_authorized'] and m['no_retry'] and len(m['groups'])==6
    assert m['SFCW']['method']=='direct' and m['fields']==['Ez','Hx','Hy']
    for g in m['groups']:
        lines=(package/g['input']).read_text('utf-8').splitlines()
        base=(package/g['baseline']/'profile.in').read_text('utf-8').splitlines()
        stripped=lambda seq:[s for s in seq if s.startswith('#') and not s.startswith(('#title:','#rx:','#snapshot:','#geometry_objects_read:'))]
        assert stripped(lines)==stripped(base)
        assert next(s for s in lines if s.startswith('#rx:'))==next(s for s in base if s.startswith('#rx:'))
        assert '#waveform: impulse 1 1 impulse' in lines and '#time_window: 2.4e-06' in lines
        snapshots=[s.split()[1:] for s in lines if s.startswith('#snapshot:')]
        indices=[int(s[9]) for s in snapshots]
        expected=[]
        if g['kind']=='local':
            for n in range(40703):
                if n*dt>550e-9:break
                if 120e-9<=n*dt<=440e-9:
                    if (n-int(np.ceil(120e-9/dt)))%9==0:expected.append(n)
                elif n%17==0:expected.append(n)
            assert g['roi']['shape']==[800,1120,1] and len(indices)==833 and len(g['probes'])==31
            bounds=[g['profile_x_m']-35,12,0,g['profile_x_m']-15,40,.025]
            stride=[.025]*3
            probe_rx=[p for p in g['probes'] if p['band']=='rx'][0]
            np.testing.assert_allclose([probe_rx['x_m'],probe_rx['y_m']],g['rx_m'][:2],rtol=0,atol=1e-12)
            assert len({a['name'] for p in g['probes'] for a in p['anchors']})==186
        else:
            expected=[n for n in range(40703) if n%17==0]
            assert g['roi']['shape']==[1050,212,1] and len(indices)==2395 and not g['probes']
            bounds=[0,0,0,210,42.4,.025];stride=[.2,.2,.025]
        assert indices==expected==g['snapshot_iterations']
        for s in snapshots:
            np.testing.assert_array_equal(np.array(s[:9],float),bounds+stride)
            assert s[10].endswith('.h5')
        ncell=np.prod(g['roi']['shape'])
        assert g['snapshot_saved_bytes']==len(indices)*ncell*3*8
        assert g['snapshot_six_buffer_bytes']==len(indices)*ncell*6*8
        rows.append(dict(id=g['id'],snapshots=len(indices),output_GiB=g['snapshot_saved_bytes']/2**30,
                         cuda_history_GiB=g['snapshot_six_buffer_bytes']/2**30,input_sha256=sha(package/g['input'])))
    return dict(status='PASS_INDEPENDENT_PHYSICS_AND_SAMPLING_INPUTS',manifest_sha256=sha(package/'manifest.json'),
                solver_called=False,groups=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();result=check(a.package)
    assert not a.out.exists();a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
