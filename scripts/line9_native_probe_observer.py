"""One immutable190m H0 with native-time, spatially collocatable TM probe receivers."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha

DL = .025
XVALUES = np.arange(154,174.001,.5)
BANDS = ['air_surface','cover_upper','cover_middle','cover_lower','mud_below','air_34','air_38']


def read(p):
    return json.loads(Path(p).read_text('utf-8'))


def save(p,j):
    Path(p).write_text(json.dumps(j,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def layout(geometry):
    with h5py.File(geometry,'r') as h:
        g=h['data'][:,:,0]
    assert g.shape==(8400,1700) and set(np.unique(g))=={0,1,2}
    rows=[]
    for col,x in enumerate(XVALUES):
        i=int(round(x/DL))
        surface=(np.flatnonzero(g[i]!=0)[-1]+1)*DL
        bottom=(np.flatnonzero(g[i]==2)[-1]+1)*DL
        ys=[surface+.25,surface-.25,(surface+bottom)/2,bottom+.25,bottom-.25,34.025,38.025]
        for band,y in zip(BANDS,ys):
            j=int(round(y/DL))
            # Three raw anchors gather Hx above/below and Hy left/right of Ez.
            points=[(i,j,['Ez','Hx','Hy']),(i,j-1,['Hx']),(i-1,j,['Hy'])]
            media=[]
            for ii,jj,fields in points:
                assert 80<ii<8320 and 80<jj<1620
                support=g[ii-1:ii+2,jj-1:jj+2]
                assert len(np.unique(support))==1, 'Collocation stencil crosses a contact'
                media.append(int(support[0,0]))
            expected=0 if band.startswith('air') else 2 if band=='mud_below' else 1
            assert media==[expected]*3
            rows.append({'id':f'{band}_c{col:02d}','band':band,'x_m':float(x),'y_m':j*DL,
                         'wanted_y_m':float(y),'material_id':expected,'surface_m':surface,'bottom_m':bottom,
                         'anchors':[{'coord':[ii,jj,0],'position_m':[ii*DL,jj*DL,0],
                                     'outputs':outputs,'name':f'p_{band}_c{col:02d}_a{k}',
                                     'receiver_index':2+3*len(rows)+k}
                                    for k,(ii,jj,outputs) in enumerate(points)]})
    assert len(rows)==287
    return rows


def prepare(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    base=a.out/'baseline';base.mkdir()
    for name in ['geometry.h5','materials.json','profile.in','native_H0.h5']:
        shutil.copyfile(a.baseline/name,base/name)
    dep=a.out/'dependency';dep.mkdir()
    for name in ['execution_contract.json','completed_verification.json']:
        shutil.copyfile(a.source/name,dep/name)
    c=read(dep/'execution_contract.json');v=read(dep/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(dep/'execution_contract.json')
    assert sha(base/'native_H0.h5')==v['groups'][0]['native_sha256']
    origin=c['groups'][0];gp=a.out/'native_probes_H0';(gp/'geometries').mkdir(parents=True)
    geometry=gp/'geometries/full2d_compact.h5';material=gp/'geometries/line9_research_materials_v1_smoothed.json'
    shutil.copyfile(base/'geometry.h5',geometry);shutil.copyfile(base/'materials.json',material)
    assert sha(geometry)==origin['geometry_sha256'] and sha(material)==origin['material_sha256']
    probes=layout(geometry)
    inp=gp/'cases/native_probes_H0/profile.in';inp.parent.mkdir(parents=True)
    original=(base/'profile.in').read_bytes();assert original.endswith(b'\n')
    lines=[]
    for probe in probes:
        for point in probe['anchors']:
            x,y,_=point['position_m']
            lines.append(f"#rx: {x:.6f} {y:.6f} 0.0125 {point['name']} "+' '.join(point['outputs']))
    inp.write_bytes(original+('\n'.join(lines)+'\n').encode('utf-8'))
    row=copy.deepcopy(origin)
    row.update(id='native_probes_H0',input=inp.relative_to(a.out).as_posix(),
               geometry=geometry.relative_to(a.out).as_posix(),material=material.relative_to(a.out).as_posix())
    for k in ['input','geometry','material']:row[k+'_sha256']=sha(a.out/row[k])
    nreceiver=1+len(lines)
    m={'status':'PREPARED_ONE_NATIVE_TIME_OBSERVER','approval_basis':'User继续研究明白; next observer-only190mH0 native-time experiment, one attempt; no stopped-line restart.',
       'groups':[row],'dt_s':row['dt_s'],'no_retry':True,'probes':probes,'receiver_count':nreceiver,
       'device_receiver_history_bytes':6*20352*nreceiver*8,'saved_receiver_history_bytes':20352*8*(1+5*len(probes)),
       'collocation':{'Ez':'raw anchor(i,j)','Hx':'mean raw anchors(i,j) and(i,j-1)',
                      'Hy':'mean raw anchors(i,j) and(i-1,j)','H_at_E_time':'mean H[n] and H[n+1], n=0..N-2',
                      'TM_invariant_axis':'z collapsed to0; x/y Yee offsets checked, not finite3D'},
       'gates_native_ns':{'first':[120,250],'late':[300,370],'full':[0,500],'later':[500,1100]},
       'SFCW':{'min_Hz':20e6,'max_Hz':170e6,'step_Hz':300000,'tones':501,'window':['Hann','Blackman'],'no_AGC':True},
       'baseline_sha256':{p.name:sha(p) for p in base.iterdir()},'dependency_sha256':{p.name:sha(p) for p in dep.iterdir()},
       'source_read_identities':read(a.runtime_sources/'source_identity.json'),
       'limits':'Only native point observations changed; total E/H are not separated rays or conserved individual-pulse energy; no whole line, unique bounce, site, finite3D or training certification.'}
    save(a.out/'manifest.json',m)
    print(json.dumps({'status':m['status'],'collocated_points':len(probes),'raw_receivers':nreceiver,
                      'receiver_device_bytes':m['device_receiver_history_bytes'],'saved_bytes':m['saved_receiver_history_bytes']}))


def audit(path,completed=False):
    import line9_v401_version_controls as native
    from audit_line9_native_probe_inputs import check
    c=read(path);package=Path(c['package']);m=read(package/'manifest.json')
    assert m==c['study_manifest'];check(package)
    result=native.audit(path,completed)
    if completed:
        row=result['groups'][0];raw=Path(row['native_path'])
        with h5py.File(raw,'r') as h,h5py.File(package/'baseline/native_H0.h5','r') as old:
            for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
                assert h[key][:].tobytes()==old[key][:].tobytes(), 'Observer changed native baseline'
            assert h.attrs['nrx']==m['receiver_count']==len(h['rxs'])
            dt=float(h.attrs['dt'])
            for probe in m['probes']:
                for anchor in probe['anchors']:
                    r=h[f"rxs/rx{anchor['receiver_index']}"]
                    assert r.attrs['Name']==anchor['name'] and set(r.keys())==set(anchor['outputs'])
                    np.testing.assert_array_equal(r.attrs['GridPosition'],anchor['coord'])
                    np.testing.assert_allclose(r.attrs['Position'],anchor['position_m'],atol=1e-12,rtol=0)
                    for field in anchor['outputs']:
                        ds=r[field];x=ds[:]
                        assert x.dtype==np.float64 and x.shape==(20352,) and np.isfinite(x).all()
                        assert ds.attrs['Quantity']==field and ds.attrs['SampleInterval']==dt
                        assert ds.attrs['TimeSampleOffset']==(0 if field=='Ez' else -.5*dt)
        row.update(observer_native_source_and_receiver_bitwise_equal=True,collocated_probe_count=287,raw_receiver_count=862)
        result['status']='PASS_NATIVE_TIME862_RX_AND_BASELINE_SOURCE_RX_BITWISE_INVARIANT'
    return result


def freeze(a):
    import gprMax
    import psutil
    import line9_v401_version_controls as native
    from audit_line9_native_probe_inputs import check
    assert check(a.package)==read(a.package/'independent_input_audit.json')
    m=read(a.package/'manifest.json');old=read(a.package/'dependency/execution_contract.json')
    assert gprMax.__version__=='4.0.1'
    installed=Path(gprMax.__file__).parent
    for name,digest in old['source_identities'].items():assert sha(installed/name)==digest
    for r in m['source_read_identities']:assert sha(installed/r['name'])==r['sha256']
    from gprMax.hash_cmds_file import get_user_objects
    card=a.package/m['groups'][0]['input']
    get_user_objects(card.read_text('utf-8').splitlines(),input_dir=card.parent)
    for p in psutil.process_iter(['pid','name','cmdline']):
        try:
            if p.pid!=__import__('os').getpid() and 'python' in (p.info['name'] or '').lower():
                assert 'line9_basal_pair_20261008_r1' not in ' '.join(p.info['cmdline'] or []), 'Prior owned solver still active'
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    native.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=read(path)
    c.update(package=str(a.package.resolve()),study_manifest=m,approval_basis=m['approval_basis'],
             max_runs=1,max_batch_wall_s=2100,gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock',
             min_free_VRAM_GiB=(8400*1700*420+2*2**30+m['device_receiver_history_bytes'])/2**30)
    for name in ['line9_native_probe_observer.py','audit_line9_native_probe_inputs.py']:
        c['code_identities'][str(Path(__file__).parent/name)]=sha(Path(__file__).parent/name)
    for f in a.package.rglob('*'):
        if f.is_file():c['file_identities'][str(f.resolve())]=sha(f)
    assert c['hardware_at_freeze']['free_VRAM_bytes']>c['min_free_VRAM_GiB']*2**30
    assert shutil.disk_usage(a.package).free>m['saved_receiver_history_bytes']+2**30
    native.save(path,c);native.save(a.out/'preflight_verification.json',audit(path))
    print(json.dumps({'status':'FROZEN_ONE_NATIVE_TIME_OBSERVER','contract_sha256':sha(path),'receiver_count':862}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    p.add_argument('--out',type=Path,required=True)
    for name in ['baseline','source','runtime-sources','package']:p.add_argument('--'+name,type=Path)
    a=p.parse_args();a.out=a.out.resolve()
    if a.action=='prepare':prepare(a)
    elif a.action=='freeze':freeze(a)
    elif a.action=='run':
        import hs4_station_grid_controls as supervisor
        supervisor.audit=audit;supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
