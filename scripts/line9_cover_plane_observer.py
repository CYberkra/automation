"""One immutable190mH0: preserve862 receivers, add369 flat-cover observers."""
import argparse,copy,json,shutil
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def read(p):return json.loads(Path(p).read_text('utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def prepare(a):
    from audit_line9_native_probe_inputs import check
    assert not a.out.exists();check(a.parent)
    old=read(a.parent/'manifest.json');proposal=read(a.proposal);cert=read(a.source/'completed_verification.json')
    assert cert['completed'] and sha(a.source/'execution_contract.json')==cert['contract_sha256']
    assert sha(a.source/'profile.h5')==cert['groups'][0]['native_sha256']== '46288bae706a5534940baf3b87748c0f68eea02ffc6e1ba6ae1ae17c6efbe68e'
    assert proposal['parent_manifest_sha256']==sha(a.parent/'manifest.json')
    assert proposal['status']=='GEOMETRY_VERIFIED_PROPOSAL_NOT_PREPARED_FROZEN_OR_RUN'
    a.out.mkdir(parents=True);shutil.copytree(a.parent,a.out/'parent')
    for n in ['execution_contract.json','completed_verification.json']:shutil.copyfile(a.source/n,a.out/('parent_'+n))
    shutil.copyfile(a.proposal,a.out/'geometry_proposal.json')
    g=copy.deepcopy(old['groups'][0]);group='native_cover_planes_H0';folder=a.out/group;gp=folder/'geometries';gp.mkdir(parents=True)
    for key in ['geometry','material']:shutil.copyfile(a.parent/g[key],gp/Path(g[key]).name)
    original=(a.parent/g['input']).read_bytes();assert original.endswith(b'\n')
    new=copy.deepcopy(proposal['new_probes'])
    with h5py.File(a.parent/g['geometry']) as h:geo=h['data'][:,:,0]
    lines=[]
    for p in new:
        i=p['anchors'][0]['coord'][0]
        p.update(band='flat_cover_'+str(round(p['y_m']/.025)),wanted_y_m=p['y_m'],
                 surface_m=(np.flatnonzero(geo[i]!=0)[-1]+1)*.025,bottom_m=(np.flatnonzero(geo[i]==2)[-1]+1)*.025)
        for q in p['anchors']:
            x,y,_=q['position_m'];lines.append(f"#rx: {x:.6f} {y:.6f} 0.0125 {q['name']} "+' '.join(q['outputs']))
    card=folder/'cases'/group/'profile.in';card.parent.mkdir(parents=True);card.write_bytes(original+('\n'.join(lines)+'\n').encode('utf-8'))
    g.update(id=group,input=card.relative_to(a.out).as_posix(),geometry=(gp/Path(g['geometry']).name).relative_to(a.out).as_posix(),
             material=(gp/Path(g['material']).name).relative_to(a.out).as_posix())
    for k in ['input','geometry','material']:g[k+'_sha256']=sha(a.out/g[k])
    m={'status':'PREPARED_ONE_FLAT_COVER_OBSERVER_NOT_RUN','approval_basis':'User持续推进直到研究明白; one newly frozen observer-only190mH0; preserve previous source/main/287 probes, no stopped-line restart.',
       'groups':[g],'dt_s':old['dt_s'],'no_retry':True,'receiver_count':1231,'probes':old['probes']+new,
       'previous_receiver_count':862,'previous_logical_points':287,'new_logical_points':123,'cover_plane_y_m':[27.025,28.525,30.025],
       'device_receiver_history_bytes':1231*20352*6*8,'saved_receiver_history_bytes':20352*8*(1+5*410),
       'parent_manifest_sha256':sha(a.out/'parent/manifest.json'),'parent_contract_sha256':sha(a.out/'parent_execution_contract.json'),
       'parent_completed_sha256':sha(a.out/'parent_completed_verification.json'),'reference_native_sha256':cert['groups'][0]['native_sha256'],
       'geometry_proposal_sha256':sha(a.out/'geometry_proposal.json'),'generator_sha256':sha(__file__),
       'source_read_identities':old['source_read_identities'],'collocation':old['collocation'],'SFCW':old['SFCW'],
       'gates_native_ns':old['gates_native_ns'],
       'analysis_plan':{'full_complex_cover_TE_up_down':True,'plane_height_steps_m':[1.5,3.],
                       'smooth_air_radiating_sector':True,'aperture_window_sensitivity':True,'native_time_gate_sensitivity':True,
                       'amplitude_phase_or_delay_fit':False,'source_main_and_previous_observers_bitwise_invariant_required':True},
       'limits':'Only observers changed; homogeneous cover support locally checked, finite20m aperture/exterior heterogeneity remain. Not unique bounce/ray/site/material/finite3D or training certification.'}
    assert len(m['probes'])==410;save(a.out/'manifest.json',m)
    print(json.dumps({'status':m['status'],'logical_points':410,'raw_receivers':1231,'device_receiver_history_bytes':m['device_receiver_history_bytes']}))


def audit(path,completed=False):
    import line9_v401_version_controls as native
    from audit_line9_cover_plane_inputs import check
    c=read(path);package=Path(c['package']);m=read(package/'manifest.json');assert m==c['study_manifest'];check(package)
    assert sha(c['reference_native'])==m['reference_native_sha256']
    result=native.audit(path,completed)
    if completed:
        row=result['groups'][0];channels=0
        with h5py.File(row['native_path']) as h,h5py.File(c['reference_native']) as old:
            assert h.attrs['nrx']==len(h['rxs'])==1231
            assert h['srcs/src1/excitation/samples'][:].tobytes()==old['srcs/src1/excitation/samples'][:].tobytes()
            for name,r in old['rxs'].items():
                nr=h['rxs/'+name]
                for attr in ['Name','GridPosition','Position']:np.testing.assert_array_equal(nr.attrs[attr],r.attrs[attr])
                for field in r:
                    assert nr[field][:].tobytes()==r[field][:].tobytes(),(name,field)
                    channels+=1
            for p in m['probes']:
                for q in p['anchors']:
                    rr=h[f"rxs/rx{q['receiver_index']}"];assert rr.attrs['Name']==q['name'] and set(rr)==set(q['outputs'])
                    np.testing.assert_array_equal(rr.attrs['GridPosition'],q['coord'])
                    np.testing.assert_allclose(rr.attrs['Position'],q['position_m'],rtol=0,atol=1e-12)
                    for field in q['outputs']:
                        d=rr[field];v=d[:];assert v.dtype==np.float64 and v.shape==(20352,) and np.isfinite(v).all()
                        assert d.attrs['Quantity']==field and d.attrs['SampleInterval']==m['dt_s']
                        assert d.attrs['TimeSampleOffset']==(0 if field=='Ez' else -.5*m['dt_s'])
        assert channels==1436
        row.update(previous_1436_receiver_channels_bitwise_equal=True,source_bitwise_equal=True,
                   logical_points=410,raw_receivers=1231)
        result['status']='PASS1231_NATIVE_RX_AND_SOURCE_PREVIOUS1436_CHANNELS_BITWISE_INVARIANT'
    return result


def freeze(a):
    import os,sys,gprMax,psutil
    import line9_v401_version_controls as native
    from audit_line9_cover_plane_inputs import check
    m=read(a.package/'manifest.json');assert check(a.package)==read(a.package/'independent_input_audit.json')
    old=read(a.package/'parent_execution_contract.json');assert gprMax.__version__=='4.0.1'
    installed=Path(gprMax.__file__).parent
    for name,digest in old['source_identities'].items():assert sha(installed/name)==digest,name
    for r in m['source_read_identities']:assert sha(installed/r['name'])==r['sha256']
    assert sha(a.reference_native)==m['reference_native_sha256']
    from gprMax.hash_cmds_file import get_user_objects
    card=a.package/m['groups'][0]['input'];get_user_objects(card.read_text('utf-8').splitlines(),input_dir=card.parent)
    # Windows venv python.exe delegates to its base Python child. Exclude only
    # ancestors with this exact script/argument vector, never arbitrary parents.
    own_launchers={p.pid for p in psutil.Process().parents()
                   if p.cmdline()[1:]==sys.argv}
    for p in psutil.process_iter(['pid','name','cmdline']):
        try:
            cmd=' '.join(p.info['cmdline'] or []).lower()
            if p.pid not in own_launchers|{os.getpid()} and 'python' in (p.info['name'] or '').lower():
                assert 'line9_basal_pair_20261008_r1' not in cmd and 'gprmax_cached_cuda_entry.py' not in cmd,'Prior scientific process active'
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    native.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    path=a.out/'execution_contract.json';c=read(path)
    c.update(package=str(a.package.resolve()),study_manifest=m,approval_basis=m['approval_basis'],reference_native=str(a.reference_native.resolve()),
             max_runs=1,max_batch_wall_s=2100,gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock',
             min_free_VRAM_GiB=(8400*1700*420+2*2**30+m['device_receiver_history_bytes'])/2**30)
    for name in ['line9_cover_plane_observer.py','audit_line9_cover_plane_inputs.py','audit_line9_native_probe_inputs.py','check_line9_cover_plane_inputs.py']:
        c['code_identities'][str(Path(__file__).parent/name)]=sha(Path(__file__).parent/name)
    for p in a.package.rglob('*'):
        if p.is_file():c['file_identities'][str(p.resolve())]=sha(p)
    c['file_identities'][str(a.reference_native.resolve())]=sha(a.reference_native)
    assert c['hardware_at_freeze']['free_VRAM_bytes']>c['min_free_VRAM_GiB']*2**30
    assert shutil.disk_usage(a.package).free>m['saved_receiver_history_bytes']+2**30
    native.save(path,c);native.save(a.out/'preflight_verification.json',audit(path))
    print(json.dumps({'status':'FROZEN_ONE_COVER_PLANE_OBSERVER','contract_sha256':sha(path),'receiver_count':1231}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify']);p.add_argument('--out',type=Path,required=True)
    for n in ['parent','proposal','source','package','reference-native']:p.add_argument('--'+n,type=Path)
    a=p.parse_args();a.out=a.out.resolve()
    if a.action=='prepare':prepare(a)
    elif a.action=='freeze':freeze(a)
    elif a.action=='run':
        import hs4_station_grid_controls as supervisor
        supervisor.audit=audit;supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
