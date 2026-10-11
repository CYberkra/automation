"""Six explicitly approved observer-only impulse wavefields, no retries."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys

import gprMax
import h5py
import numpy as np
import psutil

import line9_v401_version_controls as native
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]
DL = .025
FIELDS = ['Ez', 'Hx', 'Hy']
APPROVAL = ('2026-10-11 用户显性批准上述6项精细诊断：190/170m各H0/H1局部快照，'
            '2.5cm输出，0–550ns、120–440ns每9步其余17步、31处全2400ns时序及Yee邻点；'
            '另170m H0/H1全域20cm每17步至2400ns；现行物理/SFCW不变，'
            '先验收170m局部H1，失败停止不重试，新批3小时；配套代码、验收、图件和推送。')


def read(p):
    return json.loads(Path(p).read_text('utf-8'))


def steps(dt, kind):
    if kind == 'global':
        return list(range(0, 40703, 17))
    ids = {n for n in range(0, int(np.floor(550e-9/dt))+1, 17)
           if not 120e-9 <= n*dt <= 440e-9}
    ids.update(range(int(np.ceil(120e-9/dt)), int(np.floor(440e-9/dt))+1, 9))
    return sorted(ids)


def probe_layout(data, centre, rx):
    points = []
    for k, offset in enumerate([-12., -8., -4., 0., 4.]):
        x = centre+offset; i = round(x/DL); col = data[i]
        surface = (np.flatnonzero(col != 0)[-1]+1)*DL
        cover = (np.flatnonzero(col == 1)[0])*DL
        basal = (np.flatnonzero(col == 3)[-1]+1)*DL
        heights = [('air',surface+.3), ('cover_upper',cover+.75*(surface-cover)),
                   ('cover_middle',cover+.5*(surface-cover)), ('cover_lower',cover+.25*(surface-cover)),
                   ('mud',cover-.3), ('basal_below',basal-.3)]
        points += [dict(id=f'{band}_{k}', band=band, x_m=x, y_m=round(y/DL)*DL,
                        wanted_y_m=y, column=k, surface_y_m=surface, cover_y_m=cover, basal_y_m=basal)
                   for band,y in heights]
    points.append(dict(id='main_rx',band='rx',x_m=rx[0],y_m=rx[1],wanted_y_m=rx[1],column=-1))
    for p in points:
        i,j = round(p['x_m']/DL),round(p['y_m']/DL)
        # First three anchors collocate at Ez's native location. Positive
        # neighbours also independently reproduce the snapshot cell centre.
        stencils = [(i,j,['Ez','Hx','Hy']), (i,j-1,['Hx']), (i-1,j,['Hy']),
                    (i+1,j,['Ez','Hx']), (i,j+1,['Ez','Hy']), (i+1,j+1,['Ez'])]
        p['anchors'] = [dict(coord=[ii,jj,0], position_m=[ii*DL,jj*DL,0.], outputs=fields,
                             name=f"p_{p['id']}_a{k}",receiver_index=2+6*points.index(p)+k)
                        for k,(ii,jj,fields) in enumerate(stencils)]
    assert len(points) == 31
    return points


def prepare(a):
    assert not a.out.exists()
    parent = read(a.package/'manifest.json'); prior = read(a.prior/'manifest.json')
    assert parent['execution_authorized'] and len(parent['groups']) == 42
    assert prior['SFCW'] == parent['SFCW'] and prior['dt_s'] == parent['dt_s']
    cert = read(a.raw/'completed_verification.json')
    assert cert['completed'] and len(cert['groups']) == 42
    assert cert['contract_sha256'] == sha(a.raw/'execution_contract.json')
    receipts = {r['id']:r for r in cert['groups']}
    dt = parent['dt_s']; a.out.mkdir(parents=True)
    geometry_dir = a.out/'geometries'; geometry_dir.mkdir()
    for role in ['H0','H1']:
        old = next(g for g in parent['groups'] if g['role']==role)
        shutil.copyfile(a.package/old['geometry'], geometry_dir/f'{role}.h5')
        assert sha(geometry_dir/f'{role}.h5') == old['geometry_sha256']
    mat = geometry_dir/'line9_materials_v0_2.json'
    shutil.copyfile(a.package/parent['groups'][0]['material'],mat)
    with h5py.File(geometry_dir/'H1.h5') as h: geometry = h['data'][:,:,0]
    baselines = {}
    for x in [170.,190.]:
        for role in ['H1','H0']:
            source = parent if x == 170 else prior
            g = next(g for g in source['groups'] if g['role']==role and g.get('profile_x_m',190.)==x)
            package = a.package if x == 170 else a.prior
            raw = a.raw/(g['id']+'.h5') if x == 170 else a.reused/(g['id']+'.h5')
            expected = receipts[g['id']]['native_sha256'] if x == 170 else next(
                r for r in parent['reused_solver_outputs'] if r['role']==role)['native_sha256']
            assert sha(raw)==expected
            base = a.out/'baseline'/f'x{int(x)}_{role}'; base.mkdir(parents=True)
            shutil.copyfile(package/g['input'],base/'profile.in'); shutil.copyfile(raw,base/'native.h5')
            baselines[x,role] = dict(g, baseline=f'baseline/x{int(x)}_{role}',baseline_native_sha256=expected)
    order = [(170.,'H1','local'),(170.,'H0','local'),(190.,'H1','local'),
             (190.,'H0','local'),(170.,'H1','global'),(170.,'H0','global')]
    groups = []
    for x,role,kind in order:
        g = copy.deepcopy(baselines[x,role]); centre=x-20
        g.update(id=f'fine_x{int(x)}_{role}_{kind}', profile_x_m=x, role=role, kind=kind,
                 geometry=f'geometries/{role}.h5',material='geometries/line9_materials_v0_2.json')
        for k in ['parent_contract_sha256','native_path','native_sha256']:
            g.pop(k,None)
        g['input'] = f"cases/{g['id']}/profile.in"
        g['roi'] = dict(bounds=[centre-15,12,0,centre+5,40,.025],spacing=[.025,.025,.025],shape=[800,1120,1]) if kind=='local' else dict(
            bounds=[0,0,0,210,42.4,.025],spacing=[.2,.2,.025],shape=[1050,212,1])
        g['snapshot_iterations']=steps(dt,kind)
        g['probes']=probe_layout(geometry,centre,g['rx_m']) if kind=='local' else []
        lines = (a.out/g['baseline']/'profile.in').read_text('utf-8').splitlines()
        lines = [f'#title: {g["id"]}; explicitly approved observer-only impulse' if s.startswith('#title:') else
                 f'#geometry_objects_read: 0 0 0 ../../geometries/{role}.h5 line9_materials_v0_2 y' if s.startswith('#geometry_objects_read:') else s for s in lines]
        for p in g['probes']:
            for anchor in p['anchors']:
                px,py,_=anchor['position_m']
                lines.append(f"#rx: {px:.15g} {py:.15g} 0.0125 {anchor['name']} "+' '.join(anchor['outputs']))
        roi = g['roi']
        for n in g['snapshot_iterations']:
            lines.append('#snapshot: '+' '.join(f'{v:.15g}' for v in roi['bounds']+roi['spacing'])+f' {n} {kind}_{n:05d}.h5')
        card = a.out/g['input'];card.parent.mkdir(parents=True)
        card.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
        g['baseline_input_sha256']=sha(a.out/g['baseline']/'profile.in')
        g['receiver_count']=1+6*len(g['probes'])
        g['snapshot_saved_bytes']=len(g['snapshot_iterations'])*int(np.prod(roi['shape']))*24
        g['snapshot_six_buffer_bytes']=2*g['snapshot_saved_bytes']
        groups.append(g)
    manifest = dict(status='PREPARED_APPROVED_SIX_FINE_WAVEFIELDS',approval_basis=APPROVAL,
        execution_authorized=True, dt_s=dt, groups=groups, SFCW=parent['SFCW'], no_retry=True,
        fields=FIELDS, collocation='Raw node E; mean staggered H neighbours and H[n],H[n+1]. Snapshot cell centres independently reproduced by positive native neighbours.',
        parent_certificate_sha256=sha(a.raw/'completed_verification.json'),
        parent_manifest_sha256=sha(a.package/'manifest.json'),generator_sha256=sha(__file__),
        saved_snapshot_bytes=sum(g['snapshot_saved_bytes'] for g in groups),
        limits='Observer-only 6 cases, not new spatial line. Snapshot temporal sampling does not certify broadband impulse bandwidth. Complete native probes are used for exact501-tone checks. Total E cross H is net flux, not unique rays/bounce energies or dispersive energy closure. No 3D, mesh/time-length convergence or field calibration.')
    native.save(a.out/'manifest.json',manifest)
    native.save(a.out/'input_verification.json',audit_inputs(a.out))
    print(json.dumps(dict(status=manifest['status'],groups=6,frames=[len(g['snapshot_iterations']) for g in groups],
                         saved_GiB=manifest['saved_snapshot_bytes']/2**30)))


def audit_inputs(package):
    from gprMax.hash_cmds_file import get_user_objects
    from gprMax.user_objects.cmds_output import Snapshot
    m=read(package/'manifest.json');assert m['approval_basis']==APPROVAL and m['execution_authorized'] and m['no_retry']
    assert len(m['groups'])==6 and m['fields']==FIELDS
    expected=[(170.,'H1','local'),(170.,'H0','local'),(190.,'H1','local'),(190.,'H0','local'),(170.,'H1','global'),(170.,'H0','global')]
    assert [(g['profile_x_m'],g['role'],g['kind']) for g in m['groups']]==expected
    rows=[]
    for g in m['groups']:
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        base=package/g['baseline'];assert sha(base/'profile.in')==g['baseline_input_sha256'] and sha(base/'native.h5')==g['baseline_native_sha256']
        lines=(package/g['input']).read_text('utf-8').splitlines();old=(base/'profile.in').read_text('utf-8').splitlines()
        filtered=lambda seq:[s for s in seq if not s.startswith(('#title:','#geometry_objects_read:','#snapshot:','#rx:'))]
        assert filtered(lines)==filtered(old)
        oldrx=[s for s in old if s.startswith('#rx:')];rx=[s for s in lines if s.startswith('#rx:')]
        assert rx[:1]==oldrx and len(rx)==g['receiver_count']
        assert g['snapshot_iterations']==steps(m['dt_s'],g['kind'])
        assert len([s for s in lines if s.startswith('#snapshot:')])==len(g['snapshot_iterations'])
        objects=get_user_objects(lines,input_dir=(package/g['input']).parent)
        snaps=[o for o in objects if isinstance(o,Snapshot)]
        assert [s.iterations for s in snaps]==g['snapshot_iterations']
        for obj in snaps:
            np.testing.assert_array_equal(obj.lower_bound,g['roi']['bounds'][:3]);np.testing.assert_array_equal(obj.upper_bound,g['roi']['bounds'][3:])
            np.testing.assert_array_equal(obj.dl,g['roi']['spacing']);assert obj.file_extension=='.h5'
        if g['probes']:
            assert len(g['probes'])==31
            with h5py.File(package/g['geometry']) as h:data=h['data'][:,:,0]
            with h5py.File(package/'geometries/H1.h5') as h:target=h['data'][:,:,0]
            assert g['probes']==probe_layout(target,g['profile_x_m']-20,g['rx_m'])
            for p in g['probes']:
                for a in p['anchors']:
                    i,j,_=a['coord'];assert 80<i<8320 and 80<j<1620
                    stencil=data[i-1:i+2,j-1:j+2];assert len(np.unique(stencil))==1,'Probe interpolation crosses interface'
        rows.append(dict(id=g['id'],snapshots=len(snaps),receivers=len(rx),same_physical_input=True))
    return dict(status='PASS_OFFICIAL_PARSE_SIX_OBSERVER_ONLY_INPUTS',manifest_sha256=sha(package/'manifest.json'),groups=rows,solver_called=False)


def audit_group(c,g):
    package=Path(c['package']);raw=Path(g['input']).with_suffix('.h5');base=package/g['baseline']/'native.h5';dt=c['study_manifest']['dt_s']
    histories=[]
    with h5py.File(raw) as h,h5py.File(base) as old:
        for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
            assert h[key][:].tobytes()==old[key][:].tobytes(),'Passive observers changed main trace'
        assert len(h['rxs'])==g['receiver_count']
        for p in g['probes']:
            values=[]
            for a in p['anchors']:
                r=h[f"rxs/rx{a['receiver_index']}"]
                assert r.attrs['Name']==a['name'] and set(r.keys())==set(a['outputs'])
                np.testing.assert_array_equal(r.attrs['GridPosition'],a['coord'])
                np.testing.assert_allclose(r.attrs['Position'],a['position_m'],atol=1e-12,rtol=0)
                rr={}
                for field in a['outputs']:
                    ds=r[field];v=ds[:];assert v.dtype==np.float64 and v.shape==(40703,) and np.isfinite(v).all()
                    assert ds.attrs['SampleInterval']==dt and ds.attrs['TimeSampleOffset']==(0 if field=='Ez' else -.5*dt)
                    rr[field]=v
                values.append(rr)
            histories.append(values)
    roi=g['roi'];folder=raw.parent/'profile_snaps';records=[];maximum=0.
    for n in g['snapshot_iterations']:
        p=folder/f"{g['kind']}_{n:05d}.h5"
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['iteration']==n
            assert abs(h.attrs['time']-n*dt)<1e-20 and abs(h.attrs['magnetic_time']-(n-.5)*dt)<1e-20
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],roi['shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],roi['spacing'])
            np.testing.assert_allclose(h.attrs['origin'],roi['bounds'][:3],atol=1e-12,rtol=0)
            assert set(h)==set(FIELDS)
            for field in FIELDS:
                v=h[field][:];assert v.dtype==np.float64 and list(v.shape)==roi['shape'] and np.isfinite(v).all()
                scale=float(abs(v).max())
                if histories:
                    for probe,rr in zip(g['probes'],histories):
                        i=round((probe['x_m']-roi['bounds'][0])/DL);j=round((probe['y_m']-roi['bounds'][1])/DL)
                        if field=='Ez':expected=(rr[0]['Ez'][n]+rr[3]['Ez'][n]+rr[4]['Ez'][n]+rr[5]['Ez'][n])*.25
                        elif field=='Hx':expected=(rr[0]['Hx'][n]+rr[3]['Hx'][n])*.5
                        else:expected=(rr[0]['Hy'][n]+rr[4]['Hy'][n])*.5
                        error=abs(float(v[i,j,0])-expected)/max(abs(expected),scale,1e-300)
                        maximum=max(maximum,error);assert error<1e-12
        records.append(dict(iteration=n,file=p.name,sha256=sha(p),bytes=p.stat().st_size))
    assert sorted(p.name for p in folder.iterdir())==sorted(r['file'] for r in records)
    storage=read(raw.parent/'snapshot_storage.json')
    assert storage['native_gpu_to_host_streaming'] and storage['snapshot_count']==len(records)
    assert storage['snapshot_history_bytes']==g['snapshot_saved_bytes']
    result=dict(status='PASS_NATIVE_FP64_SNAPSHOTS_AND_OBSERVER_IDENTITY',id=g['id'],native_sha256=sha(raw),
        native_bytes=raw.stat().st_size,main_source_receiver_bitwise_equal=True,snapshot_count=len(records),
        snapshot_probe_max_normalised_abs_error=maximum,snapshots=records,snapshot_storage_sha256=sha(raw.parent/'snapshot_storage.json'))
    native.save(raw.parent/'observer_verification.json',result)
    print(json.dumps(dict(verified=g['id'],frames=len(records),max_error=maximum)),flush=True)
    return result


def audit(path,completed=False):
    c=read(path);assert c['approval_basis']==APPROVAL and c['max_runs']==6 and c['no_retry']
    assert c['solver_entrypoint']=='gprmax_fine_snapshot_entry.py' and c['max_batch_wall_s']==10800
    assert c['study_manifest']==read(Path(c['package'])/'manifest.json')
    result=native.audit(path,completed)
    if completed:
        for r,g in zip(result['groups'],c['groups']):
            p=Path(g['input']).parent/'observer_verification.json';v=read(p)
            assert v['id']==g['id'] and v['native_sha256']==r['native_sha256'] and v['main_source_receiver_bitwise_equal']
            assert v['snapshot_count']==len(g['snapshot_iterations'])
            r.update(observer_verification_sha256=sha(p),observer_verification=v)
        result['status']='PASS_SIX_OBSERVER_ONLY_NATIVE_FP64_CASES'
    return result


def freeze(a):
    inputs=audit_inputs(a.package);m=read(a.package/'manifest.json')
    assert psutil.virtual_memory().available>43*2**30
    assert shutil.disk_usage(a.package).free>m['saved_snapshot_bytes']+10*2**30
    native.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    p=a.out/'execution_contract.json';c=read(p)
    c.update(approval_basis=APPROVAL,study_manifest=m,package=str(a.package.resolve()),
        solver_entrypoint='gprmax_fine_snapshot_entry.py',max_batch_wall_s=10800,max_group_wall_s=3600,
        min_available_RAM_GiB=43,max_owned_RSS_GiB=44,min_system_available_during_run_GiB=6,
        min_free_VRAM_GiB=9,max_runs=6,no_retry=True)
    for name in ['line9_fine_wavefield.py','gprmax_fine_snapshot_entry.py','gprmax_snapshot_cuda_entry.py']:
        c['code_identities'][str(ROOT/'scripts'/name)]=sha(ROOT/'scripts'/name)
    for f in a.package.rglob('*'):
        if f.is_file():c['file_identities'][str(f.resolve())]=sha(f)
    native.save(p,c);native.save(a.out/'input_verification.json',inputs);native.save(a.out/'preflight_verification.json',audit(p))
    print(json.dumps(dict(status='FROZEN_SIX_APPROVED_FINE_WAVEFIELDS',contract_sha256=sha(p))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','freeze','run','verify','input-audit'])
    for name in ['out','package','prior','raw','reused']:parser.add_argument('--'+name,type=Path,required=name=='out')
    a=parser.parse_args();a.out=a.out.resolve()
    if a.action=='prepare':prepare(a)
    elif a.action=='freeze':freeze(a)
    elif a.action=='input-audit':print(json.dumps(audit_inputs(a.package)))
    elif a.action=='run':
        c=read(a.out/'execution_contract.json')
        native.supervisor.audit=audit
        native.supervisor.group_completed_hook=lambda g:audit_group(c,g)
        native.supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
