"""Bounded source/window factorial and two additional Line9 causal stations."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import line9_basal_pair as pair
import line9_boundary_controls as boundary
import hs4_station_grid_controls as supervisor


def excitation(kind,freq,n,dt):
    if kind=='impulse':
        s=np.zeros(n);s[0]=40.;return s
    t=(np.arange(n)+.5)*dt;u=(np.pi*freq*(t-np.sqrt(2)/freq))**2
    return 40*(1-2*u)*np.exp(-u)


def prepare(parent,reference,original,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh outputs required')
    m=json.loads((parent/'manifest.json').read_text('utf-8'));db=parent/'H0/geometries/line9_research_materials_v1_smoothed.json'
    assert pair.sha(db)==m['material_sha256']
    data={}
    for name,j in [('H1',0),('H0',1)]:
        p=parent/name/'geometries/full2d_compact.h5';assert pair.sha(p)==m['groups'][j]['geometry_sha256']
        with h5py.File(p) as h:data[name]=h['data'][:]
    far=data['H0'].copy();x=np.arange(4400)*.025;mask=(far==3)&((x>=45)&(x<60))[:,None,None];far[mask]=2
    data['far']=far;np.testing.assert_array_equal(far[~mask],data['H0'][~mask])
    with h5py.File(reference) as h:dt=float(h.attrs['dt'])
    from review_line9_result_packages import geometry_at,indices
    from diagnose_line9_layer_kinematics import primary_paths
    from audit_line9_postprocessing import FREQ,inverse,weights
    materials=json.loads(db.read_text('utf-8'))['materials'];n95=indices(materials,95e6)
    rows=[]
    for kind,time in [('ricker',800),('impulse',1600),('ricker',1600)]:rows.append((f'{kind}_{time}','c0008','H0',78.6,79.9,kind,time))
    for station,tx,rx in [('c0058',68.6,69.9),('c0108',58.6,59.9)]:
        for variant in ['H1','H0','far']:rows.append((station+'_'+variant,station,variant,tx,rx,'impulse',800))
    original_manifest=json.loads((original/'package_manifest.json').read_text('utf-8'))
    declared={g['id']:g for g in original_manifest['cases']}
    card='\n'.join(s for s in (parent/'H0/cases/full2d_c0008/profile.in').read_text('utf-8').splitlines() if not s.startswith('#snapshot:'))+'\n'
    out.mkdir(parents=True);shutil.copyfile(reference,out/'reference_H0.h5');groups=[];stations={}
    for name,station,variant,tx,rx,kind,time in rows:
        if station not in stations:
            sid='full2d_'+station;raw=original/'cases'/sid/'profile.h5'
            with h5py.File(raw) as h:
                np.testing.assert_allclose(h['srcs/src1'].attrs['Position'],[tx,35,0],rtol=0,atol=1e-12)
                np.testing.assert_allclose(h['rxs/rx1'].attrs['Position'],[rx,35,0],rtol=0,atol=1e-12)
                assert h.attrs['dt']==dt
            assert abs(declared[sid]['chainage_m']-((tx+rx)/2+100))<1e-10
            g=geometry_at(data['H1'][:,:,0],np.array([.025]*3),np.array([tx,35,0]),np.array([rx,35,0]),n95)
            p=primary_paths(g,materials)[g['boundaries'].index(g['basal_sand'])];z,t=inverse(p,FREQ,weights('hann',501));center=float(t[np.argmax(abs(z))]*1e9)
            stations[station]=dict(tx_m=[tx,35,0],rx_m=[rx,35,0],profile_midpoint_m=(tx+rx)/2+100,model_midpoint_m=(tx+rx)/2+80,geometry=g,basal_template_peak_ns=center,basal_gate_ns=[center-12,center+12],original_native_sha256=pair.sha(raw))
        gp=out/name/'geometries';gp.mkdir(parents=True);path=out/name/'cases'/('full2d_'+station)/'profile.in';path.parent.mkdir(parents=True)
        old=parent/('H1' if variant=='H1' else 'H0')/'geometries/full2d_compact.h5';shutil.copyfile(old,gp/old.name);shutil.copyfile(db,gp/db.name)
        if variant=='far':
            with h5py.File(gp/old.name,'r+') as h:h['data'][:]=far
        freq=100e6 if kind=='ricker' else 1.
        lines=[]
        for line in card.splitlines():
            if line.startswith('#title:'):line=f'#title: {name}; source/window/station causal study; local x midpoint={(tx+rx)/2:g}; no snapshots'
            elif line.startswith('#time_window:'):line=f'#time_window: {time*1e-9:.12g}'
            elif line.startswith('#waveform:'):line=f'#waveform: {kind} 40 {freq:.12g} pulse'
            elif line.startswith('#hertzian_dipole:'):line=f'#hertzian_dipole: z {tx:g} 35 0.0125 pulse'
            elif line.startswith('#rx:'):line=f'#rx: {rx:g} 35 0.0125 {name}_rx1 Ez'
            lines.append(line)
        path.write_text('\n'.join(lines)+'\n',encoding='utf-8');samples=int(np.ceil(time*1e-9/dt))+1
        source=excitation(kind,freq,samples,dt)
        spectrum=np.exp(-2j*np.pi*FREQ[:,None]*((np.arange(samples)+.5)*dt))@source
        support=float(abs(spectrum).min()/abs(spectrum).max())
        assert support>.01
        groups.append(dict(id=name,station=station,variant=variant,input=path.relative_to(out).as_posix(),input_sha256=pair.sha(path),geometry=(gp/old.name).relative_to(out).as_posix(),geometry_sha256=pair.sha(gp/old.name),
            material=(gp/db.name).relative_to(out).as_posix(),material_sha256=pair.sha(db),native_shape=[4400,1600,1],domain_m=[110,40,.025],tx_m=[tx,35.,0.],rx_m=[rx,35.,0.],
            time_window_ns=time,expected_samples=samples,source_type=kind,source_frequency_Hz=freq,source_amplitude_A=40.,source_band_min_over_max=support))
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        path=out/g['input'];get_user_objects(path.read_text('utf-8').splitlines(),input_dir=path.parent)
    manifest=dict(status='PREPARED_APPROVED_NOT_RUN',approval_basis='User persistent objective 研究明白 plus autonomous SSH simulation authorization.',parent_manifest_sha256=pair.sha(parent/'manifest.json'),parent_H0_native_sha256=pair.sha(reference),groups=groups,stations=stations,
        variables='At c0008/H0,complete2x2impulse/Ricker100MHz versus800/1600ns; reuse prior impulse800. At c0058/c0108,matchedH1,H0,farROI45-60m removal.',
        invariants='110x40m,2.5cm,FP64,HORIPML80cells,material database/interface averaging,Tx-Rx1.3m at y35m; no snapshots; exact501tones20-170MHz/.3MHz and actual-source normalization,Hann+Blackman,no tail taper,AGC,alignment or fitted phase.',
        gates='Each station basal model-template peak+-12ns frozen before solve; source/window c0008 reuses332.311-356.311ns.',
        snapshot_count=0,limits='Three informed horizontal stations do not certify entire line or finite3D/antenna/field performance. FarROI deletion adds edges and changes interactions; its role varies with station, including near-underfoot atc0108. No production removal of true geological reflections.')
    pair.save(out/'manifest.json',manifest);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',manifest)


def audit(path,completed=False):
    c=json.loads(path.read_text('utf-8'));pair.check_files(c);rows=[]
    with h5py.File(Path(c['package'])/'reference_H0.h5') as ref:dt=float(ref.attrs['dt'])
    for g in c['groups']:
        r=dict(id=g['id'],input_sha256=pair.sha(g['input']))
        if completed:
            raw=Path(g['input']).with_suffix('.h5')
            with h5py.File(raw) as h:
                assert str(h.attrs['gprMax'])=='4.0.0' and h.attrs['dt']==dt and h.attrs['Iterations']==g['expected_samples']
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-12,rtol=0)
                x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];ea=h['srcs/src1/excitation'].attrs
                assert x.dtype==s.dtype==np.float64 and np.isfinite(x).all() and np.isfinite(s).all() and x.shape==s.shape==(g['expected_samples'],)
                assert ea['WaveformType']==g['source_type'] and ea['WaveformFrequency']==g['source_frequency_Hz'] and ea['SpatialScale']==.025 and ea['Polarisation']=='z'
                np.testing.assert_allclose(s,excitation(g['source_type'],g['source_frequency_Hz'],len(s),dt),atol=1e-12,rtol=1e-12)
                r.update(native_sha256=pair.sha(raw),native_path=str(raw),dtype=str(x.dtype),dt_s=dt,samples=len(x),source_type=str(ea['WaveformType']))
        rows.append(r)
    return dict(status='PASS_NATIVE_SOURCE_GRID_IDENTITY_NOT_PHYSICAL_CLASSIFICATION',completed=completed,contract_sha256=pair.sha(path),groups=rows)


def freeze(package,out,parent):
    boundary.audit=audit;boundary.freeze(package,out,parent)
    path=out/'execution_contract.json';c=json.loads(path.read_text('utf-8'));c['code_identities'][str(Path(__file__).resolve())]=pair.sha(__file__)
    pair.save(path,c);pair.save(out/'preflight_verification.json',audit(path));print(json.dumps(dict(final_contract_sha256=pair.sha(path),groups=len(c['groups']))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for key in ['parent','reference','original','out','evidence','package']:p.add_argument('--'+key,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.parent,a.reference,a.original,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':supervisor.audit=audit;supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
