"""Three bounded v5 controls: Ricker nonflat H0 and matched planar H0/H1."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import line9_basal_pair as pair
import line9_boundary_controls as boundary
import hs4_station_grid_controls as supervisor
from line9_source_station_controls import excitation


def prepare(parent, reference, out, evidence):
    assert not out.exists() and not evidence.exists(), 'Fresh outputs required'
    old = json.loads((parent/'manifest.json').read_text('utf-8'))
    assert old['source_package']=='line9_pkg1_v5_d02m_194st'
    records={g['id']:g for g in old['groups']}
    arrays={}
    for name in ['base_H1','base_H0']:
        g=records[name]
        for key in ['input','geometry','material']:assert pair.sha(parent/g[key])==g[key+'_sha256']
        with h5py.File(parent/g['geometry']) as h:arrays[name]=h['data'][:]
    with h5py.File(reference) as h:
        assert h['rxs/rx1/Ez'].dtype==np.float64 and str(h.attrs['gprMax'])=='4.0.0'
        dt=float(h.attrs['dt']);n=int(h.attrs['Iterations'])
        assert dt==old['dt_s'] and n==old['expected_samples']==20352
        for key,pos in [('srcs/src1',records['base_H0']['tx_m']),('rxs/rx1',records['base_H0']['rx_m'])]:
            np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
    column=round(np.mean([records['base_H0']['tx_m'][0],records['base_H0']['rx_m'][0]])/.025)
    flat1=np.repeat(arrays['base_H1'][column:column+1],8400,axis=0)
    assert flat1.shape==(8400,1700,1)
    changes=np.flatnonzero(np.diff(flat1[0,:,0]))+1
    np.testing.assert_allclose(changes*.025,[17.1,24.1,31.4],rtol=0,atol=1e-12)
    flat0=flat1.copy();flat0[pair.basal_mask(flat1)]=2
    values={'nonflat_ricker_H0':arrays['base_H0'],'flat_ricker_H0':flat0,'flat_ricker_H1':flat1}
    out.mkdir(parents=True);shutil.copyfile(reference,out/'reference_impulse_H0.h5')
    original=(parent/records['base_H0']['input']).read_text('utf-8')
    assert original.count('#waveform: impulse 40 1 pulse')==1 and '#snapshot:' not in original
    groups=[];geo=parent/records['base_H1']['geometry'];db=parent/records['base_H0']['material']
    from audit_line9_postprocessing import FREQ
    source=excitation('ricker',100e6,n,dt)
    spectrum=np.exp(-2j*np.pi*FREQ[:,None]*((np.arange(n)+.5)*dt))@source
    support=float(abs(spectrum).min()/abs(spectrum).max());assert support>.01
    for name,data in values.items():
        gp=out/name/'geometries';gp.mkdir(parents=True)
        p=out/name/'cases/full2d_r0098/profile.in';p.parent.mkdir(parents=True)
        if name=='nonflat_ricker_H0':shutil.copyfile(parent/records['base_H0']['geometry'],gp/geo.name)
        else:
            with h5py.File(geo) as src,h5py.File(gp/geo.name,'x') as h:
                for key,value in src.attrs.items():h.attrs[key]=value
                for key in src:
                    if key!='data':src.copy(key,h)
                h.create_dataset('data',data=data,compression='gzip',compression_opts=4)
                for key,value in src['data'].attrs.items():h['data'].attrs[key]=value
        with h5py.File(geo) as src,h5py.File(gp/geo.name) as h:
            assert set(h)==set(src)=={'data','material_keys'}
            np.testing.assert_array_equal(h['material_keys'][:],src['material_keys'][:])
            np.testing.assert_array_equal(h['data'][:],data)
        shutil.copyfile(db,gp/db.name)
        lines=original.replace('#waveform: impulse 40 1 pulse','#waveform: ricker 40 100000000 pulse').splitlines()
        lines[0]=f'#title: {name}; Line9 v5 single station source/flat controls; no snapshots'
        p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        g=dict(records['base_H0']);g.update(id=name,input=p.relative_to(out).as_posix(),input_sha256=pair.sha(p),geometry=(gp/geo.name).relative_to(out).as_posix(),geometry_sha256=pair.sha(gp/geo.name),material=(gp/db.name).relative_to(out).as_posix(),source_type='ricker',source_frequency_Hz=100e6,source_amplitude_A=40.,source_band_min_over_max=support)
        groups.append(g)
    from gprMax.hash_cmds_file import get_user_objects
    for g in groups:
        p=out/g['input'];get_user_objects(p.read_text('utf-8').splitlines(),input_dir=p.parent)
    m=dict(status='PREPARED_APPROVED_NOT_RUN',approval_basis='Persistent user research objective and autonomous SSH simulation authorization; three bounded v5 controls.',groups=groups,parent_manifest_sha256=pair.sha(parent/'manifest.json'),reference_impulse_H0_sha256=pair.sha(reference),script_sha256=pair.sha(__file__),flat_column_index=column,flat_interfaces_y_m=[17.1,24.1,31.4],geometry_diagnostic=old['geometry_diagnostic'],dt_s=dt,expected_samples=n,basal_gate_ns=old['basal_gate_ns'],deep_gate_ns=old['deep_gate_ns'],early_gate_ns=old['early_gate_ns'],snapshot_count=0,
        invariants='210x42.5m,2.5cm,1200ns,FP64,source/rx positions and40A amplitude,complete material spectra and interface averaging,HORIPML80cells. Nonflat retains exact originalH0 voxels; flatH1 repeats midpoint column,flatH0 changes only bottom sandstone to mud.',
        processing='Exact501tones20-170MHz/.3MHz,actual-source/Yee-time/spatial normalization,Hann+Blackman,noAGC/tail taper/time alignment/amplitude or phase fit.',
        limits='Single selected v5 station. Planar geometry is a controlled different model, not corrected field truth or nonflat grid convergence. Comparison with continuum Green function still includes finite-domain,grid and interface-averaging effects. No field/3D/whole-line certification.')
    pair.save(out/'manifest.json',m);evidence.mkdir(parents=True);pair.save(evidence/'preparation.json',m)
    print(json.dumps(dict(groups=3,column=column,source_band_min_over_max=support)))


def audit(path,completed=False):
    c=json.loads(path.read_text('utf-8'));pair.check_files(c);m=c['study_manifest'];rows=[]
    for g in c['groups']:
        row=dict(id=g['id'],input_sha256=pair.sha(g['input']))
        if completed:
            raw=Path(g['input']).with_suffix('.h5')
            with h5py.File(raw) as h:
                assert str(h.attrs['gprMax'])=='4.0.0' and h.attrs['dt']==m['dt_s'] and h.attrs['Iterations']==g['expected_samples']
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
                x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];ea=h['srcs/src1/excitation'].attrs
                assert x.dtype==s.dtype==np.float64 and x.shape==s.shape==(g['expected_samples'],) and np.isfinite(x).all() and np.isfinite(s).all()
                assert ea['WaveformType']=='ricker' and ea['WaveformFrequency']==100e6 and ea['WaveformAmplitude']==40 and ea['SpatialScale']==.025 and ea['Polarisation']=='z'
                np.testing.assert_allclose(s,excitation('ricker',100e6,len(s),m['dt_s']),atol=1e-12,rtol=1e-12)
                for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
            row.update(native_path=str(raw),native_sha256=pair.sha(raw),dtype=str(x.dtype),samples=len(x),dt_s=m['dt_s'])
        rows.append(row)
    return dict(status='PASS_NATIVE_SOURCE_GRID_IDENTITY_NOT_CONTINUUM_OR_FIELD_CERTIFICATION',completed=completed,contract_sha256=pair.sha(path),groups=rows)


def freeze(package,out,parent):
    boundary.audit=audit;boundary.freeze(package,out,parent)
    p=out/'execution_contract.json';c=json.loads(p.read_text('utf-8'))
    for key in ['reused_groups','prior_execution_contract_sha256','prior_execution_log_sha256','recovery_basis']:c.pop(key,None)
    for name in ['line9_v5_source_flat_controls.py','line9_source_station_controls.py']:c['code_identities'][str(pair.ROOT/'scripts'/name)]=pair.sha(pair.ROOT/'scripts'/name)
    pair.save(p,c);pair.save(out/'preflight_verification.json',audit(p));print(json.dumps(dict(contract_sha256=pair.sha(p),groups=len(c['groups']))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify'])
    for key in ['parent','reference','out','evidence','package']:p.add_argument('--'+key,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.parent,a.reference,a.out,a.evidence)
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve(),a.parent.resolve())
    elif a.action=='run':supervisor.audit=audit;supervisor.run(a.out/'execution_contract.json')
    else:print(json.dumps(audit(a.out/'execution_contract.json',True)))
