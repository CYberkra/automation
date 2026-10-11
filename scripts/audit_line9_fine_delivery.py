"""Independent retrieved native/probe/key-frame and official-band audit."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha


def read(p):
    return json.loads(p.read_text('utf-8'))


def audit(root,package,arrays,public):
    c=read(root/'execution_contract.json');v=read(root/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(root/'execution_contract.json')
    assert len(v['groups'])==len(c['groups'])==6 and c['no_retry'] and c['max_runs']==6
    assert c['max_batch_wall_s']==10800 and c['solver_entrypoint']=='gprmax_fine_snapshot_entry.py'
    m=read(package/'manifest.json');assert c['study_manifest']==m
    events=[json.loads(s) for s in (root/'execution.jsonl').read_text('utf-8').splitlines()]
    assert not any(e['status'] in ['FAILED','GUARD_STOP'] for e in events)
    done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    assert [e['group'] for e in done]==[g['id'] for g in c['groups']]
    started=[e for e in events if e['status']=='STARTED' and 'group' in e]
    assert [e['group'] for e in started]==[g['id'] for g in c['groups']]
    assert events[-1]['traces']==6 and events[-1]['verification_sha256']==sha(root/'completed_verification.json')
    rows=[];frequencies=20e6+np.arange(501)*300000.
    for g,cert in zip(c['groups'],v['groups']):
        native=root/'native'/g['id']/'profile.h5';observer=root/'native'/g['id']/'observer_verification.json'
        ov=read(observer);assert sha(observer)==cert['observer_verification_sha256']
        assert sha(native)==cert['native_sha256']==ov['native_sha256']
        assert ov['main_source_receiver_bitwise_equal'] and ov['snapshot_probe_max_normalised_abs_error']==0
        parent=package/g['baseline']/'native.h5'
        with h5py.File(native) as h,h5py.File(parent) as old:
            for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
                np.testing.assert_array_equal(h[key][:],old[key][:]);assert h[key][:].tobytes()==old[key][:].tobytes()
            assert len(h['rxs'])==g['receiver_count']
            for p in g['probes']:
                for anchor in p['anchors']:
                    rr=h[f"rxs/rx{anchor['receiver_index']}"]
                    np.testing.assert_array_equal(rr.attrs['GridPosition'],anchor['coord'])
                    for field in anchor['outputs']:
                        ds=rr[field];assert ds.dtype==np.float64 and ds.shape==(40703,)
                        assert np.isfinite(ds[:]).all()
                        assert ds.attrs['TimeSampleOffset']==(0 if field=='Ez' else -.5*g['dt_s'])
            keys=root/'keys'/g['id'];expected={r['file']:r for r in ov['snapshots']};count=0
            for p in sorted(keys.glob('*.h5')):
                row=expected[p.name];assert sha(p)==row['sha256'];n=row['iteration']
                with h5py.File(p) as snap:
                    assert snap.attrs['iteration']==n and snap.attrs['time']==n*g['dt_s']
                    assert set(snap)=={'Ez','Hx','Hy'}
                    for field in snap:
                        z=snap[field][:];assert z.dtype==np.float64 and list(z.shape)==g['roi']['shape'] and np.isfinite(z).all()
                        for probe in g['probes']:
                            a=[h[f"rxs/rx{x['receiver_index']}"] for x in probe['anchors']]
                            i=round((probe['x_m']-g['roi']['bounds'][0])/.025);j=round((probe['y_m']-g['roi']['bounds'][1])/.025)
                            if field=='Ez':want=(a[0]['Ez'][n]+a[3]['Ez'][n]+a[4]['Ez'][n]+a[5]['Ez'][n])*.25
                            elif field=='Hx':want=(a[0]['Hx'][n]+a[3]['Hx'][n])*.5
                            else:want=(a[0]['Hy'][n]+a[4]['Hy'][n])*.5
                            assert z[i,j,0]==want
                count+=1
            assert count==(8 if g['kind']=='local' else 7)
            with h5py.File(arrays/(g['id']+'_arrays.h5')) as a:
                np.testing.assert_array_equal(a['frequency_Hz'][:],frequencies)
                assert a.attrs['source_native_sha256']==sha(native)
                source=sf.load_source(native);receiver=sf.load_receiver(native,'/rxs/rx1','Ez')
                picks=np.r_[np.arange(0,501,37),500];num=[]
                for f in frequencies[picks]:
                    num.append((np.exp(-2j*np.pi*f*receiver.times)@receiver.samples)/
                        (np.exp(-2j*np.pi*f*source.times)@source.samples)/.025)
                actual=a['main_response'][:];err=float(np.linalg.norm(actual[picks]-num)/np.linalg.norm(num));assert err<1e-9
                p=sf.reconstruct_time_response(sf.direct_frequency_response(source,receiver,frequencies,tail_taper_fraction=0),window='hann',zero_pad_factor=8,time_shift=0)
                np.testing.assert_allclose(a['hann/main_complex_bandpass'][:],p.complex_bandpass/.025,atol=1e-12,rtol=1e-12)
                if g['probes']:
                    i=len(g['probes'])-1
                    e=a['native_staggered/Ez'][:,i];hx=a['native_staggered/Hx'][:,i];hy=a['native_staggered/Hy'][:,i]
                    np.testing.assert_array_equal(a['native_aligned/Sy'][:,i],e[:-1]*(hx[:-1]+hx[1:])*.5)
                    np.testing.assert_array_equal(a['native_aligned/Sx'][:,i],-e[:-1]*(hy[:-1]+hy[1:])*.5)
                    for w in ['hann','blackman']:
                        be=a[w+'/probes/Ez'][:,i];bhx=a[w+'/probes/Hx'][:,i];bhy=a[w+'/probes/Hy'][:,i]
                        np.testing.assert_array_equal(a[w+'/probes/Sy_index'][:,i],np.real(be*np.conj(bhx)))
                        np.testing.assert_array_equal(a[w+'/probes/Sx_index'][:,i],-np.real(be*np.conj(bhy)))
        rows.append(dict(id=g['id'],native_sha256=sha(native),retrieved_key_frames=count,main_independent_DFT_relative_L2=err,
                         all_snapshots_remote_audited=ov['snapshot_count'],original_main_bitwise_equal=True))
    analysis=read(public/'analysis.json');assert analysis['contract_sha256']==sha(root/'execution_contract.json')
    for row in analysis['groups']:
        assert row['private_arrays_sha256']==sha(arrays/(row['id']+'_arrays.h5'))
    for x,gate in [(170.,[220,270]),(190.,[300,350])]:
        events=read(public/f'x{int(x)}_probe_events.json')['main_peak_directions']
        event=next(e for e in events if e['event']=='upper_return');assert event['gate_ns']==gate
        g=next(g for g in c['groups'] if g['profile_x_m']==x and g['kind']=='local' and g['role']=='H0')
        parent=package/g['baseline']/'native.h5'
        fq=sf.direct_frequency_response(sf.load_source(parent),sf.load_receiver(parent,'/rxs/rx1','Ez'),frequencies,tail_taper_fraction=0)
        profile=sf.reconstruct_time_response(fq,window='hann',zero_pad_factor=8,time_shift=0)
        t=profile.time*1e9;k=np.flatnonzero((t>=gate[0])&(t<=gate[1]));j=k[np.argmax(abs(profile.complex_bandpass[k]))]
        assert event['time_ns']==float(t[j]);np.testing.assert_allclose(event['peak_abs'],abs(profile.complex_bandpass[j])/.025,rtol=1e-12)
    from PIL import Image
    render=read(public/'snapshot_render_receipts.json')['receipts']
    reductions={r['file']:r for r in read(public/'preview_reduction_receipts.json')['files']};movies=[]
    for p in sorted(public.glob('*.gif')):
        original=next(r for r in render if r.get('movie')==p.name);row=reductions.get(p.name,original)
        if p.name in reductions:assert row['source_sha256']==original['sha256'] and row['temporal_sequence_unchanged_from_preview']
        assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
        with Image.open(p) as im:
            assert im.n_frames==row['frames'] and im.width>0 and im.height>0
            movies.append(dict(file=p.name,frames=im.n_frames,pixel_size=list(im.size),sha256=sha(p)))
    assert len(movies)==10
    return dict(status='PASS_INDEPENDENT_NATIVE_PROBE_KEY_FRAME_AND_BAND_DELIVERY',solver_runs=6,
        contract_sha256=sha(root/'execution_contract.json'),groups=rows,
        R2_events_match_original_parent=True,movies=movies,
        source_script_sha256=sha(__file__),limits='Numerical delivery audit; not mesh/PML convergence, independent ray decomposition or field validation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['root','package','arrays','public','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();result=audit(a.root,a.package,a.arrays,a.public)
    assert not a.out.exists();a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
