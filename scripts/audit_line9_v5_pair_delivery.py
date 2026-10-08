"""Independent native/geometry/event and transform checks for the v5 pair study."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.private/'execution_contract.json');v=read(a.private/'completed_verification.json');r=read(a.public/'analysis.json');m=c['study_manifest']
    events=[json.loads(s) for s in (a.private/'execution.jsonl').read_text('utf-8').splitlines()]
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e];done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    names=['base_H1','base_H0','top20_H1','top20_H0']
    groups=[*c.get('reused_groups',[]),*c['groups']];historical_starts=[]
    if c.get('reused_groups'):
        previous=read(a.private/'prior_execution_contract.json');previous_events=[json.loads(s) for s in (a.private/'prior_execution.jsonl').read_text('utf-8').splitlines()]
        assert sha(a.private/'prior_execution_contract.json')==c['prior_execution_contract_sha256'] and sha(a.private/'prior_execution.jsonl')==c['prior_execution_log_sha256']
        assert previous_events[-1]['status']=='FAILED' and previous_events[-1]['error']=='solver failed; preserve attempt'
        assert [x['group'] for x in starts]==[x['group'] for x in done]==[x['id'] for x in c['groups']]==names[1:] and events[-1]['traces']==3
        historical_starts=[x for x in previous_events if x['status']=='STARTED' and 'group' in x]
        assert [x['group'] for x in historical_starts]==['base_H1','base_H0']
        previous_done=[x for x in previous_events if x['status']=='COMPLETED' and 'group' in x];assert len(previous_done)==1 and previous_done[0]['group']=='base_H1'
        assert previous_done[0]['raw_sha256']==c['reused_groups'][0]['reused_native_sha256']
        done=previous_done+done
    else:assert [x['group'] for x in starts]==[x['group'] for x in done]==names and events[-1]['traces']==4
    assert [x['group'] for x in done]==[x['id'] for x in groups]==names
    assert events[-1]['status']=='COMPLETED' and events[-1]['verification_sha256']==r['verification_sha256']==sha(a.private/'completed_verification.json')
    assert events[0]['contract_sha256']==v['contract_sha256']==r['contract_sha256']==sha(a.private/'execution_contract.json')
    assert sha(a.prepared/'manifest.json')==c['file_identities'][c['package']+'\\manifest.json']
    for path,digest in c['code_identities'].items():assert sha(Path(__file__).parent/Path(path).name)==digest
    assert r['script_sha256']==sha(Path(__file__).parent/'analyze_line9_v5_pair_controls.py')
    geometry={};cards={};source_samples=[];direct=[]
    with h5py.File(a.private/'sfcw.h5') as h:response=h['response'][:];f=h['frequency_Hz'][:];t=h['time_s'][:]
    np.testing.assert_array_equal(f,20e6+np.arange(501)*300000)
    for j,(g,row,e) in enumerate(zip(groups,v['groups'],done)):
        original=next(x for x in m['groups'] if x['id']==g['id'])
        for key in ['input','geometry','material']:assert sha(a.prepared/original[key])==g[key+'_sha256']
        with h5py.File(a.prepared/g['geometry']) as h:geometry[g['id']]=h['data'][:]
        cards[g['id']]=(a.prepared/original['input']).read_bytes()
        path=a.private/(g['id']+'.h5');assert sha(path)==row['native_sha256']==e['raw_sha256']==r['native_sha256'][j]
        with h5py.File(path) as h:
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];dt=float(h.attrs['dt']);attrs=h['srcs/src1/excitation'].attrs
            assert str(h.attrs['gprMax'])=='4.0.0' and x.dtype==s.dtype==np.float64 and x.shape==s.shape==(20352,) and np.isfinite(x).all()
            assert dt==m['dt_s'];np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape']);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert attrs['WaveformType']=='impulse' and attrs['WaveformAmplitude']==40 and attrs['WaveformFrequency']==1 and attrs['SpatialScale']==.025
            assert s[0]==40 and np.count_nonzero(s)==1;source_samples.append(s)
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-12,rtol=0)
            k=np.array([0,29,83,151,233,387,500]);rx_time=np.arange(len(x))*dt;src_time=(np.arange(len(s))+.5)*dt
            independent=dt*(np.exp(-2j*np.pi*f[k,None]*rx_time)@x)/(dt*(np.exp(-2j*np.pi*f[k,None]*src_time)@s))/.025
            err=float(np.linalg.norm(independent-response[k,j])/np.linalg.norm(independent));assert err<1e-9;direct.append(err)
    base=geometry['base_H1'];changed=geometry['base_H0'];assert base.shape==changed.shape==(8400,1700,1)
    assert np.all(base[:,0,0]==3);stop=np.argmax(base[:,:,0]!=3,axis=1);mask=(np.arange(1700)[None,:]<stop[:,None])[:,:,None]
    basal_voxel_count=int(mask.sum())
    assert basal_voxel_count==m['changed_basal_voxels'] and np.all(changed[mask]==2);np.testing.assert_array_equal(changed[~mask],base[~mask])
    for name,values in [('top20_H1',base),('top20_H0',changed)]:
        assert geometry[name].shape==(8400,2500,1);np.testing.assert_array_equal(geometry[name][:,:1700],values);assert np.all(geometry[name][:,1700:]==0)
        assert cards[name]==cards['base_H1'].replace(b'#domain: 210 42.5 inf',b'#domain: 210 62.5 inf')
    assert cards['base_H1']==cards['base_H0']
    for s in source_samples[1:]:np.testing.assert_array_equal(s,source_samples[0])
    with h5py.File(a.prepared/'reference_fp32.h5') as h:assert sha(a.prepared/'reference_fp32.h5')==r['native_sha256'][4] and h['rxs/rx1/Ez'].dtype==np.float32
    np.testing.assert_array_equal(response[:,5],response[:,0]-response[:,1]);np.testing.assert_array_equal(response[:,6],response[:,2]-response[:,3]);np.testing.assert_array_equal(response[:,7],response[:,0]-response[:,4])
    assert sha(a.private/'sfcw.h5')==r['numerical_sha256'];inverse={}
    with h5py.File(a.private/'sfcw.h5') as h:
        for window in ['hann','blackman']:
            z=h[window+'_complex_bandpass'][:];w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();take=np.arange(7,len(t),127)
            independent=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*response)/501
            err=float(np.linalg.norm(independent-z[take])/np.linalg.norm(independent));assert err<1e-9;inverse[window]=err
            for gate,row in r['metrics'][window].items():
                lo,hi=row['gate_ns'];q=z[(t*1e9>=lo)&(t*1e9<=hi)];den=np.linalg.norm(q[:,5])
                for index,key in [(7,'FP64_minus_FP32_over_base_delta'),(8,'top20_H1_change_over_base_delta'),(9,'top20_H0_change_over_base_delta'),(10,'top20_delta_change_over_base_delta')]:assert abs(np.linalg.norm(q[:,index])/den-row[key])<1e-11
                assert abs(np.linalg.norm(q[:,1])/den-row['base_H0_over_delta'])<1e-10
                if gate=='basal' and 'planar_response' in h:
                    mask=(t*1e9>=lo)&(t*1e9<=hi);pz=(np.exp(2j*np.pi*t[mask,None]*f)@(w[:,None]*h['planar_response'][:]))/501
                    value=np.linalg.norm(q[:,5]-pz[:,5])/den
                    assert abs(value-row['planar_bottom_unfitted_relative_L2'])<1e-10
                    phase=np.angle(np.vdot(pz[:,5],q[:,5]),deg=True);assert abs(phase-row['planar_bottom_inner_phase_deg'])<1e-9
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl','prior_execution_contract.json','prior_execution.jsonl']:
        if (a.private/name).exists():shutil.copyfile(a.private/name,a.public/name)
    result=dict(status='PASS_FOUR_NATIVE_GEOMETRY_EVENTS_DIRECT_DFT_AND_INDEPENDENT_INVERSE_NOT_WHOLE_LINE_CERTIFICATION',audit_script_sha256=sha(__file__),solver_started=len(starts)+len(historical_starts),solver_completed=4,reused_completed_models=len(c.get('reused_groups',[])),preserved_failed_geometry_attempts=1 if historical_starts else 0,
        independent_DFT_relative_L2=direct,independent_inverse_relative_L2=inverse,basal_changed_voxels=basal_voxel_count,top_added_air_rows=800,
        completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30,native_sha256=e['raw_sha256']) for e in done],
        plot_sha256={p.name:sha(p) for p in a.public.glob('*.png')})
    (a.public/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],DFT=direct,inverse=inverse)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','public','prepared']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
