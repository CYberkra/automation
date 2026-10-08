"""Independent geometry/events/native transforms for completed v5 source/flat controls."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.private/'execution_contract.json');v=read(a.private/'completed_verification.json');r=read(a.public/'analysis.json');m=c['study_manifest'];old=read(a.prior/'execution_contract.json')
    assert sha(a.private/'execution_contract.json')==v['contract_sha256']==r['contract_sha256']
    assert sha(a.private/'completed_verification.json')==r['verification_sha256']
    assert sha(a.prepared/'manifest.json')==r['prepared_manifest_sha256']==c['file_identities'][c['package']+'\\manifest.json']
    assert r['script_sha256']==sha(Path(__file__).parent/'analyze_line9_v5_source_flat_controls.py')
    assert sha(a.prior/'execution_contract.json')==c['parent_execution_contract_sha256'] and sha(a.prior/'completed_verification.json')==r['prior_verification_sha256']
    for path,digest in c['code_identities'].items():assert sha(Path(__file__).parent/Path(path).name)==digest
    events=[json.loads(s) for s in (a.private/'execution.jsonl').read_text('utf-8').splitlines()]
    starts=[e for e in events if e['status']=='STARTED' and 'group' in e];done=[e for e in events if e['status']=='COMPLETED' and 'group' in e]
    names=['nonflat_ricker_H0','flat_ricker_H0','flat_ricker_H1']
    assert [e['group'] for e in starts]==[e['group'] for e in done]==[g['id'] for g in c['groups']]==names and events[-1]['status']=='COMPLETED' and events[-1]['traces']==3
    assert events[-1]['verification_sha256']==sha(a.private/'completed_verification.json') and events[0]['contract_sha256']==sha(a.private/'execution_contract.json')
    assert 'reused_groups' not in c
    arrays={};cards=[];src=[];native=[]
    with h5py.File(a.private/'sfcw.h5') as h:f=h['frequency_Hz'][:];response=h['response'][:];t=h['time_s'][:];planar=h['planar_response'][:]
    np.testing.assert_array_equal(f,20e6+np.arange(501)*300000)
    errors=[]
    for j,(g,row,event) in enumerate(zip(c['groups'],v['groups'],done)):
        prepared_group=next(x for x in m['groups'] if x['id']==g['id'])
        for key in ['input','geometry','material']:assert sha(a.prepared/prepared_group[key])==g[key+'_sha256']
        with h5py.File(a.prepared/g['geometry']) as h:
            assert set(h)=={'data','material_keys'};arrays[g['id']]=h['data'][:]
            if j==0:material_keys=h['material_keys'][:]
            else:np.testing.assert_array_equal(h['material_keys'][:],material_keys)
        cards.append((a.prepared/prepared_group['input']).read_text('utf-8').splitlines()[1:])
        raw=a.private/(g['id']+'.h5');assert sha(raw)==row['native_sha256']==event['raw_sha256']==r['native_sha256'][j+2]
        with h5py.File(raw) as h:
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];dt=float(h.attrs['dt']);ea=h['srcs/src1/excitation'].attrs
            assert str(h.attrs['gprMax'])=='4.0.0' and x.dtype==s.dtype==np.float64 and np.isfinite(x).all() and np.isfinite(s).all() and x.shape==s.shape==(20352,)
            assert dt==m['dt_s'] and ea['WaveformType']=='ricker' and ea['WaveformFrequency']==100e6 and ea['WaveformAmplitude']==40 and ea['SpatialScale']==.025 and ea['Polarisation']=='z'
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
            k=np.array([0,29,107,201,373,500]);manual=(np.exp(-2j*np.pi*f[k,None]*(np.arange(len(x))*dt))@x)/(np.exp(-2j*np.pi*f[k,None]*((np.arange(len(s))+.5)*dt))@s)/.025
            error=float(np.linalg.norm(manual-response[k,j+2])/np.linalg.norm(manual));assert error<1e-9;errors.append(error);src.append(s);native.append(x)
    assert cards[0]==cards[1]==cards[2]
    for s in src[1:]:np.testing.assert_array_equal(src[0],s)
    parent_groups={g['id']:g for g in [*old.get('reused_groups',[]),*old['groups']]}
    assert c['groups'][0]['geometry_sha256']==parent_groups['base_H0']['geometry_sha256']
    original_geometry=a.parent/'base_H1/geometries/full2d_compact.h5'
    assert sha(original_geometry)==parent_groups['base_H1']['geometry_sha256']
    with h5py.File(original_geometry) as h:original=h['data'][:];np.testing.assert_array_equal(h['material_keys'][:],material_keys)
    flat1=arrays['flat_ricker_H1'];flat0=arrays['flat_ricker_H0'];column=m['flat_column_index']
    np.testing.assert_array_equal(flat1,np.broadcast_to(original[column:column+1],flat1.shape))
    mask=flat1==3;changed_voxels=int(mask.sum());assert np.all(flat0[mask]==2);np.testing.assert_array_equal(flat0[~mask],flat1[~mask])
    assert not np.any(flat0==3) and np.all(flat1==flat1[:1]) and np.all(flat0==flat0[:1])
    np.testing.assert_allclose((np.flatnonzero(np.diff(flat1[0,:,0]))+1)*.025,[17.1,24.1,31.4],rtol=0,atol=1e-12)
    with h5py.File(a.prepared/'reference_impulse_H0.h5') as h:imp=h['rxs/rx1/Ez'][:];assert imp.dtype==np.float64
    assert sha(a.prepared/'reference_impulse_H0.h5')==m['reference_impulse_H0_sha256']==r['native_sha256'][1]
    # Independent selected-index direct convolution, separate from FFT algorithm.
    selected=np.r_[0,1,np.arange(5500,7500,61),18000,20351]
    predicted=np.array([np.dot(imp[:k+1],src[0][k::-1])/40 for k in selected]);error=float(np.linalg.norm(predicted-native[0][selected])/np.linalg.norm(native[0][selected]));assert error<1e-9
    np.testing.assert_array_equal(response[:,5],response[:,0]-response[:,1]);np.testing.assert_array_equal(response[:,6],response[:,4]-response[:,3]);np.testing.assert_array_equal(response[:,7],response[:,2]-response[:,1]);np.testing.assert_array_equal(response[:,8],response[:,3]-response[:,2])
    assert sha(a.private/'sfcw.h5')==r['numerical_sha256'];inverse={};metric_reconstruction_bounds={}
    with h5py.File(a.private/'sfcw.h5') as h:
        for window in ['hann','blackman']:
            w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();z=h[window+'_complex_bandpass'][:];take=np.arange(17,len(t),127)
            independent=np.exp(2j*np.pi*t[take,None]*f)@(w[:,None]*response)/501;err=float(np.linalg.norm(independent-z[take])/np.linalg.norm(independent));assert err<1e-9;inverse[window]=err
            padded=np.zeros((len(t),6),complex);padded[:501]=w[:,None]*planar
            reference_planar=np.fft.ifft(padded,axis=0)*(len(t)/501)*np.exp(2j*np.pi*f[0]*t[:,None])
            window_bounds={}
            for gate,row in r['metrics'][window].items():
                lo,hi=row['gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi);q=z[keep];den=np.linalg.norm(q[:,5])
                assert abs(np.linalg.norm(q[:,7])/den-row['source_change_over_nonflat_delta'])<1e-7
                pz=np.exp(2j*np.pi*t[keep,None]*f)@(w[:,None]*planar)/501
                for field,index,model,reference_model in [('flat_delta_vs_planar',6,pz[:,5],reference_planar[keep,5]),('flat_H0_vs_planar',3,pz[:,:5].sum(axis=1),reference_planar[keep,:5].sum(axis=1))]:
                    # Reverse triangle inequality bounds changes of the norm metric.
                    # Tiny late profiles amplify inverse-transform roundoff; do not
                    # relax a fixed tolerance or call this a physical error floor.
                    denominator=np.linalg.norm(q[:,index]);reconstruction=np.linalg.norm(model-reference_model)
                    bound=float(reconstruction/denominator+64*np.finfo(float).eps)
                    measured=float(np.linalg.norm(q[:,index]-model)/denominator)
                    assert abs(measured-row[field]['unfitted_relative_L2'])<=bound
                    inner=np.vdot(reference_model,q[:,index]);perturbation=reconstruction*denominator
                    phase_bound=float(2*perturbation/(abs(inner)-perturbation)+64*np.finfo(float).eps)
                    phase=float(np.angle(np.vdot(model,q[:,index]),deg=True))
                    assert abs(np.exp(1j*np.deg2rad(phase))-np.exp(1j*np.deg2rad(row[field]['inner_phase_deg'])))<=phase_bound
                    window_bounds[gate+'_'+field]=bound
            metric_reconstruction_bounds[window]=window_bounds
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:shutil.copyfile(a.private/name,a.public/name)
    result=dict(status='PASS_THREE_NATIVE_GEOMETRY_EVENTS_DIRECT_CONVOLUTION_DFT_AND_INVERSE_NOT_FIELD_CERTIFICATION',audit_script_sha256=sha(__file__),solver_started=len(starts),solver_completed=len(done),flat_basal_changed_voxels=changed_voxels,independent_DFT_relative_L2=errors,independent_direct_convolution_relative_L2=error,independent_inverse_relative_L2=inverse,model_metric_inverse_roundoff_bounds=metric_reconstruction_bounds,completed_models=[dict(id=e['group'],elapsed_s=e['elapsed_s'],peak_owned_RSS_GiB=e['peak_owned_RSS_bytes']/2**30,native_sha256=e['raw_sha256']) for e in done],plot_sha256={p.name:sha(p) for p in a.public.glob('*.png')})
    (a.public/'delivery_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','public','prepared','prior','parent']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
