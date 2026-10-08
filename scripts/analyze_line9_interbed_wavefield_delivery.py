"""Local checks and exact SFCW for the completed native Ricker wavefield pair."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from audit_line9_postprocessing import FREQ,inverse,weights,rel
from hs_capsule_identity import sha256 as sha
from analyze_line9_basal_pair import plotting,save
from trace_line9_time_origin import correlation


def main(a):
    assert not a.out.exists() and not a.numerical.exists(),'Fresh outputs required'
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.private/'execution_contract.json');v=read(a.private/'completed_verification.json');m=c['study_manifest']
    render=read(a.private/'render/render_receipt.json');events=[json.loads(s) for s in (a.private/'execution.jsonl').read_text('utf-8').splitlines()]
    assert v['completed'] and sha(a.private/'execution_contract.json')==v['contract_sha256']==render['contract_sha256']==events[0]['contract_sha256']
    assert sha(a.private/'completed_verification.json')==render['verification_sha256']==events[-1]['verification_sha256'] and events[-1]['status']=='COMPLETED'
    assert len([e for e in events if e['status']=='STARTED' and 'group' in e])==len([e for e in events if e['status']=='COMPLETED' and 'group' in e])==2
    assert sha(a.prepared/'manifest.json')==c['file_identities'][c['package']+'\\manifest.json']
    hashes=[];responses=[];errors=[];native=[]
    for g,row in zip(c['groups'],v['groups']):
        p=a.private/(g['id']+'.h5');assert sha(p)==row['native_sha256'] and row['snapshot_count']==799
        for key in ['input','geometry','material']:
            local=next(x for x in m['groups'] if x['id']==g['id'])[key]
            assert sha(a.prepared/local)==g[key+'_sha256']
        assert len(row['snapshots'])==799
        np.testing.assert_array_equal([r['iteration'] for r in row['snapshots']],m['snapshot_iterations'])
        source=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez')
        with h5py.File(p) as h:
            assert str(h.attrs['gprMax'])=='4.0.0' and h['rxs/rx1/Ez'].dtype==np.float64
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[4400,1600,1]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert h['srcs/src1/excitation'].attrs['WaveformType']=='ricker' and h['srcs/src1/excitation'].attrs['WaveformFrequency']==100e6
        value=sf.direct_frequency_response(source,rx,FREQ,tail_taper_fraction=0);assert value.source_valid.all() and source.spatial_scale==.025
        response=value.response.reshape(501)/.025;k=np.r_[0,np.arange(39,501,67),500]
        manual=rx.dt*(np.exp(-2j*np.pi*FREQ[k,None]*rx.times)@rx.samples)/(source.dt*(np.exp(-2j*np.pi*FREQ[k,None]*source.times)@source.samples))/.025
        error=rel(response[k],manual);assert error<1e-9
        hashes.append(sha(p));responses.append(response);errors.append(error);native.append(rx.samples)
    assert hashes[0]==sha(a.source_results/'ricker_800.h5'),'H0 observers changed byte identity'
    prior=read(a.localized/'completed_verification.json');row=next(r for r in prior['groups'] if r['id']=='far_interbed_removed')
    p=a.localized/'far_interbed_removed.h5';assert sha(p)==row['native_sha256']
    source=sf.load_source(p);rx=sf.load_receiver(p,'/rxs/rx1','Ez')
    old=sf.direct_frequency_response(source,rx,FREQ,tail_taper_fraction=0).response.reshape(501)/.025
    response=np.column_stack(responses+[responses[0]-responses[1],old]);a.out.mkdir(parents=True)
    plt=plotting();metrics={};profiles={};checks={}
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for k,window in enumerate(['hann','blackman']):
        w=weights(window,501);z,t=inverse(response,FREQ,w);profiles[window]=z;take=np.arange(0,len(t),97)
        manual=np.exp(2j*np.pi*t[take,None]*FREQ)@(w[:,None]*response)/501
        checks[window]=rel(z[take],manual);assert checks[window]<1e-9
        mask=(t*1e9>=332.31137724550894)&(t*1e9<=356.31137724550894)
        metrics[window]=dict(far_remaining_L2_over_H0=float(np.linalg.norm(z[mask,1])/np.linalg.norm(z[mask,0])),difference_L2_over_H0=float(np.linalg.norm(z[mask,2])/np.linalg.norm(z[mask,0])),
            far_ricker_change_over_prior_impulse=float(np.linalg.norm(z[mask,1]-z[mask,3])/np.linalg.norm(z[mask,3])),far_ricker_prior_correlation=correlation(z[mask,1],z[mask,3]),
            difference_peak_ns=float(t[np.flatnonzero(mask)[np.argmax(abs(z[mask,2]))]]*1e9))
        show=(t*1e9>=200)&(t*1e9<=410);lim=float(abs(z[show,:3].real).max());ax=axes[k,0]
        im=ax.imshow(z[show,:3].real,extent=[-.5,2.5,410,200],aspect='auto',cmap='gray',vmin=-lim,vmax=lim)
        ax.set_xticks([0,1,2],['原H0','只替换夹层','原H0−替换场'])
        ax.set(ylabel='SFCW时间 / ns',title=f'{window}：同站三列共享色标，非连续测线')
        fig.colorbar(im,ax=ax,label='双极带通信号 / (V/m)/(A·m)')
        ax=axes[k,1]
        for j,label in enumerate(['原H0','只替换夹层','原H0−替换场','历史冲激替换场']):ax.plot(t*1e9,abs(z[:,j]),lw=1,label=label)
        ax.axvspan(332.31137724550894,356.31137724550894,color='red',alpha=.08)
        zoom=(t*1e9>=300)&(t*1e9<=385)
        ax.set(xlim=(300,385),ylim=(0,1.05*float(abs(z[zoom]).max())),xlabel='SFCW时间 / ns',ylabel='复数包络 / (V/m)/(A·m)',title=f'{window}：目标窗共同纵轴，无AGC/移时/拟合');ax.legend(fontsize=9)
    fig.suptitle('密集快照批：X179.25m，AGL10.275m；精确501点源归一化SFCW，差分先复数后包络')
    fig.savefig(a.out/'wavefield_pair_sfcw_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('time_s',data=t);h.create_dataset('response',data=response)
        for window,z in profiles.items():h.create_dataset(window+'_complex_bandpass',data=z)
        h.attrs['columns']='H0,far_removed,H0_minus_far,prior_impulse_far';h.attrs['units']='(V/m)/(A*m)'
    for file in (a.private/'render').iterdir():
        if file.is_file():shutil.copyfile(file,a.out/file.name)
    assert sha(a.out/'interbed_wavefield.gif')==render['gif_sha256']
    for row in render['static_frames']:assert sha(a.out/row['file'])==row['sha256']
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl']:shutil.copyfile(a.private/name,a.out/name)
    result=dict(status='PASS_TWO_NATIVE_CASES_AND_RENDER_RECEIPT_NOT_UNIQUE_PATH_OR_FULL_LINE_CERTIFICATION',script_sha256=sha(__file__),contract_sha256=sha(a.private/'execution_contract.json'),verification_sha256=sha(a.private/'completed_verification.json'),native_sha256=hashes,
        solver_started=2,solver_completed=2,snapshot_count_each=799,H0_byte_identical_to_no_observer_Ricker=True,independent_DFT_relative_L2=errors,independent_inverse_relative_L2=checks,metrics=metrics,numerical_sha256=sha(a.numerical),
        reference_impulse_far_native_sha256=sha(p),render_script_sha256=render['script_sha256'],render_receipt_sha256=sha(a.out/'render_receipt.json'),limits=m['limits'])
    save(a.out/'analysis.json',result);print(json.dumps(dict(status=result['status'],metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['private','source-results','localized','prepared','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
