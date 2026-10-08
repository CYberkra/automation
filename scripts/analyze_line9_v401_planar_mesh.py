"""Fixed-gate domain comparison and phase-preserving mesh-convergence analysis."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import FREQ,inverse,response,relative


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    c=json.loads((a.source/'execution_contract.json').read_text('utf-8'))
    v=json.loads((a.source/'completed_verification.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    assert [g['id'] for g in c['groups']]==['planar_H0','planar_H1']
    dx=c['study_manifest']['spacing_m'];z=[];audit=[]
    for g,row in zip(c['groups'],v['groups']):
        p=a.source/(g['id']+'.h5');assert sha(p)==row['native_sha256']
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1'
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[dx]*3)
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],g['native_shape'])
            assert h['rxs/rx1/Ez'].dtype==np.float64 and len(h['rxs/rx1/Ez'])==g['expected_samples']
            assert abs(float(h.attrs['dt'])/g['dt_s']-1)<1e-13
            assert h['srcs/src1/excitation'].attrs['SpatialScale']==dx
            assert h['srcs/src1/excitation'].attrs['WaveformAmplitude']==40
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:np.testing.assert_allclose(h[key].attrs['Position'],pos,atol=1e-12,rtol=0)
        q,error=response(p,reference_scale=.025);z.append(q)
        audit.append(dict(id=g['id'],sha256=sha(p),dtype='float64',source_actual_spatial_scale_m=dx,
                          common_reference_scale_m=.025,independent_DFT_relative_L2=error))
    z=np.column_stack([*z,z[1]-z[0]])
    with h5py.File(a.reference) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ)
        wide=h['v401_response'][:,[1,2,3]]
    with h5py.File(a.continuum) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);cont=h['planar_response'][:,5]
    spectra=[wide,z];labels=['210m/2.5cm','80m/'+str(dx*100)+'cm']
    if a.coarse:
        with h5py.File(a.coarse) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);coarse=h['response'][:]
        spectra=[wide,coarse,z];labels=['210m/2.5cm','80m/2.5cm','80m/1.25cm']
    bounds=c['study_manifest']['basal_gate_ns'];metrics={};profiles={}
    for window in ['hann','blackman']:
        w=sf.spectral_window(window,501);pv=[]
        for q in spectra:
            p,t=inverse(q,w);pv.append(p)
        pcont,_=inverse(cont[:,None],w);mask=(t*1e9>=bounds[0])&(t*1e9<=bounds[1]);x=pv[-1][mask];r=pv[0][mask]
        ids=np.arange(9,len(t),113);manual=np.exp(2j*np.pi*t[ids,None]*FREQ)@(z*w[:,None])/501
        inv_error=relative(pv[-1][ids],manual);assert inv_error<1e-9
        metrics[window]=dict(basal_delta_change_vs_wide=relative(x[:,2],r[:,2]),
                             basal_H0_change_vs_wide=relative(x[:,0],r[:,0]),
                             unfitted_delta_error_over_fdtd=float(np.linalg.norm(x[:,2]-pcont[mask,0])/np.linalg.norm(x[:,2])),
                             independent_inverse_relative_L2=inv_error)
        profiles[window]=(pv,pcont[:,0])
    selected=[0,83,250,417,500];ratio=z[:,2]/cont
    full_error=relative(z[:,2],wide[:,2]);gate_limits=c['study_manifest']['domain_gate']
    gate_pass=dx==.025 and full_error<=gate_limits['full_band_delta_relative_L2_max'] and all(
        r['basal_delta_change_vs_wide']<=gate_limits['basal_delta_relative_L2_max'] and
        r['basal_H0_change_vs_wide']<=gate_limits['basal_H0_relative_L2_max'] for r in metrics.values())
    a.out.mkdir(parents=True)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z);h.create_dataset('time_s',data=t)
        h.attrs['normalisation']='E/I_line divided by SAME reference .025m for all grids'
    result=dict(status=('PASS_NARROW_PLANAR_DOMAIN_GATE' if gate_pass else 'REJECT_NARROW_PLANAR_DOMAIN_GATE') if dx==.025 else 'COMPLETED_FINE_PLANAR_MESH_DIAGNOSTIC',
                script_sha256=sha(__file__),coarse_spacing_m=.025,width_m=80,spacing_m=dx,
                contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
                reference_numerical_sha256=sha(a.reference),continuum_numerical_sha256=sha(a.continuum),
                numerical_sha256=sha(a.numerical),preceding_coarse_numerical_sha256=sha(a.coarse) if a.coarse else None,
                native=audit,metrics=metrics,full_band_delta_change_vs_wide=full_error,declared_domain_gate=gate_limits,
                frequency_comparison=dict(frequency_MHz=(FREQ[selected]/1e6).tolist(),delta_over_continuum_amplitude=abs(ratio[selected]).tolist(),delta_over_continuum_phase_deg=np.angle(ratio[selected],deg=True).tolist()),
                actual_source_spatial_scale_m=dx,common_reference_scale_m=.025,AGC=False,tail_taper=False,fit=False,
                limits='Planar selected-station numerical evidence only; finite-grid improvement is not full convergence, nonflat/whole-line or field validation.')
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,(pv,pc) in profiles.items():
        fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        mask=(t*1e9>=340)&(t*1e9<=405)
        for j,name in [(0,'H0剩余场'),(2,'H1−H0底砂差场')]:
            ax=axs[0,0 if j==0 else 1];data=np.column_stack([q[mask,j] for q in pv]).real;lim=float(abs(data).max())
            im=ax.imshow(data,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[-.5,len(pv)-.5,405,340]);ax.set_xticks(range(len(pv)),labels,rotation=12)
            ax.axhline(371.756487,color='#20a06a',ls=':',lw=1);ax.set(title=name+'：各配置共用绝对灰度',ylabel='SFCW时间 / ns');fig.colorbar(im,ax=ax,label='共同参考单位 (V/m)/(A·m)')
        for q,label in zip(pv,labels):axs[1,0].plot(t[mask]*1e9,q[mask,2].real,lw=1,label=label)
        axs[1,0].plot(t[mask]*1e9,pc[mask].real,'k--',lw=1,label='连续平层二维解析')
        axs[1,0].set(xlabel='SFCW时间 / ns',ylabel='底砂双极幅度',title='不做幅相/移时拟合');axs[1,0].legend(fontsize=8)
        for q,label in zip(spectra,labels):axs[1,1].plot(FREQ/1e6,np.angle(q[:,2]/cont,deg=True),label=label)
        axs[1,1].set(xlabel='频率 / MHz',ylabel='FDTD/连续解析相位 / 度',title='完整501点相位；线电流尺度保持一致');axs[1,1].legend(fontsize=8)
        fig.suptitle(f'平层AGL8m / {window} / 精确20–170MHz / 原生FP64\n同站配置列，非连续空间测线；PML固定2m，材料/物理界面不变')
        fig.savefig(a.out/f'v401_mesh_{window}_gray.png',dpi=140);plt.close(fig)
    print(json.dumps(dict(status=result['status'],metrics=metrics,frequency_comparison=result['frequency_comparison'],full_band_delta_change_vs_wide=full_error)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','reference','continuum','out','numerical','coarse']:p.add_argument('--'+name,type=Path,required=name!='coarse')
    main(p.parse_args())
