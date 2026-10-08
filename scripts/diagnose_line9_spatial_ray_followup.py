"""Post-hoc CPU ray hypothesis for the 146.4m residual; never modifies data."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from diagnose_line9_refracted_paths import extract_columns,ray_on_geometry,corrected_primary,checks
from review_line9_result_packages import indices,FREQ,C0
from analyze_line9_v401_version_controls import inverse
from gprMax.toolboxes.SFCW.processing import spectral_window
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(a.package/'manifest.json');pub=read(a.public/'analysis.json');aud=read(a.public/'independent_audit.json')
    assert aud['analysis_sha256']==sha(a.public/'analysis.json')
    group=next(g for g in m['groups'] if g['id']=='full2d_r0359_H1');s=m['station_design']['stations'][2];geometry=s['geometry']
    assert s['chainage_m']==146.4
    for key in ['geometry','material']:assert sha(a.package/group[key])==group[key+'_sha256']
    db=read(a.package/group['material'])['materials'];n95=indices(db,95e6)
    with h5py.File(a.package/group['geometry']) as h:data=h['data'][:,:,0]
    columns=extract_columns(data,.025);signature=tuple((b['above'],b['below']) for b in geometry['boundaries'])
    # Separate constitutive formula and transmission product; no production
    # material, SFCW or inverse functions in this independent numeric check.
    f=20e6+np.arange(501)*300000.;omega=2*np.pi*f;n=[]
    for mat in db.values():
        b=mat['base'];er=b['relative_permittivity']-1j*b['electric_conductivity_s_per_m']/(omega*8.8541878128e-12)
        for pole in mat.get('poles',[]):er+=pole['relative_permittivity_difference']/(1+1j*omega*pole['relaxation_time_s'])
        n.append(np.sqrt(er))
    n=np.column_stack(n);trans=lambda i,j:2*n[:,i]/(n[:,i]+n[:,j]);coef=trans(0,1)*trans(1,0)*trans(1,2)*trans(2,1)*(n[:,2]-n[:,3])/(n[:,2]+n[:,3])
    t=np.arange(4008)/(4008*300000.);matrix=np.exp(2j*np.pi*t[:,None]*f)/501;rows=[]
    for stride in [8,4]:
        path=ray_on_geometry(np.array(group['tx_m']),np.array(group['rx_m']),columns,signature,.025,n95,stride,data)
        assert path['geometry_check']=='NO_OUTSIDE_BAND_MISMATCH'
        local,ray=corrected_primary(geometry,db,path)
        points=np.array(path['nodes_xy_m']);media=np.array(path['segment_materials']);lengths=np.linalg.norm(np.diff(points,axis=0),axis=1)
        independent_lengths=np.bincount(media,weights=lengths,minlength=4)
        np.testing.assert_allclose(independent_lengths,path['path_length_by_material_m'],atol=1e-13,rtol=0)
        phase=float(np.dot(independent_lengths,n[250].real)/C0*1e9);assert abs(phase-path['phase_time95_ns'])<1e-10
        expected=coef*np.exp(-1j*omega*(n@independent_lengths)/C0)
        spectrum_error=float(np.linalg.norm(expected-ray)/np.linalg.norm(expected));assert spectrum_error<1e-10
        windows={}
        for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
            w=w/w.mean();profiles,_=inverse(np.column_stack([local,ray]),spectral_window(window,501))
            direct=matrix@(expected*w);ray_peak=float(t[np.argmax(abs(profiles[:,1]))]*1e9);direct_peak=float(t[np.argmax(abs(direct))]*1e9)
            assert abs(ray_peak-direct_peak)<1e-9
            actual=next(r for r in pub['metrics'][window] if r['chainage_m']==146.4)['gates']['basal']['delta_peak_ns']
            windows[window]=dict(local_peak_ns=float(t[np.argmax(abs(profiles[:,0]))]*1e9),ray_peak_ns=ray_peak,native_delta_peak_ns=actual,
                actual_minus_local_ns=float(actual-t[np.argmax(abs(profiles[:,0]))]*1e9),actual_minus_ray_ns=float(actual-ray_peak),independent_direct_inverse_peak_error_ns=abs(ray_peak-direct_peak))
        rows.append(dict(stride=stride,path=path,independent_spectrum_relative_L2=spectrum_error,windows=windows))
    result=dict(status='POSTHOC_EXPLORATORY_1464_OFF_NADIR_PATH_NOT_UNIQUE_PHYSICAL_ATTRIBUTION',script_sha256=sha(__file__),
        ray_helper_sha256=sha(Path(__file__).with_name('diagnose_line9_refracted_paths.py')),analysis_sha256=sha(a.public/'analysis.json'),
        package_manifest_sha256=sha(a.package/'manifest.json'),selected_after_observing_residual=True,new_solves=0,fit_to_field_data=False,
        flat_and_inclined_ray_checks=checks(),results=rows,
        limits='Uses existing five-start phase95 Fermat optimizer and piecewise-linear sampled interfaces. Not independent optimizer/global minimum proof. Multi-start spread, unresolved two-cell contact band, oblique coefficients, spreading, diffraction and full discrete interfaces remain. Only mathematical spectra/length/inverse consistency audited; does not independently certify actual peak path or field model. Original FDTD and predeclared metrics unchanged.')
    assert sha(a.numerical)==pub['numerical_sha256']
    with h5py.File(a.numerical) as h:z=h['response'][:]
    profiles,_=inverse(z,spectral_window('hann',501));keep=(t*1e9>=230)&(t*1e9<=290)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap,BoundaryNorm
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,2,figsize=(15,6),layout='constrained')
    lo,hi=137.,157.;ix0=round((lo-20)/.025);ix1=round((hi-20)/.025)
    cmap=ListedColormap(['#e8f4fa','#e0b992','#9c8c87','#e6cf70'])
    axes[0].imshow(data[ix0:ix1:4,::4].T,origin='lower',extent=[lo,hi,0,42.5],aspect='auto',interpolation='nearest',cmap=cmap,norm=BoundaryNorm(np.arange(5)-.5,4))
    for row,color in zip(rows,['#246bb2','#8e44ad']):
        points=np.array(row['path']['nodes_xy_m']);axes[0].plot(points[:,0]+20,points[:,1],'.-',c=color,lw=1,label=f'界面采样{row["stride"]*.025:g}m候选射线')
    axes[0].axvline(146.4,color='green',ls=':',label='天线中点正下方')
    axes[0].set(xlim=(lo,hi),ylim=(16,40),title='几何决定的候选离轴路径，非波场快照',xlabel='剖面里程 / m',ylabel='模型y / m');axes[0].legend(fontsize=8)
    axes[1].plot(t[keep]*1e9,abs(profiles[keep,10]),label='原生H1总场',lw=1.5)
    axes[1].plot(t[keep]*1e9,abs(profiles[keep,11]),label='原生H1−H0差场',ls='--',lw=1.2)
    q=rows[0]['windows']['hann']
    for value,color,label in [(q['local_peak_ns'],'green','局部正下方预测'),(q['ray_peak_ns'],'#246bb2','离轴射线预测'),(q['native_delta_peak_ns'],'#c0392b','实际差场峰')]:axes[1].axvline(value,color=color,ls=':',label=f'{label} {value:.2f}ns')
    axes[1].set(title='Hann绝对包络：不做延时或幅相拟合',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');axes[1].legend(fontsize=8)
    fig.suptitle('146.4m事后路径探索：侧偏约3.1–3.3m，候选峰259.48ns，实际260.31ns\n仍有多初值局部极小/界面两格未分辨带；接近不等于唯一物理归因，原数据与事前指标不改')
    figure=a.public/'posthoc_off_nadir_geometry_and_arrival.png';assert not figure.exists();fig.savefig(figure,dpi=140);plt.close(fig)
    result['figure_sha256']=sha(figure);result['numerical_sha256']=sha(a.numerical)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','public','numerical','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
