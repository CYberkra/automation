"""Unfitted component propagation predictions for planar H0; never filter data."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    m=json.loads((a.prepared/'manifest.json').read_text('utf-8'))
    db=a.prepared/'flat_ricker_H1/geometries/line9_research_materials_v1_smoothed.json'
    mat=json.loads(db.read_text('utf-8'))['materials']
    with h5py.File(a.continuum) as h:f=h['frequency_Hz'][:];parts=h['planar_response'][:];t=h['time_s'][:]
    with h5py.File(a.wide) as h:wide=h['v401_response'][:,[1,2,3]]
    np.testing.assert_array_equal(f,20e6+np.arange(501)*300000)
    media=[]
    for name in ['material_000_air','material_001_cover','material_002_mudstone','material_003_sandstone']:
        value=mat[name];base=value['base'];omega=2*np.pi*f
        er=base['relative_permittivity']+base['electric_conductivity_s_per_m']/(1j*omega*8.8541878128e-12)
        for pole in value.get('poles',[]):er=er+pole['relative_permittivity_difference']/(1+1j*omega*pole['relaxation_time_s'])
        media.append(np.sqrt(er*base['relative_permeability']))
    n=np.column_stack(media);kc=2*np.pi*f[:,None]*n/299792458
    g=m['geometry_diagnostic'];tx,rx=m['groups'][0]['tx_m'],m['groups'][0]['rx_m']
    air=tx[1]+rx[1]-2*g['surface_y_m'];cover=2*g['cover_base']['depth_m'];mud=2*(g['basal_sand']['depth_m']-g['cover_base']['depth_m'])
    # The unresolved sum of third/higher cover trips is deliberately left as
    # its original continuum term. No invented single path or fitted phase.
    paths=np.array([[np.hypot(rx[0]-tx[0],rx[1]-tx[1]),0,0,0],
                    [air,0,0,0],[air,cover,0,0],[air,2*cover,0,0],
                    [0,0,0,0],[air,cover,mud,0]])
    cases=[('210m_2.5cm',.025,wide)]
    if a.fine:
        with h5py.File(a.fine) as h:np.testing.assert_array_equal(h['frequency_Hz'][:],f);fine=h['response'][:]
        cases.append(('80m_1.25cm',.0125,fine))
    rows={};curves={}
    mask=(t*1e9>=m['basal_gate_ns'][0])&(t*1e9<=m['basal_gate_ns'][1])
    for name,dl,actual in cases:
        dt=dl/(299792458*np.sqrt(2));kg=2/dl*np.arcsin(n*dl/(299792458*dt)*np.sin(np.pi*f*dt)[:,None])
        factors=np.exp(-1j*((kg-kc)@paths.T));pred=parts*factors;rr={}
        for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
            w=w/w.mean();response=np.column_stack([actual[:,0],parts[:,:5].sum(axis=1),pred[:,:5].sum(axis=1),pred[:,3],pred[:,[0,1,2,4]].sum(axis=1),actual[:,2],pred[:,5]])
            # Direct full inverse; independent of the production CZT/IFFT path.
            z=np.exp(2j*np.pi*t[:,None]*f)@(response*w[:,None])/501;q=z[mask]
            norm=lambda x:np.linalg.norm(x)
            rr[window]=dict(H0_continuum_unfitted_error=float(norm(q[:,0]-q[:,1])/norm(q[:,0])),
                            H0_propagation_prediction_unfitted_error=float(norm(q[:,0]-q[:,2])/norm(q[:,0])),
                            cover_twice_vs_H0_correlation=float(abs(np.vdot(q[:,3],q[:,0]))/(norm(q[:,3])*norm(q[:,0]))),
                            predicted_other_terms_over_H0=float(norm(q[:,4])/norm(q[:,0])),
                            H0_over_bottom_delta=float(norm(q[:,0])/norm(q[:,5])),
                            H0_vs_bottom_inner_phase_deg=float(np.angle(np.vdot(q[:,5],q[:,0]),deg=True)),
                            bottom_prediction_unfitted_error=float(norm(q[:,5]-q[:,6])/norm(q[:,5])))
            curves[name,window]=z
        selected=[0,83,250,417,500]
        rr['bottom_phase_prediction']=dict(frequency_MHz=(f[selected]/1e6).tolist(),
            observed_phase_deg=np.angle((actual[:,2]/parts[:,5])[selected],deg=True).tolist(),
            predicted_phase_deg=np.angle(factors[selected,5],deg=True).tolist(),
            full_band_phase_residual_max_deg=float(abs(np.angle((actual[:,2]/parts[:,5])/factors[:,5],deg=True)).max()))
        rows[name]=rr
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(len(cases),2,figsize=(13,4.5*len(cases)),squeeze=False,layout='constrained')
    show=(t*1e9>=340)&(t*1e9<=405)
    for row,(name,dl,_) in enumerate(cases):
        for col,window in enumerate(['hann','blackman']):
            z=curves[name,window];ax=axes[row,col]
            for j,label,style in [(0,'实际H0','-'),(1,'连续解析完整H0',':'),(2,'含网格传播的完整H0预测','--'),(3,'覆盖层第二次往返预测','-'),(5,'实际底砂差场','-')]:ax.plot(t[show]*1e9,abs(z[show,j]),style,lw=1,label=label)
            ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.07);ax.set(title=name+' / '+window,xlabel='SFCW时间 / ns',ylabel='包络 / 共同参考 (V/m)/(A·m)');ax.legend(fontsize=8)
    fig.suptitle('平层强剩余回波：解析路径与网格传播预测\n事后机制诊断，零拟合参数；仅改变解析预测，原始FDTD数据不改')
    fig.savefig(a.out/'v401_planar_overlap_prediction.png',dpi=140);plt.close(fig)
    result=dict(status='UNFITTED_PLANAR_OVERLAP_DIAGNOSTIC_NOT_NONFLAT_PATH_CERTIFICATION',script_sha256=sha(__file__),
                material_sha256=sha(db),prepared_manifest_sha256=sha(a.prepared/'manifest.json'),
                continuum_sha256=sha(a.continuum),wide_sha256=sha(a.wide),fine_sha256=sha(a.fine) if a.fine else None,
                paths_m=paths.tolist(),metrics=rows,fit=False,production_correction=False,
                limits='Continuum Debye spectra in normal bulk Yee relation; not exact discrete interfaces/angular spectrum. Third/higher cover sum remains unmodified. Norm ratios are not energy fractions; flat model only, not a unique path attribution of nonflat H0 or an instruction to delete real reflections.')
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print(json.dumps(rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prepared','continuum','wide','out','fine']:p.add_argument('--'+key,type=Path,required=key!='fine')
    main(p.parse_args())
