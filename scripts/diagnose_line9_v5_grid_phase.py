"""Unfitted propagation-only Yee-dispersion prediction for the planar v5 basal pair."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import newton
from review_line9_result_packages import indices,C0
from audit_line9_postprocessing import inverse,weights
from hs_capsule_identity import sha256 as sha
from analyze_line9_basal_pair import plotting,save


def main(a):
    assert not a.out.exists() and not a.numerical.exists(),'Fresh outputs required'
    c=json.loads((a.prepared/'manifest.json').read_text('utf-8'));dt=c['dt_s'];dl=.025
    db=a.prepared/'flat_ricker_H1/geometries/line9_research_materials_v1_smoothed.json'
    assert sha(db)==c['groups'][2]['material_sha256'];materials=json.loads(db.read_text('utf-8'))['materials']
    with h5py.File(a.source) as h:f=h['frequency_Hz'][:];actual=h['response'][:,6];planar=h['planar_response'][:,5]
    np.testing.assert_array_equal(f,20e6+np.arange(501)*300000)
    n=np.array([indices(materials,x) for x in f]);k=2*np.pi*f[:,None]*n/C0
    argument=n*dl/(C0*dt)*np.sin(np.pi*f*dt)[:,None]
    numerical_k=2/dl*np.arcsin(argument)
    # Independent nonlinear root solve of the original sine relation, no arcsin.
    check=[]
    for j in [0,83,250,417,500]:
        for medium in [0,1,2]:
            value=newton(lambda z:np.sin(z*dl/2)-argument[j,medium],k[j,medium],fprime=lambda z:dl/2*np.cos(z*dl/2),tol=1e-12,maxiter=30)
            check.append(abs(value-numerical_k[j,medium])/abs(value))
    assert max(check)<1e-12
    # Lossless 1D magic time step is exact; refinement must reduce bulk error.
    ideal=2/dl*np.arcsin(np.sin(np.pi*f*dl/C0));assert np.max(abs(ideal-k[:,0])/abs(k[:,0]))<1e-12
    fine=4/dl*np.arcsin(n*dl/(C0*dt)*np.sin(np.pi*f*dt/2)[:,None])
    assert np.all(abs(fine-k)<abs(numerical_k-k))
    g=c['geometry_diagnostic'];hs=(c['groups'][0]['tx_m'][1]+c['groups'][0]['rx_m'][1])/2-g['surface_y_m']
    cover=g['cover_base']['depth_m'];mud=g['basal_sand']['depth_m']-cover
    path=np.array([2*hs,2*cover,2*mud,0.])
    delta=-np.sum((numerical_k-k)*path,axis=1);factor=np.exp(1j*delta)
    fine_delta=-np.sum((fine-k)*path,axis=1)
    observed=actual/planar;phase_residual=np.angle(observed/factor,deg=True);predicted=planar*factor
    metrics={};plt=plotting();fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    axes[0,0].plot(f/1e6,np.angle(observed,deg=True),label='平层FDTD/连续解析的实际相位');axes[0,0].plot(f/1e6,np.degrees(delta.real),'--',label='仅由网格、dt、材料和路程预测')
    axes[0,0].plot(f/1e6,np.degrees(fine_delta.real),':',label='网格减半的解析预测，尚未正演')
    axes[0,0].set(xlabel='频率 / MHz',ylabel='相位差 / °',title='不拟合相位、幅度或时延');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(f/1e6,abs(observed),label='实际幅度比');axes[0,1].plot(f/1e6,abs(factor),'--',label='复数传播常数预测幅度比')
    axes[0,1].set(xlabel='频率 / MHz',ylabel='幅度比',title='只预测体传播，不含完整离散界面/角谱');axes[0,1].legend(fontsize=8)
    for ax,window in zip(axes[1],['hann','blackman']):
        z,t=inverse(np.column_stack([actual,planar,predicted]),f,weights(window,501));mask=(t*1e9>=c['basal_gate_ns'][0])&(t*1e9<=c['basal_gate_ns'][1]);show=(t*1e9>=340)&(t*1e9<=405)
        for j,label in enumerate(['平层FDTD底砂差场','连续平层二维解析','加入预计算网格传播因子的解析预测']):ax.plot(t[show]*1e9,abs(z[show,j]),lw=1,label=label)
        metrics[window]=dict(original_unfitted_relative_L2=float(np.linalg.norm(z[mask,0]-z[mask,1])/np.linalg.norm(z[mask,0])),predicted_unfitted_relative_L2=float(np.linalg.norm(z[mask,0]-z[mask,2])/np.linalg.norm(z[mask,0])))
        ax.axvspan(*c['basal_gate_ns'],color='red',alpha=.08);ax.set(xlabel='SFCW时间 / ns',ylabel='包络 / (V/m)/(A·m)',title=window+'：只改解析预测，原正演数据不改');ax.legend(fontsize=8)
    fig.suptitle('v5平层底砂：Yee网格传播误差的无拟合诊断\n连续Debye谱代入均匀介质关系，仅法向体传播近似；不是完整离散层栈或网格收敛证明')
    a.out.mkdir(parents=True);fig.savefig(a.out/'v5_grid_phase_prediction.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        for key,value in [('frequency_Hz',f),('continuum_k',k),('numerical_k',numerical_k),('bulk_factor',factor),('observed_ratio',observed),('predicted_response',predicted)]:h.create_dataset(key,data=value)
    j=np.array([0,83,250,417,500])
    save(a.out/'analysis.json',dict(status='UNFITTED_BULK_YEE_PROPAGATION_DIAGNOSTIC_NOT_MESH_CONVERGENCE',script_sha256=sha(__file__),source_numerical_sha256=sha(a.source),prepared_manifest_sha256=sha(a.prepared/'manifest.json'),material_sha256=sha(db),numerical_sha256=sha(a.numerical),independent_complex_root_relative_error_max=float(max(check)),magic_timestep_check=True,half_mesh_bulk_error_decreases=True,phase_residual_max_deg=float(abs(phase_residual).max()),complex_ratio_residual_relative_L2=float(np.linalg.norm(observed-factor)/np.linalg.norm(observed)),metrics=metrics,
        selected=dict(frequency_MHz=(f[j]/1e6).tolist(),observed_phase_deg=np.angle(observed[j],deg=True).tolist(),predicted_phase_deg=np.degrees(delta[j].real).tolist(),observed_amplitude_ratio=abs(observed[j]).tolist(),predicted_amplitude_ratio=abs(factor[j]).tolist(),half_mesh_predicted_phase_deg=np.degrees(fine_delta[j].real).tolist()),
        method='k_grid=2/dl*asin(n(omega)*dl/(c*dt)*sin(omega*dt/2)); factor=exp(-i*sum[(k_grid-k_cont)*two-way-layer-length]). Normal-incidence propagation approximation, zero fitted coefficients. Input n is continuum Debye+DC conductivity, not the complete discrete constitutive model.',
        source_reading='Schneider Understanding FDTD, sections7.3-7.4, printedpp163-166, Eqs7.19-7.40; lossless homogeneous derivation. Complex continuum-index substitution and layered phase accumulation are our diagnostic extension, not claimed exact discrete Debye/slab theory.',sources=['https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf','https://docs.gprmax.com/en/latest/gprmodelling.html'],
        limits='Posthoc mechanism diagnostic with prior geometric/material parameters; no data-fit correction or production output change. Planar finite-domain2D pair only, no nonflat phase correction, mesh-refined FDTD or field/3D certification. Small remaining error does not certify an FDTD error floor.'))
    print(json.dumps(dict(phase_residual_max_deg=float(abs(phase_residual).max()),metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','prepared','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
