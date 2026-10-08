"""Recompute dielectric spectra, nonlinear roots and direct inverse for grid prediction."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import newton
from hs_capsule_identity import sha256 as sha


def main(a):
    r=json.loads((a.public/'analysis.json').read_text('utf-8'));m=json.loads((a.prepared/'manifest.json').read_text('utf-8'))
    assert r['script_sha256']==sha(Path(__file__).parent/'diagnose_line9_v5_grid_phase.py')
    assert r['source_numerical_sha256']==sha(a.source) and r['numerical_sha256']==sha(a.numerical) and r['prepared_manifest_sha256']==sha(a.prepared/'manifest.json')
    db=a.prepared/'flat_ricker_H1/geometries/line9_research_materials_v1_smoothed.json';assert sha(db)==r['material_sha256']
    materials=json.loads(db.read_text('utf-8'))['materials']
    with h5py.File(a.source) as h:f=h['frequency_Hz'][:];actual=h['response'][:,6];planar=h['planar_response'][:,5];t=h['time_s'][:]
    with h5py.File(a.numerical) as h:factor=h['bulk_factor'][:];numerical=h['numerical_k'][:];continuous=h['continuum_k'][:]
    dt=m['dt_s'];dl=.025;eps0=8.8541878128e-12;c0=299792458.;independent=[];root_errors=[]
    for j in range(501):
        ks=[]
        for k,name in enumerate(['material_000_air','material_001_cover','material_002_mudstone']):
            mat=materials[name];base=mat['base'];omega=2*np.pi*f[j]
            er=base['relative_permittivity']+base['electric_conductivity_s_per_m']/(1j*omega*eps0)
            er+=sum(p['relative_permittivity_difference']/(1+1j*omega*p['relaxation_time_s']) for p in mat.get('poles',[]))
            kc=omega/c0*np.sqrt(er);assert abs(kc-continuous[j,k])<1e-12
            value=newton(lambda z:np.sin(z*dl/2)-np.sqrt(er)*dl/(c0*dt)*np.sin(omega*dt/2),kc,fprime=lambda z:dl/2*np.cos(z*dl/2),tol=1e-12,maxiter=30)
            root_errors.append(abs(value-numerical[j,k])/abs(value));ks.append(value-kc)
        # All distances are independently read from pre-solve geometry and positions.
        g=m['geometry_diagnostic'];air=(m['groups'][0]['tx_m'][1]+m['groups'][0]['rx_m'][1])-2*g['surface_y_m'];cover=2*g['cover_base']['depth_m'];mud=2*(g['basal_sand']['depth_m']-g['cover_base']['depth_m'])
        independent.append(np.exp(-1j*(air*ks[0]+cover*ks[1]+mud*ks[2])))
    independent=np.array(independent);error=float(np.linalg.norm(independent-factor)/np.linalg.norm(factor));assert error<1e-11
    residual=float(abs(np.angle((actual/planar)/independent,deg=True)).max());assert abs(residual-r['phase_residual_max_deg'])<1e-9
    metrics={}
    for window in ['hann','blackman']:
        w=np.hanning(501) if window=='hann' else np.blackman(501);w/=w.mean();lo,hi=m['basal_gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi)
        z=np.exp(2j*np.pi*t[keep,None]*f)@(w[:,None]*np.column_stack([actual,planar,planar*independent]))/501
        value=float(np.linalg.norm(z[:,0]-z[:,2])/np.linalg.norm(z[:,0]));assert abs(value-r['metrics'][window]['predicted_unfitted_relative_L2'])<1e-10;metrics[window]=value
    result=dict(status='PASS_INDEPENDENT_CONSTITUTIVE_SPECTRA_ROOTS_AND_DIRECT_INVERSE_NOT_FDTD_REFINEMENT',audit_script_sha256=sha(__file__),independent_root_relative_error_max=float(max(root_errors)),independent_bulk_factor_relative_L2=error,phase_residual_max_deg=residual,direct_inverse_predicted_unfitted_relative_L2=metrics,plot_sha256=sha(a.public/'v5_grid_phase_prediction.png'),limits='Only unfitted normal-propagation approximation audited. No full discrete-Debye/interface solution, physical calibration, production correction or refined-grid solver proof.')
    (a.public/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['public','source','numerical','prepared']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
