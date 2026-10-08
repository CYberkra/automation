"""Independent sine-equation roots and direct inverse of the overlap prediction."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import newton


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(a):
    r=json.loads((a.public/'analysis.json').read_text('utf-8'))
    assert sha(a.continuum)==r['continuum_sha256'] and sha(a.wide)==r['wide_sha256'] and sha(a.fine)==r['fine_sha256']
    m=json.loads((a.prepared/'manifest.json').read_text('utf-8'));g=m['geometry_diagnostic']
    path=a.prepared/'flat_ricker_H1/geometries/line9_research_materials_v1_smoothed.json';assert sha(path)==r['material_sha256']
    db=json.loads(path.read_text('utf-8'))['materials']
    with h5py.File(a.continuum) as h:f=h['frequency_Hz'][:];parts=h['planar_response'][:];t=h['time_s'][:]
    with h5py.File(a.wide) as h:coarse=h['v401_response'][:,[1,2,3]]
    with h5py.File(a.fine) as h:fine=h['response'][:]
    airm=m['groups'][0]['tx_m'][1]+m['groups'][0]['rx_m'][1]-2*g['surface_y_m']
    cover=2*g['cover_base']['depth_m'];mud=2*(g['basal_sand']['depth_m']-g['cover_base']['depth_m'])
    horizontal=abs(m['groups'][0]['rx_m'][0]-m['groups'][0]['tx_m'][0]);vertical=m['groups'][0]['rx_m'][1]-m['groups'][0]['tx_m'][1]
    paths=np.array([[np.hypot(horizontal,vertical),0,0],[airm,0,0],[airm,cover,0],[airm,2*cover,0],[0,0,0],[airm,cover,mud]])
    metrics={};max_root_residual=0
    for name,spacing,actual in [('210m_2.5cm',.025,coarse),('80m_1.25cm',.0125,fine)]:
        timestep=spacing/(299792458*np.sqrt(2));differences=np.zeros((501,3),complex)
        for j,frequency in enumerate(f):
            for column,key in enumerate(['material_000_air','material_001_cover','material_002_mudstone']):
                v=db[key];base=v['base'];omega=2*np.pi*frequency
                eps=complex(base['relative_permittivity'],-base['electric_conductivity_s_per_m']/(omega*8.8541878128e-12))
                eps+=sum(p['relative_permittivity_difference']/(1+1j*omega*p['relaxation_time_s']) for p in v.get('poles',[]))
                n=np.sqrt(eps);continuous=omega*n/299792458
                rhs=n*spacing/(299792458*timestep)*np.sin(omega*timestep/2)
                root=newton(lambda z:np.sin(z*spacing/2)-rhs,continuous,
                            fprime=lambda z:spacing/2*np.cos(z*spacing/2),tol=1e-12,maxiter=30)
                max_root_residual=max(max_root_residual,float(abs(np.sin(root*spacing/2)-rhs)))
                differences[j,column]=root-continuous
        factors=np.exp(-1j*(differences@paths.T));pred=parts*factors
        mask=(t*1e9>=m['basal_gate_ns'][0])&(t*1e9<=m['basal_gate_ns'][1]);e=np.exp(2j*np.pi*t[mask,None]*f);rows={}
        for window,weights in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
            weights=weights/weights.mean();z=e@(weights[:,None]*np.column_stack([actual[:,0],pred[:,:5].sum(axis=1),actual[:,2],pred[:,5]]))/501
            rows[window]=dict(H0_prediction_error=float(np.linalg.norm(z[:,0]-z[:,1])/np.linalg.norm(z[:,0])),
                              bottom_prediction_error=float(np.linalg.norm(z[:,2]-z[:,3])/np.linalg.norm(z[:,2])))
            expected=r['metrics'][name][window]
            assert abs(rows[window]['H0_prediction_error']-expected['H0_propagation_prediction_unfitted_error'])<1e-10
            assert abs(rows[window]['bottom_prediction_error']-expected['bottom_prediction_unfitted_error'])<1e-10
        indices=[0,83,250,417,500];predicted=np.angle(factors[indices,5],deg=True)
        np.testing.assert_allclose(predicted,r['metrics'][name]['bottom_phase_prediction']['predicted_phase_deg'],atol=1e-9,rtol=0)
        rows['predicted_bottom_phase_deg']=predicted.tolist();metrics[name]=rows
    assert max_root_residual<1e-12
    out=dict(status='PASS_INDEPENDENT_ROOTS_AND_INVERSE',audit_script_sha256=sha(__file__),
             analysis_sha256=sha(a.public/'analysis.json'),max_sine_equation_residual=max_root_residual,metrics=metrics,
             limits='Audits the normal-propagation approximation; not exact discrete Debye/angle/interface theory or field evidence.')
    (a.public/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['public','prepared','continuum','wide','fine']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
