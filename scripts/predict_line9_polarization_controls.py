"""Unfitted local planar 2D prediction for the constitutive factorial, no FDTD."""
import argparse
import copy
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.special import hankel2
from diagnose_line9_v5_planar_green import integral,reflect_parts,checks,MU0
from review_line9_result_packages import indices,C0
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    original=a.parent/'baseline_H1/geometries/line9_research_materials_v1_smoothed.json'
    baseline=json.loads(original.read_text('utf-8'));f=20e6+np.arange(501)*300000.
    names=['span3_dc003','span3_dc0003','span03_dc003','span03_dc0003'];groups={g['id']:g for g in m['groups']}
    diagnostics=checks();values=[];rows=[]
    geo=m['geometry_diagnostic'];tx,rx=m['groups'][0]['tx_m'],m['groups'][0]['rx_m']
    height=tx[1]+rx[1]-2*geo['surface_y_m'];dx=abs(tx[0]-rx[0])
    d1=geo['cover_base']['depth_m'];d2=geo['basal_sand']['depth_m']-d1
    assert abs(d1-7.3)<1e-12 and abs(d2-7)<1e-12
    for name in names:
        db=copy.deepcopy(baseline) if name.startswith('span3_') else json.loads((a.package/groups[name+'_H1']['material']).read_text('utf-8'))
        if name=='span3_dc0003':db['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m']=.0003
        n=np.array([indices(db['materials'],q) for q in f]);k=2*np.pi*f[:,None]/C0*n
        kernel=lambda q,ky:reflect_parts(q,ky,k[:,1],k[:,2],k[:,3],d1,d2)
        prefactor=-f*MU0/.025
        coarse=integral(f,height,dx,kernel,order=256)*prefactor[:,None]
        fine=integral(f,height,dx,kernel,order=512)*prefactor[:,None]
        error=np.linalg.norm(fine-coarse,axis=0)/np.linalg.norm(fine,axis=0)
        assert max(error)<1e-6
        direct=-2*np.pi*f*MU0/(4*.025)*hankel2(0,k[:,0]*np.hypot(dx,tx[1]-rx[1]))
        h0=direct+fine[:,:4].sum(axis=1);delta=fine[:,4]
        values.extend([h0,h0+delta,delta])
        j=250;nm=n[j,2];nc=n[j,1];ns=n[j,3]
        rt=(1-((nc-nm)/(nc+nm))**2)*((nm-ns)/(nm+ns))
        rows.append(dict(id=name,quadrature_component_relative_L2=error.tolist(),
            normal_incidence95=dict(absorption_7m_roundtrip_dB=float(20*np.log10(abs(np.exp(-2j*k[j,2]*d2)))),
                cover_mud_transmission_times_basal_reflection_amplitude=float(abs(rt)),
                cover_mud_transmission_times_basal_reflection_phase_deg=float(np.angle(rt,deg=True))),
            material_content=db))
    z=np.column_stack(values);t=np.arange(4008)/(4008*300000.);metrics={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();metrics[window]={}
        for gate in ['basal','deep']:
            lo,hi=m[gate+'_gate_ns'];keep=(t*1e9>=lo)&(t*1e9<=hi)
            x=np.exp(2j*np.pi*t[keep,None]*f)@(z*w[:,None])/501
            norms=[np.linalg.norm(x[:,3*j+2]) for j in range(4)]
            metrics[window][gate]={name:dict(delta_over_original=float(norms[j]/norms[0]),
                H0_over_delta=float(np.linalg.norm(x[:,3*j])/norms[j]),
                H1_vs_delta_correlation=float(abs(np.vdot(x[:,3*j+1],x[:,3*j+2]))/(np.linalg.norm(x[:,3*j+1])*norms[j]))) for j,name in enumerate(names)}
    a.out.mkdir(parents=True)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=f);h.create_dataset('response',data=z)
    result=dict(status='UNFITTED_CONTINUUM_PLANAR_PREDICTION_NOT_NONFLAT_PATH_CERTIFICATION',script_sha256=sha(__file__),
        planar_kernel_sha256=sha(Path(__file__).with_name('diagnose_line9_v5_planar_green.py')),
        package_manifest_sha256=sha(a.package/'manifest.json'),parent_material_sha256=sha(original),
        numerical_sha256=sha(a.numerical),checks=diagnostics,source_height_sum_m=height,horizontal_separation_m=dx,
        cover_thickness_m=d1,mud_thickness_m=d2,metrics=metrics,cases=rows,fit=False,
        limits='Infinite flat local column, continuous medium, exact scalar2D angular spectrum. No phase/grid/data fit. Nonflat lateral structure and discrete FDTD absent; not a correction to original data or site material validation.')
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(metrics=metrics,normal_incidence95=[r['normal_incidence95'] for r in rows])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','parent','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
