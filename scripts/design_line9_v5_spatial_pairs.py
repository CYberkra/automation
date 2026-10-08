"""Geometry-only next-station design, not a solver execution contract."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from review_line9_result_packages import geometry_at,indices
from diagnose_line9_layer_kinematics import primary_paths
from audit_line9_postprocessing import inverse,weights,FREQ


def main(a):
    assert not a.out.exists()
    m=json.loads((a.parent/'manifest.json').read_text('utf-8'));groups={g['id']:g for g in m['groups']}
    base=groups['span03_dc0003_H1'];geometry=a.parent/base['geometry'];material=a.parent/base['material']
    for p,key in [(geometry,'geometry'),(material,'material')]:assert sha(p)==base[key+'_sha256']
    db=json.loads(material.read_text('utf-8'))['materials'];n95=indices(db,95e6)
    plan=json.loads(a.stations.read_text('utf-8'));assert plan['x_offset_model_eq_chainage_minus']==20.65
    with h5py.File(geometry) as h:data=h['data'][:,:,0];dl=h.attrs['dx_dy_dz']
    rows=[]
    for position in [190.,180.,146.4]:
        candidates=[g for g in plan['cases'] if abs(g['chainage_m']-position)<1e-10];assert len(candidates)==1
        case=candidates[0];tx=np.array(case['tx_m']);rx=np.array(case['rx_m'])
        assert abs(tx[0]-(position-20.65))<1e-10
        assert abs((tx[0]+rx[0])/2-(position-20))<1e-10
        assert abs(rx[0]-tx[0]-1.3)<1e-10
        g=geometry_at(data,dl,tx,rx,n95);assert abs(g['midpoint_agl_m']-8)<.025
        for p in [tx,rx]:assert data[round(p[0]/.025),round(p[1]/.025)]==0
        template=primary_paths(g,db)[g['boundaries'].index(g['basal_sand'])]
        x,t=inverse(template,FREQ,weights('hann',501));center=float(t[np.argmax(abs(x))]*1e9)
        side=float(min(tx[0]-2,rx[0]-2,208-tx[0],208-rx[0]))
        top=float(40.5-max(tx[1],rx[1]));assert side>29 and top>1.075-1e-12
        rows.append(dict(id=case['id'],chainage_m=position,tx_m=tx.tolist(),rx_m=rx.tolist(),
            cover_thickness_m=g['cover_base']['depth_m'],mud_thickness_m=g['basal_sand']['depth_m']-g['cover_base']['depth_m'],
            basal_depth_m=g['basal_sand']['depth_m'],geometry=g,basal_template_peak_ns=center,
            proposed_basal_gate_ns=[center-12,center+12],proposed_wide_gate_ns=[center-60,center+60],
            min_side_PML_clearance_m=side,min_top_PML_clearance_m=top))
    result=dict(status='GEOMETRY_ONLY_NEXT_STATION_DESIGN_NOT_EXECUTION_CONTRACT',script_sha256=sha(__file__),
        geometry_sha256=sha(geometry),material_sha256=sha(material),planned_station_manifest_sha256=sha(a.stations),
        parent_manifest_sha256=sha(a.parent/'manifest.json'),stations=rows,new_solves_started=0,proposed_new_pairs=3,
        proposed_reused_station_chainage_m=198.6,proposed_new_native=6,solver_version='4.0.1',precision='native_float64',
        source='40A/100MHz Ricker',domain_m=[210,42.5,.025],spacing_m=[.025]*3,time_window_ns=1200,
        source_antenna_model='2D z-electric line-source field proxy, not actual 3D antenna or port S21',
        material_status='Reuse span0.3/DC0.0003 counterfactual, no new fitting or production replacement',
        processing='Exact501tones20-170MHz/.3MHz,actual-source/Yee offsets,.025m reference,Hann+Blackman,complex H1-H0,noAGC/fit',
        prerequisite='Review current low-loss boundary outputs; generate inputs, independently audit and freeze new bounded execution before solving. This document cannot launch a solver.',
        limits='Three geometry-chosen sites plus existing reference. Original planned points are not completed simulations. Local template windows are model-informed, not an observed-event truth or separability criterion. New sites have farther PML clearance, which is not a convergence proof.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],stations=[{k:r[k] for k in ['chainage_m','cover_thickness_m','mud_thickness_m','basal_depth_m','basal_template_peak_ns','min_side_PML_clearance_m','min_top_PML_clearance_m']} for r in rows])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['parent','stations','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
