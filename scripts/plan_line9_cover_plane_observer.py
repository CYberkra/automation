"""Geometry-checked proposal only: native parallel planes inside the cover."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    m=json.loads((a.package/'manifest.json').read_text('utf-8'));geometry=a.package/'baseline/geometry.h5'
    assert sha(geometry)==m['baseline_sha256']['geometry.h5']
    planes=[27.025,28.525,30.025];xs=np.arange(154,174.001,.5);rows=[];receiver=m['receiver_count']+1
    with h5py.File(geometry) as h:
        g=h['data'][:,:,0]
        for yy in planes:
            j=round(yy/.025)
            for xx in xs:
                i=round(xx/.025);points=[]
                for k,(di,dj,fields) in enumerate([(0,0,['Ez','Hx','Hy']),(0,-1,['Hx']),(-1,0,['Hy'])]):
                    ii=i+di;jj=j+dj
                    assert 80<ii<8320 and 80<jj<1620 and (g[ii-1:ii+2,jj-1:jj+2]==1).all()
                    points.append({'coord':[ii,jj,0],'position_m':[ii*.025,jj*.025,0],
                                   'receiver_index':receiver,'name':f'flatcover_y{j:04d}_x{i:04d}_a{k}','outputs':fields})
                    receiver+=1
                rows.append({'id':f'flatcover_y{j:04d}_x{i:04d}','x_m':float(xx),'y_m':j*.025,'material_id':1,'anchors':points})
    original=m['receiver_count'];n=receiver-1;assert len(rows)==123 and n==original+369==1231
    out={'status':'GEOMETRY_VERIFIED_PROPOSAL_NOT_PREPARED_FROZEN_OR_RUN','script_sha256':sha(__file__),
         'parent_manifest_sha256':sha(a.package/'manifest.json'),'geometry_sha256':sha(geometry),
         'existing_native_observer_contract_sha256':'8b538b98c839ccc0b74af7a00dc89b5e372fec48ce1439a9418dd0fc11612820',
         'cover_plane_y_m':planes,'x_m':xs.tolist(),'new_probes':rows,'new_logical_points':123,'new_raw_receivers':369,
         'original_raw_receivers':original,'proposed_total_receivers':n,'proposed_total_logical_points':410,
         'receiver_device_history_bytes_lower_bound':6*20352*n*8,
         'new_saved_history_bytes_lower_bound':5*20352*123*8,
         'planned_source_and_main_rx_unchanged':True,'planned_existing_287_observers_retained':True,
         'planned_new_solves':1,'planned_max_attempts':1,'input_card_created':False,'execution_contract_frozen':False,'solver_started':False,
         'analysis_plan':'Full complex cover angular TE split on 1.5m and3m source-free uniform-cover slabs; no fitted amplitude/phase/time. Compare smooth air-radiating sectors and native time-window sensitivity. Require existing source/main/287 fields bitwise invariant, then trace upgoing cover packet without forcing a bounce count.',
         'required_before_run':['Live runtime/hash/GPU/host-capacity and owned-process preflight','Fresh immutable card and independent field/stencil guards','New single-attempt contract, stop/heartbeat ownership and raw completion audit'],
         'limits':'Geometry of observed support only; finite20m aperture and exterior heterogeneity persist. Budget is receiver-history lower bound, not total GPU need. No site/material/unique-ray certification.'}
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:out[k] for k in ['status','new_logical_points','proposed_total_receivers','receiver_device_history_bytes_lower_bound','new_saved_history_bytes_lower_bound']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
