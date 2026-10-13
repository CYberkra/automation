"""Prepare the approved102 matched-coordinate runs and analyze fixed operators."""
from __future__ import annotations
import argparse
import copy
import math
from pathlib import Path
import shutil
import sys
import h5py
import numpy as np

from study_line9_v31_smoothing_domain import DB, load, save, sha, write_geometry, vertical_events


def raster_at(contacts,selection,origin,h):
    n=round(50/h);xs=origin+(np.arange(n)+.5)*h;ys=410+(np.arange(n)+.5)*h
    lines={k:np.interp(xs,np.array(v)[:,0],np.array(v)[:,1])[:,None] for k,v in contacts.items()}
    surface,clay,upper,lt,lb=[lines[k] for k in ['surface','clay_base','upper_mudstone_base','lower_mudstone_top','lower_mudstone_base']]
    lo,hi=selection['upper_mudstone_extent_m'];u=((xs>=lo)&(xs<hi))[:,None];y=ys[None,:]
    data=np.where(y>=surface,0,np.where(y>=clay,1,np.where(u&(y>=upper)&(y<clay),2,
        np.where((y>=lb)&(y<np.minimum(lt,clay)),2,3)))).astype(np.int16)
    independent=np.full((n,n),3,dtype=np.int16)
    independent[(y>=lb)&(y<np.minimum(lt,clay))]=2
    independent[u&(y>=upper)&(y<clay)]=2
    independent[(y>=clay)&(y<surface)]=1;independent[y>=surface]=0
    assert np.array_equal(data,independent)
    return data[:,:,None]


def profile_text(case):
    h=case['grid_m'];tx,rx=case['tx_m'],case['rx_m'];pc=case['pml_cells']
    return '\n'.join([f"#title: {case['id']}; matched-coordinate mesh comparison",
        '#domain_mode: TM','#domain: 50 50 inf',f'#dx_dy_dz: {h:.17g} {h:.17g} {h:.17g}',
        '#time_window: 1.1e-6','#time_step_stability_factor: 0.95',
        '#pml_cells: '+' '.join(str(v) for v in pc),'#pml_formulation: HORIPML',
        f"#geometry_objects_read: 0 0 0 ../{case['geometry']} {DB} y",
        '#waveform: impulse 1 1 line9_impulse',
        f'#hertzian_dipole: z {tx[0]:.17g} {tx[1]:.17g} inf line9_impulse',
        f'#rx: {rx[0]:.17g} {rx[1]:.17g} inf line9_rx Ez',''])


def prepare(a):
    root=a.root.resolve();root.mkdir(parents=True,exist_ok=False);geo=root/'geometries';geo.mkdir();(root/'profiles').mkdir()
    audit=load(a.compatibility)['stations'];assert len(audit)==51
    cs=load(a.package/'model/contacts_m.json');selpath=a.package/'config/selection_pdf_reconstructed.json';sel=load(selpath)
    assert cs['selection_sha256']==sha(selpath) and cs['geometry_model_version']=='V3.1_VERTICAL_PDF_SCALE_FIXED'
    contacts=cs['contact_points_x_elevation_m'];prior=load(a.seed/'contract.json')
    db=a.package/'model'/f'{DB}.json';assert sha(db)==prior['material_database_sha256'];shutil.copy2(db,geo/db.name)
    cases=[];proofs=[]
    for row in audit:
        s=row['station_m'];origin=s-25.05;ym=row['proposed_common_y_m'];pair={}
        oldpath=a.old_geometries/f'line9_x{s:08.3f}.h5'
        with h5py.File(oldpath) as old:olddata=old['data'][:];attrs=dict(old.attrs)
        for label,h in [('fine',.025),('coarse',1/30)]:
            data=raster_at(contacts,sel,origin,h)
            if label=='fine':assert np.array_equal(data[2:],olddata[:-2])
            cid=f'{label}_x{s:.1f}';path=geo/f'{cid}.h5';write_geometry(path,data,attrs,h,origin)
            tx=[24.4,float(ym[0]),0.];rx=[25.7,float(ym[1]),0.]
            for pos in [tx,rx]:
                assert np.allclose(np.array(pos[:2])/h,np.round(np.array(pos[:2])/h),atol=1e-8,rtol=0)
                assert data[round(pos[0]/h),round(pos[1]/h),0]==0
                assert 50-3-pos[1]>=2 and min(pos[0]-3,50-3-pos[0])>=10
            surface_y=np.max(np.flatnonzero(data[data.shape[0]//2,:,0]!=0))*h+h
            patch=[[math.floor(23/h)*h,math.floor((surface_y-.3)/h)*h,0.],
                   [math.ceil(27.1/h)*h,math.ceil((surface_y+.3)/h)*h,h]]
            # Include actual surface variation over the entire export interval.
            xidx=np.arange(round(patch[0][0]/h),round(patch[1][0]/h))
            heights=np.array([(np.flatnonzero(data[ix,:,0]!=0)[-1]+1)*h for ix in xidx])
            patch[0][1]=math.floor((heights.min()-.15)/h)*h;patch[1][1]=math.ceil((heights.max()+.15)/h)*h
            dt=.95*h/(299792458*np.sqrt(2));pml=round(3/h)
            case=dict(id=cid,stage='A',family='matched_continuous_bscan',station_m=s,grid_m=h,
                geometry=path.relative_to(root).as_posix(),geometry_sha256=sha(path),averaging='y',domain_m=[50.,50.],
                tx_m=tx,rx_m=rx,tx_global_x_elevation_m=[s-.65,410+tx[1]],rx_global_x_elevation_m=[s+.65,410+rx[1]],
                origin_global_x_elevation_m=[origin,410.],pml_cells=[pml,pml,0,pml,pml,0],
                dt_s=float(dt),expected_samples=math.ceil(1100e-9/dt)+1,geometry_only=False,
                surface_coefficient_export_bounds=patch,old_fine_source_y_m=row['old_fine_y_m'],
                source_y_adjustment_m=row['proposed_y_change_m'])
            pair[label]=case
            if label=='fine':events=vertical_events(data,tx,rx,h)
        for label in ['coarse','fine']:
            case=pair[label];case['vertical_reference_events']=events
            (root/'profiles'/f"{case['id']}.in").write_text(profile_text(case),encoding='utf-8')
            cases.append(case)
        for key in ['tx_global_x_elevation_m','rx_global_x_elevation_m','domain_m','origin_global_x_elevation_m']:
            assert pair['coarse'][key]==pair['fine'][key]
        proofs.append(dict(station_m=s,same_actual_tx_rx=True,same_domain_edges=True,
            fine_retained_old_voxels_identical=True,common_y_m=ym,max_height_adjustment_m=max(abs(v) for v in row['proposed_y_change_m'])))
    c=copy.deepcopy(prior)
    for key in ['grid_fraction_m','horizontal_domain_shift_m','reference_checks','stage_C_condition','quantization_policy']:
        c.pop(key,None)
    c['official_sfcw'].update(f_stop_hz=110_000_000,steps=301)
    c.update(id='line9_matched_bscan_20261013_r1',schema='line9-matched-bscan/1',
        authorization='User approved strict same-position pair,51 traces each102 new solves, common0.1m source height and station-25.05m origin',
        cases=cases,anchors_m=[r['station_m'] for r in audit],trace_count_per_grid=51,station_step_m=.6,
        max_new_solves=102,max_wall_seconds=3600,reuse_references=[],references=[],
        no_retry_after_case_failure=True,source_contacts_sha256=sha(a.package/'model/contacts_m.json'),
        source_selection_sha256=sha(selpath),compatibility_audit_sha256=sha(a.compatibility),
        actual_pose_rule='nearest common0.1m y node; x origin=station-25.05m; global TX/RX=station±0.65m',
        coarse_grid_fraction_m=[1,30],fine_grid_m=.025,reference_checks=proofs,
        source_or_sfcw_changed=True,parent_contract_sha256=sha(a.seed/'contract.json'),prepare_script_sha256=sha(__file__),
        processing_input_ns=[0,1100],evaluation_windows_ns=[[100,400],[400,1000]],
        background_methods=[dict(name='identity',config_id='B0_G1_BG'),dict(name='mean',config_id='B3_G1_BG'),
            *[dict(name=f'sliding{w}',width=w) for w in [5,11,21]],
            *[dict(name=f'svd{k}',config_id=f'B{k+3}_G1_BG') for k in [1,2,3]]],
        current_model_modified=True,automatic_expansion=False,
        comparison_scope='Matched actual source/receiver,domain,materials; finite2.5cm reference, not clean truth; shared0.1m acquisition height steps retained')
    assert len(cases)==102 and max(v['max_height_adjustment_m'] for v in proofs)<.050000001
    save(root/'contract.json',c);save(root/'preparation_audit.json',dict(status='PASS',paired_checks=proofs,case_count=102))
    shutil.copy2(Path(__file__).with_name('study_line9_v31_smoothing_domain.py'),root/'study.py')
    save(root/'input_manifest.json',dict(files=[dict(path=f.relative_to(root).as_posix(),sha256=sha(f),bytes=f.stat().st_size)
        for f in sorted(root.rglob('*')) if f.is_file()]))
    print('Prepared102 new inputs; all51 pairs actual pose and domain identical; fine geometry overlaps verified')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare']);p.add_argument('--root',type=Path,required=True)
    for name in ['package','compatibility','seed','old-geometries']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a)
