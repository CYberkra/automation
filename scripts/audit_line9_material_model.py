"""Independent identity, acquisition, extrusion and borehole audit; no solver."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from build_pdf_profile_geometry import digest, save_json


def audit(out):
    target=out/'independent_verification.json'
    if target.exists(): raise ValueError('Refuse to overwrite verification')
    m=json.loads((out/'manifest.json').read_text('utf-8'))
    src=json.loads((out/'contacts_m.json').read_text('utf-8'))
    db=out/'geometries/line9_research_materials_v1.json'
    if digest(db)!=m['material_database_sha256']: raise ValueError('Material database changed')
    materials=json.loads(db.read_text('utf-8'))['materials']
    for entry in materials.values():
        if any(v is None for v in entry['base'].values()): raise ValueError('Unassigned constitutive property')
    geo_sha={}
    for name,g in m['geometries'].items():
        path=out/'geometries'/g['file']
        if digest(path)!=g['sha256']: raise ValueError('Geometry identity mismatch')
        geo_sha[name]=g['sha256']
        with h5py.File(path,'r') as h:
            assert tuple(h['data'].shape)==tuple(g['shape_nxyz']) and h['data'].dtype==np.int16
            assert set(h['material_keys'].asstr()[:])==set(materials) and len(materials)==4
            assert str(h.attrs['MaterialDatabase'])=='line9_research_materials_v1'
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],g['grid_xyz_m'])
            assert sum(g['material_counts'])==g['cells']
    max_baseline_error=0.; max_height_error=0.
    for c in m['cases']:
        path=out/c['input']; assert digest(path)==c['input_sha256']
        lines={s.split(':',1)[0]:s.split(':',1)[1].strip().split() for s in path.read_text('utf-8').splitlines()}
        tx=np.array(list(map(float,lines['#hertzian_dipole'][1:4])))
        rx=np.array(list(map(float,lines['#rx'][:3])))
        np.testing.assert_allclose(tx,c['tx_m'],rtol=0,atol=1e-9)
        np.testing.assert_allclose(rx,c['rx_m'],rtol=0,atol=1e-9)
        err=abs(np.linalg.norm(rx-tx)-1.3);max_baseline_error=max(max_baseline_error,err);assert err<1e-9
        expected=[0,0,1.3] if c['baseline']=='cross_track' else [1.3,0,0]
        np.testing.assert_allclose(rx-tx,expected,rtol=0,atol=1e-9)
        assert abs(c['acquisition_s_m']-(220-c['profile_x_m']))<1e-9
        g=next(g for g in m['geometries'].values() if g['file']==c['geometry'])
        surf=np.array(src['points_x_elevation_m']['surface'])
        agl=tx[1]+g['elevation_range_m'][0]-np.interp(c['profile_x_m'],surf[:,0],surf[:,1])
        max_height_error=max(max_height_error,abs(agl-15));assert abs(agl-15)<=g['grid_xyz_m'][0]/2+1e-8
        with h5py.File(out/'geometries'/g['file'],'r') as h:
            for point in [tx,rx]:
                idx=np.floor(point/np.array(g['grid_xyz_m'])).astype(int)
                idx[2]=min(idx[2],g['shape_nxyz'][2]-1)
                assert h['data'][tuple(idx)]==0,'Source/Rx not in native air voxel'
    stations=[c for c in m['cases'] if c['id'].startswith('full2d_s')]
    assert len(stations)==391
    np.testing.assert_allclose([c['profile_x_m'] for c in stations],220-np.arange(391)*.5,atol=0,rtol=0)
    for i in range(3):
        a=m['geometries'][f'local2d_{i}'];b=m['geometries'][f'local3d_{i}']
        with h5py.File(out/'geometries'/a['file'],'r') as h2,h5py.File(out/'geometries'/b['file'],'r') as h3:
            plane=h2['data'][:,:,0]
            for z in [0,b['shape_nxyz'][2]//2,b['shape_nxyz'][2]-1]:
                np.testing.assert_array_equal(plane,h3['data'][:,:,z])
    errors=[]
    bx=src['borehole_x_m'];collar=src['calibration']['borehole_collar_elevation_m']
    for name in ['full2d','local2d_1','local3d_1']:
        g=m['geometries'][name];dl=g['grid_xyz_m'][0]
        with h5py.File(out/'geometries'/g['file'],'r') as h:
            col=h['data'][int((bx-g['profile_x_range_m'][0])/dl),:,g['shape_nxyz'][2]//2]
        jumps=np.flatnonzero(col[1:]!=col[:-1])+1
        assert col[np.r_[0,jumps]][::-1].tolist()==[0,1,2,3,2,3]
        depths=collar-(g['elevation_range_m'][0]+jumps[:-1][::-1]*dl)
        error=float(np.max(abs(depths-[7,9.8,12.6,14.3])));assert error<.08;errors.append(error)
    result=dict(status='PASS_INDEPENDENT_MODEL_AUDIT_NOT_FDTD_VALIDATION',calls_solver=False,calls_training=False,
        manifest_sha256=digest(out/'manifest.json'),material_database_sha256=digest(db),script_sha256=digest(Path(__file__)),
        native_geometry_sha256=geo_sha,input_count=len(m['cases']),full2d_station_count=391,geometry_count=len(geo_sha),
        all_four_materials_assigned=True,maximum_baseline_absolute_error_m=max_baseline_error,
        maximum_midpoint_AGL_quantisation_error_m=max_height_error,maximum_borehole_interface_absolute_error_m=max(errors),
        matched_local2d_and_three_local3d_cross_sections_equal=True,
        checks=['geometry/material/input identities','all407 native air Tx/Rx pairs','195m reverse acquisition','along/cross1.3m baselines','borehole ordering/depths','matched local2D/3D extrusions'],
        resource_preflight_sha256=digest(out/'resource_preflight.json'))
    save_json(target,result);print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    audit(p.parse_args().out)
