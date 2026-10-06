"""Build a fresh private Line9 material/input package without running FDTD."""
import argparse
import copy
import json
import math
import subprocess
from pathlib import Path

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT/'configs/research/line9_model_plan_v0_1.json'
DATABASE_ID = 'line9_research_materials_v1'


def classify(x, elevation, contacts):
    curves = {k:np.interp(x,np.asarray(v)[:,0],np.asarray(v)[:,1])
              for k,v in contacts.items()}
    e = elevation[None,:]
    result = np.full((len(x),len(elevation)),3,dtype=np.int16)
    result[(e>=curves['lower_mudstone_base'][:,None]) &
           (e<curves['lower_mudstone_top'][:,None])] = 2
    upper = x>=contacts['upper_mudstone_base'][0][0]
    result[upper[:,None] & (e>=curves['upper_mudstone_base'][:,None]) &
           (e<curves['clay_base'][:,None])] = 2
    result[e>=curves['clay_base'][:,None]]=1
    result[e>=curves['surface'][:,None]]=0
    return result


def make_geometry(path, xbounds, zbounds, dl, yrange, contacts, source_sha):
    shape = tuple(int(round(v/dl)) for v in [xbounds[1]-xbounds[0],yrange[1]-yrange[0],zbounds[1]-zbounds[0]])
    assert np.allclose(np.array(shape)*dl,[xbounds[1]-xbounds[0],yrange[1]-yrange[0],zbounds[1]-zbounds[0]],rtol=0,atol=1e-9)
    nx,ny,nz=shape
    elevation=yrange[0]+(np.arange(ny)+.5)*dl
    counts=np.zeros(4,dtype=np.int64)
    with h5py.File(path,'x') as h:
        data=h.create_dataset('data',shape=shape,dtype='int16',chunks=(min(nx,8),min(ny,256),min(nz,32)),compression='gzip',compression_opts=1)
        for first in range(0,nx,8):
            last=min(nx,first+8)
            x=xbounds[0]+(np.arange(first,last)+.5)*dl
            plane=classify(x,elevation,contacts)
            data[first:last]=np.broadcast_to(plane[:,:,None],(last-first,ny,nz))
            counts+=np.bincount(plane.ravel(),minlength=4)*nz
        h.create_dataset('material_keys',data=np.asarray(['material_000_air','material_001_cover','material_002_mudstone','material_003_sandstone'],dtype='S'))
        h.attrs['dx_dy_dz']=[dl]*3
        h.attrs['origin_xyz']=[0.,0.,0.]
        h.attrs['shape_nxyz']=shape
        h.attrs['MaterialDatabase']=DATABASE_ID
        h.attrs['MaterialDatabaseSchemaVersion']=1
        h.attrs['profile_x_offset_m']=xbounds[0]
        h.attrs['elevation_offset_m']=yrange[0]
        h.attrs['cross_track_offset_m']=zbounds[0]
        h.attrs['SourcePDFSHA256']=source_sha
        h.attrs['CandidateStatus']='COMPLETE_RESEARCH_CONSTITUTIVE_MODEL_NOT_FDTD_VALIDATED'
    # Check every voxel; use separate scalar-column intervals for the reference.
    mismatches=0
    with h5py.File(path,'r') as h:
        for first in range(0,nx,8):
            last=min(nx,first+8)
            expected=np.full((last-first,ny),3,dtype=np.int16)
            for j,i in enumerate(range(first,last)):
                px=xbounds[0]+(i+.5)*dl
                height=lambda name:float(np.interp(px,np.asarray(contacts[name])[:,0],np.asarray(contacts[name])[:,1]))
                expected[j,(elevation>=height('lower_mudstone_base')) & (elevation<height('lower_mudstone_top'))]=2
                if px>=contacts['upper_mudstone_base'][0][0]:
                    expected[j,(elevation>=height('upper_mudstone_base')) & (elevation<height('clay_base'))]=2
                expected[j,elevation>=height('clay_base')]=1
                expected[j,elevation>=height('surface')]=0
            actual=h['data'][first:last]
            mismatches+=int(np.count_nonzero(actual!=expected[:,:,None]))
    if mismatches or counts.sum()!=math.prod(shape):
        raise ValueError('Full voxel readback or count failure')
    nodes=math.prod([v+1 for v in shape])
    cells=math.prod(shape)
    # Exact primary array payloads from FDTDGrid initialisers; all-Debye histories are real.
    host_lower=22*cells+(6*8+6*4+3*8)*nodes
    gpu_lower=(6*8+6*4+3*8)*nodes
    return dict(file=path.name,sha256=digest(path),shape_nxyz=shape,grid_xyz_m=[dl]*3,
                profile_x_range_m=xbounds,elevation_range_m=yrange,cross_track_range_m=zbounds,
                material_counts=counts.tolist(),full_readback_mismatches=mismatches,cells=cells,
                host_primary_array_lower_bytes=host_lower,device_primary_array_lower_bytes=gpu_lower,
                logical_int16_geometry_bytes=2*cells,file_bytes=path.stat().st_size)


def material_database(out, design_paths, source):
    from gprMax.material_database import validate_material_database
    old=json.loads((source/'line9_pdf_materials.json').read_text('utf-8'))
    db=copy.deepcopy(old)
    db['database'].update(id=DATABASE_ID,name='Line9 fixed research constitutive assignments',description='Complete four-material research database; not site measured spectra')
    mapping={}
    for path in design_paths:
        design=json.loads(path.read_text('utf-8'))
        mapping.update({name:(p,path) for name,p in design['materials'].items()})
    for key,name in [('material_001_cover','cover'),('material_002_mudstone','mudstone'),('material_003_sandstone','rock')]:
        p,path=mapping[name]
        entry=db['materials'][key]
        entry['model']='debye'
        entry['base']=dict(relative_permittivity=p['epsilon_infinity'],electric_conductivity_s_per_m=p['sigma_DC_S_m'],relative_permeability=p['mu_r'],magnetic_conductivity_s_per_m=p['sigma_magnetic'])
        entry['poles']=[dict(relative_permittivity_difference=q['delta_epsilon'],relaxation_time_s=q['tau_s']) for q in p['debye']]
        entry['metadata'].update(parameter_status='FIXED_RESEARCH_NOT_SITE_CALIBRATION',parameter_design_sha256=digest(path))
        entry['averagable']=False
    target=out/'geometries'/f'{DATABASE_ID}.json'
    save_json(target,db)
    catalogue=validate_material_database(DATABASE_ID,search_directory=target.parent)
    if len(catalogue)!=4 or any(s.relative_permittivity is None for _,s in catalogue):
        raise ValueError('Incomplete V4 material catalogue')
    return target


def write_case(out, name, geo, px, baseline, p, contacts):
    from gprMax.hash_cmds_file import get_user_objects
    dl=geo['grid_xyz_m'][0]
    width=geo['profile_x_range_m'][1]-geo['profile_x_range_m'][0]
    height=p['elevation_range_m'][1]-p['elevation_range_m'][0]
    cross=geo['cross_track_range_m'][1]-geo['cross_track_range_m'][0]
    y=float(np.interp(px,np.asarray(contacts['surface'])[:,0],np.asarray(contacts['surface'])[:,1]))+p['flight_height_midpoint_AGL_m']-p['elevation_range_m'][0]
    y=round(y/dl)*dl
    center=np.array([px-geo['profile_x_range_m'][0],y,cross/2])
    offset=np.array([.65,0,0]) if baseline=='along_track' else np.array([0,0,.65])
    tx,rx=center-offset,center+offset
    dimensions=np.array([width,height,cross])
    active_axes=[0,1] if geo['shape_nxyz'][2]==1 else [0,1,2]
    for point in [tx,rx]:
        if any(not (2<point[i]<dimensions[i]-2) for i in active_axes):
            raise ValueError('Station inside PML')
        if classify(np.array([point[0]+geo['profile_x_range_m'][0]]),np.array([point[1]+p['elevation_range_m'][0]]),contacts)[0,0]!=0:
            raise ValueError('Antenna point outside air')
        if not np.allclose(point[active_axes]/dl,np.rint(point[active_axes]/dl),atol=1e-8,rtol=0):
            raise ValueError('Active source/receiver coordinates not grid aligned')
    n=int(round(2/dl))
    pml=[n,n,0 if geo['shape_nxyz'][2]==1 else n]*2
    fmt=lambda values:' '.join(f'{v:.12g}' for v in values)
    text='\n'.join([f'#title: {name}; Line9 research hypothesis;15m AGL midpoint; not calibrated S21',
        f'#domain: {fmt(dimensions)}',f'#dx_dy_dz: {fmt([dl]*3)}',
        f"#time_window: {p['time_window_s']:.12g}",'#omp_threads: 2',
        f'#pml_cells: {fmt(pml)}','#pml_formulation: HORIPML',
        f"#waveform: ricker {1/dl:.12g} {p['source']['frequency_hz']:.12g} pulse",
        f'#hertzian_dipole: z {fmt(tx)} pulse',f'#rx: {fmt(rx)} rx1 Ez',
        f"#geometry_objects_read: 0 0 0 ../../geometries/{geo['file']} {DATABASE_ID} n",''])
    case=out/'cases'/name
    case.mkdir(parents=True)
    path=case/'profile.in'
    path.write_text(text,encoding='utf-8')
    objects=get_user_objects(text.splitlines(),input_dir=case)
    if len(objects)<9:
        raise ValueError('V4 command parsing incomplete')
    return dict(id=name,input=str(path.relative_to(out)).replace('\\','/'),input_sha256=digest(path),geometry=geo['file'],profile_x_m=px,acquisition_s_m=220-px,baseline=baseline,tx_m=tx.tolist(),rx_m=rx.tolist(),nominal_current_moment_A_m=1,command_parser_object_count=len(objects),dt_CFL_initial_s=dl/(299792458*math.sqrt(len(active_axes))))


def prepare(source,out):
    import gprMax
    from gprMax.grid import fdtd_grid
    from gprMax import config
    if out.exists():
        raise ValueError('Refuse to overwrite model package')
    if gprMax.__version__!='4.0.0':
        raise ValueError('Use audited gprMax4.0.0 environment')
    p=json.loads(PLAN.read_text('utf-8'))
    contacts_path=source/'contacts_m.json'
    source_geo=source/'line9_pdf_geometry.h5'
    src=json.loads(contacts_path.read_text('utf-8'))
    contacts=src['points_x_elevation_m']
    tracked_inputs=[PLAN,ROOT/p['acquisition_selection'],*[ROOT/v for v in p['material_designs']]]
    input_hashes={str(path):digest(path) for path in [*tracked_inputs,contacts_path,source_geo,source/'line9_pdf_materials.json']}
    if digest(source_geo)!='2fd42ce39dd9270ab252b4e4fa5efa24ea31f1a7ff14a5b0d4e0a700a6ce48d0':
        raise ValueError('Reviewed source geometry hash changed')
    out.mkdir(parents=True)
    (out/'geometries').mkdir()
    database=material_database(out,[ROOT/v for v in p['material_designs']],source)
    geos={}
    yrange=p['elevation_range_m']
    geos['full2d']=make_geometry(out/'geometries/full2d.h5',p['full2d']['profile_x_domain_m'],[0,.025],.025,yrange,contacts,src['source_pdf_sha256'])
    print('Full2D geometry: all voxels and complete V4 materials PASS',flush=True)
    cases=[]
    stations=220-np.arange(p['full2d']['station_count'])*p['full2d']['station_spacing_m']
    if stations[-1]!=25 or not np.all(np.diff(stations)<0):
        raise ValueError('Acquisition endpoints/order mismatch')
    for i,px in enumerate(stations):
        cases.append(write_case(out,f'full2d_s{i:04d}',geos['full2d'],float(px),'along_track',p,contacts))
    pilot=[]
    for i,px in enumerate(p['full2d']['pilot_profile_x_m']):
        name=f'full2d_pilot_{i}'
        cases.append(write_case(out,name,geos['full2d'],px,'along_track',p,contacts)); pilot.append(name)
    for i,px in enumerate(p['local3d']['pilot_profile_x_m']):
        for dim,cross in [('local2d',.05),('local3d',24.)]:
            key=f'{dim}_{i}'
            geos[key]=make_geometry(out/f'geometries/{key}.h5',[px-12,px+12],[0,cross],.05,yrange,contacts,src['source_pdf_sha256'])
            print(f'{key} geometry full readback PASS',flush=True)
            for baseline in (['along_track'] if dim=='local2d' else ['along_track','cross_track']):
                name=f'{key}_{baseline}'
                cases.append(write_case(out,name,geos[key],px,baseline,p,contacts))
        print(f'Local station{i}: matched2D, along/cross3D inputs parsed',flush=True)
    px=p['local3d']['width_controls_at_profile_x_m']
    for i,(xw,zw) in enumerate(p['local3d']['width_control_x_z_m']):
        key=f'local3d_wide_{i}'
        geos[key]=make_geometry(out/f'geometries/{key}.h5',[px-xw/2,px+xw/2],[0,zw],.05,yrange,contacts,src['source_pdf_sha256'])
        cases.append(write_case(out,key,geos[key],px,'cross_track',p,contacts))
        print(f'{key} geometry full readback PASS',flush=True)
    if any(digest(Path(path))!=expected for path,expected in input_hashes.items()):
        raise ValueError('Sources changed during preparation')
    manifest=dict(status='COMPLETE_MODEL_INPUT_PACKAGE_NOT_EXECUTED',calls_solver=False,calls_training=False,
        gprMax_version=gprMax.__version__,plan_sha256=digest(PLAN),source_inputs_sha256=input_hashes,
        script_sha256=digest(Path(__file__)),material_database_sha256=digest(database),
        installed_native_sources_sha256={str(path):digest(path) for path in [Path(fdtd_grid.__file__),Path(config.__file__)]},
        geometries=geos,cases=cases,full2d_pilot_ids=pilot,
        private_scope='Contains PDF-derived geology and station elevations; keep this whole package local',
        coordinate_convention=p['coordinates'],limits=p['limits'])
    save_json(out/'manifest.json',manifest)
    save_json(out/'plan_snapshot.json',p)
    save_json(out/'contacts_m.json',src)
    (out/'README.txt').write_text('Prepared V4 model/input package, not executed.\nStart with resource_preflight.json. Never run391 stations before pilots, resource checks and a new machine-specific execution contract.\nStrict2D is along-track surrogate; true cross-track1.3m baseline is local3D.\nUse gprMax4.0.0 -gpu_precision double; raw dtype must be checked.\nNo physical antenna-port S21 or site electrical calibration.\n',encoding='utf-8')
    print(f'Model preparation PASS: {len(cases)} parsed inputs; {len(geos)} audited geometries',flush=True)


def preflight(out):
    import psutil
    path=out/'resource_preflight.json'
    if path.exists():
        raise ValueError('Refuse to overwrite resource snapshot')
    m=json.loads((out/'manifest.json').read_text('utf-8'))
    p=json.loads((out/'plan_snapshot.json').read_text('utf-8'))
    v=psutil.virtual_memory()
    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.free','--format=csv,noheader,nounits'],text=True).strip().splitlines()[0].split(',')
    free_vram=int(gpu[2].strip())*2**20
    records={}
    factor=p['runtime_policy']['allocation_margin_factor']
    for name,g in m['geometries'].items():
        need_ram=math.ceil(g['host_primary_array_lower_bytes']*factor)+int(1.5*2**30)
        need_vram=math.ceil(g['device_primary_array_lower_bytes']*factor)+int(.5*2**30)
        records[name]=dict(required_available_RAM_bytes=need_ram,required_free_VRAM_bytes=need_vram,
            passes_available_RAM=v.available>=need_ram,passes_free_VRAM=free_vram>=need_vram)
    full_ready=all(records['full2d'][k] for k in ['passes_available_RAM','passes_free_VRAM'])
    local_ready=all(r[k] for name,r in records.items() if name.startswith('local3d') for k in ['passes_available_RAM','passes_free_VRAM'])
    save_json(path,dict(status='PASS_CAPACITY_ONLY_NOT_EXECUTION_AUTHORIZATION' if full_ready and local_ready else 'PREPARED_NOT_RUN_RESOURCE_PREFLIGHT_REJECTED',
        manifest_sha256=digest(out/'manifest.json'),calls_solver=False,calls_training=False,
        observed=dict(available_RAM_bytes=v.available,total_RAM_bytes=v.total,gpu_name=gpu[0].strip(),total_VRAM_bytes=int(gpu[1].strip())*2**20,free_VRAM_bytes=free_vram),
        records=records,estimation_scope=p['runtime_policy']['memory_estimate'],full2d_capacity_pass=full_ready,all_local3d_capacity_pass=local_ready,
        note='Local2D fit does not authorize replacing the full model with a cropped domain. Recheck target machine and freeze execution separately. No process was terminated.'))
    print(f'Capacity full2D={full_ready}; all local3D={local_ready}; no solver called')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='stage',required=True)
    prep=sub.add_parser('prepare')
    prep.add_argument('--source',type=Path,required=True)
    prep.add_argument('--out',type=Path,required=True)
    pre=sub.add_parser('preflight')
    pre.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.stage=='prepare': prepare(args.source,args.out)
    else: preflight(args.out)
