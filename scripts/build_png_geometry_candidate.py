"""Write a bounded-memory V4 geometry candidate from confirmed four-colour IDs.

Run with the project's V4 Python. No solver is called. Extent is an explicit
design assumption, not measured site geometry. Mudstone remains unassigned;
air and the approved cover/sandstone design populate the editable database.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import gprMax
from gprMax.toolboxes.GeometryImport.common import write_null_material_database


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(4*1024*1024):
            value.update(block)
    return value.hexdigest()


def build(labels_path, design_path, out, extent, spacing, zcells):
    if out.exists():
        raise ValueError('Refuse to overwrite a geometry candidate')
    if gprMax.__version__ != '4.0.0':
        raise ValueError('Use the verified project V4.0.0 environment')
    if (not np.isfinite(extent).all() or not np.isfinite(spacing).all() or
            min(extent) <= 0 or min(spacing) <= 0 or zcells < 1):
        raise ValueError('Invalid positive dimensions, spacing, or zcells')
    ratio = np.asarray(extent)/np.asarray(spacing[:2])
    nx, ny = np.rint(ratio).astype(int)
    if not np.allclose(ratio, [nx, ny], rtol=0, atol=1e-8):
        raise ValueError('Extent must be a whole number of cells')
    labels_hash, design_hash = digest(labels_path), digest(design_path)
    with h5py.File(labels_path, 'r') as h5:
        labels = h5['nearest_colour_suggestion'][:]
        palette = h5['palette_rgb'][:]
        source_png_hash = h5.attrs['source_sha256']
    expected_palette = np.array([[217,217,217],[244,177,131],
                                 [191,144,0],[127,96,0]], dtype=np.uint8)
    if not np.array_equal(palette, expected_palette) or labels.min() < 0 or labels.max() > 3:
        raise ValueError('The confirmed four-colour project legend is required')
    x_pixels = np.floor((np.arange(nx)+.5)*labels.shape[1]/nx).astype(int)
    y_pixels = labels.shape[0]-1-np.floor((np.arange(ny)+.5)*labels.shape[0]/ny).astype(int)
    # PNG rows increase downward; model Y increases upward. X stays left-right.
    out.mkdir(parents=True)
    database_id = 'line9_candidate_materials'
    database_path = out/(database_id+'.json')
    names = ('air', 'cover', 'mudstone', 'sandstone')
    keys = write_null_material_database(
        database_path, database_id, names,
        source='User-confirmed PNG legend; candidate scale inherited from legacy H5',
        metadata=[dict(selected_colour=colour.tolist(),
                       source_png_sha256=source_png_hash,
                       geological_role=role)
                  for colour, role in zip(palette, ('air','silty clay','mudstone','sandstone'))])
    database = json.loads(database_path.read_text('utf-8'))
    database['database']['description'] = 'Geometry candidate: mudstone properties pending; not runnable as supplied.'
    materials = database['materials']
    materials[keys[0]]['base'] = dict(relative_permittivity=1.,
                                    electric_conductivity_s_per_m=0.,
                                    relative_permeability=1., magnetic_conductivity_s_per_m=0.)
    design = json.loads(design_path.read_text('utf-8'))
    for index, design_name in ((1,'cover'), (3,'rock')):
        item = design['materials'][design_name]
        entry = materials[keys[index]]
        entry['model'] = 'debye'
        entry['base'] = dict(relative_permittivity=item['epsilon_infinity'],
                            electric_conductivity_s_per_m=item['sigma_DC_S_m'],
                            relative_permeability=item['mu_r'],
                            magnetic_conductivity_s_per_m=item['sigma_magnetic'])
        entry['poles'] = [dict(relative_permittivity_difference=p['delta_epsilon'],
                               relaxation_time_s=p['tau_s']) for p in item['debye']]
        entry['metadata']['parameter_source'] = 'approved research material design; not measured site properties'
        entry['metadata']['parameter_design_sha256'] = design_hash
    materials[keys[2]]['metadata']['parameter_status'] = 'UNASSIGNED; legacy high-loss recipe not automatically adopted'
    database_path.write_text(json.dumps(database, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    geometry = out/'line9_geometry_candidate.h5'
    with h5py.File(geometry, 'x') as h5:
        data = h5.create_dataset('data', shape=(nx,ny,zcells), dtype=np.int16,
                                 chunks=(min(128,nx), min(512,ny), zcells), compression='gzip')
        for start in range(0,nx,256):
            block = labels[y_pixels[None,:], x_pixels[start:start+256,None]]
            data[start:start+256,:,:] = np.repeat(block[:,:,None], zcells, axis=2)
        h5.create_dataset('material_keys', data=np.asarray(keys, dtype='S'))
        h5.attrs['dx_dy_dz'] = spacing
        h5.attrs['shape_nxyz'] = (nx,ny,zcells)
        h5.attrs['origin_xyz'] = (0.,0.,0.)
        h5.attrs['MaterialDatabase'] = database_id
        h5.attrs['MaterialDatabaseSchemaVersion'] = 1
        h5.attrs['CandidateStatus'] = 'USER_LEGEND_CONFIRMED_SCALE_AND_MUDSTONE_PENDING'
        h5.attrs['SourcePNGSHA256'] = source_png_hash
    # Full chunk readback, plus independent weighted source-pixel conservation.
    counts = np.zeros(4, dtype=np.int64)
    with h5py.File(geometry, 'r') as h5:
        assert h5['data'].shape == (nx,ny,zcells)
        assert [k.decode() for k in h5['material_keys'][:]] == list(keys)
        for start in range(0,nx,256):
            block = h5['data'][start:start+256,:,:]
            expected = labels[y_pixels[None,:], x_pixels[start:start+256,None]]
            for z in range(zcells):
                if not np.array_equal(block[:,:,z], expected):
                    raise ValueError('Geometry readback or image orientation differs')
            counts += np.bincount(block.ravel(), minlength=4)
    wx = np.bincount(x_pixels, minlength=labels.shape[1])
    wy = np.bincount(y_pixels, minlength=labels.shape[0])
    expected_counts = np.array([np.sum((labels==i)*wy[:,None]*wx[None,:])*zcells for i in range(4)])
    if not np.array_equal(counts, expected_counts):
        raise ValueError('Source-pixel material conservation failed')
    if digest(labels_path) != labels_hash or digest(design_path) != design_hash:
        raise ValueError('Input changed during geometry preparation')
    summary = dict(status='PASS_GEOMETRY_READBACK_NOT_SOLVER_READY',
        source_png_sha256=source_png_hash, labels_sha256=labels_hash,
        material_design_sha256=design_hash, script_sha256=digest(Path(__file__)),
        geometry_sha256=digest(geometry), materials_sha256=digest(database_path),
        dimensions_xyz_m=[float(extent[0]),float(extent[1]),spacing[2]*zcells],
        spacing_xyz_m=list(spacing), shape_nxyz=[int(nx),int(ny),zcells],
        axis_order='XYZ; image top maps to high Y; image left maps to low X',
        scale_status='ASSUMED_LEGACY_EXTENT_NOT_MEASURED_PNG_SCALE',
        counts_by_material=counts.tolist(), material_keys=list(keys),
        user_legend='top to bottom: air, silty clay, mudstone, sandstone',
        segmentation_status='nearest-colour candidate; pixel changes and ties retained in source preparation',
        logical_geometry_bytes=int(nx*ny*zcells*2),
        missing_material_properties=['mudstone'],
        scope='Line9 site adaptation; ineligible for independent Line9 validation or generic development',
        is_strict_2d=(zcells==1), calls_solver=False, calls_training=False,
        originals_unchanged=True)
    (out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--labels', type=Path, required=True)
    parser.add_argument('--design', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--extent', type=float, nargs=2, required=True)
    parser.add_argument('--spacing', type=float, nargs=3, required=True)
    parser.add_argument('--zcells', type=int, required=True)
    args = parser.parse_args()
    build(args.labels,args.design,args.out,args.extent,args.spacing,args.zcells)
