"""Extract reviewed PDF polylines, build a private V4 geometry, and preview it.

Three separate CLI stages keep PDF/VTK dependencies outside the V4 environment.
The selection file is private site data: explicit paths, calibration and crop.
No solver, antenna, acquisition model, or geological interpretation is inferred.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(4*1024*1024):
            value.update(block)
    return value.hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def extract(pdf, out):
    import pypdfium2 as pdfium
    from pypdfium2 import raw
    from ctypes import c_float, c_uint, byref
    if out.exists():
        raise ValueError('Refuse to overwrite vector extraction')
    source_hash = digest(pdf)
    doc = pdfium.PdfDocument(pdf)
    if len(doc) != 1:
        raise ValueError('Reviewed profile must have exactly one page')
    page = doc[0]
    paths = []
    for obj in page.get_objects():
        if obj.type != raw.FPDF_PAGEOBJ_PATH:
            continue
        bounds = obj.get_bounds()
        if bounds[2]-bounds[0] < 250:
            continue
        matrix = obj.get_matrix()
        segments = []
        for i in range(raw.FPDFPath_CountSegments(obj.raw)):
            seg = raw.FPDFPath_GetPathSegment(obj.raw, i)
            x, y = c_float(), c_float()
            raw.FPDFPathSegment_GetPoint(seg, byref(x), byref(y))
            segments.append(dict(type=raw.FPDFPathSegment_GetType(seg),
                                 point=matrix.on_point(x.value, y.value)))
        rgba = [c_uint() for _ in range(4)]
        raw.FPDFPageObj_GetStrokeColor(obj.raw, *[byref(v) for v in rgba])
        paths.append(dict(bounds=bounds, stroke_rgba=[v.value for v in rgba],
                          segments=segments))
    if digest(pdf) != source_hash:
        raise ValueError('PDF changed during extraction')
    save_json(out, dict(source_pdf_sha256=source_hash,
                       extraction='Pdfium transformed paths; width >=250 page points',
                       page=1, page_size_pt=page.get_size(), paths=paths))
    print(f'Extracted {len(paths)} wide paths without changing PDF')


def ordered_points(paths, indices):
    blocks = []
    for index in indices:
        path = paths[index]
        if path['stroke_rgba'] != [0, 0, 0, 255]:
            raise ValueError('Reviewed contacts must be black paths')
        if any(s['type'] not in (0, 2) for s in path['segments']):
            raise ValueError('Only reviewed straight polylines are supported')
        points = np.array([s['point'] for s in path['segments']], dtype=float)
        if np.all(np.diff(points[:, 0]) < 0):
            points = points[::-1]
        if not np.all(np.diff(points[:, 0]) > 0):
            raise ValueError('Contact X must increase strictly')
        blocks.append(points)
    blocks.sort(key=lambda a: a[0, 0])
    # At a tiny CAD overlap, retain the preceding segment up to the join.
    combined = blocks[0]
    for block in blocks[1:]:
        gap = block[0, 0]-combined[-1, 0]
        if abs(gap) > 1.0 or abs(block[0, 1]-combined[-1, 1]) > 1.0:
            raise ValueError('Contact pieces do not meet within one PDF point')
        combined = np.vstack([combined[combined[:, 0] < block[0, 0]], block])
    return combined


def write_materials(out, design_path, source_hash):
    from gprMax.toolboxes.GeometryImport.common import write_null_material_database
    database_id = 'line9_pdf_materials'
    path = out/(database_id+'.json')
    keys = write_null_material_database(path, database_id,
        ('air', 'cover', 'mudstone', 'sandstone'),
        source='Reviewed geological PDF contacts; electrical properties are research assumptions',
        metadata=[dict(source_pdf_sha256=source_hash, geological_role=role)
                  for role in ('air', 'silty clay', 'mudstone', 'sandstone')])
    database = json.loads(path.read_text('utf-8'))
    database['database']['description'] = 'PDF geometry only; mudstone unassigned; no solver execution.'
    materials = database['materials']
    materials[keys[0]]['base'] = dict(relative_permittivity=1., electric_conductivity_s_per_m=0.,
        relative_permeability=1., magnetic_conductivity_s_per_m=0.)
    design = json.loads(design_path.read_text('utf-8'))
    for index, name in ((1, 'cover'), (3, 'rock')):
        item = design['materials'][name]
        entry = materials[keys[index]]
        entry['model'] = 'debye'
        entry['base'] = dict(relative_permittivity=item['epsilon_infinity'],
            electric_conductivity_s_per_m=item['sigma_DC_S_m'],
            relative_permeability=item['mu_r'], magnetic_conductivity_s_per_m=item['sigma_magnetic'])
        entry['poles'] = [dict(relative_permittivity_difference=p['delta_epsilon'],
                             relaxation_time_s=p['tau_s']) for p in item['debye']]
        entry['metadata']['parameter_design_sha256'] = digest(design_path)
        entry['metadata']['parameter_status'] = 'RESEARCH_DESIGN_NOT_SITE_MEASUREMENT'
    materials[keys[2]]['metadata']['parameter_status'] = 'UNASSIGNED'
    save_json(path, database)
    return keys, database_id, path


def build(vectors_path, selection_path, design_path, out):
    import gprMax
    if out.exists():
        raise ValueError('Refuse to overwrite an earlier geometry package')
    if gprMax.__version__ != '4.0.0':
        raise ValueError('Use the verified V4.0.0 Python environment')
    inputs = {str(p): digest(p) for p in (vectors_path, selection_path, design_path)}
    vectors = json.loads(vectors_path.read_text('utf-8'))
    selection = json.loads(selection_path.read_text('utf-8'))
    if vectors['source_pdf_sha256'] != selection['source_pdf_sha256']:
        raise ValueError('Selection must be reviewed against this exact PDF hash')
    calibration = selection['calibration']
    width = float(calibration['profile_width_m'])
    source = {name: ordered_points(vectors['paths'], ids)
              for name, ids in selection['contacts'].items()}
    surface = source['surface']
    scale = (surface[-1, 0]-surface[0, 0])/width
    anchor_x = float(calibration['borehole_pdf_x_pt'])
    anchor_y = np.interp(anchor_x, surface[:, 0], surface[:, 1])
    collar = float(calibration['borehole_collar_elevation_m'])
    contacts = {name: np.column_stack([(p[:, 0]-surface[0, 0])/scale,
                                       (p[:, 1]-anchor_y)/scale+collar])
                for name, p in source.items()}
    # Equal horizontal/vertical scale is explicitly verified in the source.
    borehole_x = (anchor_x-surface[0, 0])/scale
    depths = {name: collar-np.interp(borehole_x, points[:, 0], points[:, 1])
              for name, points in contacts.items() if name != 'surface'}
    for name, expected in calibration['borehole_contact_depths_m'].items():
        if abs(depths[name]-expected) > .02:
            raise ValueError(f'Borehole depth mismatch for {name}: {depths[name]}')
    bottom, top = map(float, selection['crop_elevation_m'])
    spacing = np.asarray(selection['spacing_xyz_m'], dtype=float)
    zcells = int(selection['zcells'])
    if (not np.isfinite(spacing).all() or min(spacing) <= 0 or zcells < 1 or top <= bottom):
        raise ValueError('Invalid crop or grid')
    ratios = np.array([width, top-bottom])/spacing[:2]
    nx, ny = np.rint(ratios).astype(int)
    if not np.allclose(ratios, [nx, ny], rtol=0, atol=1e-8):
        raise ValueError('Crop must contain an integer number of cells')
    x = (np.arange(nx)+.5)*spacing[0]
    y = bottom+(np.arange(ny)+.5)*spacing[1]
    curves = {name: np.interp(x, p[:, 0], p[:, 1]) for name, p in contacts.items()}
    upper_exists = x >= contacts['upper_mudstone_base'][0, 0]
    # Finite upper bed terminates at the mapped clay contact, not by extension leftward.
    if (np.any(curves['surface'] <= curves['clay_base']) or
            np.any(curves['lower_mudstone_top'] <= curves['lower_mudstone_base']) or
            np.any(curves['upper_mudstone_base'][upper_exists] < curves['lower_mudstone_top'][upper_exists])):
        raise ValueError('Contact ordering is inconsistent')
    out.mkdir(parents=True)
    keys, database_id, materials_path = write_materials(out, design_path, vectors['source_pdf_sha256'])
    geometry_path = out/'line9_pdf_geometry.h5'
    counts = np.zeros(4, dtype=np.int64)
    with h5py.File(geometry_path, 'x') as handle:
        data = handle.create_dataset('data', shape=(nx, ny, zcells), dtype='int16',
            chunks=(min(nx, 128), min(ny, 512), zcells), compression='gzip')
        for start in range(0, nx, 128):
            stop = min(nx, start+128)
            z = y[None, :]
            block = np.full((stop-start, ny), 3, dtype=np.int16)
            mud_lower = ((z >= curves['lower_mudstone_base'][start:stop, None]) &
                         (z < curves['lower_mudstone_top'][start:stop, None]))
            mud_upper = ((z >= curves['upper_mudstone_base'][start:stop, None]) &
                         (z < curves['clay_base'][start:stop, None]) & upper_exists[start:stop, None])
            block[mud_lower | mud_upper] = 2
            block[z >= curves['clay_base'][start:stop, None]] = 1
            block[z >= curves['surface'][start:stop, None]] = 0
            data[start:stop] = np.repeat(block[:, :, None], zcells, axis=2)
            counts += np.bincount(block.ravel(), minlength=4)*zcells
        handle.create_dataset('material_keys', data=np.asarray(keys, dtype='S'))
        handle.attrs['dx_dy_dz'] = spacing
        handle.attrs['shape_nxyz'] = (nx, ny, zcells)
        handle.attrs['origin_xyz'] = (0., 0., 0.)
        handle.attrs['elevation_offset_m'] = bottom
        handle.attrs['MaterialDatabase'] = database_id
        handle.attrs['MaterialDatabaseSchemaVersion'] = 1
        handle.attrs['SourcePDFSHA256'] = vectors['source_pdf_sha256']
        handle.attrs['CandidateStatus'] = 'PDF_CONTACTS_GRID_VIEW_ONLY_MUDSTONE_UNASSIGNED'
    # Independently classify every read-back column with interval priority.
    read_counts = np.zeros(4, dtype=np.int64)
    failures = 0
    with h5py.File(geometry_path, 'r') as handle:
        for i in range(nx):
            expected = np.full(ny, 3, dtype=np.int16)
            for lo, hi in [(curves['lower_mudstone_base'][i], curves['lower_mudstone_top'][i])]:
                expected[(y >= lo) & (y < hi)] = 2
            if upper_exists[i]:
                expected[(y >= curves['upper_mudstone_base'][i]) & (y < curves['clay_base'][i])] = 2
            expected[y >= curves['clay_base'][i]] = 1
            expected[y >= curves['surface'][i]] = 0
            actual = handle['data'][i, :, :]
            failures += np.count_nonzero(actual != expected[:, None])
            read_counts += np.bincount(actual.ravel(), minlength=4)
        if failures or not np.array_equal(counts, read_counts):
            raise ValueError('Full geometry read-back failed')
    if any(digest(Path(path)) != value for path, value in inputs.items()):
        raise ValueError('Inputs changed during build')
    save_json(out/'contacts_m.json', dict(source_pdf_sha256=vectors['source_pdf_sha256'],
        points_x_elevation_m={k: p.tolist() for k, p in contacts.items()},
        calibration=calibration, pdf_points_per_m=scale,
        borehole_x_m=borehole_x, borehole_depths_m=depths))
    save_json(out/'summary.json', dict(status='PASS_GEOMETRY_NOT_SOLVER_READY',
        source_pdf_sha256=vectors['source_pdf_sha256'], input_sha256=inputs,
        script_sha256=digest(Path(__file__)), geometry_sha256=digest(geometry_path),
        materials_sha256=digest(materials_path), shape_nxyz=[int(nx), int(ny), zcells],
        spacing_xyz_m=spacing.tolist(), dimensions_xyz_m=[width, top-bottom, zcells*spacing[2]],
        axis_convention='XYZ; X=profile distance, Y=elevation minus crop bottom; Z=invariant thin direction',
        crop_elevation_m=[bottom, top], borehole_x_m=borehole_x, borehole_depths_m=depths,
        material_counts=counts.tolist(), full_readback_mismatches=int(failures),
        missing_material_properties=['mudstone'], geometry_mesh_status='viewing discretisation, not FDTD convergence',
        source_contact_selection=selection['contacts'], source_reconciliation=selection['source_reconciliation'],
        below_borehole='Mapped lowest sandstone continues to chosen viewing crop; not borehole-confirmed below 20.2m',
        is_strict_2d=zcells == 1, scope='Line9 site adaptation; no independent Line9 validation',
        calls_solver=False, calls_training=False))
    print(f'Geometry readback PASS: {nx} x {ny} x {zcells}, ZK08 depths {depths}')


def preview(out):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    summary = json.loads((out/'summary.json').read_text('utf-8'))
    contacts = json.loads((out/'contacts_m.json').read_text('utf-8'))
    vti = out/'line9_pdf_geometry.vti'
    png = out/'line9_pdf_geometry.png'
    if vti.exists() or png.exists():
        raise ValueError('Refuse to overwrite a preview')
    with h5py.File(out/'line9_pdf_geometry.h5', 'r') as handle:
        data = handle['data'][:]
    spacing = summary['spacing_xyz_m']
    bottom, top = summary['crop_elevation_m']
    width = summary['dimensions_xyz_m'][0]
    image = vtk.vtkImageData()
    image.SetDimensions(*(np.array(data.shape)+1))
    image.SetSpacing(*spacing)
    image.SetOrigin(0., bottom, 0.)
    values = numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_SHORT)
    values.SetName('MaterialID')
    image.GetCellData().SetScalars(values)
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(str(vti))
    writer.SetInputData(image)
    if writer.Write() != 1:
        raise ValueError('VTK write failed')
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti))
    reader.Update()
    actual = reader.GetOutput()
    expected_bounds = [0, width, bottom, top, 0, data.shape[2]*spacing[2]]
    if (not np.allclose(actual.GetBounds(), expected_bounds, rtol=0, atol=1e-10) or
            not np.array_equal(vtk_to_numpy(actual.GetCellData().GetScalars()), data.ravel(order='F'))):
        raise ValueError('VTK bounds or cell ordering failed')
    colours = ['#d9d9d9', '#f4b183', '#bf9000', '#7f6000']
    names = ['空气', '粉质黏土', '泥岩', '砂岩']
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(2, 1, figsize=(18, 7), layout='constrained')
    for ax in axes:
        ax.imshow(data[:, :, 0].T, origin='lower', extent=[0, width, bottom, top],
            aspect='equal', interpolation='nearest', cmap=ListedColormap(colours),
            norm=BoundaryNorm(np.arange(-.5, 4.5), 4))
        for name, points in contacts['points_x_elevation_m'].items():
            points = np.array(points)
            ax.plot(points[:, 0], points[:, 1], color='#403d35', lw=.7)
        x = contacts['borehole_x_m']
        collar = contacts['calibration']['borehole_collar_elevation_m']
        ax.plot([x, x], [collar-20.2, collar], color='#234488', lw=1.5)
        ax.text(x+1, collar+.8, 'ZK08', color='#234488', fontsize=10)
        ax.set(xlabel='剖面水平距离 (m)', ylabel='高程 (m)')
        ax.grid(alpha=.15)
    axes[0].set_title('九号线 PDF 折线重建：横纵等比例，水平范围 0–350 m', fontsize=15)
    axes[0].set_xlim(0, width)
    axes[0].set_ylim(bottom, top)
    axes[1].set_title('ZK08 附近局部放大：仍保持横纵等比例', fontsize=13)
    axes[1].set_xlim(145, 250)
    axes[1].set_ylim(420, 450)
    fig.legend(handles=[Patch(facecolor=c, label=n) for c, n in zip(colours, names)],
        loc='outside lower center', ncol=4)
    fig.savefig(png, dpi=150)
    plt.close(fig)
    save_json(out/'preview_validation.json', dict(status='PASS', vtk_bounds=actual.GetBounds(),
        vtk_cells_equal_h5=True, vtk_sha256=digest(vti), png_sha256=digest(png),
        grid_spacing_xyz_m=spacing, vertical_exaggeration=1, calls_solver=False))
    print('Equal-scale PNG and native VTK cell order/bounds PASS')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    p = sub.add_parser('extract')
    p.add_argument('--pdf', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('build')
    for name in ('vectors', 'selection', 'design', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p = sub.add_parser('preview')
    p.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'extract':
        extract(args.pdf, args.out)
    elif args.stage == 'build':
        build(args.vectors, args.selection, args.design, args.out)
    else:
        preview(args.out)
