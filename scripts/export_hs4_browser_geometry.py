"""Export archived HS4 input geometry for third-party browser inspection.

No gprMax invocation, geometry-only run, smoothing or vertical exaggeration.
Meshes represent input boxes, not a freshly verified solved Yee material grid.
Only positions are read from H5. Receiver time series are not accessed.
"""
import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import h5py
import numpy as np

from hs_capsule_identity import read_manifest, sha256, verify_file

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Mesh:
    def __init__(self):
        self.vertices, self.triangles, self.lookup = [], [], {}

    def vertex(self, point):
        key = tuple(float(v) for v in point)
        if key not in self.lookup:
            self.lookup[key] = len(self.vertices)
            self.vertices.append(key)
        return self.lookup[key]

    def quad(self, points):
        a, b, c, d = [self.vertex(p) for p in points]
        self.triangles.extend([(a, b, c), (a, c, d)])

    def clone(self, reverse=False):
        other = Mesh()
        other.vertices = list(self.vertices)
        other.lookup = dict(self.lookup)
        other.triangles = [(c, b, a) if reverse else (a, b, c)
                           for a, b, c in self.triangles]
        return other

    def volume(self):
        triangles = np.asarray(self.vertices)[np.asarray(self.triangles)]
        return float(np.einsum('ij,ij->i', triangles[:, 0],
                               np.cross(triangles[:, 1], triangles[:, 2])).sum() / 6)


def interface_mesh(z, step):
    mesh = Mesh()
    for i, j in np.ndindex(z.shape):
        x0, x1, y0, y1 = i * step, (i + 1) * step, j * step, (j + 1) * step
        a = float(z[i, j])
        mesh.quad([(x0, y0, a), (x1, y0, a), (x1, y1, a), (x0, y1, a)])
        if i + 1 < z.shape[0] and a != z[i + 1, j]:
            b = float(z[i + 1, j]); lo, hi = sorted((a, b))
            points = [(x1, y0, lo), (x1, y1, lo), (x1, y1, hi), (x1, y0, hi)]
            mesh.quad(points if a > b else points[::-1])
        if j + 1 < z.shape[1] and a != z[i, j + 1]:
            b = float(z[i, j + 1]); lo, hi = sorted((a, b))
            points = [(x0, y1, lo), (x0, y1, hi), (x1, y1, hi), (x1, y1, lo)]
            mesh.quad(points if a > b else points[::-1])
    return mesh


def solid_layer(interface, z, step, top, is_cover):
    mesh = interface.clone(reverse=is_cover)
    nx, ny = z.shape; width, length = nx * step, ny * step
    for j in range(ny):
        y0, y1 = j * step, (j + 1) * step
        left, right = float(z[0, j]), float(z[-1, j])
        a, b = (left, top) if is_cover else (0, left)
        mesh.quad([(0, y0, a), (0, y0, b), (0, y1, b), (0, y1, a)])
        a, b = (right, top) if is_cover else (0, right)
        mesh.quad([(width, y0, a), (width, y1, a), (width, y1, b), (width, y0, b)])
    for i in range(nx):
        x0, x1 = i * step, (i + 1) * step
        front, back = float(z[i, 0]), float(z[i, -1])
        a, b = (front, top) if is_cover else (0, front)
        mesh.quad([(x0, 0, a), (x1, 0, a), (x1, 0, b), (x0, 0, b)])
        a, b = (back, top) if is_cover else (0, back)
        mesh.quad([(x0, length, a), (x0, length, b), (x1, length, b), (x1, length, a)])
    if is_cover:
        mesh.quad([(0, 0, top), (width, 0, top), (width, length, top), (0, length, top)])
    else:
        mesh.quad([(0, 0, 0), (0, length, 0), (width, length, 0), (width, 0, 0)])
    return mesh


def markers(positions):
    mesh = Mesh()
    # Symbolic 0.12 m octahedra centred at H5 positions; not antenna dimensions.
    for position in positions:
        points = np.asarray(position) + .06 * np.asarray(
            [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])
        ids = [mesh.vertex(p) for p in points]
        for triangle in ((0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
                         (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)):
            mesh.triangles.append(tuple(ids[k] for k in triangle))
    return mesh


def write_obj(path, objects, material=False):
    with path.open('x', encoding='ascii', newline='\n') as f:
        f.write('# HS4 archived input geometry. Coordinates in metres; Z up.\n# No smoothing or exaggeration.\n')
        if material:
            f.write('mtllib hs4_scene.mtl\n')
        offset = 1
        for name, mesh in objects:
            f.write(f'o {name}\ns off\n')
            if material:
                f.write(f'usemtl {name}\n')
            for point in mesh.vertices:
                f.write('v ' + ' '.join(f'{v:.9g}' for v in point) + '\n')
            for triangle in mesh.triangles:
                f.write('f ' + ' '.join(str(v + offset) for v in triangle) + '\n')
            offset += len(mesh.vertices)


def write_vtp(path, mesh):
    root = ET.Element('VTKFile', type='PolyData', version='0.1', byte_order='LittleEndian')
    piece = ET.SubElement(ET.SubElement(root, 'PolyData'), 'Piece',
        NumberOfPoints=str(len(mesh.vertices)), NumberOfPolys=str(len(mesh.triangles)),
        NumberOfVerts='0', NumberOfLines='0', NumberOfStrips='0')
    def array(parent, name, dtype, values, components='1'):
        node = ET.SubElement(parent, 'DataArray', type=dtype, Name=name,
                             NumberOfComponents=components, format='ascii')
        node.text = ' '.join(str(v) for v in values)
    array(ET.SubElement(piece, 'PointData', Scalars='z_m'), 'z_m', 'Float64',
          [v[2] for v in mesh.vertices])
    array(ET.SubElement(piece, 'Points'), 'Points', 'Float64',
          [v for point in mesh.vertices for v in point], '3')
    polys = ET.SubElement(piece, 'Polys')
    array(polys, 'connectivity', 'Int32', [v for triangle in mesh.triangles for v in triangle])
    array(polys, 'offsets', 'Int32', range(3, len(mesh.triangles) * 3 + 1, 3))
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    require(not args.out.exists(), 'output directory must not exist')
    require(not args.out.resolve().is_relative_to(CAPSULE.resolve()), 'output inside raw capsule')
    records = read_manifest(CAPSULE)
    used, geometry = [], None
    z = np.full((48, 48), np.nan); step, top = .25, 12.0
    with np.load(CAPSULE / 'hs4_interface_binned_table.npz') as table:
        expected = table['z0_m'].copy()
        require(float(table['bin_m']) == step, 'unexpected bin size')
    used.append(verify_file(CAPSULE, records, 'hs4_interface_binned_table.npz'))
    positions = []
    for k in range(1, 26):
        prefix = f'hs4_col6_t{k:02d}'
        used.extend(verify_file(CAPSULE, records, prefix + suffix) for suffix in ('.in', '.h5'))
        lines = (CAPSULE / (prefix + '.in')).read_text('utf-8').splitlines()
        boxes = [line for line in lines if line.startswith('#box:')]
        if geometry is None:
            geometry = boxes
            rock = boxes[0].split()[1:]
            require(rock == ['0', '0', '0', '12', '12', '12', 'rock'], 'unexpected base rock')
            for line in boxes[1:]:
                x0, y0, bottom, x1, y1, upper = map(float, line.split()[1:7])
                i, j = int(round(x0 / step)), int(round(y0 / step))
                require(line.split()[7] == 'cover' and upper == top
                        and abs(x1 - x0 - step) < 1e-12 and abs(y1 - y0 - step) < 1e-12
                        and np.isnan(z[i, j]), 'unexpected/duplicate cover box')
                z[i, j] = bottom
        require(boxes == geometry, 'geometry differs between COL6 inputs')
        with h5py.File(CAPSULE / (prefix + '.h5'), 'r') as h:
            tx = np.asarray(h['srcs/src1'].attrs['Position'], dtype=float)
            rx = np.asarray(h['rxs/rx1'].attrs['Position'], dtype=float)
        require(np.array_equal(tx[[0, 2]], [1.6, 27])
                and np.array_equal(rx[[0, 2]], [1.6, 27]), 'unexpected solved positions')
        positions.append({'trace': k, 'tx_m': tx.tolist(), 'rx_m': rx.tolist()})
    require(np.array_equal(z, expected) and np.isfinite(z).all(), 'input/table geometry mismatch')
    interface = interface_mesh(z, step)
    rock = solid_layer(interface, z, step, top, False)
    cover = solid_layer(interface, z, step, top, True)
    expected_rock_volume = float(z.sum() * step ** 2)
    require(abs(rock.volume() - expected_rock_volume) < 1e-9, 'rock mesh volume mismatch')
    require(abs(cover.volume() - (12 * 12 * 12 - expected_rock_volume)) < 1e-9,
            'cover mesh volume mismatch')
    args.out.mkdir(parents=True)
    write_obj(args.out / 'hs4_interface.obj', [('interface', interface)])
    objects = [('rock', rock), ('cover', cover),
               ('Tx_COL6_markers', markers([p['tx_m'] for p in positions])),
               ('Rx_COL6_markers', markers([p['rx_m'] for p in positions]))]
    write_obj(args.out / 'hs4_scene.obj', objects, True)
    (args.out / 'hs4_scene.mtl').write_text(''.join(
        f'newmtl {name}\nKd {color}\nd 1\nillum 1\n\n' for name, color in (
            ('rock', '.45 .55 .65'), ('cover', '.75 .58 .32'),
            ('Tx_COL6_markers', '.9 .15 .15'), ('Rx_COL6_markers', '.1 .7 .2'))), encoding='ascii')
    write_vtp(args.out / 'hs4_interface.vtp', interface)
    with (args.out / 'col6_positions_and_profile.csv').open('x', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['trace', 'tx_x_m', 'tx_y_m', 'rx_x_m', 'rx_y_m', 'height_z_m',
                         'midpoint_y_m', 'interface_z_at_midpoint_m', 'cover_depth_at_midpoint_m'])
        for p in positions:
            tx, rx = np.asarray(p['tx_m']), np.asarray(p['rx_m'])
            mid = (tx + rx) / 2; i, j = np.floor(mid[:2] / step).astype(int)
            writer.writerow([p['trace'], *tx[:2], *rx[:2], tx[2], mid[1], z[i, j], top - z[i, j]])
    for entry in used:
        verify_file(CAPSULE, records, entry['file'])
    (args.out / 'export_provenance.json').write_text(json.dumps({
        'scope': 'archived input-box surfaces plus H5 position metadata, NOT solved material-grid export',
        'units': 'metres', 'z_up': True, 'vertical_exaggeration': 1, 'smoothing': False,
        'source_code_sha256': sha256(__file__), 'inputs': used,
        'cover_box_count': int(z.size), 'interface_z_range_m': [float(z.min()), float(z.max())],
        'surface_z_m': top, 'air_height_m': 15, 'marker_diameter_m': .12,
        'markers_are_symbolic_not_antenna_geometry': True,
        'volumes_m3': {'rock': rock.volume(), 'cover': cover.volume()},
        'positions': positions,
        'meshes': {name: {'vertices': len(mesh.vertices), 'triangles': len(mesh.triangles)}
                   for name, mesh in [('interface', interface), *objects]}},
        ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    manifest = [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha256(p)}
                for p in sorted(args.out.iterdir())]
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print('Verified input geometry exported:', args.out)


if __name__ == '__main__':
    main()
