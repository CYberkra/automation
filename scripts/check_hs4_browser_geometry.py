"""Reparse browser exports independently and compare against archived HS4 inputs."""
import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from hs_capsule_identity import read_manifest, sha256, verify_file

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_obj(path):
    vertices, objects, name = [], {}, None
    for line in path.read_text('ascii').splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == 'o':
            name = words[1]
            require(name not in objects, 'duplicate object')
            objects[name] = []
        elif words[0] == 'v':
            vertices.append([float(v) for v in words[1:]])
        elif words[0] == 'f':
            require(len(words) == 4 and name is not None, 'non-triangle or missing object')
            objects[name].append([int(v) - 1 for v in words[1:]])
    vertices = np.asarray(vertices)
    require(np.isfinite(vertices).all(), 'nonfinite vertex')
    for faces in objects.values():
        require(np.min(faces) >= 0 and np.max(faces) < len(vertices), 'invalid vertex index')
    return vertices, {k: np.asarray(v) for k, v in objects.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--directory', type=Path, required=True)
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    records = json.loads((args.directory / 'manifest.json').read_text('utf-8'))
    for entry in records:
        path = args.directory / entry['file']
        require(path.stat().st_size == entry['bytes'] and sha256(path) == entry['sha256'],
                'export manifest mismatch: ' + entry['file'])
    provenance = json.loads((args.directory / 'export_provenance.json').read_text('utf-8'))
    require(sha256(ROOT / 'scripts/export_hs4_browser_geometry.py') == provenance['source_code_sha256'],
            'export script identity mismatch')
    original = read_manifest(CAPSULE)
    for entry in provenance['inputs']:
        require(verify_file(CAPSULE, original, entry['file']) == entry, 'input identity mismatch')
    with np.load(CAPSULE / 'hs4_interface_binned_table.npz') as table:
        z = table['z0_m'].copy()
    vertices, objects = read_obj(args.directory / 'hs4_interface.obj')
    require(set(objects) == {'interface'}, 'unexpected interface objects')
    triangles = vertices[objects['interface']]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    require((np.linalg.norm(cross, axis=1) > 0).all(), 'degenerate surface')
    flat = np.ptp(triangles[:, :, 2], axis=1) == 0
    require(int(flat.sum()) == 2 * z.size and (cross[flat, 2] > 0).all(), 'plateau count/orientation')
    centers = triangles[flat].mean(axis=1)
    bins = np.floor(centers[:, :2] / .25).astype(int)
    require(np.array_equal(centers[:, 2], z[bins[:, 0], bins[:, 1]]), 'plateau heights mismatch')
    area = np.zeros(z.shape)
    np.add.at(area, (bins[:, 0], bins[:, 1]), cross[flat, 2] / 2)
    require(np.array_equal(area, np.full(z.shape, .25 ** 2)), 'plateau area mismatch')
    jump_area = float(np.linalg.norm(cross[~flat], axis=1).sum() / 2)
    expected_jump = float(.25 * (np.abs(np.diff(z, axis=0)).sum() + np.abs(np.diff(z, axis=1)).sum()))
    require(abs(jump_area - expected_jump) < 1e-10, 'vertical step area mismatch')
    piece = ET.parse(args.directory / 'hs4_interface.vtp').find('./PolyData/Piece')
    vtp_vertices = np.fromstring(piece.find('./Points/DataArray').text, sep=' ').reshape(-1, 3)
    connectivity = np.fromstring(piece.find('./Polys/DataArray[@Name="connectivity"]').text,
                                 sep=' ', dtype=int).reshape(-1, 3)
    offsets = np.fromstring(piece.find('./Polys/DataArray[@Name="offsets"]').text, sep=' ', dtype=int)
    scalar = np.fromstring(piece.find('./PointData/DataArray[@Name="z_m"]').text, sep=' ')
    require(np.array_equal(vtp_vertices, vertices) and np.array_equal(connectivity, objects['interface'])
            and np.array_equal(offsets, np.arange(3, len(connectivity) * 3 + 1, 3))
            and np.array_equal(scalar, vertices[:, 2]), 'VTP/OBJ/scalar mismatch')
    require(int(piece.get('NumberOfPoints')) == len(vertices)
            and int(piece.get('NumberOfPolys')) == len(connectivity), 'VTP declared count mismatch')
    scene_vertices, scene_objects = read_obj(args.directory / 'hs4_scene.obj')
    require(set(scene_objects) == {'rock', 'cover', 'Tx_COL6_markers', 'Rx_COL6_markers'}, 'scene objects')
    volumes = {}
    for name, expected in [('rock', float(z.sum() * .25 ** 2)),
                           ('cover', float(12 ** 3 - z.sum() * .25 ** 2))]:
        t = scene_vertices[scene_objects[name]]
        volumes[name] = float(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6)
        require(abs(volumes[name] - expected) < 1e-9, 'scene volume mismatch')
    for name, key in [('Tx_COL6_markers', 'tx_m'), ('Rx_COL6_markers', 'rx_m')]:
        t = scene_vertices[scene_objects[name]]
        require(len(t) == 25 * 8, 'marker face count')
        for k, p in enumerate(provenance['positions']):
            points = np.unique(t[k * 8:(k + 1) * 8].reshape(-1, 3), axis=0)
            require(len(points) == 6 and np.allclose(points.mean(axis=0), p[key], rtol=0, atol=1e-12),
                    'marker center mismatch')
    rows = list(csv.DictReader((args.directory / 'col6_positions_and_profile.csv').open(encoding='utf-8')))
    require(len(rows) == 25, 'CSV trace count')
    depths = []
    for row, p in zip(rows, provenance['positions']):
        tx, rx = np.asarray(p['tx_m']), np.asarray(p['rx_m'])
        mid = (tx + rx) / 2
        i, j = np.floor(mid[:2] / .25).astype(int)
        expected = [p['trace'], *tx[:2], *rx[:2], tx[2], mid[1], z[i, j], 12 - z[i, j]]
        require(np.array_equal([float(v) for v in row.values()], expected), 'CSV profile mismatch')
        depths.append(expected[-1])
    result = {'status': 'PASS', 'input_hashes_verified': len(provenance['inputs']),
              'manifest_entries_verified': len(records), 'plateaus_verified': int(z.size),
              'interface_vertices': len(vertices), 'interface_triangles': len(triangles),
              'vertical_step_area_m2': jump_area, 'volumes_m3': volumes,
              'col6_midpoint_depth_range_m': [min(depths), max(depths)],
              'scene_vertices': len(scene_vertices),
              'scene_triangles': sum(len(v) for v in scene_objects.values()),
              'checker_sha256': sha256(__file__),
              'scope': 'input geometry and export consistency; not solver-grid or waveform validation'}
    if args.out:
        with args.out.open('x', encoding='utf-8') as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
