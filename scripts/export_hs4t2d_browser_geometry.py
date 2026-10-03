"""Export the executed HS4T2D arch/step input cross-sections without solving."""
import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from export_hs4_browser_geometry import Mesh, write_obj
from hs_capsule_identity import read_manifest, verify_file, sha256

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_hs4t2d_transect'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists() or args.out.resolve().is_relative_to(CAPSULE.resolve()):
        raise ValueError('use a new output directory outside the raw capsule')
    records = read_manifest(CAPSULE)
    used, results = [], {}
    used.append(verify_file(CAPSULE, records, 'hs4t2d_transect_table.npz'))
    with np.load(CAPSULE / 'hs4t2d_transect_table.npz') as t:
        expected_arch = t['prof_m'].copy()
    args.out.mkdir(parents=True)
    for prefix, name in [('hs4t2d', 'arch'), ('hs4t2ds', 'step')]:
        boxes, positions = None, []
        for k in range(1, 122):
            tag = f'{prefix}_t{k:02d}'
            used.extend(verify_file(CAPSULE, records, tag + suffix) for suffix in ('.in', '.h5'))
            lines = (CAPSULE / (tag + '.in')).read_text('utf-8').splitlines()
            current = [line for line in lines if line.startswith('#box:')]
            if boxes is None:
                boxes = current
            if current != boxes or '#domain: 12 0.05 33' not in lines:
                raise ValueError('geometry/domain differs across archived inputs')
            with h5py.File(CAPSULE / (tag + '.h5'), 'r') as h:
                tx = np.asarray(h['srcs/src1'].attrs['Position'], dtype=float)
                rx = np.asarray(h['rxs/rx1'].attrs['Position'], dtype=float)
            if not np.array_equal(tx[[1, 2]], [0, 27]) or not np.array_equal(rx[[1, 2]], [0, 27]):
                raise ValueError('unexpected solved positions')
            positions.append({'trace': k, 'tx_m': tx.tolist(), 'rx_m': rx.tolist()})
        if boxes[0].split()[1:] != ['0', '0', '0', '12', '0.05', '12', 'rock']:
            raise ValueError('unexpected base geometry')
        rock, cover, bins, cursor = Mesh(), Mesh(), [], 0.0
        for line in boxes[1:]:
            words = line.split()
            x0, y0, z0, x1, y1, top = map(float, words[1:7])
            if x0 != cursor or not (x1 > x0 and y0 == 0 and y1 == .05 and top == 12
                                    and 0 < z0 < top and words[7] == 'cover'):
                raise ValueError('unsupported cover geometry')
            # A planar section at Y=0.025; invariant Y is not a finite 3D target.
            rock.quad([(x0, .025, 0), (x1, .025, 0), (x1, .025, z0), (x0, .025, z0)])
            cover.quad([(x0, .025, z0), (x1, .025, z0), (x1, .025, 12), (x0, .025, 12)])
            bins.append([x0, x1, z0, 12 - z0])
            cursor = x1
        if cursor != 12:
            raise ValueError('incomplete cross-section')
        if name == 'arch' and not np.array_equal([b[2] for b in bins], expected_arch):
            raise ValueError('input/table profile differs')
        if name == 'step' and bins != [[0, 6, 9.4, 12 - 9.4], [6, 12, 8.6, 12 - 8.6]]:
            raise ValueError('unexpected archived step')
        write_obj(args.out / f'hs4t2d_{name}.obj', [('rock', rock), ('cover', cover)], True)
        with (args.out / f'{name}_profile.csv').open('x', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['x0_m', 'x1_m', 'interface_z_m', 'cover_depth_m'])
            w.writerows(bins)
        results[name] = {'bins': len(bins), 'interface_z_range_m': [min(b[2] for b in bins), max(b[2] for b in bins)],
                         'positions': positions, 'vertices': len(rock.vertices) + len(cover.vertices),
                         'triangles': len(rock.triangles) + len(cover.triangles)}
    (args.out / 'hs4_scene.mtl').write_text(
        'newmtl rock\nKd .45 .55 .65\nd 1\nillum 1\n\n'
        'newmtl cover\nKd .75 .58 .32\nd 1\nillum 1\n', encoding='ascii')
    for entry in used:
        verify_file(CAPSULE, records, entry['file'])
    provenance = {'scope': 'executed input XZ cross-sections; not solved material-grid export',
                  'units': 'metres', 'z_up': True, 'vertical_exaggeration': 1, 'section_y_m': .025,
                  'invariant_axis': 'y', 'domain_m': [12, .05, 33], 'ground_z_m': 12,
                  'air_domain_omitted_from_geology_mesh': True,
                  'source_code_sha256': sha256(__file__),
                  'mesh_writer_sha256': sha256(ROOT / 'scripts/export_hs4_browser_geometry.py'),
                  'inputs': used, 'models': results}
    (args.out / 'export_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    files = [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in sorted(args.out.iterdir())]
    (args.out / 'manifest.json').write_text(json.dumps(files, indent=2) + '\n', encoding='utf-8')
    print('Exported arch and step sections; verified input identities:', len(used))


if __name__ == '__main__':
    main()
