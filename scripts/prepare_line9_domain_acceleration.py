"""Prepare exact HDF5 slices for a Line9 domain-width study; never run a solver."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import numpy as np

from build_pdf_profile_geometry import digest, save_json


def prepare(package, out, evidence):
    if out.exists() or evidence.exists():
        raise ValueError('Fresh private package and evidence required')
    manifest_path = package/'manifest.json'
    manifest_sha = digest(manifest_path)
    manifest = json.loads(manifest_path.read_text('utf-8'))
    reference = manifest['geometries']['full2d']
    original = package/'geometries'/reference['file']
    if digest(original) != reference['sha256']:
        raise ValueError('Full-domain geometry changed')
    database = package/'geometries/line9_research_materials_v1.json'
    if digest(database) != manifest['material_database_sha256']:
        raise ValueError('Material database changed')
    cases = {c['id']: c for c in manifest['cases']}
    pilots = [cases[n] for n in manifest['full2d_pilot_ids']]
    dl = .025
    original_left, original_right = reference['profile_x_range_m']
    out.mkdir(parents=True)
    (out/'geometries').mkdir()
    shutil.copyfile(database, out/'geometries'/database.name)
    rows = []
    with h5py.File(original, 'r') as src:
        np.testing.assert_array_equal(src.attrs['dx_dy_dz'], [dl]*3)
        if src['data'].shape != (16000, 3000, 1):
            raise ValueError('Unexpected original grid')
        for width in [200, 160]:
            nx = round(width/dl)
            for pilot in pilots:
                source_input = package/pilot['input']
                if digest(source_input) != pilot['input_sha256']:
                    raise ValueError('Original input changed')
                # Clamp inside the original domain: no extrapolated geology.
                left = max(original_left, min(pilot['profile_x_m']-width/2,
                                              original_right-width))
                offset = (left-original_left)/dl
                if not np.isclose(offset, round(offset), atol=1e-9, rtol=0):
                    raise ValueError('Crop does not preserve the original Yee lattice')
                first = round(offset)
                name = f'w{width}_{pilot["id"]}'
                geometry = out/'geometries'/f'{name}.h5'
                expected_digest = hashlib.sha256()
                with h5py.File(geometry, 'x') as target:
                    for k, v in src.attrs.items():
                        target.attrs[k] = v
                    target.attrs['shape_nxyz'] = [nx, 3000, 1]
                    target.attrs['profile_x_offset_m'] = left
                    target.attrs['CandidateStatus'] = 'EXACT_CROP_NOT_DOMAIN_EQUIVALENCE_ACCEPTED'
                    src.copy('material_keys', target)
                    data = target.create_dataset('data', shape=(nx, 3000, 1), dtype='int16',
                        chunks=(8, 256, 1), compression='gzip', compression_opts=1)
                    for j in range(0, nx, 64):
                        block = src['data'][first+j:first+min(j+64, nx), :, :]
                        data[j:min(j+64, nx)] = block
                        expected_digest.update(block.tobytes(order='C'))
                observed_digest = hashlib.sha256()
                with h5py.File(geometry, 'r') as target:
                    np.testing.assert_array_equal(target['material_keys'][:], src['material_keys'][:])
                    for j in range(0, nx, 64):
                        block = target['data'][j:min(j+64, nx), :, :]
                        np.testing.assert_array_equal(block, src['data'][first+j:first+min(j+64, nx), :, :])
                        observed_digest.update(block.tobytes(order='C'))
                    translated = []
                    for key in ['tx_m', 'rx_m']:
                        position = np.array(pilot[key], dtype=float)
                        position[0] -= left-original_left
                        indices = np.rint(position[:2]/dl).astype(int)
                        if not np.allclose(position[:2]/dl, indices, atol=1e-9, rtol=0):
                            raise ValueError('Source/receiver grid phase changed')
                        if target['data'][indices[0], indices[1], 0] != 0:
                            raise ValueError('Source/receiver is not air')
                        if position[0] <= 2 or position[0] >= width-2:
                            raise ValueError('Source/receiver enters PML')
                        translated.append(position.tolist())
                if expected_digest.digest() != observed_digest.digest():
                    raise ValueError('All-voxel crop identity failed')
                tx, rx = translated
                fmt = lambda p: ' '.join(f'{v:.12g}' for v in p)
                lines = []
                for line in source_input.read_text('utf-8').splitlines():
                    if line.startswith('#domain:'):
                        line = f'#domain: {width} 75 0.025'
                    elif line.startswith('#hertzian_dipole:'):
                        line = f'#hertzian_dipole: z {fmt(tx)} pulse'
                    elif line.startswith('#rx:'):
                        line = f'#rx: {fmt(rx)} rx1 Ez'
                    elif line.startswith('#geometry_objects_read:'):
                        line = f'#geometry_objects_read: 0 0 0 ../../geometries/{name}.h5 line9_research_materials_v1 n'
                    lines.append(line)
                input_path = out/'cases'/name/'profile.in'
                input_path.parent.mkdir(parents=True)
                input_path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
                commands = lambda text: {s.split(':', 1)[0]: s.split(':', 1)[1] for s in text.splitlines() if s.startswith('#')}
                before, after = commands(source_input.read_text('utf-8')), commands(input_path.read_text('utf-8'))
                allowed = {'#domain', '#hertzian_dipole', '#rx', '#geometry_objects_read'}
                if set(before) != set(after) or any(before[k] != after[k] for k in before.keys()-allowed):
                    raise ValueError('Unapproved input factor changed')
                rows.append(dict(id=name, reference_id=pilot['id'], width_m=width,
                    source_input_sha256=digest(source_input), input=str(input_path.relative_to(out)),
                    input_sha256=digest(input_path), geometry=str(geometry.relative_to(out)),
                    geometry_sha256=digest(geometry), profile_x_range_m=[left, left+width],
                    tx_m=tx, rx_m=rx, all_voxel_reference_match=True,
                    cells=nx*3000, work_fraction_of_reference=width/400,
                    voxel_content_sha256=observed_digest.hexdigest()))
    if digest(original) != reference['sha256'] or digest(manifest_path) != manifest_sha or digest(database) != manifest['material_database_sha256']:
        raise ValueError('Original package changed during preparation')
    save_json(out/'manifest.json', dict(status='PREPARED_NOT_RUN_NO_QUALITY_ACCEPTANCE',
        calls_solver=False, calls_training=False, reference_manifest_sha256=manifest_sha,
        reference_geometry_sha256=digest(original), material_database_sha256=digest(database),
        cases=rows, invariant='Exact original voxels and Yee positions; same height/grid/PML/material/source/time/dtype/SFCW. Only horizontal domain and translated local coordinates change.',
        limitations='No FDTD equivalence demonstrated. This package is not an execution contract; do not reuse the full-domain runner/auditor for these shapes.'))
    evidence.mkdir(parents=True)
    save_json(evidence/'preparation_summary.json', dict(status='EXACT_CROPS_PREPARED_NOT_SIMULATED',
        calls_solver=False, case_count=len(rows), widths_m=[200, 160],
        exact_voxel_checks_passed=len(rows), source_and_receiver_grid_phase_preserved=True,
        vertical_domain_m=75, grid_m=dl, time_window_s=1.2e-6,
        primary_array_work_fractions=[.5, .4],
        reference_geometry_sha256=reference['sha256'], private_manifest_sha256=digest(out/'manifest.json'),
        script_sha256=digest(Path(__file__)), quality_equivalence='NOT_TESTED',
        estimated_speedups_are_workload_ratios_not_measured=True))
    print(f'Prepared {len(rows)} exact cropped inputs; no solver or quality acceptance')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.package.resolve(), args.out.resolve(), args.evidence.resolve())
