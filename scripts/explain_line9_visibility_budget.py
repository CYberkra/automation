"""Local plane-wave primary-path budgets for audited Line9 packages; no solver."""
import argparse
import json
from pathlib import Path

import numpy as np

from review_line9_result_packages import C0, indices, save, sha


def main(root, review, out):
    if out.exists():
        raise ValueError('Use a new output directory')
    out.mkdir(parents=True)
    results = []
    for audit_path in sorted(review.glob('*_audit.json')):
        audit = json.loads(audit_path.read_text(encoding='utf-8'))
        material_path = next((root / audit['package'] / 'geometries').glob('*.json'))
        if sha(material_path) != audit['materials_sha256']:
            raise ValueError('Material identity changed')
        materials = json.loads(material_path.read_text(encoding='utf-8'))['materials']
        rows = audit['records']
        frequencies = []
        for frequency in [20e6, 40e6, 95e6, 170e6]:
            n = indices(materials, frequency)
            alpha = -2 * np.pi * frequency / C0 * n.imag
            reflection = lambda a, b: (n[a] - n[b]) / (n[a] + n[b])
            absorption = []
            relative_to_surface = []
            for row in rows:
                geometry = row['geometry']
                target = geometry['basal_sand']
                depth_previous = 0.
                transmission = 1. + 0j
                thickness = np.zeros(len(n))
                for boundary in geometry['boundaries']:
                    a, b = boundary['above'], boundary['below']
                    thickness[a] += boundary['depth_m'] - depth_previous
                    depth_previous = boundary['depth_m']
                    if boundary == target:
                        break
                    # Down/up electric-field transmission at each overlying interface.
                    transmission *= 4 * n[a] * n[b] / (n[a] + n[b]) ** 2
                else:
                    raise ValueError('Basal boundary not found')
                loss_db = -40 / np.log(10) * float(thickness @ alpha)
                absorption.append(loss_db)
                relative_to_surface.append(loss_db + 20 * np.log10(abs(
                    transmission * reflection(target['above'], target['below']) / reflection(0, 1))))
            frequencies.append(dict(
                frequency_MHz=frequency / 1e6,
                two_way_field_absorption_dB_per_m=(-40 / np.log(10) * alpha).tolist(),
                local_interface_reflection_dB={
                    'air_cover': float(20 * np.log10(abs(reflection(0, 1)))),
                    'cover_mudstone': float(20 * np.log10(abs(reflection(1, 2)))),
                    'mudstone_sandstone': float(20 * np.log10(abs(reflection(2, 3))))},
                basal_two_way_absorption_dB=[float(min(absorption)), float(max(absorption))],
                basal_primary_over_surface_field_dB=[float(min(relative_to_surface)), float(max(relative_to_surface))]))
        def limits(function):
            values = [function(row['geometry']) for row in rows]
            return [float(min(values)), float(max(values))]
        results.append(dict(package=audit['package'], native_traces=len(rows),
            audit_sha256=sha(audit_path), materials_sha256=sha(material_path),
            surface_y_m=limits(lambda g: g['surface_y_m']),
            cover_base_y_m=limits(lambda g: g['cover_base']['y_m']),
            cover_thickness_m=limits(lambda g: g['cover_base']['depth_m']),
            cover_time95_ns=limits(lambda g: g['cover_base']['time95_ns']),
            basal_time95_ns=limits(lambda g: g['basal_sand']['time95_ns']),
            frequencies=frequencies))
    if len(results) != 3:
        raise ValueError('Expected three package audits')
    save(out / 'visibility_budget.json', dict(calls_solver=False, calls_training=False,
        script_sha256=sha(__file__), packages=results,
        assumptions='Local horizontal nonmagnetic plane-wave primary paths; complex-index normal-incidence Fresnel coefficients and material absorption only. No geometrical spreading, antenna directivity, direct coupling, multiples, or 2D lateral scattering. Per-frequency field ratios are not measured broadband SNR or verified detections.'))
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'out']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.review, args.out)
