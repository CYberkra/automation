"""Independent official-material, dimensional-CFL and timing limit checks."""
import argparse
import json
from pathlib import Path

import numpy as np
import gprMax.materials as official_materials
from gprMax import config
from gprMax.materials import DispersiveMaterial
from hs4t2d_physics_diagnostic import (Material, read_model, profile_at_x,
    vertical_roundtrip_ns, analyse_capsule, C0, EPS0, DEFAULT_CAPSULE)
from hs_capsule_identity import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--capsule', type=Path, default=DEFAULT_CAPSULE)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new check output required')
    model = read_model(args.capsule / 'hs4t2d_t01.in')
    cover = model['materials']['cover']
    official = DispersiveMaterial(4, 'cover')
    official.er, official.se, official.mr = cover.epsilon_inf, cover.conductivity, cover.mu_r
    official.type, official.poles = 'debye', len(cover.poles)
    official.deltaer, official.tau = [p[0] for p in cover.poles], [p[1] for p in cover.poles]
    frequencies = np.linspace(20e6, 170e6, 501)
    reference = np.array([official.calculate_er(f) for f in frequencies])
    error = float(np.max(np.abs(reference - cover.epsilon(frequencies)) / np.abs(reference)))
    # Difference is bounded by the two documented epsilon0 constants, not a
    # physical acceptance tolerance or a rounded material amplitude criterion.
    constant_error = np.max(np.abs(cover.conductivity / (2j*np.pi*frequencies) * (1/EPS0 - 1/config.e0)) / np.abs(reference))
    if error > constant_error + 20*np.finfo(float).eps:
        raise ValueError('official complex material implementation mismatch')
    if cover.static_dielectric_permittivity != 25.895:
        raise ValueError('static dielectric permittivity differs')
    water = Material(4.9, 0, 1, ((75.2, 9.231e-12),))
    if not np.isclose(water.static_dielectric_permittivity, 80.1):
        raise ValueError('official water limit differs')
    constant = Material(9, 0, 1)
    if not np.allclose(constant.index(frequencies).real, 3) or not np.allclose(constant.group_index(frequencies), 3):
        raise ValueError('nondispersive phase/group limit differs')
    if not np.isclose(vertical_roundtrip_ns(15, 3, 9, 3, 3), 2*(15 + 3*3 + 9*3)/C0*1e9):
        raise ValueError('roundtrip missed a leg')
    boxes = model['cover_boxes']
    j = next(k for k in range(1, len(boxes)) if boxes[k][2] != boxes[k-1][2])
    edge = boxes[j][0]
    if not np.array_equal(profile_at_x(model, [edge-1e-6, edge+1e-6]), [boxes[j-1][2], boxes[j][2]]):
        raise ValueError('staircase replaced with interpolated geometry')
    for bad in (0, -1, np.nan):
        try:
            cover.epsilon(bad)
        except ValueError:
            pass
        else:
            raise ValueError('invalid frequency accepted')
    report = analyse_capsule(args.capsule)
    if report['numerics']['raw_to_CFL_ratio'] != 1 or report['numerics']['active_axes'] != ['x', 'z']:
        raise ValueError('wrong dimensional CFL')
    if not report['event_183_8ns_verdict'].startswith('UNRESOLVED'):
        raise ValueError('unsupported identity from correlation')
    result = {'status': 'PASS', 'checks_passed': 9, 'calls_solver': False,
              'official_complex_epsilon_max_relative_difference': error,
              'epsilon0_constant_difference_bound': float(constant_error),
              'official_materials_py_sha256': sha256(official_materials.__file__),
              'code_sha256': sha256(__file__),
              'scope': 'material formula vs official implementation, analytic limits, actual2D CFL and geometry; not physical acceptance'}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
