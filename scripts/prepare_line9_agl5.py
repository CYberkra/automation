"""Derive Line9 AGL-5m standard inputs from the frozen AGL-8m baselines; never solve.

User approval 2026-10-11: flight height default changed 8 m -> 5 m, no validation
batch. Only #title/#hertzian_dipole/#rx lines may change; antenna y is recomputed
from the frozen geometry surface, not hardcoded.
"""
import argparse
import json
from pathlib import Path

import hashlib
import h5py
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

ROOT = Path(__file__).resolve().parents[1]
DL = 0.025
AGL_NEW = 5.0
AGL_OLD_EXPECTED = 8.0


def surface_y(data, x_m):
    col = data[round(x_m / DL)]
    return (np.flatnonzero(col != 0)[-1] + 1) * DL


def shift_antennas(text, surf_tx, surf_rx):
    lines = text.splitlines(keepends=True)
    dipoles = [i for i, s in enumerate(lines) if s.startswith('#hertzian_dipole:')]
    rxs = [i for i, s in enumerate(lines) if s.startswith('#rx: ')]
    titles = [i for i, s in enumerate(lines) if s.startswith('#title:')]
    assert len(dipoles) == 1 and len(rxs) == 1 and len(titles) == 1
    d = lines[dipoles[0]].split()
    r = lines[rxs[0]].split()
    old_tx_y, old_rx_y = float(d[3]), float(r[2])
    assert abs((old_tx_y - surf_tx) - AGL_OLD_EXPECTED) < 1e-9, 'tx baseline AGL != 8 m'
    assert abs((old_rx_y - surf_rx) - AGL_OLD_EXPECTED) < 1e-9, 'rx baseline AGL != 8 m'
    new_tx_y, new_rx_y = surf_tx + AGL_NEW, surf_rx + AGL_NEW
    assert abs(new_tx_y / DL - round(new_tx_y / DL)) < 1e-9, 'new tx y off grid'
    assert abs(new_rx_y / DL - round(new_rx_y / DL)) < 1e-9, 'new rx y off grid'
    d[3] = f'{new_tx_y:.15g}'
    r[2] = f'{new_rx_y:.15g}'
    for i, parts in [(dipoles[0], d), (rxs[0], r)]:
        end = '\r\n' if lines[i].endswith('\r\n') else '\n'
        lines[i] = ' '.join(parts) + end
    t = lines[titles[0]]
    end = '\r\n' if t.endswith('\r\n') else '\n'
    lines[titles[0]] = t.rstrip('\r\n') + '; AGL 5 m (user-approved 2026-10-11, was 8 m)' + end
    return ''.join(lines), dict(tx_y_m=new_tx_y, rx_y_m=new_rx_y)


def prepare(baseline_dir, geometry_dir, out):
    if out.exists():
        raise ValueError('Fresh package required; historical inputs stay read-only')
    out.mkdir(parents=True)
    groups = []
    for role in ['H0', 'H1']:
        with h5py.File(geometry_dir / f'{role}.h5') as h:
            data = h['data'][:, :, 0]
        card = baseline_dir / f'x190_{role}' / 'profile.in'
        original = card.read_bytes().decode('utf-8')
        d = next(s.split() for s in original.splitlines() if s.startswith('#hertzian_dipole:'))
        r = next(s.split() for s in original.splitlines() if s.startswith('#rx: '))
        tx_x, rx_x = float(d[2]), float(r[1])
        changed, pos = shift_antennas(original, surface_y(data, tx_x), surface_y(data, rx_x))
        dst = out / f'x190_{role}' / 'profile.in'
        dst.parent.mkdir(parents=True)
        dst.write_bytes(changed.encode('utf-8'))
        before, after = original.splitlines(), changed.splitlines()
        assert [s for s in before if not s.startswith(('#title:', '#hertzian_dipole:', '#rx: '))] == \
               [s for s in after if not s.startswith(('#title:', '#hertzian_dipole:', '#rx: '))]
        groups.append(dict(id=f'x190_{role}_agl5', role=role, profile_x_m=190.0,
                           tx_m=[tx_x, pos['tx_y_m'], 0.0125], rx_m=[rx_x, pos['rx_y_m'], 0.0125],
                           input=str(dst.relative_to(out)), input_sha256=sha(dst),
                           parent_input_sha256=sha(card)))
    manifest = dict(status='PREPARED_AGL5_NOT_RUN',
                    approval_basis='2026-10-11 用户显性批准：航高默认 8 m 改为 5 m，不跑验证批',
                    generator_sha256=sha(__file__), agl_m=AGL_NEW, groups=groups,
                    no_retry=True, calls_solver=False)
    (out / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n',
                                       encoding='utf-8')
    print(json.dumps(dict(status=manifest['status'], inputs=len(groups), new_solver_runs=0)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--geometries', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    prepare(a.baseline, a.geometries, a.out)
