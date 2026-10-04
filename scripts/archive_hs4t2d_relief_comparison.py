"""Archive completed raw traces, one geometry/group and failed-attempt evidence."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np

from hs_capsule_identity import sha256


def copy_verified(source, target, digest=None):
    before = digest or sha256(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha256(source) != before or sha256(target) != before:
        raise ValueError('copy identity failed')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new archive required')
    c = json.loads(args.contract.read_text('utf-8'))
    run = Path(c['run_directory'])
    events = [json.loads(s) for s in (run/'series_execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('traces') != 363 or events[-1].get('status') != 'COMPLETED':
        raise ValueError('incomplete run')
    if events[0]['contract_sha256'] != sha256(args.contract):
        raise ValueError('contract identity changed')
    copy_verified(run/'series_execution.jsonl', args.out/'series_execution.jsonl')
    grids = []
    for g in c['groups']:
        folder = run/g['id']
        audit = json.loads((folder/'audit.json').read_text('utf-8'))
        if len(audit) != 121 or sha256(folder/'profile.in') != g['sha256']:
            raise ValueError('input/audit differs')
        for row in audit:
            copy_verified(folder/row['file'], args.out/g['id']/row['file'], row['sha256'])
        for name in ['profile.in', 'audit.json', 'stdout.log', 'stderr.log', 'cuda_cache.jsonl']:
            copy_verified(folder/name, args.out/g['id']/name)
        copy_verified(folder/'hs4t2d_geom1.vtkhdf', args.out/g['id']/'hs4t2d_geom1.vtkhdf')
        reference = None
        for k in range(1,122):
            p = folder/f'hs4t2d_geom{k}.vtkhdf'
            with h5py.File(p,'r') as h:
                material = h['VTKHDF/CellData/Material'][:]
                if reference is None:
                    reference = material
                elif not np.array_equal(material, reference):
                    raise ValueError('geometry varies across geometry-fixed traces')
            grids.append({'group':g['id'], 'file':p.name, 'sha256':sha256(p),
                          'same_material_array_as_trace1':True})
    attempts = ['2026-10-03_hs4t2d_relief08_local', '2026-10-03_hs4t2d_relief08_cuda_fresh',
                '2026-10-03_hs4t2d_relief08_cuda_ascii_att2', '2026-10-03_hs4t2d_relief08_cuda_ascii_att3']
    for name in attempts:
        folder = run.parent/name
        for p in sorted(folder.rglob('*')):
            if p.is_file() and p.suffix in ('.h5','.json','.jsonl','.log','.in'):
                copy_verified(p, args.out/'prior_attempts'/name/p.relative_to(folder))
    (args.out/'geometry_exports.json').write_text(json.dumps(grids,indent=2)+'\n',encoding='utf-8')
    report = {'contract_sha256':sha256(args.contract), 'raw_h5':363,
              'geometry_exports_verified':363, 'representative_geometry_files_saved':3,
              'geometry_scope':'all material arrays within each fixed-geometry group equal trace1; full file hashes retained',
              'archive_code_sha256':sha256(__file__), 'prior_attempts':attempts}
    (args.out/'archive_provenance.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    rows = [{'file':p.relative_to(args.out).as_posix(), 'bytes':p.stat().st_size, 'sha256':sha256(p)}
            for p in sorted(args.out.rglob('*')) if p.is_file()]
    (args.out/'manifest.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
