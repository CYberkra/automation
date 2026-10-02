"""Fail-closed acceptance and provenance for development-only carrier reruns.

No solver imports. Historical results are read-only. H5 manifests for t3/B2
must be supplied explicitly; S1/S3 already has a frozen per-trace manifest.
Missing source hashes are never inferred from the output being accepted.
"""

import hashlib
import json
import platform
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def file_record(path):
    path = Path(path).resolve()
    try:
        name = path.relative_to(ROOT).as_posix()
    except ValueError:
        name = str(path)
    return {'path': name, 'sha256': sha256(path)}


def input_manifest(path, expected_ids):
    """Require an independently supplied manifest for exactly this dev run."""
    require(path is not None, 'a pre-established --input-manifest is required; '
            'do not generate historical H5 identities from the rerun itself')
    doc = json.loads(Path(path).read_text(encoding='utf-8'))
    rows = doc['inputs']
    by_id = {r['run_id']: r for r in rows}
    require(len(by_id) == len(rows), 'duplicate input run_id')
    require(set(by_id) == set(expected_ids), 'input manifest must cover exactly the expected dev run_ids')
    for row in rows:
        digest = row.get('h5_sha256', '')
        require(isinstance(digest, str) and len(digest) == 64
                and all(c in '0123456789abcdef' for c in digest), 'missing/invalid H5 SHA-256')
        source_dir = row.get('source_dir', '')
        require(bool(source_dir) and Path(source_dir).name == source_dir
                and '/' not in source_dir and '\\' not in source_dir
                and source_dir not in ('.', '..'), 'source_dir must be a single directory name')
    return by_id


def verify_input(path, expected):
    path = Path(path)
    require(path.name == expected['run_id'] + '.h5', 'H5 run_id/path mismatch')
    require(path.parent.name == expected['source_dir'], 'H5 recovery/source directory mismatch')
    actual = file_record(path)
    require(actual['sha256'] == expected['h5_sha256'], 'H5 SHA-256 mismatch: ' + str(path))
    return dict(expected, **actual)


def provenance(script, inputs, sources):
    import gprMax.toolboxes.SFCW.processing as processing
    import sfcw_official_loader_v0_2 as loader
    loader.verify_official_runtime()
    return {'schema': 'carrierfix-source-audit/1', 'python': platform.python_version(),
            'numpy': np.__version__, 'official_processing': file_record(processing.__file__),
            'script': file_record(script), 'loader': file_record(loader.__file__),
            'source_files': [file_record(p) for p in sources], 'h5_inputs': inputs}


def accept_t3(records, archived, keys):
    def index(rows):
        result = {(r['family'], r['damage'], r['level']): r for r in rows}
        require(len(result) == len(rows), 'duplicate t3 damage row')
        return result
    old = index(archived['records'])
    by_rep = {}
    for rep in ('legacy95', 'official20'):
        by_rep[rep] = index([r for r in records if r['representation'] == rep])
        require(set(by_rep[rep]) == set(old), 't3 missing/extra damage rows')
    require(len(records) == 2 * len(old), 'unknown t3 representation')
    for key, row in old.items():
        for metric in keys:
            require(by_rep['legacy95'][key][metric] == row[metric],
                    f't3 legacy/archive mismatch: {key} {metric}')
    return {'legacy_archive_equal': True, 'metrics_compared': len(old) * len(keys)}


def accept_b2(results, archived):
    old = {r['family']: r for r in archived['results']}
    require(len(old) == len(archived['results']), 'duplicate archived B2 family')
    require(len(results) == len(old) and {r['family'] for r in results} == set(old),
            'B2 missing/extra/duplicate families')
    rank_changes = []
    for result in results:
        fam = result['family']
        baseline = old[fam]
        legacy = result['representations']['legacy95']
        old_rows = {r['config']: r for r in baseline['rows']}
        new_rows = {r['config']: r for r in legacy['rows']}
        require(len(new_rows) == len(legacy['rows']) and len(old_rows) == len(baseline['rows'])
                and set(old_rows) == set(new_rows), 'B2 configuration coverage mismatch')
        for cid, row in new_rows.items():
            # preservation_db is absent from the new serialization; all fields
            # that it does report must reproduce the historical result exactly.
            require(set(row) == set(old_rows[cid]) - {'preservation_db'}, 'B2 metric field coverage mismatch')
            for name, value in row.items():
                require(value == old_rows[cid][name], f'B2 legacy/archive mismatch: {fam} {cid} {name}')
        for name in ('selection', 'selection_R_bg_db', 'ranking'):
            require(legacy[name] == baseline[name], f'B2 legacy/archive mismatch: {fam} {name}')
        official = result['representations']['official20']
        changed = legacy['ranking'] != official['ranking']
        require(result['ranking_changed'] == changed, 'B2 ranking flag mismatch')
        if changed:
            rank_changes.append(fam)
    return {'legacy_archive_equal': True, 'ranking_changed_families': rank_changes}


def accept_s1s3(rows, frozen, constants):
    require(len(rows) == len(frozen['t_measured_ns']) == len(frozen['t_fermat_er18_ns']),
            'S1/S3 trace coverage mismatch')
    require([r['trace'] for r in rows] == list(range(1, len(rows) + 1)), 'S1/S3 trace order mismatch')
    for i, row in enumerate(rows):
        require(row['t_measured_legacy_ns'] == frozen['t_measured_ns'][i]
                == row['t_measured_frozen_ns'], 'S1/S3 legacy/frozen peak mismatch')
    tf = np.asarray(frozen['t_fermat_er18_ns'], dtype=float)
    peaks = np.asarray([r['t_measured_official_ns'] for r in rows], dtype=float)
    half = constants['ev_half_ns']
    nc0 = constants['nc_window_ns'][0]
    delta = float(np.max(np.abs(peaks - tf)))
    require(np.isfinite(peaks).all() and delta < half, 'S1/S3 official peak outside Fermat window')
    require(float(tf.max()) + half < nc0, 'S1/S3 Fermat window overlaps NC')
    clearance = float(nc0 - peaks.max() - half)
    require(clearance > 0, 'S1/S3 measured window overlaps NC')
    return {'legacy_archive_equal': True, 'peak_window_and_nc_pass': True,
            'max_abs_official_minus_fermat_ns': delta,
            'measured_window_nc_clearance_ns': clearance}


def write_new_pair(path1, path2, doc1, doc2):
    """Only accepted, repeat-identical results can create new evidence files."""
    text1 = json.dumps(doc1, ensure_ascii=False, indent=1, allow_nan=False) + '\n'
    text2 = json.dumps(doc2, ensure_ascii=False, indent=1, allow_nan=False) + '\n'
    require(text1 == text2, 'r1/r2 mismatch')
    path1, path2 = Path(path1), Path(path2)
    require(path1.resolve() != path2.resolve(), 'r1/r2 paths must differ')
    require(not path1.exists() and not path2.exists(), 'refusing to overwrite archived evidence')
    path1.parent.mkdir(parents=True, exist_ok=True)
    path2.parent.mkdir(parents=True, exist_ok=True)
    with path1.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(text1)
    with path2.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(text2)
