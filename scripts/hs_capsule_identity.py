"""Read historical dict/current list manifests without rewriting capsules."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import re


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON key: {key}')
        result[key] = value
    return result


def read_manifest(base):
    manifest = json.loads((Path(base) / 'manifest.json').read_text(encoding='utf-8'),
                          object_pairs_hook=_unique_object)
    if isinstance(manifest, dict):
        records = [{'file': key, 'sha256': value} for key, value in manifest.items()]
    elif isinstance(manifest, list):
        records = manifest
    else:
        raise ValueError('manifest must be a filename/hash object or a record list')
    result = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError('manifest record must be an object')
        name, digest = record.get('file'), record.get('sha256')
        if (not isinstance(name, str) or not name or '\\' in name or ':' in name
                or PurePosixPath(name).is_absolute()
                or any(p in ('', '.', '..') for p in name.split('/'))):
            raise ValueError(f'unsafe manifest path: {name!r}')
        if name in result:
            raise ValueError(f'duplicate manifest file: {name}')
        if not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest):
            raise ValueError(f'invalid SHA-256 for {name}')
        size = record.get('bytes')
        if size is not None and (type(size) is not int or size < 0):
            raise ValueError(f'invalid byte count for {name}')
        result[name] = dict(record)
    return result


def verify_file(base, records, name):
    if name not in records:
        raise ValueError(f'input absent from manifest: {name}')
    root = Path(base).resolve()
    path = (root / name).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise ValueError(f'missing/unsafe input: {name}')
    record = records[name]
    actual = sha256(path)
    if actual != record['sha256'] or ('bytes' in record and path.stat().st_size != record['bytes']):
        raise ValueError(f'input identity mismatch: {name}')
    return {'file': name, 'sha256': actual, 'bytes': path.stat().st_size}
