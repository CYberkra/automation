"""Archive a private prepared model and publish only aggregate/hash evidence."""
import argparse
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from build_pdf_profile_geometry import digest, save_json

ROOT = Path(__file__).resolve().parents[1]


def archive(out, target, evidence):
    out, target, evidence = out.resolve(), target.resolve(), evidence.resolve()
    if target.exists() or evidence.exists() or (out/'reproduction_tools').exists():
        raise ValueError('Refuse to overwrite archive, evidence or source snapshots')
    if target.is_relative_to(out):
        raise ValueError('Archive must be outside the packaged directory')
    manifest = json.loads((out/'manifest.json').read_text('utf-8'))
    audit = json.loads((out/'independent_verification.json').read_text('utf-8'))
    capacity = json.loads((out/'resource_preflight.json').read_text('utf-8'))
    preview = json.loads((out/'preview_validation.json').read_text('utf-8'))
    assert digest(out/'manifest.json') == audit['manifest_sha256'] == capacity['manifest_sha256'] == preview['manifest_sha256']
    assert manifest['script_sha256'] == digest(ROOT/'scripts/prepare_line9_material_model.py')
    assert audit['script_sha256'] == digest(ROOT/'scripts/audit_line9_material_model.py')
    assert not manifest['calls_solver'] and not audit['calls_solver'] and not capacity['calls_solver']
    for filename, expected in manifest['source_inputs_sha256'].items():
        assert digest(Path(filename)) == expected, 'Original source changed'
    for g in manifest['geometries'].values():
        assert digest(out/'geometries'/g['file']) == g['sha256']
    for c in manifest['cases']:
        assert digest(out/c['input']) == c['input_sha256']
    scripts = ['prepare_line9_material_model.py', 'audit_line9_material_model.py',
               'preview_line9_material_model.py', 'view_line9_material_model_paraview.py',
               'build_pdf_profile_geometry.py', 'archive_line9_material_model.py']
    for name in scripts:
        ast.parse((ROOT/'scripts'/name).read_text('utf-8'))
    tools_dir = out/'reproduction_tools'
    tools_dir.mkdir()
    for name in scripts:
        shutil.copyfile(ROOT/'scripts'/name, tools_dir/name)
    files = sorted(p for p in out.rglob('*') if p.is_file())
    source_hashes = {p.relative_to(out).as_posix(): digest(p) for p in files}
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=1) as z:
        for p in files:
            z.write(p, arcname=p.relative_to(out).as_posix())
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        assert set(z.namelist()) == set(source_hashes)
        for name, expected in source_hashes.items():
            h = hashlib.sha256()
            with z.open(name) as stream:
                while block := stream.read(4*1024*1024):
                    h.update(block)
            assert h.hexdigest() == expected, 'Archive byte mismatch'
    duplicate_checks = {}
    v4python = Path('D:/gprmax_v4_gpu_env/Scripts/python.exe')
    checks = [([str(v4python), str(ROOT/'scripts/prepare_line9_material_model.py'),
                'prepare', '--source', str(ROOT/'artifacts/local_checks/2026-10-06_line9_pdf_model_r1'),
                '--out', str(out)], 'Refuse to overwrite model package'),
              ([sys.executable, str(ROOT/'scripts/audit_line9_material_model.py'),
                '--out', str(out)], 'Refuse to overwrite verification')]
    for command, expected in checks:
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode != 0 and expected in result.stderr
        duplicate_checks[expected] = True
    assert source_hashes == {p.relative_to(out).as_posix(): digest(p) for p in files}
    evidence.mkdir(parents=True)
    resource = {}
    for name, g in manifest['geometries'].items():
        resource[name] = {k: g[k] for k in ['shape_nxyz', 'grid_xyz_m', 'cells', 'file_bytes',
            'host_primary_array_lower_bytes', 'device_primary_array_lower_bytes']}
        resource[name].update(capacity['records'][name])
    summary = dict(status='PREPARED_VERIFIED_ARCHIVED_NOT_SIMULATED', calls_solver=False,
        calls_training=False, independent_verification=audit, resource_preflight=capacity,
        geometry_resources=resource, preview_validation=preview,
        installed_native_sources_sha256=manifest['installed_native_sources_sha256'],
        original_source_files_verified_unchanged=True, overwrite_refusal_checks=duplicate_checks,
        archive=dict(filename=target.name, sha256=digest(target), bytes=target.stat().st_size,
                     file_count=len(files), CRC_and_every_entry_SHA256_verified=True,
                     scope='Private local archive; contains PDF-derived geology; do not upload'),
        preview_files_sha256={p.name:digest(p) for p in out.iterdir()
                              if p.suffix in ['.png','.vti','.pvsm']},
        reproducibility_script_sha256={n:digest(ROOT/'scripts'/n) for n in scripts},
        public_scope='Aggregate counts, hardware requirements and identities only; native model and figures stay private',
        limitations=['No native FDTD results or B-scan', 'No mesh/boundary/physical validation',
                     'Target-machine capacity and execution contract must be checked separately'])
    save_json(evidence/'summary.json', summary)
    print(json.dumps(dict(status=summary['status'], archive_MiB=target.stat().st_size/2**20,
                          archived_files=len(files), all_entry_hashes_match=True,
                          source_files_unchanged=True), ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True)
    args = p.parse_args()
    archive(args.out, args.archive, args.evidence)
