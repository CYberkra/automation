"""Check the latest demotion and reproduce a clean-checkout verification failure.

No solver, no historical evaluator, no test-family data. All subprocess output
and the default report are written to a fresh ignored local_checks directory.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/research/2026-10-02_large_legacy_demotion_manifest.json'
ADAPTER = Path('artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    work = Path(tempfile.mkdtemp(prefix='demotion_review_', dir=ROOT / 'artifacts/local_checks'))
    demoted = json.loads(MANIFEST.read_text(encoding='utf-8'))
    tracked = set(subprocess.check_output(
        ['git', 'ls-tree', '-r', '--name-only', 'cbe203e'], cwd=ROOT).decode('utf-8').splitlines())
    local_matches = sum((ROOT / name).is_file() and sha(ROOT / name) == rec['sha256']
                        for name, rec in demoted.items())
    adapter_manifest = json.loads((ROOT / ADAPTER / 'manifest.json').read_text(encoding='utf-8'))
    contract = json.loads((ROOT / 'configs/research/sfcw_c3_pair_diagnostic_v1.json').read_text(encoding='utf-8'))
    array = adapter_manifest['array_archive']
    mirror = work / 'clean_checkout'
    dest = mirror / ADAPTER
    dest.mkdir(parents=True)
    (mirror / 'START_HERE.md').write_text('Review-only empty root marker\n', encoding='utf-8')
    for name in ('root_verify.py', 'manifest.json'):
        shutil.copy2(ROOT / ADAPTER / name, dest / name)
    run = subprocess.run([sys.executable, str(dest / 'root_verify.py')],
                         cwd=mirror, capture_output=True, text=True, encoding='utf-8')
    (work / 'stderr.txt').write_text(run.stderr, encoding='utf-8')
    if run.returncode == 0 or 'FileNotFoundError' not in run.stderr or array['path'] not in run.stderr:
        raise RuntimeError('Expected clean-checkout missing-array failure did not occur: ' + run.stderr)
    report = {
        'schema': 'storage_demotion_review/1', 'reviewed_commit': 'cbe203e',
        'manifest_sha256': sha(MANIFEST), 'demoted_count': len(demoted),
        'demoted_files_absent_from_current_git_tree': sum(n not in tracked for n in demoted),
        'preserved_local_files_matching_original_sha256': local_matches,
        'frozen_contract_manifest_sha_matches': contract['input_manifest_sha256'] == sha(ROOT / ADAPTER / 'manifest.json'),
        'transitively_required_array': str(ADAPTER / array['path']).replace('\\', '/'),
        'array_manifest_sha_matches_demotion': array['sha256'] == demoted[(ADAPTER / array['path']).as_posix()]['sha256'],
        'clean_checkout_verifier_exit_code': run.returncode,
        'clean_checkout_error': run.stderr.splitlines()[-1],
        'recovery': 'Git history still retains original blobs; restore them explicitly, without rerunning scientific evaluation.',
        'solver_invoked': False, 'test_family_read': False,
        'source_sha256': sha(Path(__file__)),
    }
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else work / 'review.json'
    with output.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps({'report': str(output), 'demoted_count': len(demoted),
                      'clean_checkout_exit_code': run.returncode}, indent=2))


if __name__ == '__main__':
    main()
