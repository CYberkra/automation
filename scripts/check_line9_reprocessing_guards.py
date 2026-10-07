"""Read-only negative checks and post-processing native identity verification."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

from reprocess_line9_result_packages import main as reprocess, process_package
from review_line9_result_packages import sha, save


def main(root, review, out):
    if out.exists():
        raise ValueError('Fresh guard-report directory required')
    checks = {}
    def reject(name, operation):
        try:
            operation()
        except ValueError as error:
            checks[name] = str(error)
        else:
            raise AssertionError(name+' did not reject')
    never_created = out/'must_not_be_created'
    reject('existing_input_directory_as_output', lambda: reprocess(root, review, root, never_created))
    reject('same_output_alias', lambda: reprocess(root, review, never_created, never_created))
    audits = [json.loads(p.read_text('utf-8')) for p in sorted(review.glob('*_audit.json'))]
    if len(audits) != 3:
        raise ValueError('Expected three audits')
    a = audits[0]
    b = deepcopy(a)
    b['records'][0]['native_sha256'] = '0'*64
    reject('changed_native_hash', lambda: process_package(root/a['package'], b))
    b = deepcopy(a)
    b['materials_sha256'] = '0'*64
    reject('changed_material_hash', lambda: process_package(root/a['package'], b))
    assert not never_created.exists()
    total = 0
    for a in audits:
        assert len(a['records']) == a['raw_count']
        for r in a['records']:
            assert sha(root/a['package']/'cases'/r['id']/'profile.h5') == r['native_sha256']
            total += 1
    assert total == sum(a['raw_count'] for a in audits)
    out.mkdir(parents=True)
    save(out/'guard_checks.json', dict(calls_solver=False, calls_training=False,
        script_sha256=sha(__file__), processing_script_sha256=sha(Path(__file__).with_name('reprocess_line9_result_packages.py')),
        rejections=checks, native_hashes_unchanged_after_processing=total))
    print(f'{len(checks)} rejection checks passed; {total} native hashes unchanged')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.review, a.out)
