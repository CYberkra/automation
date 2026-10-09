"""Verify preflight errors cancel the owned pipeline before any SSH or solver."""
import argparse
import json
from pathlib import Path

import run_line9_bounded_continuation_pipeline as pipeline


def main(out):
    assert not out.exists()
    out.mkdir(parents=True)
    connection = out / 'fixture_connection.json'
    connection.write_text(json.dumps(dict(host='fixture.invalid', options=[])))
    checks = []
    original = pipeline.subprocess.run
    calls = []
    def forbidden(*args, **kwargs):
        calls.append(args)
        raise AssertionError('No subprocess permitted in preflight rejection checks')
    pipeline.subprocess.run = forbidden
    try:
        for name in ['changed_source', 'cancelled']:
            cancel = out / (name + '_USER_STOP')
            if name == 'cancelled':
                cancel.write_text('USER_STOP')
            config = out / (name + '.json')
            config.write_text(json.dumps(dict(work=str(out / name), connection=str(connection),
                cancel=str(cancel), source_identities={str(Path(pipeline.__file__)): '0' * 64},
                continuation_root='E:/fixture/unreachable')))
            try:
                pipeline.main(config)
            except AssertionError:
                pass
            else:
                raise AssertionError('Invalid preflight accepted')
            events = [json.loads(s) for s in (out / name / 'pipeline.jsonl').read_text('utf-8').splitlines()]
            assert events[-1]['status'] == 'PIPELINE_FAILED_PRESERVE_ATTEMPTS_NO_RETRY'
            assert cancel.exists() and not calls
            checks.append(name + ':failed/cancelled before any SSH/subprocess')
    finally:
        pipeline.subprocess.run = original
    (out / 'checks.json').write_text(json.dumps(dict(status='PASS_PIPELINE_PREFLIGHT_FAILSTOP_NO_SSH',
        checks=checks, fixture=True, calls_solver=False, calls_SSH=False,
        pipeline_sha256=pipeline.sha(pipeline.__file__), checker_sha256=pipeline.sha(__file__),
        limits='Only source-mismatch and pre-existing cancellation preflight;full future pipeline not yet executed.'), indent=2) + '\n')
    print(json.dumps(dict(status='PASS', checks=len(checks))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args().out.resolve())
