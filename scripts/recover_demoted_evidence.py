"""Recover demoted large legacy evidence from pinned git-history blobs.

Context: commit cbe203e demoted 33 large legacy files (459.1 MB) from git
tracking to T2-local-in-place (see
docs/research/2026-10-02_large_legacy_demotion_manifest.json). The blobs
remain in git history at the pinned commit below. A fresh clone that needs
these files (e.g. to re-run a frozen-contract verification entrypoint such as
the sfcw_c3_pair_diagnostic_v1 chain) can restore them with byte identity
verified against the manifest SHA-256.

This restores *file content* from history; it does NOT re-run old evaluation
procedures, and it does not alter the archive freeze (restored bytes are
hash-identical to the archived bytes).

Usage (repo root, any Python):
  python scripts/recover_demoted_evidence.py            # verify + restore missing
  python scripts/recover_demoted_evidence.py --check    # verify only, no writes
Exit 0 = all present & hash-verified (or restored & verified); 1 = failure.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/research/2026-10-02_large_legacy_demotion_manifest.json'
PINNED_COMMIT = 'c921881'  # last commit in which all 33 files were tracked


def git_blob_bytes(path):
    r = subprocess.run(['git', 'show', f'{PINNED_COMMIT}:{path}'],
                       capture_output=True, cwd=ROOT)
    if r.returncode != 0:
        raise RuntimeError(f'git show failed for {path}: {r.stderr.decode(errors="replace")[:200]}')
    return r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help='verify only; never write')
    args = ap.parse_args()
    man = json.loads(MANIFEST.read_text(encoding='utf-8'))
    ok = restored = 0
    failures = []
    for rel, meta in man.items():
        p = ROOT / rel
        want = meta['sha256']
        if p.exists():
            got = hashlib.sha256(p.read_bytes()).hexdigest()
            if got == want:
                ok += 1
                continue
            failures.append(f'hash mismatch (existing file differs from archived bytes): {rel}')
            continue
        if args.check:
            failures.append(f'missing (check mode, not restored): {rel}')
            continue
        blob = git_blob_bytes(rel)
        got = hashlib.sha256(blob).hexdigest()
        if got != want:
            failures.append(f'history blob hash mismatch: {rel} got={got[:16]}… want={want[:16]}…')
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(blob)
        restored += 1
    print(f'manifest entries: {len(man)} | present+verified: {ok} | restored from history+verified: {restored} | failures: {len(failures)}')
    for f in failures:
        print('FAIL:', f)
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
