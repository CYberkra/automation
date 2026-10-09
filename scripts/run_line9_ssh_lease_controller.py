"""Bounded local SSH owner renews a specific frozen remote lease; no credentials saved."""
import argparse
import base64
import json
from pathlib import Path
import subprocess
import threading
import time


def probe_script(root, digest, stop=False):
    # Literal PS strings; never interpolate arbitrary shell syntax.
    root = root.replace("'", "''")
    assert len(digest) == 64 and all(c in '0123456789abcdef' for c in digest)
    action = "Set-Content -LiteralPath ($r+'\\USER_STOP') -Value 'Bounded local controller cancellation'" if stop else "Set-Content -LiteralPath ($r+'\\session_heartbeat') -Value 'Bounded local SSH owner active'"
    return f"$ErrorActionPreference='Stop';$ProgressPreference='SilentlyContinue';$r='{root}';if((Get-FileHash -LiteralPath ($r+'\\execution_contract.json') -Algorithm SHA256).Hash.ToLower() -ne '{digest}'){{throw 'Contract changed'}};{action};$e=@(Get-Content -LiteralPath ($r+'\\execution.jsonl') -ErrorAction SilentlyContinue | ForEach-Object{{ConvertFrom-Json $_}});[pscustomobject]@{{completed=@($e|Where-Object{{$_.status -eq 'COMPLETED' -and $_.group}}).Count;last=$e[-1].status}}|ConvertTo-Json -Compress"


def main(a):
    assert 60 <= a.max_seconds <= 7200 and not a.log.exists()
    config = json.loads(a.connection.read_text('utf-8'))
    ssh = ['ssh', *config['options'], config['host']]
    end = time.monotonic() + a.max_seconds
    halted = threading.Event()
    failures = []

    def event(value):
        value['unix_s'] = time.time()
        with a.log.open('a', encoding='utf-8') as f:
            f.write(json.dumps(value) + '\n')

    def probe(stop=False):
        code = probe_script(a.remote_execution, a.contract_sha256, stop)
        encoded = base64.b64encode(code.encode('utf-16-le')).decode()
        r = subprocess.run([*ssh, 'powershell', '-NoProfile', '-EncodedCommand', encoded], capture_output=True, text=True, timeout=35)
        assert r.returncode == 0, f'SSH status {r.returncode}'
        return r.stdout.strip()

    def keeper():
        while not halted.wait(45):
            try:
                if a.cancel.exists() or time.monotonic() >= end:
                    event(dict(status='CANCEL_SENT', remote=probe(True)))
                    failures.append('controller deadline or USER_STOP')
                    return
                event(dict(status='LEASE_RENEWED', remote=probe()))
            except Exception as exc:
                event(dict(status='LEASE_CONNECTION_FAILED_NO_FURTHER_RENEWAL', error=str(exc)))
                failures.append(str(exc))
                return

    assert not a.cancel.exists()
    event(dict(status='CONTROLLER_STARTED', remote_execution=a.remote_execution,
        contract_sha256=a.contract_sha256, maximum_seconds=a.max_seconds,
        lease_authority='This bounded local process,not desktop sampling. Loss of this process/SSH stops renewal;remote600s expiry remains.', remote=probe()))
    thread = threading.Thread(target=keeper, daemon=True)
    thread.start()
    process = subprocess.Popen([*ssh, 'cmd', '/c', a.remote_run_cmd])
    try:
        code = process.wait(timeout=a.max_seconds + 660)
        assert code == 0 and not failures, (code, failures)
        event(dict(status='CONTROLLER_FINISHED_SSH_EXIT0_VERIFY_NATIVE_SEPARATELY'))
    except BaseException:
        try:
            probe(True)
        except Exception:
            pass
        if process.poll() is None:
            process.terminate()  # Only our local SSH;remote owned solver stops via USER_STOP/lease.
        event(dict(status='CONTROLLER_FAILED_OR_CANCELLED_NATIVE_STATE_MUST_BE_READ'))
        raise
    finally:
        halted.set()
        thread.join(timeout=40)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['connection', 'log', 'cancel']:
        p.add_argument('--' + key, type=Path, required=True)
    for key in ['remote-execution', 'contract-sha256', 'remote-run-cmd']:
        p.add_argument('--' + key, required=True)
    p.add_argument('--max-seconds', type=int, default=7200)
    main(p.parse_args())
