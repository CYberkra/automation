"""Run the approved Line9 V3.1 scripts once on CPU/single.

The supplied archive stays outside Git. This runner adds receipts and bounded
process supervision. The separately approved 1100 ns task changes only the
copied time-window token; solver calls and SFCW settings remain original.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def run(root):
    import gprMax
    import h5py
    import numpy as np
    import psutil

    contract = json.loads((root / "contract.json").read_text(encoding="utf-8"))
    assert contract["authorization"] in ("严格原包CPU单精度试跑", "批准1100ns单道对照")
    assert gprMax.__version__ == "4.0.1"
    assert sha(root / "original.zip") == contract["archive_sha256"]
    attempt = root / "attempt.json"
    with attempt.open("x", encoding="utf-8") as f:
        json.dump({"started_unix": time.time(), "pid": os.getpid(), "no_retry": True}, f)
    package = root / "package"
    assert not package.exists()
    package.mkdir()
    with zipfile.ZipFile(root / "original.zip") as z:
        for info in z.infolist():
            rel = Path(info.filename)
            assert not rel.is_absolute() and ".." not in rel.parts
            target = package / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(z.read(info))
    for name, expected in contract["original_file_sha256"].items():
        assert sha(package / name) == expected, name
    expected_files = dict(contract["original_file_sha256"])
    if contract["authorization"] == "批准1100ns单道对照":
        assert contract["record_ns"] == 1100
        config_path = package / "config/pilot_ascan.json"
        original = config_path.read_bytes()
        token = b'"time_window_ns": 800.0'
        assert original.count(token) == 1
        changed = original.replace(token, b'"time_window_ns": 1100.0', 1)
        before, after = json.loads(original), json.loads(changed)
        after["time_window_ns"] = before["time_window_ns"]
        assert before == after, "Only the approved time window may change"
        assert hashlib.sha256(changed).hexdigest() == contract["modified_pilot_sha256"]
        config_path.write_bytes(changed)
        expected_files["config/pilot_ascan.json"] = contract["modified_pilot_sha256"]
    logs = root / "logs"
    logs.mkdir()
    batch_start = time.monotonic()
    receipts = []
    environment = dict(os.environ, PYTHONUTF8="1", PYTHONUNBUFFERED="1")

    def execute(label, args):
        command = [sys.executable, *args]
        start = time.monotonic()
        peak_rss = 0
        reason = None
        with (logs / f"{label}.txt").open("wb") as log:
            child = subprocess.Popen(command, cwd=package, env=environment, stdout=log, stderr=subprocess.STDOUT)
            try:
                while child.poll() is None:
                    if (root / "USER_STOP").exists():
                        reason = "USER_STOP"
                    elif time.monotonic() - batch_start > contract["max_batch_seconds"]:
                        reason = "batch_time_limit"
                    elif time.time() - (root / "heartbeat").stat().st_mtime > 600:
                        reason = "session_lease_expired"
                    try:
                        processes = [psutil.Process(child.pid)]
                        processes += processes[0].children(recursive=True)
                        rss = sum(p.memory_info().rss for p in processes if p.is_running())
                        peak_rss = max(peak_rss, rss)
                        if rss > 8 * 2**30:
                            reason = "owned_memory_limit"
                    except psutil.NoSuchProcess:
                        pass
                    if psutil.virtual_memory().available < 1 * 2**30:
                        reason = "system_memory_limit"
                    save(root / "progress.json", {"stage": label, "seconds": time.monotonic() - start,
                                                 "peak_rss_bytes": peak_rss, "stop_reason": reason})
                    if reason:
                        raise RuntimeError(reason)
                    time.sleep(2)
                code = child.wait()
            except BaseException:
                if child.poll() is None:
                    try:
                        owned = psutil.Process(child.pid)
                        descendants = owned.children(recursive=True)
                        for p in reversed(descendants):
                            p.kill()
                        owned.kill()
                    except psutil.NoSuchProcess:
                        pass
                    child.wait()
                raise
        receipt = {"stage": label, "command": command, "exit_code": code,
                   "wall_seconds": time.monotonic() - start, "peak_rss_bytes": peak_rss,
                   "log_sha256": sha(logs / f"{label}.txt")}
        receipts.append(receipt)
        save(root / "receipts.json", receipts)
        if code != 0:
            raise RuntimeError(f"{label} exited {code}; no retry")

    try:
        execute("01_build", ["01_build_geometry.py"])
        execute("02_preflight", ["02_run_ascan.py"])
        execute("03_solve_cpu_single", ["02_run_ascan.py", "--solve"])
        raw = package / "results" / "line9_pilot_impulse.h5"
        assert raw.is_file(), "Expected original output filename absent; do not guess"
        with h5py.File(raw, "r") as h:
            rx = h["rxs/rx1/Ez"][:]
            assert rx.dtype == np.dtype("float32")
            assert rx.ndim == 1 and len(rx) == int(h.attrs["Iterations"])
            assert np.isfinite(rx).all() and np.any(rx != 0)
            assert np.isclose(float(h.attrs["dt"]), contract["dt_s"], rtol=1e-12)
            tail = float(np.max(np.abs(rx[-int(np.ceil(len(rx) * .05)):])) / np.max(np.abs(rx)))
        execute("04_sfcw_inspect", ["03_run_sfcw.py", "--raw", str(raw), "--inspect"])
        execute("05_sfcw_process", ["03_run_sfcw.py", "--raw", str(raw), "--process"])
        for name, expected in expected_files.items():
            assert sha(package / name) == expected, f"Original input changed: {name}"
        outputs = [{"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha(p)}
                   for p in sorted(package.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
        save(root / "completion.json", {"status": "COMPLETE_CPU_FLOAT32_NOT_CONVERGENCE_CERTIFIED",
                                       "authorization": contract["authorization"], "record_ns": contract["record_ns"],
                                       "version": gprMax.__version__, "python": sys.executable,
                                       "gprMax_path": gprMax.__file__, "gprMax_init_sha256": sha(Path(gprMax.__file__)),
                                       "native_dtype": "float32", "native_samples": len(rx),
                                       "tail_last5percent_db": float(20 * np.log10(max(tail, np.finfo(float).tiny))),
                                       "receipts": receipts, "files": outputs,
                                       "batch_wall_seconds": time.monotonic() - batch_start})
        save(root / "progress.json", {"stage": "COMPLETE"})
    except BaseException as e:
        save(root / "failure.json", {"error": repr(e), "receipts": receipts, "no_retry": True})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    run(parser.parse_args().root.resolve())
