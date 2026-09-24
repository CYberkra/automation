"""Windows Job Object supervisor. Stdlib only; contains no solver invocation.

Commit-memory cap applies to the job, not physical RSS. File budget is polled,
so bounded overshoot is possible between checks; do not call it a disk quota.
"""
import ctypes as C
from ctypes import wintypes as W
import json
import os
from pathlib import Path
import subprocess
import time


class IO(C.Structure):
    _fields_ = [(x,C.c_ulonglong) for x in ('read_ops','write_ops','other_ops','read_bytes','write_bytes','other_bytes')]


class Basic(C.Structure):
    _fields_ = [('process_time',C.c_longlong),('job_time',C.c_longlong),('flags',W.DWORD),
                ('min_ws',C.c_size_t),('max_ws',C.c_size_t),('active_limit',W.DWORD),
                ('affinity',C.c_size_t),('priority',W.DWORD),('scheduling',W.DWORD)]


class Extended(C.Structure):
    _fields_ = [('basic',Basic),('io',IO),('process_memory',C.c_size_t),
                ('job_memory',C.c_size_t),('peak_process_memory',C.c_size_t),('peak_job_memory',C.c_size_t)]


def supervise(command, directory, *, wall_s, memory_bytes, output_bytes, poll_s=.1):
    if os.name != 'nt':
        raise RuntimeError('This supervisor requires Windows')
    if min(wall_s,memory_bytes,output_bytes,poll_s)<=0:
        raise ValueError('Limits must be positive')
    directory=Path(directory).resolve()
    directory.mkdir(parents=True,exist_ok=False)
    k=C.WinDLL('kernel32',use_last_error=True)
    k.CreateJobObjectW.argtypes=[C.c_void_p,W.LPCWSTR]; k.CreateJobObjectW.restype=W.HANDLE
    k.SetInformationJobObject.argtypes=[W.HANDLE,C.c_int,C.c_void_p,W.DWORD]; k.SetInformationJobObject.restype=W.BOOL
    k.QueryInformationJobObject.argtypes=[W.HANDLE,C.c_int,C.c_void_p,W.DWORD,C.c_void_p]; k.QueryInformationJobObject.restype=W.BOOL
    k.AssignProcessToJobObject.argtypes=[W.HANDLE,W.HANDLE]; k.AssignProcessToJobObject.restype=W.BOOL
    k.TerminateJobObject.argtypes=[W.HANDLE,W.UINT]; k.TerminateJobObject.restype=W.BOOL
    k.CloseHandle.argtypes=[W.HANDLE]; k.CloseHandle.restype=W.BOOL
    job=k.CreateJobObjectW(None,None)
    if not job: raise C.WinError(C.get_last_error())
    proc=None
    started=time.monotonic()
    reason='not_started'
    peak=0
    size=0
    try:
        limits=Extended()
        limits.basic.flags=0x2000|0x200  # KILL_ON_JOB_CLOSE | JOB_MEMORY
        limits.job_memory=memory_bytes
        if not k.SetInformationJobObject(job,9,C.byref(limits),C.sizeof(limits)):
            raise C.WinError(C.get_last_error())
        with (directory/'stdout.log').open('wb') as out,(directory/'stderr.log').open('wb') as err:
            proc=subprocess.Popen(command,cwd=directory,stdout=out,stderr=err,
                                  creationflags=0x4|0x08000000)  # suspended, no window
            if not k.AssignProcessToJobObject(job,W.HANDLE(int(proc._handle))):
                raise C.WinError(C.get_last_error())
            nt=C.WinDLL('ntdll')
            nt.NtResumeProcess.argtypes=[W.HANDLE]; nt.NtResumeProcess.restype=C.c_long
            status=nt.NtResumeProcess(W.HANDLE(int(proc._handle)))
            if status!=0: raise RuntimeError(f'NtResumeProcess failed: {status}')
            reason='running'
            while True:
                usage=Extended()
                if not k.QueryInformationJobObject(job,9,C.byref(usage),C.sizeof(usage),None):
                    raise C.WinError(C.get_last_error())
                peak=max(peak,usage.peak_job_memory)
                size=sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())
                if size>output_bytes: reason='output_limit'; break
                if time.monotonic()-started>wall_s: reason='wall_limit'; break
                if proc.poll() is not None:
                    reason='completed' if proc.returncode==0 else 'nonzero_exit'
                    break
                time.sleep(poll_s)
            if reason!='completed':
                if not k.TerminateJobObject(job,124): raise C.WinError(C.get_last_error())
            proc.wait(timeout=10)
    except BaseException:
        # Includes failed assignment while the child is still suspended.
        if proc is not None and proc.poll() is None:
            k.TerminateJobObject(job,125)
            proc.kill()
            proc.wait(timeout=10)
        raise
    finally:
        k.CloseHandle(job)  # also kills any child surviving its root process
    result=dict(reason=reason,exit_code=proc.returncode,wall_s=time.monotonic()-started,
                peak_job_commit_bytes=peak,observed_output_bytes=size,
                memory_cap_bytes=memory_bytes,wall_cap_s=wall_s,output_cap_bytes=output_bytes,
                memory_semantics='Windows job committed memory, not RSS',
                output_semantics='polled file bytes; overshoot possible',
                command=command)
    (directory/'supervision.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
