"""Exercise the job supervisor with harmless Python children; never calls FDTD."""
import argparse
import ctypes as C
from ctypes import wintypes as W
import hashlib
import json
from pathlib import Path
import sys

from bounded_windows_process import supervise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    cases={
        'normal':("print('ok')",5,128*2**20,2**20),
        'nonzero':('raise SystemExit(7)',5,128*2**20,2**20),
        'wall_child':("import subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)']); print(p.pid,flush=True); time.sleep(60)",1,256*2**20,2**20),
        'memory':("try:\n x=bytearray(256*2**20)\nexcept MemoryError:\n raise SystemExit(86)",5,96*2**20,2**20),
        'output':("from pathlib import Path; import time; Path('payload.bin').write_bytes(b'x'*2**20); time.sleep(5)",5,128*2**20,65536),
    }
    results={name:supervise([sys.executable,'-c',code],a.output/name,
                            wall_s=wall,memory_bytes=memory,output_bytes=output)
             for name,(code,wall,memory,output) in cases.items()}
    pid=int((a.output/'wall_child/stdout.log').read_text().strip())
    k=C.WinDLL('kernel32',use_last_error=True)
    k.OpenProcess.argtypes=[W.DWORD,W.BOOL,W.DWORD]; k.OpenProcess.restype=W.HANDLE
    k.WaitForSingleObject.argtypes=[W.HANDLE,W.DWORD]; k.WaitForSingleObject.restype=W.DWORD
    k.CloseHandle.argtypes=[W.HANDLE]
    handle=k.OpenProcess(0x100000,False,pid)
    if handle:
        gone=k.WaitForSingleObject(handle,2000)==0
        k.CloseHandle(handle)
    else:
        gone=C.get_last_error()==87  # no such PID; access denied is not success
    checks={'normal':results['normal']['reason']=='completed',
            'nonzero':results['nonzero']['exit_code']==7,
            'timeout':results['wall_child']['reason']=='wall_limit',
            'descendant_terminated':gone,
            'memory_allocation_denied':results['memory']['exit_code']==86,
            'output_stopped':results['output']['reason']=='output_limit'}
    result={'solver_executed':False,'checks':checks,'passed':all(checks.values()),
            'supervisor_sha256':hashlib.sha256(Path(__file__).with_name('bounded_windows_process.py').read_bytes()).hexdigest(),
            'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'results':results}
    (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':result['passed'],'checks':checks}))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
