"""Sample device-wide free/total VRAM to a JSONL file until terminated.

Read-only monitoring via a dedicated CUDA context; performs no FDTD and no
kernel compilation. The context itself occupies a small amount of VRAM, so
samples are recorded as device-wide values, not per-process quotas.
"""
import json
import sys
import time
from pathlib import Path


def main():
    output = Path(sys.argv[1])
    interval = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
    import pycuda.driver as cuda
    cuda.init()
    context = cuda.Device(0).make_context()
    started = time.monotonic()
    try:
        with output.open('x', encoding='utf-8', newline='\n') as f:
            while True:
                free, total = cuda.mem_get_info()
                f.write(json.dumps({'t_s': time.monotonic() - started,
                                    'free_bytes': free, 'total_bytes': total}) + '\n')
                f.flush()
                time.sleep(interval)
    finally:
        context.pop()
        context.detach()


if __name__ == '__main__':
    main()
