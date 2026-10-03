import sys
try:
    import pycuda.driver as cuda
    cuda.init()
    n = cuda.Device.count()
    print(f"CUDA init OK, devices: {n}")
    for i in range(n):
        d = cuda.Device(i)
        print(f"  [{i}] {d.name()}")
    ctx = cuda.Device(0).make_context()
    free, total = cuda.mem_get_info()
    print(f"context OK, free={free/2**20:.0f} MiB total={total/2**20:.0f} MiB")
    ctx.pop()
    print("PROBE PASS")
except Exception as e:
    print(f"PROBE FAIL: {type(e).__name__}: {e}")
    sys.exit(1)
