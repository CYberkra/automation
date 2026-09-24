"""Compile and execute a tiny double-precision CUDA kernel; no FDTD."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise SystemExit('Refusing overwrite')
    import numpy as np
    import pycuda.driver as cuda
    from pycuda.compiler import SourceModule
    import pycuda._driver as native
    cuda.init()
    device=cuda.Device(0)
    context=device.make_context()
    try:
        free,total=cuda.mem_get_info()
        values=np.arange(4096,dtype=np.float64)
        allocation=cuda.mem_alloc(values.nbytes)
        cuda.memcpy_htod(allocation,values)
        module=SourceModule('''__global__ void check_double(double *x) {
            const int i=blockIdx.x*blockDim.x+threadIdx.x;
            if(i<4096) x[i]=x[i]*x[i]+1.0;
        }''')
        kernel=module.get_function('check_double')
        start,end=cuda.Event(),cuda.Event()
        start.record();kernel(allocation,block=(128,1,1),grid=(32,1));end.record();end.synchronize()
        elapsed=start.time_till(end)
        output=np.empty_like(values);cuda.memcpy_dtoh(output,allocation)
        error=float(np.max(abs(output-(values*values+1))))
        if error!=0: raise ValueError('CUDA double kernel check failed')
        allocation.free()
        result=dict(executable=sys.executable,device_id=0,device_name=device.name(),
            pci_bus_id=device.pci_bus_id(),compute_capability=device.compute_capability(),
            total_device_bytes=total,available_device_bytes=free,
            driver_version=cuda.get_driver_version(),compiled_cuda_version=cuda.get_version(),
            pycuda_version=importlib.metadata.version('pycuda'),
            native_module_path=native.__file__,native_module_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
            kernel_float64=True,kernel_elements=4096,max_absolute_error=error,
            synchronized_kernel_ms=elapsed,kernel_time_is_not_fdtd_benchmark=True,
            solver_executed=False,passed=True)
    finally:
        context.pop();context.detach()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':main()
