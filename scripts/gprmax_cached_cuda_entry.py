"""Process-local CUBIN cache for identical compiler calls; no solver edits.

Every key contains verbatim CUDA source, resolved architecture and compiler
arguments. Cached bytes are loaded into each fresh context by PyCUDA normally.
This avoids repeated nvcc work; it does not reuse field arrays or GPU modules.
"""
import hashlib
import inspect
import json
import os
from pathlib import Path
import runpy
import sys

import pycuda.compiler as compiler
import pycuda.driver as driver

original_compile = compiler.compile
signature = inspect.signature(original_compile)
cache = {}
log = Path(os.environ['HS4_CUDA_CACHE_LOG'])


def cached_compile(*args, **kwargs):
    bound = signature.bind(*args, **kwargs)
    bound.apply_defaults()
    arguments = dict(bound.arguments)
    # SourceModule's default arch depends on the current device/context.
    if arguments['arch'] is None:
        cc = driver.Context.get_device().compute_capability()
        arguments['arch'] = f'sm_{cc[0]}{cc[1]}'
    if arguments['options'] is None:
        arguments['options'] = list(compiler.DEFAULT_NVCC_FLAGS)
    key = hashlib.sha256(json.dumps(arguments,sort_keys=True).encode('utf-8')).hexdigest()
    hit = key in cache
    if not hit:
        cache[key] = original_compile(**arguments)
    cubin = cache[key]
    with log.open('a',encoding='utf-8') as f:
        f.write(json.dumps({'key':key,'cache_hit':hit,'cubin_sha256':hashlib.sha256(cubin).hexdigest(),
                            'bytes':len(cubin),'arch':arguments['arch']})+'\n')
    return cubin


compiler.compile = cached_compile
sys.argv = ['gprMax', *sys.argv[1:]]
runpy.run_module('gprMax',run_name='__main__')
