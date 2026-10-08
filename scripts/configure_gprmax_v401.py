"""Install an isolated 4.0.1 solver, reusing a verified read-only dependency tree.

The old environment is never changed. The new environment's site-packages takes
precedence over the fallback dependency directory; audit the imported solver path.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import venv


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--venv', type=Path, required=True)
    p.add_argument('--wheel', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--dependency-site', type=Path,
                   help='Optional ASCII alias for an existing dependency directory (NVCC on Windows).')
    a = p.parse_args()
    assert not a.venv.exists(), 'Fresh environment required; preserve prior attempts'
    # Preserve an ASCII junction path: resolving it can reintroduce Unicode into
    # PyCUDA's NVCC include argument even though the files themselves are valid.
    old_site = (a.dependency_site or Path(sys.prefix)/'Lib/site-packages').absolute()
    old_package = old_site / 'gprMax'
    before = {str(x.relative_to(old_package)): digest(x) for x in old_package.rglob('*')
              if x.is_file() and x.suffix in ('.py', '.pyd', '.tmpl')}
    venv.EnvBuilder(with_pip=True).create(a.venv)
    site = a.venv / 'Lib/site-packages'
    (site / 'verified_cuda_dependencies.pth').write_text(str(old_site) + '\n', encoding='utf-8')
    python = a.venv / 'Scripts/python.exe'
    subprocess.run([str(python), '-m', 'pip', 'install', '--no-deps', '--ignore-installed',
                    str(a.wheel.resolve())], check=True)
    probe = '''import gprMax,sys,json,importlib.metadata as m
from pathlib import Path
from gprMax.cython import fields_updates_dispersive, fields_updates_normal, virtual_waveguide
import pycuda.driver as cuda
cuda.init()
print(json.dumps(dict(version=gprMax.__version__,package=str(Path(gprMax.__file__).parent),
python=sys.executable,dependencies={k:m.version(k) for k in ['numpy','h5py','pycuda','scipy','psutil']},
gpu=[dict(name=cuda.Device(i).name(),capability=cuda.Device(i).compute_capability()) for i in range(cuda.Device.count())])))
'''
    r = subprocess.run([str(python), '-c', probe], check=True, capture_output=True, text=True)
    audit = json.loads(r.stdout.strip().splitlines()[-1])
    assert audit['version'] == '4.0.1'
    assert Path(audit['package']).resolve() == (site / 'gprMax').resolve()
    after = {str(x.relative_to(old_package)): digest(x) for x in old_package.rglob('*')
             if x.is_file() and x.suffix in ('.py', '.pyd', '.tmpl')}
    assert before == after, 'Old solver changed'
    audit.update(status='PASS_INSTALL_IMPORT_GPU_ENUMERATION_NOT_SOLVER_VALIDATION',
                 wheel_sha256=digest(a.wheel), wheel=a.wheel.name,
                 shared_dependency_site=str(old_site), old_solver_unchanged=True,
                 source_identities={str(x.relative_to(site / 'gprMax')): digest(x)
                                    for x in (site / 'gprMax').rglob('*')
                                    if x.is_file() and x.suffix in ('.py', '.pyd', '.tmpl')})
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k != 'source_identities'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
