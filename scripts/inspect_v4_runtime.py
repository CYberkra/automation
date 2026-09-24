"""Import-only V4 identity check; no model parsing, allocation or solve.

Run with the isolated V4 Python, supplying a NEW output JSON path.
"""
import argparse
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise SystemExit('Refusing overwrite')
    import gprMax
    import numpy
    import Cython.Compiler.Code as cython_code
    from gprMax.toolboxes.SFCW.cli import build_parser
    assert gprMax.__version__=='4.0.0',gprMax.__version__
    assert (3,11)<=sys.version_info[:2]<(3,14)
    package_path=Path(gprMax.__file__).resolve()
    assert package_path.is_relative_to(Path(sys.prefix).resolve()),package_path
    modules={}
    for name in ('fields_updates_normal','geometry_primitives',
                 'pml_updates_electric_HORIPML','pml_updates_magnetic_HORIPML'):
        module=importlib.import_module('gprMax.cython.'+name)
        path=Path(module.__file__)
        modules[name]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    # Parser only; the named HDF5 need not exist and is never opened.
    args=build_parser().parse_args(['process','not_executed.h5','--receiver','name:measurement',
        '--component','Ex','--f-start','20e6','--f-stop','170e6','--steps','501','--method','direct'])
    assert args.steps==501 and args.method=='direct'
    result={'schema':'gprmax-runtime-identity/1','python':sys.version,'executable':sys.executable,
        'prefix':sys.prefix,'platform':platform.platform(),'gprmax_version':gprMax.__version__,
        'gprmax_module_path':str(package_path),'numpy_version':numpy.__version__,
        'cython_compiler_path':cython_code.__file__, 'compiled_modules_imported':modules,
        'official_sfcw_parser_verified':True,'solver_executed':False,'model_parsed':False,
        'packages':dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions())),
        'inspection_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'version':gprMax.__version__,'compiled_imports':len(modules),
                      'sfcw_parser_verified':True,'solver_executed':False}))


if __name__=='__main__':main()
