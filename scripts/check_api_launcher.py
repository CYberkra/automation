"""Exercise reviewed API copy/execution and changed-input rejection; no solver."""
import os,tempfile
from pathlib import Path
from run_approved_airborne import solver_code

with tempfile.TemporaryDirectory() as td:
    source=Path(td)/'input.py'
    source.write_text("from pathlib import Path\nPath('marker').write_text(__name__)\n")
    code=solver_code(source,'probe',0);cwd=Path.cwd()
    try:
        os.chdir(td);exec(code)
        assert Path('marker').read_text()=='__main__'
        assert Path('probe.py').read_bytes()==source.read_bytes()
        source.write_text('raise RuntimeError("must not execute")')
        try:exec(code)
        except AssertionError:pass
        else:raise AssertionError('Changed model was not rejected')
    finally:os.chdir(cwd)
print('API launcher execution/copy/hash checks passed; no solver invoked')
