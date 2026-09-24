"""Run metadata audit, official SFCW CLI, and independent diagnostics; no FDTD.

Input must be a completed supervised M00 run. Output goes to a fresh directory.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from bounded_windows_process import supervise

ROOT = Path(__file__).resolve().parents[1]


def worker(raw, output):
    def command(label, args):
        cmd = [sys.executable, *map(str, args)]
        with (output/f'{label}.stdout.log').open('wb') as out, (output/f'{label}.stderr.log').open('wb') as err:
            result = subprocess.run(cmd, stdout=out, stderr=err, timeout=540,
                                    creationflags=0x08000000)
        records.append(dict(label=label, command=cmd, returncode=result.returncode))
        (output/'commands.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
        if result.returncode:
            raise RuntimeError(f'{label} failed; inspect preserved logs')
    records = []
    command('raw_audit', [ROOT/'scripts/audit_airborne_output.py', raw, '--output', output/'raw_audit.json'])
    audit = json.loads((output/'raw_audit.json').read_text(encoding='utf-8'))
    command('official_inspect', ['-m','gprMax.toolboxes.SFCW','inspect',raw])
    command('official_direct', ['-m','gprMax.toolboxes.SFCW','process',raw,
        '--source',audit['source_path'],'--receiver','name:measurement','--component','Ex',
        '--f-start','20e6','--f-stop','170e6','--steps','501','--method','direct',
        '--window','rectangular','--zero-pad','1','--time-shift','0','--tail-taper','0',
        '--source-floor-db','-100','--output',output/'official_direct.h5'])
    command('sfcw_audit', [ROOT/'scripts/audit_official_sfcw.py',output/'official_direct.h5',
        '--raw-audit',output/'raw_audit.json','--output-dir',output/'diagnostics'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory',type=Path)
    parser.add_argument('--output-directory',type=Path,required=True)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    run=args.run_directory.resolve()
    output=args.output_directory.resolve()
    status=json.loads((run/'supervision.json').read_text(encoding='utf-8'))
    if status['reason']!='completed' or status['exit_code']!=0:
        raise SystemExit('Refusing to process an incomplete/failed supervised run')
    raw=run/'M00_x_3d.h5'
    if not raw.is_file(): raise SystemExit('Missing expected raw HDF5')
    if args.worker:
        worker(raw,output)
    else:
        result=supervise([sys.executable,str(Path(__file__).resolve()),str(run),
            '--output-directory',str(output),'--worker'],output,
            wall_s=600,memory_bytes=4*1024**3,output_bytes=1024**3)
        print(json.dumps(result))
        if result['reason']!='completed': raise SystemExit(1)


if __name__=='__main__':main()
