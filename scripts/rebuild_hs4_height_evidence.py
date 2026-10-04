"""CPU-only fresh rebuild and byte comparison of the height attribution unit."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from hs_capsule_identity import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--full-local',action='store_true',help='Also requires local dense snapshots and duplicate geometries.')
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new rebuild directory required')
    a.out.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]; evidence=root/'artifacts/research_checks'
    def study(name): return evidence/('2026-10-04_'+name)
    jobs=[('patch','analyze_hs4_local_patch_controls.py',['--study',study('hs4_local_patch_controls')],study('hs4_local_patch_results')),
        ('refinement','analyze_hs4_patch_refinement.py',['--study',study('hs4_patch_refinement'),'--patch-results',study('hs4_local_patch_results')],study('hs4_patch_refinement_results')),
        ('finest','analyze_hs4_patch_finest.py',['--study',study('hs4_patch_finest'),'--refinement',study('hs4_patch_refinement'),'--patch-results',study('hs4_local_patch_results')],study('hs4_patch_finest_results'))]
    if a.full_local:
        jobs.extend([('fields','analyze_hs4_height_wavefield.py',['--study',study('hs4_height_wavefield_continuation')],study('hs4_height_wavefield_results')),
            ('visual','plot_hs4_wavefield_directions.py',['--study',study('hs4_height_wavefield_continuation'),'--results',a.out/'fields'],study('hs4_height_wavefield_visual'))])
    records=[]
    for name,script,arguments,original in jobs:
        fresh=a.out/name
        command=[sys.executable,str(root/'scripts'/script),*map(str,arguments),'--out',str(fresh)]
        completed=subprocess.run(command,cwd=root,capture_output=True,text=True)
        (a.out/(name+'.stdout.log')).write_text(completed.stdout,encoding='utf-8')
        (a.out/(name+'.stderr.log')).write_text(completed.stderr,encoding='utf-8')
        if completed.returncode: raise RuntimeError('rebuild failed:'+name)
        for file in sorted(original.rglob('*')):
            if file.is_file():
                rebuilt=fresh/file.relative_to(original)
                if not rebuilt.is_file() or sha256(file)!=sha256(rebuilt): raise ValueError('byte mismatch:'+str(file))
                records.append({'original':file.relative_to(root).as_posix(),'sha256':sha256(file)})
        print('Exact CPU rebuild:'+name,flush=True)
    if a.full_local:
        for script,args,original,target in [
            ('check_hs4_height_wavefield.py',['--study',study('hs4_height_wavefield_continuation')],study('hs4_height_wavefield_continuation')/'independent_verification.json',a.out/'independent_verification.json'),
            ('audit_hs4_height_native_archive.py',[],study('hs4_height_archive_verification')/'native_archive.json',a.out/'native_archive.json')]:
            subprocess.run([sys.executable,str(root/'scripts'/script),*map(str,args),'--out',str(target)],cwd=root,check=True)
            if sha256(original)!=sha256(target): raise ValueError('audit differs:'+script)
            records.append({'original':original.relative_to(root).as_posix(),'sha256':sha256(original)})
    report={'status':'PASS','calls_solver':False,'code_sha256':sha256(__file__),'full_local':a.full_local,
        'exactly_rebuilt_files':len(records),'files':records,
        'scope':'Exact execution reproducibility, not absolute FDTD convergence/physical acceptance. Full-local mode needs ignored dense snapshot/field and duplicate material files. Default rebuild needs archived native H5 only.'}
    (a.out/'rebuild_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS exact files:'+str(len(records)),flush=True)


if __name__=='__main__': main()
