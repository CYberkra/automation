"""Fresh CPU rebuild: exact products, explicitly bounded scalar reductions."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from hs_capsule_identity import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--full-local',action='store_true')
    p.add_argument('--reuse-fresh',type=Path,help='Audit existing independent fresh CPU products, without overwriting them.')
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new verification directory required')
    a.out.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]; evidence=root/'artifacts/research_checks'
    def study(name): return evidence/('2026-10-04_'+name)
    freshroot=a.reuse_fresh or a.out
    jobs=[('patch','analyze_hs4_local_patch_controls.py',['--study',study('hs4_local_patch_controls')],study('hs4_local_patch_results')),
        ('refinement','analyze_hs4_patch_refinement.py',['--study',study('hs4_patch_refinement'),'--patch-results',study('hs4_local_patch_results')],study('hs4_patch_refinement_results')),
        ('finest','analyze_hs4_patch_finest.py',['--study',study('hs4_patch_finest'),'--refinement',study('hs4_patch_refinement'),'--patch-results',study('hs4_local_patch_results')],study('hs4_patch_finest_results'))]
    if a.full_local:
        jobs.extend([('fields','analyze_hs4_height_wavefield.py',['--study',study('hs4_height_wavefield_continuation')],study('hs4_height_wavefield_results')),
            ('visual','plot_hs4_wavefield_directions.py',['--study',study('hs4_height_wavefield_continuation'),'--results',freshroot/'fields'],study('hs4_height_wavefield_visual'))])
    jobs.append(('relation','plot_hs4_model_response_relation.py',[],study('hs4_model_response_relation')))
    records=[]; reductions=[]
    # Two complex norms plus division: conservative gamma(8*N) bound on the
    # positive square-sum reductions. All fields/products are still exact.
    n=501*36; u=np.finfo(np.float64).eps; gamma=(8*n*u)/(1-8*n*u)
    allowed={f'/field_sampling_checks/{group}/decimation/{hop}/hann_complex_relative_L2'
        for group in ('high_rough','high_halfspace','low_rough','low_halfspace') for hop in ('5','10')}
    def compare_json(x,y,path=''):
        if isinstance(x,dict):
            if x.keys()!=y.keys(): raise ValueError('JSON keys differ')
            for k in x: compare_json(x[k],y[k],path+'/'+k)
        elif isinstance(x,list):
            if len(x)!=len(y): raise ValueError('JSON length differs')
            for i,(xx,yy) in enumerate(zip(x,y)): compare_json(xx,yy,path+'/'+str(i))
        elif x!=y:
            if path not in allowed or not isinstance(x,float) or not isinstance(y,float) or not np.isfinite([x,y]).all():
                raise ValueError('unapproved scalar difference:'+path)
            bound=gamma*max(abs(x),abs(y))
            if abs(x-y)>bound: raise ValueError('norm reduction roundoff exceeded')
            reductions.append({'path':path,'original':x,'rebuilt':y,'absolute_difference':abs(x-y),'derived_reduction_bound':bound})
    for name,script,arguments,original in jobs:
        fresh=freshroot/name
        if not (a.reuse_fresh and fresh.is_dir()):
            command=[sys.executable,str(root/'scripts'/script),*map(str,arguments),'--out',str(fresh)]
            completed=subprocess.run(command,cwd=root,capture_output=True,text=True)
            (a.out/(name+'.stdout.log')).write_text(completed.stdout,encoding='utf-8')
            (a.out/(name+'.stderr.log')).write_text(completed.stderr,encoding='utf-8')
            if completed.returncode: raise RuntimeError('rebuild failed:'+name)
        for file in sorted(original.rglob('*')):
            if not file.is_file(): continue
            rebuilt=fresh/file.relative_to(original)
            if not rebuilt.is_file(): raise ValueError('fresh file missing')
            mode='BYTE_IDENTICAL'
            if sha256(file)!=sha256(rebuilt):
                if name=='fields' and file.name=='summary.json':
                    compare_json(json.loads(file.read_text('utf-8')),json.loads(rebuilt.read_text('utf-8')))
                    mode='ONLY_DECLARED_REDUCTION_SCALARS_DIFFER_WITHIN_DERIVED_BOUND'
                elif name=='visual' and file.name=='summary.json':
                    x=json.loads(file.read_text('utf-8')); y=json.loads(rebuilt.read_text('utf-8'))
                    if x['input_summary_sha256']!=sha256(study('hs4_height_wavefield_results')/'summary.json') or y['input_summary_sha256']!=sha256(freshroot/'fields/summary.json'):
                        raise ValueError('visual provenance differs from its actual input')
                    x.pop('input_summary_sha256'); y.pop('input_summary_sha256')
                    if x!=y: raise ValueError('visual products/metrics differ')
                    mode='IDENTICAL_EXCEPT_HASH_OF_VERIFIED_FRESH_INPUT_SUMMARY'
                else: raise ValueError('byte mismatch:'+str(file))
            records.append({'original':file.relative_to(root).as_posix(),'original_sha256':sha256(file),'rebuilt_sha256':sha256(rebuilt),'mode':mode})
        print('CPU rebuild checked:'+name,flush=True)
    if a.full_local:
        for script,args,original,target in [
            ('check_hs4_height_wavefield.py',['--study',study('hs4_height_wavefield_continuation')],study('hs4_height_wavefield_continuation')/'independent_verification.json',a.out/'independent_verification.json'),
            ('audit_hs4_height_native_archive.py',[],study('hs4_height_archive_verification')/'native_archive.json',a.out/'native_archive.json')]:
            subprocess.run([sys.executable,str(root/'scripts'/script),*map(str,args),'--out',str(target)],cwd=root,check=True)
            if sha256(original)!=sha256(target): raise ValueError('independent audit differs:'+script)
            records.append({'original':original.relative_to(root).as_posix(),'original_sha256':sha256(original),'rebuilt_sha256':sha256(target),'mode':'BYTE_IDENTICAL'})
    report={'status':'PASS','calls_solver':False,'code_sha256':sha256(__file__),'full_local':a.full_local,
        'fresh_products_reused':a.reuse_fresh is not None,'checked_files':len(records),
        'byte_identical_files':sum(r['mode']=='BYTE_IDENTICAL' for r in records),
        'bounded_scalar_reduction_differences':reductions,'norm_bound_derivation':{'positive_squared_terms_N':n,'binary64_eps':u,'gamma_8N':gamma},
        'files':records,'scope':'Execution reproducibility only, not FDTD absolute error/physics. Only eight named metadata reductions may differ under derived binary64 bound. All saved scientific arrays/images require exact bytes. Full mode needs local ignored files; default archived receivers suffice.'}
    (a.out/'rebuild_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('PASS files:'+str(len(records))+' exact:'+str(report['byte_identical_files']),flush=True)


if __name__=='__main__': main()
