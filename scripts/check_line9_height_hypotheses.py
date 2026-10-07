"""Independent artifact/paired-field checks for the local height hypotheses."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_line9_postprocessing import FREQ,inverse,weights,sha,save
from diagnose_line9_layer_kinematics import primary_paths


def main(root,review,cache,pair,ray1,ray2,stack,out):
    if out.exists():raise ValueError('Fresh output directory required')
    pair_report=json.loads(pair.read_text('utf-8'));stack_report=json.loads(stack.read_text('utf-8'))
    ray_reports=[json.loads(p.read_text('utf-8')) for p in [ray1,ray2]]
    reports=[pair_report,stack_report,*ray_reports]
    sources=['diagnose_line9_height_pair.py','diagnose_line9_layer_stack.py',
             'diagnose_line9_refracted_paths.py','diagnose_line9_refracted_paths.py']
    for report,source in zip(reports,sources):
        assert report['script_sha256']==sha(Path(__file__).parent/source)
        assert report['calls_solver'] is False and report['calls_training'] is False
    assert pair_report['ray_json_sha256']==sha(ray2)
    assert 60 < min(pair_report['air_advance_ns']) < max(pair_report['air_advance_ns']) < 70
    native_count=0;audits=[]
    for ap in sorted(review.glob('*_audit.json')):
        a=json.loads(ap.read_text('utf-8'));audits.append(a)
        for row in a['records']:
            assert sha(root/a['package']/'cases'/row['id']/'profile.h5')==row['native_sha256']
            native_count+=1
        cache_hash=sha(cache/(a['package']+'.npz'))
        for rr in ray_reports:
            p=next(p for p in rr['packages'] if p['package']==a['package'])
            assert p['audit_sha256']==sha(ap) and p['private_cache_sha256']==cache_hash
        p=next(p for p in stack_report['packages'] if p['package']==a['package'])
        assert p['audit_sha256']==sha(ap) and p['cache_sha256']==cache_hash
    assert native_count==406
    a,b=audits[1:];count=70
    mf=next((root/b['package']/'geometries').glob('*.json'))
    materials=json.loads(mf.read_text('utf-8'))['materials']
    high=np.load(cache/(a['package']+'.npz'))['response'][:,:count]
    low=np.load(cache/(b['package']+'.npz'))['response']
    # Derive air delays again from audited geometries, independently of report values.
    delays=np.array([2*(ra['geometry']['midpoint_agl_m']-rb['geometry']['midpoint_agl_m'])/299792458.
                     for ra,rb in zip(a['records'][:count],b['records'])])
    np.testing.assert_allclose(delays*1e9,pair_report['air_advance_ns'],atol=1e-12,rtol=0)
    high*=np.cos(2*np.pi*FREQ[:,None]*delays)+1j*np.sin(2*np.pi*FREQ[:,None]*delays)
    vh,t=inverse(high,FREQ,weights('hann',501));vl,_=inverse(low,FREQ,weights('hann',501))
    cross=0j;aa=bb=prednorm=0.;predcross=0j;difference=0.;individual=[]
    for j,row in enumerate(b['records']):
        g=row['geometry'];paths=primary_paths(g,materials);centers={}
        for key in ['cover_base','basal_sand']:
            hp=paths[g['boundaries'].index(g[key])]
            vp,_=inverse(hp,FREQ,weights('hann',501))
            centers[key]=t[np.argmax(abs(vp))]
        shallow=abs(t-centers['cover_base'])<=12e-9
        deep=abs(t-centers['basal_sand'])<=12e-9
        lo=vl[deep,j];hi=vh[deep,j]
        cij=np.sum(lo.conj()*hi);al=float(np.sum(abs(lo)**2));bh=float(np.sum(abs(hi)**2))
        cross+=cij;aa+=al;bb+=bh;individual.append(abs(cij)/np.sqrt(al*bh))
        coefficient=np.sum(vl[shallow,j].conj()*vh[shallow,j])/np.sum(abs(vl[shallow,j])**2)
        prediction=lo*coefficient
        prednorm+=float(np.sum(abs(prediction)**2));predcross+=np.sum(prediction.conj()*hi)
        difference+=float(np.sum(abs(prediction-hi)**2))
    corr=float(abs(cross)/np.sqrt(aa*bb));relative_L2=float(np.sqrt(difference/bb))
    predcorr=float(abs(predcross)/np.sqrt(prednorm*bb))
    declared=pair_report['variants']['air_delay_only_hann']
    np.testing.assert_allclose(corr,declared['basal_sand']['global_complex_correlation'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(individual,declared['basal_sand']['per_trace_complex_correlation'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(relative_L2,declared['cover_fit_tested_on_basal']['relative_L2'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(predcorr,declared['cover_fit_tested_on_basal']['global_complex_correlation'],atol=1e-12,rtol=0)
    sensitivity=[]
    for p,q in zip(ray_reports[0]['packages'],ray_reports[1]['packages']):
        assert [r['id'] for r in p['records']]==[r['id'] for r in q['records']]
        delta=np.array([r['phase_time95_ns'] for r in p['records']])-np.array([r['phase_time95_ns'] for r in q['records']])
        sensitivity.append(dict(package=p['package'],number_of_paths=len(delta),
            max_time_change_ns=float(max(abs(delta))),
            largest_converged_start_time_spread_ns=float(max(-r['best_minus_worst_converged_time_ns'] for r in q['records']))))
    out.mkdir(parents=True)
    save(out/'independent_checks.json',dict(calls_solver=False,calls_training=False,
        script_sha256=sha(__file__),reports_sha256={str(p):sha(p) for p in [pair,ray1,ray2,stack]},
        native_hashes_unchanged_count=native_count,source_hashes_match=True,
        independent_hann_basal_pair_global_correlation=corr,
        independent_shallow_fit_deep_relative_L2=relative_L2,
        independent_shallow_fit_deep_correlation=predcorr,ray_curve_sensitivity=sensitivity,
        limits='Artifact/math verification only, not FDTD convergence, global ray optimality or independent physical attribution.'))
    print(json.dumps(dict(native_hashes_unchanged=native_count,independent_pair_correlation=corr,
                         independent_shallow_fit_deep_relative_L2=relative_L2,sensitivity=sensitivity),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','review','cache','pair','ray1','ray2','stack','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.root,a.review,a.cache,a.pair,a.ray1,a.ray2,a.stack,a.out)
