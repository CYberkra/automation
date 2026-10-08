"""Fresh bounded 190->185m paired-loss line; historical capsules stay read-only."""
import argparse, copy, json, shutil
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner

ROOT = Path(__file__).resolve().parents[1]
AUTHORITY = ('User ongoing autonomous goal: 持续推进直到研究明白; implements the proposed '
             '5m/0.25m continuous paired-loss diagnostic,46 new attempts,2 verified reuse. '
             'No production material replacement,full195m scan,training or snapshots.')

def read(p): return json.loads(p.read_text('utf-8'))

def prepare(a):
    from review_line9_result_packages import geometry_at, indices
    from diagnose_line9_layer_kinematics import primary_paths
    from analyze_line9_v401_version_controls import inverse
    from gprMax.toolboxes.SFCW.processing import spectral_window
    assert not a.out.exists()
    old = read(a.parent/'manifest.json')
    evidence = ROOT/'artifacts/research_checks/2026-10-08_v401_line_prebatch_results_r1'
    assert read(evidence/'completed_verification.json')['completed']
    assert read(evidence/'independent_audit.json')['status'].startswith('PASS')
    a.out.mkdir(parents=True)
    checks = a.out/'prebatch_evidence'; checks.mkdir()
    for name in ['execution_contract.json','completed_verification.json','analysis.json','independent_audit.json','completion_summary.json']:
        shutil.copyfile(evidence/name, checks/name)
    bases = {}; dbs = {}
    for name, b in old['baseline_sources'].items():
        parent = Path(b['parent']); g = b['group']; folder = a.out/'baselines'/name
        folder.mkdir(parents=True); r = copy.deepcopy(g)
        for key in ['input','geometry','material']:
            src = parent/g[key]; assert sha(src)==g[key+'_sha256']
            dst = folder/Path(g[key]).name; shutil.copyfile(src,dst)
            r[key] = dst.relative_to(a.out).as_posix()
        bases[name] = r
        dbs[name.split('_')[0]] = read(a.out/r['material'])['materials']
    with h5py.File(a.out/bases['low_H1']['geometry']) as h:
        data = h['data'][:,:,0]; spacing = h.attrs['dx_dy_dz']
    positions = [190.-i*.25 for i in range(21)]
    anchors = [190.,187.5,185.]; groups=[]; reused=[]; stations=[]
    oldreuse = {(r['chainage_m'],r['variant'],r['role']):r for r in old['reused']}
    for x in positions:
        tx=np.array([x-20.65,0.,.0125]); rx=np.array([x-19.35,0.,.0125])
        for q in [tx,rx]: q[1]=(np.flatnonzero(data[round(q[0]/.025)]!=0)[-1]+1)*.025+8
        templates={}; geometry={}
        for variant in ['low','high']:
            db=dbs[variant]; geo=geometry_at(data,spacing,tx,rx,indices(db,95e6)); geometry[variant]=geo
            spectra=primary_paths(geo,db); j=geo['boundaries'].index(geo['basal_sand'])
            z,t=inverse(spectra[j][:,None],spectral_window('hann',501)); peak=float(t[np.argmax(abs(z[:,0]))]*1e9)
            templates[variant]=dict(peak_ns=peak,basal_gate_ns=[peak-12,peak+12],wide_gate_ns=[peak-60,peak+60])
        s=dict(chainage_m=x,tx_m=tx.tolist(),rx_m=rx.tolist(),geometry=geometry,templates=templates)
        stations.append(s)
        # Consecutive paired H1 columns; anchor H0 before its H1 gives early usable pairs.
        requests=([(v,'H0') for v in ['low','high']] if x in anchors else [])+[(v,'H1') for v in ['low','high']]
        for variant,role in requests:
            sid=f'{variant}_x{round(x*100):05d}_{role}'; base=bases[variant+'_'+role]
            if (x,variant,role) in oldreuse:
                prior=oldreuse[(x,variant,role)]; src=a.parent/prior['path']; assert sha(src)==prior['native_sha256']
                with h5py.File(src) as h:
                    for node,q in [('srcs/src1',tx),('rxs/rx1',rx)]: np.testing.assert_allclose(h[node].attrs['Position'],[*q[:2],0.],rtol=0,atol=1e-11)
                dst=a.out/'reused'/(sid+'.h5'); dst.parent.mkdir(exist_ok=True); shutil.copyfile(src,dst)
                reused.append(dict(id=sid,variant=variant,role=role,chainage_m=x,path=dst.relative_to(a.out).as_posix(),native_sha256=sha(dst),tx_m=[*tx[:2],0.],rx_m=[*rx[:2],0.]))
                continue
            folder=a.out/sid; gp=folder/'geometries'; gp.mkdir(parents=True)
            for key in ['geometry','material']: shutil.copyfile(a.out/base[key],gp/Path(base[key]).name)
            lines=[]
            for text in (a.out/base['input']).read_text('utf-8').splitlines():
                if text.startswith('#title:'): text=f'#title: {sid}; dense paired-loss diagnostic; no snapshots'
                elif text.startswith('#hertzian_dipole:'): text='#hertzian_dipole: z '+' '.join(f'{q:.12g}' for q in tx)+' pulse'
                elif text.startswith('#rx:'): text='#rx: '+' '.join(f'{q:.12g}' for q in rx)+' '+sid+'_rx1 Ez'
                lines.append(text)
            inp=folder/'cases'/sid/'profile.in'; inp.parent.mkdir(parents=True); inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            g=copy.deepcopy(base); g.update(id=sid,variant=variant,role=role,chainage_m=x,input=inp.relative_to(a.out).as_posix(),geometry=(gp/Path(base['geometry']).name).relative_to(a.out).as_posix(),material=(gp/Path(base['material']).name).relative_to(a.out).as_posix(),tx_m=[*tx[:2],0.],rx_m=[*rx[:2],0.],input_tx_m=tx.tolist(),input_rx_m=rx.tolist(),**templates[variant])
            for key in ['input','geometry','material']: g[key+'_sha256']=sha(a.out/g[key])
            groups.append(g)
    assert len(groups)==46 and len(reused)==2 and len(stations)==21
    m=dict(status='PREPARED_AUTONOMOUS_DENSE_DIAGNOSTIC_NOT_RUN',approval_basis=AUTHORITY,groups=groups,reused=reused,stations=stations,chainage_m=positions,anchors_m=anchors,baseline_copies=bases,dt_s=old['dt_s'],generator_sha256=sha(__file__),locator_manifest_sha256=sha(a.parent/'manifest.json'),
        evidence_hashes={p.name:sha(p) for p in checks.iterdir()},
        processing=dict(frequency_start_Hz=20e6,frequency_step_Hz=300000.,tones=501,source_reference_m=.025,tail_taper_fraction=0,windows=['hann','blackman'],complex_subtraction_before_magnitude=True),
        limits='5m paired-loss shape diagnostic,not complete195m/global spatial bandlimit/field fit. Low mud is a counterfactual; both real dispersion and DC change. H0 only3anchors,not continuous clean reference. Main-window boundary sensitivity tested at220m,not all-station/all-window convergence. No retries,snapshots,training,AGC,taper,delay/amplitude/phase fitting.')
    runner.save(a.out/'manifest.json',m)
    print(json.dumps(dict(status=m['status'],new=46,reused=2,stations=21)))

def freeze(a):
    m=read(a.package/'manifest.json'); v=read(a.package/'independent_input_audit.json')
    assert runner.gprMax.__version__=='4.0.1'
    assert v['manifest_sha256']==sha(a.package/'manifest.json') and v['status'].startswith('PASS') and v['bulk_execution_ready'] is True
    assert v['auditor_sha256']==sha(Path(__file__).with_name('audit_line9_dense_loss_inputs.py'))
    for name,digest in m['evidence_hashes'].items(): assert sha(a.package/'prebatch_evidence'/name)==digest
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    p=a.out/'execution_contract.json'; c=read(p)
    c.update(study_manifest=m,approval_basis=AUTHORITY,max_batch_wall_s=14400)
    for name in ['line9_v401_dense_loss_line.py','audit_line9_dense_loss_inputs.py']:
        path=Path(__file__).with_name(name);c['code_identities'][str(path.resolve())]=sha(path)
    for path in [a.package/'manifest.json',a.package/'independent_input_audit.json']+[a.package/r['path'] for r in m['reused']]+list((a.package/'prebatch_evidence').iterdir()):
        c['file_identities'][str(path.resolve())]=sha(path)
    runner.save(p,c);runner.save(a.out/'preflight_verification.json',runner.audit(p))
    print(json.dumps(dict(status='FROZEN_46_NEW_2_REUSED_DENSE_DIAGNOSTIC',contract_sha256=sha(p))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze'])
    for key in ['parent','package','out']:p.add_argument('--'+key,type=Path,required=key=='out')
    a=p.parse_args();a.out=a.out.resolve()
    prepare(a) if a.action=='prepare' else freeze(a)
