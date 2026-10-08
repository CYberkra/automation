"""Bounded 220->140m low-loss diagnostic line, with reused and paired controls."""
import argparse,copy,json,shutil
from pathlib import Path
import numpy as np
import h5py
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner

ROOT=Path(__file__).resolve().parents[1]
APPROVAL='User explicitly authorizes the proposed 220->140m,2m diagnostic line plus few same-station controls. Single attempts; no full195m line, snapshots, training or production material change.'

def read(p):return json.loads(p.read_text('utf-8'))

def prepare(out):
    from review_line9_result_packages import geometry_at,indices
    from diagnose_line9_layer_kinematics import primary_paths
    from analyze_line9_v401_version_controls import inverse
    from gprMax.toolboxes.SFCW.processing import spectral_window
    assert not out.exists()
    low=ROOT/'artifacts/local_checks/2026-10-08_v401_spatial_prepared_r1'
    lowresult=ROOT/'artifacts/local_checks/2026-10-08_v401_spatial_results_r1'
    high=ROOT/'artifacts/local_checks/2026-10-08_v401_nonflat_material_prepared_r1'
    highresult=ROOT/'artifacts/local_checks/2026-10-08_v401_nonflat_material_results_r1'
    plans=ROOT/'artifacts/local_checks/2026-10-07_line9_three_packages_review_r1/line9_pkg1_v5_d02m_194st/package_manifest.json'
    lm=read(low/'manifest.json');lv=read(lowresult/'completed_verification.json');hm=read(high/'manifest.json');hv=read(highresult/'completed_verification.json')
    assert lv['completed'] and hv['completed']
    assert lv['contract_sha256']==sha(lowresult/'execution_contract.json') and hv['contract_sha256']==sha(highresult/'execution_contract.json')
    bases={('low',r): (low,next(g for g in lm['groups'] if g['id']=='full2d_r0141_'+r)) for r in ['H0','H1']}
    bases[('high','H1')]=(high,next(g for g in hm['groups'] if g['id']=='baseline_H1'))
    # The reused high H0's exact original input/map is from the v5 source controls.
    source=ROOT/'artifacts/local_checks/2026-10-08_line9_v5_source_flat_prepared_r1'
    bases[('high','H0')]=(source,next(g for g in read(source/'manifest.json')['groups'] if g['id']=='nonflat_ricker_H0'))
    out.mkdir(parents=True);reused=[];groups=[];stations=[]
    existing={}
    for g in lm['groups']:
        record=next(r for r in lv['groups'] if r['id']==g['id']);existing[('low',g['chainage_m'],g['id'][-2:])]=(lowresult/(g['id']+'.h5'),record['native_sha256'])
    for r in lm['reused']:existing[('low',198.6,r['id'][-2:])]=(low/r['path'],r['native_sha256'])
    for role in ['H0','H1']:
        if role=='H1':record=next(r for r in hv['groups'] if r['id']=='baseline_H1');raw=highresult/'baseline_H1.h5'
        else:record=hm['reused_H0'];raw=high/'reused_nonflat_H0.h5'
        existing[('high',198.6,role)]=(raw,record['native_sha256'])
    for (variant,role),(parent,b) in bases.items():
        for key in ['input','geometry','material']:assert sha(parent/b[key])==b[key+'_sha256'],(variant,role,key,parent/b[key])
    parent,b=bases[('low','H1')];db=read(parent/b['material'])['materials'];n95=indices(db,95e6)
    with h5py.File(parent/b['geometry']) as h:data=h['data'][:,:,0];dl=h.attrs['dx_dy_dz']
    plan=read(plans);cases={round(g['chainage_m'],6):g for g in plan['cases']}
    line=list(np.arange(220.,139.,-2));positions=sorted(set(line+[198.6,146.4]),reverse=True)
    for x in positions:
        if x in cases:case=cases[x];tx=np.array(case['tx_m']);rx=np.array(case['rx_m'])
        else:
            assert x==220.;tx=np.array([x-20.65,0.,.0125]);rx=np.array([x-19.35,0.,.0125])
            for q in [tx,rx]:q[1]=(np.flatnonzero(data[round(q[0]/.025)]!=0)[-1]+1)*.025+8
        geometry=geometry_at(data,dl,tx,rx,n95);assert abs(geometry['midpoint_agl_m']-8)<.025
        template=primary_paths(geometry,db)[geometry['boundaries'].index(geometry['basal_sand'])]
        profile,t=inverse(template[:,None],spectral_window('hann',501));peak=float(t[np.argmax(abs(profile[:,0]))]*1e9)
        station=dict(chainage_m=x,tx_m=tx.tolist(),rx_m=rx.tolist(),geometry=geometry,template_peak_ns=peak,basal_gate_ns=[peak-12,peak+12],wide_gate_ns=[peak-60,peak+60],regular_line_station=x in line)
        stations.append(station)
        requested=[('low','H1')]
        if x in [210.,200.,198.6,190.,180.,160.,146.4]:requested.append(('low','H0'))
        if x in [198.6,190.,180.,146.4]:requested.extend([('high','H0'),('high','H1')])
        for variant,role in requested:
            sid=f'{variant}_x{round(x*10):04d}_{role}';parent,b=bases[(variant,role)]
            if (variant,x,role) in existing:
                raw,digest=existing[(variant,x,role)];assert sha(raw)==digest
                dst=out/'reused'/(sid+'.h5');dst.parent.mkdir(exist_ok=True);shutil.copyfile(raw,dst)
                reused.append(dict(id=sid,variant=variant,role=role,chainage_m=x,path=dst.relative_to(out).as_posix(),native_sha256=digest,tx_m=tx[:2].tolist()+[0.],rx_m=rx[:2].tolist()+[0.]))
                continue
            folder=out/sid;gp=folder/'geometries';gp.mkdir(parents=True)
            for key in ['geometry','material']:shutil.copyfile(parent/b[key],gp/Path(b[key]).name)
            inp=folder/'cases'/sid/'profile.in';inp.parent.mkdir(parents=True)
            lines=[]
            for text in (parent/b['input']).read_text('utf-8').splitlines():
                if text.startswith('#title:'):text=f'#title: {sid}; authorized diagnostic line; no snapshots'
                elif text.startswith('#hertzian_dipole:'):text='#hertzian_dipole: z '+' '.join(f'{v:.12g}' for v in tx)+' pulse'
                elif text.startswith('#rx:'):text='#rx: '+' '.join(f'{v:.12g}' for v in rx)+' '+sid+'_rx1 Ez'
                lines.append(text)
            inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            g=copy.deepcopy(b);g.update(id=sid,variant=variant,role=role,chainage_m=x,input=inp.relative_to(out).as_posix(),geometry=(gp/Path(b['geometry']).name).relative_to(out).as_posix(),material=(gp/Path(b['material']).name).relative_to(out).as_posix(),tx_m=tx[:2].tolist()+[0.],rx_m=rx[:2].tolist()+[0.],input_tx_m=tx.tolist(),input_rx_m=rx.tolist(),basal_gate_ns=station['basal_gate_ns'],wide_gate_ns=station['wide_gate_ns'])
            for key in ['input','geometry','material']:g[key+'_sha256']=sha(out/g[key])
            groups.append(g)
    # Low total fields first; then remaining counterfactual controls.
    groups.sort(key=lambda g:(0 if (g['variant'],g['role'])==('low','H1') else 1,-g['chainage_m'],g['id']))
    assert len(groups)==48 and len(reused)==10 and len(stations)==43
    manifest=dict(status='PREPARED_APPROVED_NOT_RUN',approval_basis=APPROVAL,groups=groups,reused=reused,stations=stations,regular_line_chainage_m=line,
        dt_s=lm['dt_s'],generator_sha256=sha(__file__),parent_manifest_sha256=sha(low/'manifest.json'),high_manifest_sha256=sha(high/'manifest.json'),station_plan_sha256=sha(plans),
        baseline_sources={variant+'_'+role:dict(parent=str(p),group=b) for (variant,role),(p,b) in bases.items()},
        invariants='210x42.5m/2.5cm/1200ns/AGL8m/40A100MHz Ricker/HORIPML80/averaging/FP64, exact501tones20-170MHz/.3MHz, .025m reference, no taper/AGC/fits',
        limits='2m coarse diagnostic sampling, not spatial Nyquist certification or complete195m line. Low-loss mudstone counterfactual only, no site calibration. High controls sparse; original194 traces are4.0.0 FP32. H1-H0 includes all bottom-replacement interactions. No snapshots/training/production-material change.')
    runner.save(out/'manifest.json',manifest);print(json.dumps(dict(status=manifest['status'],new=len(groups),reused=len(reused),stations=len(stations))))

def freeze(out,package):
    assert runner.gprMax.__version__=='4.0.1';m=read(package/'manifest.json');audit=read(package/'independent_input_audit.json')
    assert audit['manifest_sha256']==sha(package/'manifest.json') and audit['status'].startswith('PASS')
    assert audit.get('bulk_execution_ready') is True, 'Unresolved prebatch controls; bulk freeze prohibited'
    runner.prepare(argparse.Namespace(action='prepare',out=out,package=package))
    path=out/'execution_contract.json';c=read(path);c.update(study_manifest=m,approval_basis=APPROVAL,max_batch_wall_s=14400)
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    for p in [package/'manifest.json',package/'independent_input_audit.json']+[package/r['path'] for r in m['reused']]:c['file_identities'][str(p.resolve())]=sha(p)
    runner.save(path,c);runner.save(out/'preflight_verification.json',runner.audit(path));print(json.dumps(dict(status='FROZEN_48_NEW_10_REUSED',contract_sha256=sha(path))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze']);p.add_argument('--out',type=Path,required=True);p.add_argument('--package',type=Path)
    a=p.parse_args();a.out=a.out.resolve();prepare(a.out) if a.action=='prepare' else freeze(a.out,a.package.resolve())
