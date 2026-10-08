"""Independent dense-line cards/voxels/materials/reuse/local sampling audit."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha

def read(p):return json.loads(p.read_text('utf-8'))

def fields(p):
    result={}
    for line in p.read_text('utf-8').splitlines():
        key,*val=line.split();assert key not in result;result[key]=val
    return result

def audit(package):
    m=read(package/'manifest.json');bases=m['baseline_copies'];maps={};cards={};dbs={}
    expected=[190-i/4 for i in range(21)]
    assert m['chainage_m']==expected and m['anchors_m']==[190.,187.5,185.]
    assert [s['chainage_m'] for s in m['stations']]==expected
    for name,g in bases.items():
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        cards[name]=fields(package/g['input']);dbs[name]=read(package/g['material'])['materials']
        with h5py.File(package/g['geometry']) as h:
            assert h['data'].shape==(8400,1700,1);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            maps[name]=h['data'][:,:,0]
            assert list(dbs[name])==[k.decode() for k in h['material_keys'][:]]
    for role in ['H0','H1']:np.testing.assert_array_equal(maps['low_'+role],maps['high_'+role]);assert dbs['low_'+role]==dbs['low_H1'] and dbs['high_'+role]==dbs['high_H1']
    changed=maps['low_H1']!=maps['low_H0'];assert changed.sum()==5606675
    assert np.all(maps['low_H1'][changed]==3) and np.all(maps['low_H0'][changed]==2)
    low,high=dbs['low_H1'],dbs['high_H1'];mud='material_002_mudstone'
    for key in low:
        if key!=mud:assert low[key]==high[key]
    for variant,db,eps,delta,sigma in [('low',low,11.973951024930017,.6200368247790599,.0003),('high',high,11.739510249300162,6.200368247790598,.003)]:
        b=db[mud]['base'];p=db[mud]['poles'];assert b['relative_permittivity']==eps and b['electric_conductivity_s_per_m']==sigma
        assert len(p)==1 and p[0]['relative_permittivity_difference']==delta and p[0]['relaxation_time_s']==8e-9
        assert abs(eps+delta/(1+(2*np.pi*95e6*8e-9)**2)-12)<1e-12
    rows=[];stations={s['chainage_m']:s for s in m['stations']};allowed={'#title:','#hertzian_dipole:','#rx:'}
    for s in m['stations']:
        tx=np.array(s['tx_m']);rx=np.array(s['rx_m']);assert abs((tx[0]+rx[0])/2+20-s['chainage_m'])<1e-10 and abs(rx[0]-tx[0]-1.3)<1e-10
        for q in [tx,rx]:
            np.testing.assert_allclose(q[:2]/.025,np.rint(q[:2]/.025),rtol=0,atol=1e-9);assert q[2]==.0125
            col=maps['low_H1'][round(q[0]/.025)];surface=(np.flatnonzero(col!=0)[-1]+1)*.025
            assert abs(q[1]-surface-8)<1e-10 and maps['low_H1'][round(q[0]/.025),round(q[1]/.025)]==0
        rows.append(dict(chainage_m=s['chainage_m'],top_PML_entry_clearance_m=float(40.5-max(tx[1],rx[1])),side_PML_entry_clearance_m=float(min(tx[0]-2,208-rx[0]))))
    keys=[]
    for g in m['groups']:
        name=g['variant']+'_'+g['role'];b=bases[name];s=stations[g['chainage_m']]
        for key in ['input','geometry','material']:assert sha(package/g[key])==g[key+'_sha256']
        for key in ['geometry','material']:assert g[key+'_sha256']==b[key+'_sha256']
        c=fields(package/g['input']);base=cards[name]
        assert {k:v for k,v in c.items() if k not in allowed}=={k:v for k,v in base.items() if k not in allowed}
        for key,start,role in [('#hertzian_dipole:',1,'tx'),('#rx:',0,'rx')]:
            q=np.array(c[key][start:start+3],float);np.testing.assert_allclose(q,s[role+'_m'],rtol=0,atol=1e-11)
            np.testing.assert_allclose([*q[:2],0.],g[role+'_m'],rtol=0,atol=1e-11)
        assert c['#hertzian_dipole:'][0]=='z' and c['#rx:'][-1]=='Ez'
        keys.append((g['variant'],g['chainage_m'],g['role']))
    for r in m['reused']:
        assert (r['variant'],r['chainage_m'])==('low',190.) and sha(package/r['path'])==r['native_sha256']
        with h5py.File(package/r['path']) as h:
            assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64 and h['rxs/rx1/Ez'].shape==(20352,)
            assert abs(h.attrs['dt']/m['dt_s']-1)<1e-14
            for node,role in [('srcs/src1','tx'),('rxs/rx1','rx')]:np.testing.assert_allclose(h[node].attrs['Position'],r[role+'_m'],rtol=0,atol=1e-11)
        keys.append((r['variant'],r['chainage_m'],r['role']))
    wanted={(v,x,'H1') for v in ['low','high'] for x in expected}|{(v,x,'H0') for v in ['low','high'] for x in [190.,187.5,185.]}
    assert len(keys)==len(set(keys))==48 and set(keys)==wanted and len(m['groups'])==46 and len(m['reused'])==2
    checks=package/'prebatch_evidence'
    for name,digest in m['evidence_hashes'].items():assert sha(checks/name)==digest
    c=read(checks/'execution_contract.json');v=read(checks/'completed_verification.json');a=read(checks/'analysis.json');r=read(checks/'independent_audit.json')
    assert v['completed'] and len(v['groups'])==6 and v['contract_sha256']==a['contract_sha256']==r['contract_sha256']==sha(checks/'execution_contract.json')
    assert r['status'].startswith('PASS') and r['analysis_sha256']==sha(checks/'analysis.json')
    sampling={}
    for variant in ['low','high']:
        times=np.array([s['geometry'][variant]['basal_sand']['time95_ns'] for s in m['stations']]);phase=np.abs(np.diff(times))*2*np.pi*.095
        assert np.all(phase<np.pi)
        sampling[variant]=dict(max_local_phase95_increment_rad=float(phase.max()),exceeds_pi_intervals=int(sum(phase>np.pi)),limits='Local95MHz propagation proxy only; not complete501-tone/nonlocal spatial bandlimit certificate.')
    return dict(status='PASS_DENSE_INPUT_VOXEL_MATERIAL_REUSE_PREBATCH_IDENTITY',auditor_sha256=sha(__file__),manifest_sha256=sha(package/'manifest.json'),new=46,reused=2,stations=21,positions=rows,sampling=sampling,bulk_execution_ready=True,authorization_scope='Bounded5m paired-loss diagnostic only; original80m coarse batch remains gated.',boundary_interpretation='Predeclared basal window finite extensions show<0.3% change at220m; no formal physical threshold or all-station/full-window convergence certificate.',limits=m['limits'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert not a.out.exists()
    r=audit(a.package);a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps({k:r[k] for k in ['status','new','reused','stations','sampling','bulk_execution_ready']}))
