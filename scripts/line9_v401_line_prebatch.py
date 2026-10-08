"""Six single-attempt worst-clearance controls before the coarse diagnostic line."""
import argparse,copy,json,shutil
from pathlib import Path
import h5py,numpy as np
from hs_capsule_identity import sha256 as sha
import line9_v401_version_controls as runner

def read(p):return json.loads(p.read_text('utf-8'))

def prepare(a):
    assert not a.out.exists();m=read(a.parent/'manifest.json');audit=read(a.parent/'independent_input_audit.json')
    assert audit['manifest_sha256']==sha(a.parent/'manifest.json') and not audit['bulk_execution_ready']
    station=next(s for s in m['stations'] if s['chainage_m']==220.)
    a.out.mkdir(parents=True);groups=[]
    for factor in ['base','top20','right40']:
        for role in ['H0','H1']:
            ref=m['baseline_sources']['low_'+role];parent=Path(ref['parent']);b=ref['group'];name=factor+'_'+role;gp=a.out/name/'geometries';gp.mkdir(parents=True)
            for key in ['input','geometry','material']:assert sha(parent/b[key])==b[key+'_sha256']
            geom=gp/Path(b['geometry']).name
            if factor=='base':shutil.copyfile(parent/b['geometry'],geom);shape=b['native_shape']
            else:
                with h5py.File(parent/b['geometry']) as src,h5py.File(geom,'x') as dst:
                    old=src['data'][:]
                    data=np.pad(old,((0,0),(0,800),(0,0)),mode='constant',constant_values=0) if factor=='top20' else np.pad(old,((0,1600),(0,0),(0,0)),mode='edge')
                    for key,value in src.attrs.items():dst.attrs[key]=value
                    for key in src:
                        if key!='data':src.copy(key,dst)
                    ds=dst.create_dataset('data',data=data,compression='gzip',compression_opts=4)
                    for key,value in src['data'].attrs.items():ds.attrs[key]=value
                    shape=list(data.shape);dst.attrs['shape_nxyz']=shape
                    dst.attrs['BoundaryControl']='Original8400x1700 interior exact; top20 pure air or right40 edge-column continuation; no coordinate translation'
            mat=gp/Path(b['material']).name;shutil.copyfile(parent/b['material'],mat)
            lines=[]
            for text in (parent/b['input']).read_text('utf-8').splitlines():
                if text.startswith('#title:'):text=f'#title: {name}; X220 low-loss prebatch boundary control; no snapshots'
                elif text.startswith('#domain:'):text=f'#domain: {shape[0]*.025:g} {shape[1]*.025:g} inf'
                elif text.startswith('#hertzian_dipole:'):text='#hertzian_dipole: z '+' '.join(f'{q:.12g}' for q in station['tx_m'])+' pulse'
                elif text.startswith('#rx:'):text='#rx: '+' '.join(f'{q:.12g}' for q in station['rx_m'])+' '+name+'_rx1 Ez'
                lines.append(text)
            inp=a.out/name/'cases'/name/'profile.in';inp.parent.mkdir(parents=True);inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            g=copy.deepcopy(b);g.update(id=name,factor=factor,role=role,chainage_m=220.,input=inp.relative_to(a.out).as_posix(),geometry=geom.relative_to(a.out).as_posix(),material=mat.relative_to(a.out).as_posix(),native_shape=shape,tx_m=station['tx_m'][:2]+[0.],rx_m=station['rx_m'][:2]+[0.],domain_height_m=shape[1]*.025)
            for key in ['input','geometry','material']:g[key+'_sha256']=sha(a.out/g[key])
            groups.append(g)
    out=dict(status='PREPARED_APPROVED_PREBATCH_ONLY',groups=groups,reused=[],station=station,dt_s=m['dt_s'],parent_manifest_sha256=sha(a.parent/'manifest.json'),baseline_sources=m['baseline_sources'],generator_sha256=sha(__file__),
        approval_basis='User authorizes diagnostic line then asks resolve outstanding issues before bulk execution. Six bounded220m low-loss top/right boundary pairs; bulk remains gated.',
        gates_ns=dict(basal=station['basal_gate_ns'],deep=[180.,450.],early=[0.,120.],late=[450.,1100.]),
        limits='Worst-clearance single-site finite extension sensitivity, not global convergence/PML-only attribution; right extension also adds exterior geology. No joint extension/finite3D/field certification. Original48 bulk attempts not started.')
    runner.save(a.out/'manifest.json',out);print(json.dumps(dict(status=out['status'],groups=len(groups))))

def freeze(a):
    assert runner.gprMax.__version__=='4.0.1';m=read(a.package/'manifest.json');r=read(a.package/'independent_input_audit.json')
    assert r['manifest_sha256']==sha(a.package/'manifest.json') and r['status'].startswith('PASS')
    runner.prepare(argparse.Namespace(action='prepare',out=a.out,package=a.package))
    p=a.out/'execution_contract.json';c=read(p);c.update(study_manifest=m,approval_basis=m['approval_basis'])
    c['code_identities'][str(Path(__file__).resolve())]=sha(__file__)
    for path in [a.package/'manifest.json',a.package/'independent_input_audit.json']:c['file_identities'][str(path.resolve())]=sha(path)
    runner.save(p,c);runner.save(a.out/'preflight_verification.json',runner.audit(p));print(json.dumps(dict(status='FROZEN_SIX_PREBATCH_NOT_BULK',contract_sha256=sha(p))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze'])
    for key in ['parent','package','out']:p.add_argument('--'+key,type=Path,required=key=='out')
    a=p.parse_args();a.out=a.out.resolve();prepare(a) if a.action=='prepare' else freeze(a)
