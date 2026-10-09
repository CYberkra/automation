"""Rehashed invalid cover observer inputs must reject without a solver."""
import argparse,json,shutil,tempfile
from pathlib import Path
import h5py
from hs_capsule_identity import sha256 as sha
from audit_line9_cover_plane_inputs import check


def main(a):
    positive=check(a.package);rejected=[];a.work.mkdir(parents=True,exist_ok=True)
    cases=['source_changed_rehashed','new_Hx_wrong_neighbor','new_Hy_wrong_neighbor','old_receiver_changed','retry_enabled',
           'wrong_time_collocation','new_plane_in_air','geometry_changed_rehashed','main_receiver_changed','duplicate_receiver_name','SFCW_grid_changed']
    for name in cases:
        with tempfile.TemporaryDirectory(prefix=name+'_',dir=a.work) as folder:
            package=Path(folder)/'package';shutil.copytree(a.package,package)
            path=package/'manifest.json';m=json.loads(path.read_text('utf-8'));g=m['groups'][0];card=package/g['input']
            if name=='source_changed_rehashed':card.write_bytes(card.read_bytes().replace(b'ricker 40 ',b'ricker 41 '))
            elif name in ['new_Hx_wrong_neighbor','new_Hy_wrong_neighbor']:
                q=m['probes'][287]['anchors'][1 if 'Hx' in name else 2];q['coord'][1 if 'Hx' in name else 0]+=2
            elif name=='old_receiver_changed':card.write_bytes(card.read_bytes().replace(b'p_air_surface_c00_a0',b'p_changed_c00_a0'))
            elif name=='retry_enabled':m['no_retry']=False
            elif name=='wrong_time_collocation':m['collocation']['H_at_E_time']='unaligned H[n]'
            elif name=='new_plane_in_air':
                p=m['probes'][287];shift=280;p['y_m']+=shift*.025;p['wanted_y_m']=p['y_m']
                for q in p['anchors']:q['coord'][1]+=shift;q['position_m'][1]+=shift*.025
            elif name=='geometry_changed_rehashed':
                p=m['probes'][287]['anchors'][0]['coord']
                with h5py.File(package/g['geometry'],'r+') as h:h['data'][p[0],p[1],0]=0
                g['geometry_sha256']=sha(package/g['geometry'])
            elif name=='main_receiver_changed':card.write_bytes(card.read_bytes().replace(b'#rx: 170.65 ',b'#rx: 170.675 ',1))
            elif name=='duplicate_receiver_name':
                q=m['probes'][287]['anchors'][0];before=q['name'];q['name']=m['probes'][0]['anchors'][0]['name'];card.write_bytes(card.read_bytes().replace(before.encode(),q['name'].encode()))
            else:m['SFCW']['step_Hz']=500000
            g['input_sha256']=sha(card);path.write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
            try:check(package)
            except (AssertionError,ValueError):rejected.append(name)
            else:raise AssertionError('Invalid observer accepted: '+name)
    result={'status':'PASS_ELEVEN_REHASHED_COVER_OBSERVER_REJECTIONS','positive':positive,'rejected':rejected,'script_sha256':sha(__file__),'new_solver_runs':0}
    assert not a.out.exists();a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['package','work','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
