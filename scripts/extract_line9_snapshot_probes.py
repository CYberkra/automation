"""Extract five fixed Ez probe strips from completed local-ROI snapshots, no solver."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.execution/'execution_contract.json');v=read(a.execution/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(a.execution/'execution_contract.json')
    row=v['groups'][0];assert row['snapshot_count']==500 and row['observer_receiver_and_source_bitwise_equal']
    m=c['study_manifest'];roi=m['ROIs']['local_view'];assert roi['shape']==[350,315,1]
    geometry=Path(c['package'])/m['groups'][0]['geometry']
    assert sha(geometry)==m['groups'][0]['geometry_sha256']
    with h5py.File(geometry) as h:data=h['data'][:,:,0]
    xx=145+(np.arange(350)+.5)*.1;iy=10+(np.arange(315)+.5)*.1
    ix=np.floor(xx/.025).astype(int)
    surface=np.max(np.where(data!=0,np.arange(1700)[None,:]+1,0),axis=1)*.025
    mudtop=np.max(np.where(data==2,np.arange(1700)[None,:]+1,0),axis=1)*.025
    wanted=np.array([np.full(350,c['groups'][0]['rx_m'][1]),surface[ix]+.15,surface[ix]-.15,mudtop[ix]+.15,mudtop[ix]-.15])
    nearest=abs(iy[None,None,:]-wanted[:,:,None]).argmin(axis=2)
    actual=iy[nearest];assert abs(actual-wanted).max()<=.05+1e-12
    records=[r for r in row['snapshots'] if r['roi']=='local_view'];assert len(records)==250
    assert [r['iteration'] for r in records]==list(range(0,8500,34))
    fields=[];times=[]
    for r in records:
        p=Path(r['file']);assert sha(p)==r['sha256']
        with h5py.File(p) as h:
            assert h.attrs['iteration']==r['iteration'] and abs(h.attrs['time']-r['iteration']*m['dt_s'])<1e-20
            x=h['Ez'][:,:,0];assert x.dtype==np.float64 and x.shape==(350,315) and np.isfinite(x).all()
            fields.append(x[np.arange(350)[None,:],nearest]);times.append(float(h.attrs['time']))
    a.out.mkdir(parents=True)
    np.savez_compressed(a.out/'probes.npz',time_s=np.array(times),x_m=xx,desired_y_m=wanted,actual_y_m=actual,
        iy=nearest,Ez_V_m=np.array(fields),iterations=np.array([r['iteration'] for r in records]))
    names=['接收器高度附近','地表空气侧+0.15m','地表覆盖层侧−0.15m','覆盖层底上方+0.15m','覆盖层底泥岩侧−0.15m']
    result=dict(status='PASS250_EXISTING_SNAPSHOT_EZ_PROBES_NO_SOLVER',script_sha256=sha(__file__),
        contract_sha256=sha(a.execution/'execution_contract.json'),verification_sha256=sha(a.execution/'completed_verification.json'),
        probes_sha256=sha(a.out/'probes.npz'),geometry_sha256=sha(geometry),probe_names=names,snapshot_records=records,
        sampling='Native snapshot coarse-cell centres; nearest centre to desired y, no new interpolation or per-frame normalization.',
        max_y_rounding_m=float(abs(actual-wanted).max()),receiver_x_m=c['groups'][0]['rx_m'][0],
        receiver_y_m=c['groups'][0]['rx_m'][1],solver_runs=0,
        limits='Native Ricker Ez observations; temporal2.005ns and spatial0.1m sampling, probe not exact Yee receiver; no E/H energy or unique-bounce attribution.')
    (a.out/'extraction.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],files=250,probe_shape=[250,5,350],probes_sha256=result['probes_sha256'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['execution','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
