"""Independently read 14 cached native frames and render fixed-scale probe strips."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    r=read(a.source/'extraction.json');c=read(a.parent/'execution_contract.json');v=read(a.parent/'completed_verification.json')
    assert r['contract_sha256']==v['contract_sha256']==sha(a.parent/'execution_contract.json')
    assert r['verification_sha256']==sha(a.parent/'completed_verification.json') and v['completed']
    assert r['script_sha256']==sha(Path(__file__).with_name('extract_line9_snapshot_probes.py'))
    assert r['probes_sha256']==sha(a.source/'probes.npz') and r['solver_runs']==0
    assert sha(a.geometry)==r['geometry_sha256']==c['study_manifest']['groups'][0]['geometry_sha256']
    records=[x for x in v['groups'][0]['snapshots'] if x['roi']=='local_view']
    assert records==r['snapshot_records'] and len(records)==250
    with np.load(a.source/'probes.npz',allow_pickle=False) as z:arr={k:z[k] for k in z.files}
    field=arr['Ez_V_m'];assert field.shape==(250,5,350) and field.dtype==np.float64 and np.isfinite(field).all()
    np.testing.assert_array_equal(arr['iterations'],[x['iteration'] for x in records])
    np.testing.assert_allclose(arr['time_s'],arr['iterations']*c['study_manifest']['dt_s'],rtol=0,atol=1e-20)
    xx=145+.1*(np.arange(350)+.5);yy=10+.1*(np.arange(315)+.5)
    np.testing.assert_array_equal(arr['x_m'],xx)
    with h5py.File(a.geometry) as h:g=h['data'][:,:,0]
    indices=np.floor(xx/.025).astype(int)
    surface=np.array([np.flatnonzero(g[x]!=0)[-1]+1 for x in indices])*.025
    mudtop=np.array([np.flatnonzero(g[x]==2)[-1]+1 for x in indices])*.025
    wanted=np.array([np.full(350,r['receiver_y_m']),surface+.15,surface-.15,mudtop+.15,mudtop-.15])
    np.testing.assert_array_equal(wanted,arr['desired_y_m'])
    assert arr['iy'].shape==(5,350) and np.all((arr['iy']>=0)&(arr['iy']<315))
    np.testing.assert_array_equal(arr['actual_y_m'],yy[arr['iy']])
    assert np.max(abs(arr['actual_y_m']-wanted))<=.05+1e-12
    assert np.all(abs(arr['actual_y_m']-wanted)<=np.min(abs(yy[None,None,:]-wanted[:,:,None]),axis=2)+1e-12)
    selected=[]
    for p in sorted(a.snapshots.glob('local_view_*.h5')):
        step=int(p.stem.rsplit('_',1)[1]);j=np.flatnonzero(arr['iterations']==step).item()
        assert sha(p)==records[j]['sha256']
        with h5py.File(p) as h:
            assert int(h.attrs['iteration'])==step and h['Ez'].dtype==np.float64
            expected=np.empty((5,350))
            for k in range(5):
                for i in range(350):expected[k,i]=h['Ez'][i,int(arr['iy'][k,i]),0]
            np.testing.assert_array_equal(expected,field[j])
            assert expected.tobytes()==field[j].tobytes()
        selected.append(dict(iteration=step,sha256=sha(p)))
    assert len(selected)==14
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    ns=arr['time_s']*1e9;extent=[145,180,ns[-1]+(ns[1]-ns[0])/2,-(ns[1]-ns[0])/2]
    for tag,title,norm,ylim in [
        ('native_log','原生Ez时空探针：所有面板共用对数灰度',SymLogNorm(1e-4,vmin=-4000,vmax=4000),[500,0]),
        ('late_linear','晚到原生Ez：所有面板固定±0.05V/m灰度（超限截色）',None,[400,200])]:
        fig,axes=plt.subplots(5,1,figsize=(13,15),sharex=True,layout='constrained')
        for k,ax in enumerate(axes):
            opts=dict(norm=norm) if norm else dict(vmin=-.05,vmax=.05)
            im=ax.imshow(field[:,k,:],aspect='auto',extent=extent,cmap='gray',interpolation='nearest',**opts)
            ax.set_ylim(*ylim);ax.axvline(r['receiver_x_m'],color='tab:red',ls='--',lw=1)
            ax.set_title(r['probe_names'][k]+'；最近快照单元中心，非精确Yee接收点',fontsize=12)
            ax.set_ylabel('原生时间 / ns')
        axes[-1].set_xlabel('模型x / m（测线里程=x+20m）；红线为Rx横坐标')
        fig.colorbar(im,ax=axes,label='Ez / (V/m)，带符号共用尺度',shrink=.65)
        fig.suptitle(title+'\n190m站高损耗H0；已有250帧，无新求解；Ricker含源时延，不能直接当SFCW时标',fontsize=14)
        fig.savefig(a.out/(tag+'.png'),dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(13,5),layout='constrained')
    ax.plot(xx,surface,'k',label='真实体素地表');ax.plot(xx,mudtop,'brown',label='真实覆盖层底')
    for k in range(5):ax.plot(xx,arr['actual_y_m'][k],lw=1,label=r['probe_names'][k])
    ax.axvline(r['receiver_x_m'],ls='--',color='red');ax.set(xlabel='模型x / m',ylabel='模型y / m',title='五条实际取样线；每点取最近已有0.1m快照单元中心')
    ax.legend(ncol=2);fig.savefig(a.out/'probe_geometry.png',dpi=150);plt.close(fig)
    proof=dict(status='PASS14_INDEPENDENT_NATIVE_FRAMES_BITWISE_PROBES',audit_script_sha256=sha(__file__),
        extraction_sha256=sha(a.source/'extraction.json'),probes_sha256=r['probes_sha256'],
        contract_sha256=r['contract_sha256'],verification_sha256=r['verification_sha256'],
        remote_verified_frames=250,locally_independently_read_frames=14,selected=selected,solver_runs=0,
        probe_shape=list(field.shape),max_y_rounding_m=float(abs(arr['actual_y_m']-wanted).max()),
        display=dict(native_log=dict(linthresh=1e-4,vmin=-4000,vmax=4000),late_linear=dict(vmin=-.05,vmax=.05)),
        limits='Fourteen frames independently reread locally; other 236 source hashes checked remotely. Scalar total Ez and coarse sampling do not identify unique path, reflection count or power fractions.')
    (a.out/'independent_audit.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (a.out/'extraction.json').write_bytes((a.source/'extraction.json').read_bytes())
    print(json.dumps({k:proof[k] for k in ['status','remote_verified_frames','locally_independently_read_frames','solver_runs']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','parent','geometry','snapshots','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
